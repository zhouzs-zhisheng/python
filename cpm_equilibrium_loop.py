#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
cpm_equilibrium_loop.py

恒电势 (CPM) 体系平衡检测循环脚本 (0V 优先 + 多电压点版本)。

流程：
    Phase 0: 检测 0V 是否已完成 (有 new_equilibrium_result.log)
    Phase 1: 若 0V 未完成，准备 0V (从 fine 取 nvt.gro + index.ndx)
             并循环跑 0V 直到电荷收敛
    Phase 2: 0V 完成后，用 0V 的 nve.gro 作为 1V/2V/3V/4V 的输入结构，
             准备各电压点 (ln -s allMatrixA.bin + cp CPM_ControlFile.dat)
    Phase 3: 循环跑 1V/2V/3V/4V 直到全部电荷收敛

续跑机制：
    - 首轮: grompp 生成 nve.tpr (50ns)，mdrun -deffnm nve
    - 续跑: mv nve.tpr → nve_loop{N-1}.tpr，
            convert-tpr -extend 10000 (延长 10ns) 生成新 nve.tpr，
            mdrun -deffnm nve -cpi nve.cpt -append (轨迹追加到原 xtc)

收敛判据：连续两个 5ns 窗口的平均电极电荷相对变化 < 0.5%。

用法：
    # 完整流程 (0V 优先 + 所有电压点)
    python cpm_equilibrium_loop.py <system_root> [--gmx GMX] [--max-loops 10]

    # 单电压点模式 (只跑指定电压点)
    python cpm_equilibrium_loop.py <system_root> --voltage 1V [--gmx GMX]

    # 强制重新开始 (清空既有 nve.* / density.log，从首轮 50ns 重跑)
    python cpm_equilibrium_loop.py <system_root> --mode restart [--voltage 2V]

    # 指定续跑 (用 tpr 归档名 + density.log 双重判定已完成轮次，
    #           读取轨迹总时长作密度窗口起点，并回填缺失的 density.log)
    python cpm_equilibrium_loop.py <system_root> --mode continue [--voltage 2V]

    system_root : ACN 目录 (其下有 fine/, system_summary.json)
    --voltage   : 只跑指定电压点 (0V/1V/2V/3V/4V)，不进入其他电压点
    --gmx       : gmx 可执行文件路径，默认 gmx
    --max-loops : 每个电压点的最大循环轮次，默认 10
    --mode      : auto=自动判定(默认) / restart=强制重新开始 / continue=指定续跑

目录结构要求：
    qmof-xxx/
    ├── ACN/                      # = system_root
    │   ├── fine/
    │   │   ├── first/nvt.gro     # 0V 输入结构来源
    │   │   ├── index.ndx         # 所有电压点共用
    │   │   └── topol.top         # 拓扑 (所有电压点共用)
    │   ├── system_summary.json
    │   ├── grompp.mdp            # 10ns NVE 参数
    │   ├── grompp_50ns.mdp       # 50ns NVE 参数 (首轮)
    │   ├── 0V/                   # 工作目录 (脚本自动创建)
    │   ├── 1V/
    │   ├── 2V/
    │   ├── 3V/
    │   └── 4V/
    ├── 0V/                       # CPM_ControlFile.dat 源
    │   └── CPM_ControlFile.dat
    ├── 1V/
    │   └── CPM_ControlFile.dat
    ├── ...
    └── allMatrixA.bin
"""

# ============================================================
# 更新记录
#   2026-08-31 : 首轮 NVE 时长改回 50ns，使用 grompp_50ns.mdp
#                (FIRST_NVE_PS=50000)；density 窗口 begin 随之调整。
#                注：此前临时改为 30ns 已作废。
# ============================================================

import argparse
import json
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path


# ============================================================
# 常量
# ============================================================

VOLTAGE_DIRS = ["0V", "1V", "2V", "3V", "4V"]
ZERO_V = "0V"

CHARGE_FILE = "CPM_electrodeCharge.dat"
CONTROL_FILE = "CPM_ControlFile.dat"
MATRIX_FILE = "allMatrixA.bin"

DENSITY_GROUP = 6
CHARGE_INTERVAL_STEPS = 5000
CHARGE_CONVERGENCE_THRESHOLD = 0.005
# 0V 收敛采用绝对判据：连续两电荷窗口平均值的绝对差 < 此值(e) 即收敛。
# 因为 0V 的电荷均值≈0，相对判据分母→0 会失效。
CHARGE_ABS_CONV_TOL = 0.5

DEFAULT_MAX_LOOPS = 10

# 运行模式开关
#   auto      : 自动判定——优先用 tpr 归档名 + density.log 的磁盘事实推断轮次，
#               据此决定续跑还是首轮 (推荐默认，可修复旧代码崩溃导致的误判)。
#   restart   : 强制重新开始。清空本电压点已有 nve.*/nve_loop*/density.log 产物，
#               从首轮 (grompp 50ns) 重新跑。
#   continue  : 因特殊原因强制续跑。用 tpr 归档名 + density.log 双重判定已完成轮次，
#               并读取 nve.xtc 实测总时长作为密度窗口起点，回填缺失的 density.log。
MODE_AUTO = "auto"
MODE_RESTART = "restart"
MODE_CONTINUE = "continue"

DENSITY_LOG = "density.log"
EQUILIBRIUM_LOG = "new_equilibrium_result.log"

# NVE 相关文件名 (从 NVT 改为 NVE)
START_GRO = "start.gro"
NVE_GRO = "nve.gro"
NVE_TPR = "nve.tpr"
NVE_CPT = "nve.cpt"
NVE_XTC = "nve.xtc"
NVE_EDR = "nve.edr"
NVE_LOG = "nve.log"

GROMPP_MDP = "grompp.mdp"
GROMPP_50NS_MDP = "grompp_50ns.mdp"
INDEX_NDX = "index.ndx"
TOPOL_TOP = "topol.top"
DENSITY_XVG = "density.xvg"
SYSTEM_SUMMARY = "system_summary.json"

# 续跑延长的时间 (ps)
EXTEND_PS = 10000  # 10 ns = 10000 ps

# ---------- 密度统计参数 ----------
# 首轮 NVE 时长 (ps)：默认 50 ns。若回退到 grompp.mdp (10ns) 使用，
# 请相应改为 10000。
FIRST_NVE_PS = 50000  # 50 ns
# density 固定切片数（不再按盒子高度动态计算，避免多余计算）
DENSITY_SL = 1500
# 每轮密度只统计整条累计轨迹的“最后这段窗口”（ns）
DENSITY_WINDOW_NS = 5
DENSITY_WINDOW_PS = DENSITY_WINDOW_NS * 1000


# ============================================================
# 工具函数
# ============================================================

def fail(message, code=1):
    print(f"错误：{message}")
    sys.exit(code)


def warn(message):
    print(f"警告：{message}")


def build_gmx_env(gmx_path):
    """
    根据 gmx 可执行文件的位置，构造 subprocess 运行时使用的环境。

    背景：用户交互式 shell 中 (通过 bashrc / module load / GMXRC) 配好了
    LD_LIBRARY_PATH 以找到 libomp.so / libgromacs.so 等动态库；但 Python
    subprocess 启动的是非登录非交互 shell，不一定继承这些设置。

    策略：
      1. 传入的 gmx_path 如果含目录 (非裸 "gmx")，则：
           标准布局    <prefix>/bin/gmx    → 查 <prefix>/{lib,lib64}
           嵌套布局    <prefix>/bin_all/bin/gmx
                                         → 查 <prefix>/bin_all/{lib,lib64}
                                           并 查 <prefix>/{lib,lib64}
         将找到的目录注入子进程 LD_LIBRARY_PATH 最前端。
      2. 继承 os.environ，这样用户在父 shell 里已设的 LD_LIBRARY_PATH
         依然有效，不会丢失。
    """
    env = os.environ.copy()
    if gmx_path is None:
        return env
    gmx_p = Path(gmx_path)
    # 含目录的路径 (相对或绝对)：推断 lib 目录
    if "/" in str(gmx_p) or "\\" in str(gmx_p):
        try:
            gmx_abs = gmx_p.resolve()
        except (OSError, RuntimeError):
            return env
        # <prefix>/bin/gmx        → bin_dir = <prefix>/bin
        # <prefix>/bin_all/bin/gmx → bin_dir = <prefix>/bin_all/bin
        bin_dir = gmx_abs.parent
        parents_to_try = [bin_dir.parent]               # <prefix> 或 <prefix>/bin_all
        parents_to_try.append(parents_to_try[0].parent) # <prefix>/.. 或 <prefix>
        extra_libs = []
        seen = set()
        for p in parents_to_try:
            for libname in ("lib", "lib64"):
                candidate = p / libname
                if candidate.is_dir() and str(candidate) not in seen:
                    extra_libs.append(str(candidate))
                    seen.add(str(candidate))
        if extra_libs:
            existing = env.get("LD_LIBRARY_PATH", "")
            prepend = ":".join(extra_libs)
            if existing:
                env["LD_LIBRARY_PATH"] = f"{prepend}:{existing}"
            else:
                env["LD_LIBRARY_PATH"] = prepend
    return env


def run_command(args, cwd=None, stdin_text=None, env=None,
                tolerate_artifacts=None):
    """
    运行外部命令。

    env=None  时：subprocess 自动继承当前 Python 进程环境 (最常见情况)。
    env=<dict>时：直接使用该字典作为子进程完整环境。

    tolerate_artifacts: 可选路径列表 (相对 cwd 或绝对路径)。当命令返回非零时，
        若这些产物文件全部已生成(存在且非空)，则视为“运行已完成、仅是收尾
        报错”(如 GROMACS 的 CUDA teardown)，忽略返回码放行，而不是硬失败；
        若任一产物缺失则仍按失败处理。

    用法：所有 gmx 调用都应传入 build_gmx_env() 返回的 env，
         确保 LD_LIBRARY_PATH 被正确注入以便找到 libomp.so 等。
    """
    print(f"[CMD] {' '.join(args)}")
    if env and env.get("LD_LIBRARY_PATH"):
        injected_prefix = env["LD_LIBRARY_PATH"].split(os.pathsep)[0]
        print(f"[ENV] LD_LIBRARY_PATH prepend: {injected_prefix}")

    stdin_pipe = subprocess.PIPE if stdin_text is not None else None
    try:
        proc = subprocess.Popen(
            args,
            cwd=cwd,
            stdin=stdin_pipe,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            bufsize=1,   # 行缓冲，配合逐行回显实现实时输出
            env=env,
        )
    except FileNotFoundError:
        fail(f"找不到外部命令：{args[0]}。请确认 gmx 路径正确（使用 --gmx 绝对路径推荐）。")

    # 先喂 stdin（gmx density 需传入组选择行），再关闭，避免管道写满死锁
    if stdin_text is not None:
        proc.stdin.write(stdin_text)
        proc.stdin.close()

    buffered = []
    # 逐行实时回显到终端，同时累加留作失败诊断（长任务如 mdrun -v 进度实时可见）
    for line in proc.stdout:
        buffered.append(line)
        print(f"  {line}", end="", flush=True)

    proc.wait()
    if proc.returncode != 0:
        # 软失败：若指定了「完成后产物清单」且全部已生成，说明模拟本体已完成，
        # 仅是非零收尾报错(如 GROMACS CUDA teardown 的 cudaFreeHost)，放行继续。
        if tolerate_artifacts:
            base = Path(cwd) if cwd else Path.cwd()
            all_paths = [
                (name if isinstance(name, Path) and name.is_absolute()
                 else base / name)
                for name in tolerate_artifacts
            ]
            present = [p for p in all_paths
                       if p.is_file() and p.stat().st_size > 0]
            if len(present) == len(all_paths):
                warn(
                    "命令返回码非零，但指定产物均已生成，视为完成后的收尾"
                    "报错（如 CUDA teardown），放行继续"
                )
                return proc
        detail = "".join(buffered)[-4000:]
        fail(
            f"命令执行失败 (返回码 {proc.returncode})：\n"
            f"  {' '.join(args)}\n{detail}"
        )
    return proc


def enable_line_buffering():
    """
    非 TTY 重定向下（任务系统/nohup/管道）让 print 逐行落盘。

    Python 在 stdout 非终端时改用 ~8KB 块缓冲：未带 flush=True 的 print
    会攒到缓冲满或进程退出才一次性写出，导致 job 日志里看不到程序自身的
    实时输出。这里把 stdout 切到行缓冲，使"每遇到一个换行即落盘"。
    调用发生在 __main__ 启动处，仅作用于以脚本运行时。
    """
    try:
        sys.stdout.reconfigure(line_buffering=True)
        if sys.stderr is not None and hasattr(sys.stderr, "reconfigure"):
            sys.stderr.reconfigure(line_buffering=True)
    except Exception:
        # stdout/stderr 可能被替换为不支持 reconfigure 的自定义流，忽略保持默认
        pass


def load_system_summary(system_root):
    summary_path = Path(system_root) / SYSTEM_SUMMARY
    if not summary_path.is_file():
        fail(f"找不到 {summary_path}")
    with open(summary_path, "r", encoding="utf-8") as f:
        return json.load(f)


def derive_density_region(summary):
    """
    从 system_summary.json 推导体相密度 z 区间。
    取 middle_vacuum 中间 4 nm 作为体相区间。
    """
    regions = summary.get("z_structure_regions_nm", [])
    vacuum_regions = [r for r in regions if r.get("type") == "vacuum"]
    if len(vacuum_regions) < 1:
        fail("system_summary.json 中找不到 vacuum 区域")

    if len(vacuum_regions) >= 3:
        middle = vacuum_regions[len(vacuum_regions) // 2]
    else:
        middle = vacuum_regions[0]

    z_low = float(middle["z_low"])
    z_high = float(middle["z_high"])
    thickness = z_high - z_low

    bulk_half = 2.0
    if thickness >= 4.0:
        center = (z_low + z_high) / 2.0
        return (center - bulk_half, center + bulk_half)
    else:
        return (z_low, z_high)


def process_electrode_charge(charge_file, interval):
    """
    读取电极电荷数据文件，每 interval 行求一次平均。
    文件格式：每行两列为 (负极电荷, 正极电荷)，开头可有带 # 的注释行。
    当前 0V 收敛只统计“负极”列，故取每行第 1 列；# 开头和空行跳过。
    """
    data = []
    with open(charge_file, "r") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            try:
                first_value = float(line.split()[0])
                data.append(first_value)
            except (ValueError, IndexError):
                continue

    if not data:
        return []

    means = []
    for i in range(0, len(data), interval):
        chunk = data[i:i + interval]
        if len(chunk) == 0:
            continue
        means.append(sum(chunk) / len(chunk))
    return means


def _charge_windows_converged(avg_charges, is_zero_voltage):
    """
    判断连续两个电荷窗口的平均值是否收敛。
    返回 (converged, delta)：
      - 0V  ：delta = c_cur - c_prev (绝对差)，
              |delta| < CHARGE_ABS_CONV_TOL (0.5) 视为收敛。
      - 其他：delta = (c_cur - c_prev) / c_prev (相对差)，
              |delta| < CHARGE_CONVERGENCE_THRESHOLD 视为收敛。
    窗口不足或(非0V)分母为0时返回 (False, None)。
    """
    if len(avg_charges) < 2:
        return False, None
    c_prev, c_cur = avg_charges[-2], avg_charges[-1]
    if is_zero_voltage:
        delta = c_cur - c_prev
        return abs(delta) < CHARGE_ABS_CONV_TOL, delta
    if c_prev == 0:
        return False, None
    delta = (c_cur - c_prev) / c_prev
    return abs(delta) < CHARGE_CONVERGENCE_THRESHOLD, delta


def calc_average_density(xvg_path, z_low, z_high):
    if not Path(xvg_path).is_file():
        warn(f"密度文件不存在：{xvg_path}")
        return 0.0

    z_coords = []
    densities = []
    with open(xvg_path, "r") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith(("#", "@", "&")):
                continue
            parts = line.split()
            if len(parts) < 2:
                continue
            try:
                z = float(parts[0])
                d = float(parts[1])
            except ValueError:
                continue
            z_coords.append(z)
            densities.append(d)

    if not z_coords:
        warn(f"密度文件无有效数据：{xvg_path}")
        return 0.0

    total = 0.0
    count = 0
    for z, d in zip(z_coords, densities):
        if z_low <= z <= z_high:
            total += d
            count += 1

    if count == 0:
        warn(
            f"区间 [{z_low}, {z_high}] 内无密度数据点 "
            f"(z 范围: {z_coords[0]:.3f}~{z_coords[-1]:.3f})"
        )
        return 0.0

    return total / count


def read_last_loop(log_file):
    log_path = Path(log_file)
    if not log_path.is_file():
        return 1, True

    last_loop = 0
    has_data = False
    with open(log_path, "r") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) < 1:
                continue
            try:
                last_loop = int(parts[0])
                has_data = True
            except ValueError:
                continue

    if not has_data:
        return 1, True
    return last_loop + 1, False


def append_density_log(log_file, loop, density):
    log_path = Path(log_file)
    if not log_path.is_file():
        with open(log_path, "w") as f:
            f.write("# Loop   bulk_density\n")
    with open(log_path, "a") as f:
        f.write(f"{loop:5d}    {density:10.4f}\n")


def existing_loop_numbers(log_file):
    """返回 density.log 里已记录的 loop 编号集合。"""
    log_path = Path(log_file)
    loops = set()
    if not log_path.is_file():
        return loops
    with open(log_path, "r") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if not parts:
                continue
            try:
                loops.add(int(parts[0]))
            except ValueError:
                continue
    return loops


def count_archived_loops(voltage_dir):
    """
    从 nve_loop{N}.tpr 归档文件名反推整条模拟已经完成到第几轮 (磁盘事实)。

    规则：首轮 (round1) 结束时 nve.tpr 常驻；第 2 轮起把上一轮 tpr 归档为
    nve_loop1.tpr、nve_loop2.tpr ……。因此：
      - 存在 nve.tpr   -> 至少完成 1 轮
      - 归档最大序号 N -> 又额外完成 N 轮
    已完成轮数 = max_arch + (nve.tpr 存在 ? 1 : 0)。

    该值与 density.log (记录流口径) 互为印证，用于双重判定。当旧代码崩溃导致
    density.log 缺失但仍留有 nve.tpr 时，此函数能据磁盘事实判断「应续跑而非首轮」。
    """
    voltage_dir = Path(voltage_dir)
    max_arch = 0
    for p in voltage_dir.glob("nve_loop*.tpr"):
        try:
            n = int(p.stem[len("nve_loop"):])
        except ValueError:
            continue
        max_arch = max(max_arch, n)
    if max_arch == 0 and not (voltage_dir / NVE_TPR).is_file():
        return 0
    return max_arch + 1 if (voltage_dir / NVE_TPR).is_file() else max_arch


def resolve_loop_state(voltage_dir, log_file, mode):
    """
    双重判定「已完成轮次 completed / 本轮 loop / 是否首轮」。

    - restart : 强制 completed=0，从首轮开始。
    - auto/continue : completed = max(磁盘事实 arch, 记录流 log)，
      ——取较高者，避免把已跑过的体系（如旧代码崩溃留下 nve.tpr 但无 density.log）
      误判为首轮而重新 grompp 覆盖掉既有轨迹。
    返回 (completed, next_loop, is_first)。
    """
    voltage_dir = Path(voltage_dir)
    arch_completed = count_archived_loops(voltage_dir)
    log_next, _ = read_last_loop(log_file)
    log_completed = max(log_next - 1, 0)  # 下一轮-1

    if mode == MODE_RESTART:
        completed = 0
    else:
        completed = max(arch_completed, log_completed)

    return completed, completed + 1, completed == 0


def _run_capture(args, cwd, env, stdin_text=None):
    """后台静默运行并捕获全部输出 (不回显)，用于读取 gmx check 等返回信息。失败返回 None。"""
    try:
        proc = subprocess.run(
            args, cwd=cwd, input=stdin_text,
            capture_output=True, text=True, env=env,
        )
    except FileNotFoundError:
        return None
    if proc.returncode != 0:
        return None
    return (proc.stdout or "") + (proc.stderr or "")


def _last_time_ps(output):
    """从 gmx check 输出解析最后一个 'time <X> ps' 作为轨迹最后时刻 (ps)。"""
    if not output:
        return None
    matches = re.findall(
        r"[Tt]ime\s+([0-9]+(?:\.[0-9]+)?(?:[eE][-+]?[0-9]+)?)\s+ps", output
    )
    if not matches:
        return None
    try:
        return float(matches[-1])
    except ValueError:
        return None


def read_traj_total_time(voltage_dir, gmx, gmx_env):
    """
    读取累计轨迹 nve.xtc 的实测总时长 (ps)。

    用 `gmx check -f nve.xtc` 解析 "Last frame ... time <X> ps"。
    轨迹为续跑累积 (+append)，因此这是判断"这次模拟真正开始了多久"的最本质依据。
    读取失败返回 None，调用方回退到按轮次外推的时间口径。
    """
    xtc = Path(voltage_dir) / NVE_XTC
    if not xtc.is_file():
        return None
    out = _run_capture([gmx, "check", "-f", NVE_XTC],
                       cwd=str(voltage_dir), env=gmx_env)
    return _last_time_ps(out)


def reset_round(voltage_dir):
    """
    restart 模式下清空本电压点既有模拟产物，从首轮重新开始。
    保留 start.gro / index.ndx / topol.top / 矩阵 等输入准备产物。
    """
    voltage_dir = Path(voltage_dir)
    patterns = (
        "nve.gro", "nve.xtc", "nve.cpt", "nve.tpr", "nve.edr", "nve.log",
        "nve_loop*.tpr", DENSITY_LOG, DENSITY_XVG,
    )
    for pat in patterns:
        for p in voltage_dir.glob(pat):
            p.unlink()
            print(f"  [restart] 移除 {p.name}")
    elog = voltage_dir / EQUILIBRIUM_LOG
    if elog.is_file():
        elog.unlink()
        print(f"  [restart] 移除 {elog.name}")


def backfill_density_log(voltage_dir, up_to_loop, gmx, gmx_env, params):
    """
    回填 density.log 中缺失的轮次行 (1..up_to_loop 中未记录的 loop)。

    对缺失轮次 k：用 gmx density -b <begin> -e <begin+窗口> 限定窗口补算，
    再追加写入 density.log，保证 loop 编号连续，read_last_loop 下次能读到正确
    的下一轮编号（否则崩溃后日志缺行会让下一轮误判为首轮）。
    """
    voltage_dir = Path(voltage_dir)
    if up_to_loop < 1:
        return
    log_file = voltage_dir / DENSITY_LOG
    done = existing_loop_numbers(log_file)
    missing = [k for k in range(1, up_to_loop + 1) if k not in done]
    if not missing:
        return

    xtc = voltage_dir / NVE_XTC
    if not xtc.is_file():
        warn(
            f"density.log 缺 {len(missing)} 行，但缺少 {NVE_XTC}，"
            f"无法回填，跳过"
        )
        return

    for k in missing:
        begin = FIRST_NVE_PS + (k - 1) * EXTEND_PS - DENSITY_WINDOW_PS
        end = begin + DENSITY_WINDOW_PS
        print(f"  [回填] 补算密度 loop={k}  (window {begin:.0f}~{end:.0f} ps)")
        run_command(
            [gmx, "density", "-f", NVE_XTC, "-s", NVE_TPR, "-n", INDEX_NDX,
             "-sl", str(DENSITY_SL), "-d", "Z",
             "-b", str(begin), "-e", str(end),
             "-o", DENSITY_XVG],
            cwd=str(voltage_dir),
            stdin_text=f"{DENSITY_GROUP}\n",
            env=gmx_env,
        )
        density = calc_average_density(
            voltage_dir / DENSITY_XVG,
            params["bulk_z_low"], params["bulk_z_high"]
        )
        append_density_log(log_file, k, density)
        print(f"  [回填] loop={k}  bulk_density={density:.4f}")


def write_equilibrium_log(voltage_dir, loop, avg_charge, density):
    log_path = Path(voltage_dir) / EQUILIBRIUM_LOG
    with open(log_path, "w") as f:
        f.write("reached equilibrium\n")
        f.write(
            f"loop={loop}  avg_charge={avg_charge:.4f}  "
            f"bulk_density={density:.4f}\n"
        )
    print(f"  -> {Path(voltage_dir).name} 平衡记录已写入 {log_path}")


def is_voltage_converged(voltage_dir):
    """检查电压点是否已收敛 (有 new_equilibrium_result.log)。"""
    return (Path(voltage_dir) / EQUILIBRIUM_LOG).is_file()


# ============================================================
# 准备工作：符号链接 + 控制文件 + 输入结构
# ============================================================

def prepare_voltage_files(voltage_dir, voltage_name, voltage_root):
    """
    在电压点工作目录内准备 allMatrixA.bin 符号链接和 CPM_ControlFile.dat。

    使用绝对路径，不依赖 ../../ 相对路径：
      - allMatrixA.bin 从 voltage_root/allMatrixA.bin 链接
      - CPM_ControlFile.dat 从 voltage_root.parent/<电压值>/CPM_ControlFile.dat 复制
        (对应原命令 cp ../../<电压值>/CPM_ControlFile.dat .)
    """
    voltage_dir = Path(voltage_dir).resolve()

    # 1. allMatrixA.bin 符号链接
    matrix_src = Path(voltage_root).resolve() / MATRIX_FILE
    link_path = voltage_dir / MATRIX_FILE

    if not matrix_src.is_file():
        fail(f"allMatrixA.bin 不存在：{matrix_src}")

    if not link_path.exists():
        os.symlink(str(matrix_src), str(link_path))
        print(f"  创建符号链接: {link_path} -> {matrix_src}")
    else:
        if link_path.is_symlink():
            target = os.readlink(str(link_path))
            if not Path(target).resolve().samefile(matrix_src):
                link_path.unlink()
                os.symlink(str(matrix_src), str(link_path))
                print(f"  更新符号链接: {link_path} -> {matrix_src}")
            else:
                print(f"  allMatrixA.bin 符号链接已存在，跳过")
        else:
            print(f"  allMatrixA.bin 已存在 (非符号链接)，跳过")

    # 2. CPM_ControlFile.dat
    #    源在 voltage_root/<电压值>/CPM_ControlFile.dat
    #    (对应原命令 cp ../../<电压值>/CPM_ControlFile.dat .
    #     其中工作目录为 ACN/<电压值>，所以 voltage_root = <mof>/)
    ctrl_src = Path(voltage_root).resolve() / voltage_name / CONTROL_FILE
    ctrl_dst = voltage_dir / CONTROL_FILE

    if ctrl_dst.is_file():
        print(f"  {CONTROL_FILE} 已存在，跳过复制")
    elif ctrl_src.is_file():
        shutil.copy(str(ctrl_src), str(ctrl_dst))
        print(f"  复制控制文件: {ctrl_src} -> {ctrl_dst}")
    else:
        fail(
            f"找不到 {CONTROL_FILE} 源文件：{ctrl_src}\n"
            f"请确认 {voltage_root.parent}/{voltage_name}/ 下有 {CONTROL_FILE}"
        )


def prepare_input_structure(voltage_dir, source_gro, source_ndx, source_topol):
    """
    准备电压点的输入结构 start.gro、index.ndx 和 topol.top。

    - 0V: 从 fine/first/nvt.gro 复制为 start.gro
    - 1V/2V/3V/4V: 从 0V/nve.gro 复制为 start.gro
    - index.ndx: 从 fine/index.ndx 复制 (所有电压点共用)
    - topol.top: 从 fine/topol.top 复制 (所有电压点共用，拓扑不随电压变化)
    """
    voltage_dir = Path(voltage_dir).resolve()

    # start.gro (输入结构)
    dst_gro = voltage_dir / START_GRO
    if not dst_gro.is_file():
        if not Path(source_gro).is_file():
            fail(f"输入结构源文件不存在：{source_gro}")
        shutil.copy(str(source_gro), str(dst_gro))
        print(f"  复制输入结构: {source_gro} -> {dst_gro}")
    else:
        print(f"  {START_GRO} 已存在，跳过复制")

    # index.ndx
    dst_ndx = voltage_dir / INDEX_NDX
    if not dst_ndx.is_file():
        if not Path(source_ndx).is_file():
            fail(f"index.ndx 源文件不存在：{source_ndx}")
        shutil.copy(str(source_ndx), str(dst_ndx))
        print(f"  复制 index.ndx: {source_ndx} -> {dst_ndx}")
    else:
        print(f"  {INDEX_NDX} 已存在，跳过复制")

    # topol.top (grompp 默认 -p topol.top，必须在 cwd 中存在)
    dst_topol = voltage_dir / TOPOL_TOP
    if not dst_topol.is_file():
        if not Path(source_topol).is_file():
            fail(f"topol.top 源文件不存在：{source_topol}")
        shutil.copy(str(source_topol), str(dst_topol))
        print(f"  复制 topol.top: {source_topol} -> {dst_topol}")
    else:
        print(f"  {TOPOL_TOP} 已存在，跳过复制")
# ============================================================
# 核心：单电压点单轮运行 (支持 convert-tpr 续跑)
# ============================================================

def run_one_voltage(voltage_dir, gmx, params, shared_files_dir,
                    use_gpu=False, gmx_env=None):
    """
    对单个电压点执行一轮平衡检测。
    返回 True 表示已收敛，False 表示未收敛需继续循环。

    续跑逻辑：
        Loop 1 (首次):
            grompp -c start.gro -o nve.tpr (50ns)
            mdrun -deffnm nve (无 cpt)
        Loop 2+ (续跑):
            mv nve.tpr → nve_loop{loop-1}.tpr
            convert-tpr -s nve_loop{loop-1}.tpr -extend 10000 -o nve.tpr
            mdrun -deffnm nve -cpi nve.cpt -append (轨迹追加)

    use_gpu=True 时，mdrun 加 -nb gpu -pme gpu -pmefft gpu。
    """
    voltage_dir = Path(voltage_dir).resolve()
    print(f"\n{'=' * 60}")
    print(f"处理电压点：{voltage_dir.name} ({voltage_dir})")
    print(f"{'=' * 60}")

    if use_gpu:
        print("  mdrun 模式: GPU 加速 (-nb gpu -pme gpu -pmefft gpu)")
    else:
        print("  mdrun 模式: 纯 CPU")

    # 1. 双重判定当前轮次 (或按 restart 清空从首轮开始)
    log_file = voltage_dir / DENSITY_LOG
    completed, loop, is_first = resolve_loop_state(
        voltage_dir, log_file, args_mode
    )
    if args_mode == MODE_RESTART:
        reset_round(voltage_dir)
        completed, loop, is_first = 0, 1, True
    print(f"  运行模式   : {args_mode}")
    print(f"  已完成轮次 : {completed}, 本轮为 loop={loop}, is_first={is_first}")

    # 2. 续跑前先检查是否已收敛
    #    注意：首轮 / 文件不存在时跳过，mdrun 会首次生成电荷文件；
    #    只有续跑 (非首轮) 时缺少电荷文件才视为异常。
    charge_file = voltage_dir / CHARGE_FILE
    if charge_file.is_file():
        avg_charges = process_electrode_charge(
            charge_file, CHARGE_INTERVAL_STEPS
        )
        is_zero_voltage = voltage_dir.name == ZERO_V
        converged, delta = _charge_windows_converged(
            avg_charges, is_zero_voltage
        )
        if delta is not None:
            if is_zero_voltage:
                print(
                    f"  [续跑前] 0V 电荷窗口数={len(avg_charges)}, "
                    f"最近两窗口: {avg_charges[-2]:.4f} -> "
                    f"{avg_charges[-1]:.4f}, |Δcharge|={abs(delta):.4f} e "
                    f"(阈<{CHARGE_ABS_CONV_TOL})"
                )
            else:
                print(
                    f"  [续跑前] 电荷窗口数={len(avg_charges)}, "
                    f"最近两窗口: {avg_charges[-2]:.4f} -> "
                    f"{avg_charges[-1]:.4f}, delta={delta*100:.4f}%"
                )
            if converged:
                density = calc_average_density(
                    voltage_dir / DENSITY_XVG,
                    params["bulk_z_low"], params["bulk_z_high"]
                ) if (voltage_dir / DENSITY_XVG).is_file() else 0.0
                write_equilibrium_log(
                    voltage_dir, loop, avg_charges[-1], density
                )
                return True
        else:
            print(f"  [续跑前] 电荷窗口数={len(avg_charges)} (<2)，需继续跑")
    elif is_first:
        print(f"  [续跑前] 首轮，电荷文件尚未生成，将在 mdrun 结束后检查")
    else:
        fail(
            f"{voltage_dir} 缺少电荷文件：{charge_file}\n"
            f"(非首轮模拟要求存在电荷文件以续跑)"
        )

    # 4. 执行 NVE 模拟 (首次或续跑)
    start_gro = voltage_dir / START_GRO
    nve_gro = voltage_dir / NVE_GRO
    nve_tpr = voltage_dir / NVE_TPR
    nve_cpt = voltage_dir / NVE_CPT

    if is_first:
        # ---- 首轮: grompp + mdrun ----
        print(f"\n  --- 首轮 NVE (50ns) ---")

        if not start_gro.is_file():
            fail(
                f"{voltage_dir} 缺少 {START_GRO}，"
                f"无法启动首轮 NVE"
            )

        # 选择 mdp (优先 50ns，回退 10ns)
        mdp_50ns = Path(shared_files_dir) / GROMPP_50NS_MDP
        mdp_10ns = Path(shared_files_dir) / GROMPP_MDP

        if mdp_50ns.is_file():
            use_mdp = mdp_50ns
            print(f"  首轮使用 50ns mdp：{use_mdp}")
        elif mdp_10ns.is_file():
            use_mdp = mdp_10ns
            warn(f"50ns mdp 不存在，首轮使用 10ns mdp：{use_mdp}")
        else:
            fail(
                f"找不到 mdp 文件：{mdp_10ns} 或 {mdp_50ns}\n"
                f"请在 {shared_files_dir} 下准备 {GROMPP_MDP} 和 {GROMPP_50NS_MDP}"
            )

        # 保持原文件名复制到工作目录 (不重命名为 grompp.mdp)。
        # 这样 grompp 命令行的 -f 参数直接反映是 50ns 还是 10ns，
        # 避免用户从日志看到 -f grompp.mdp 误以为用的是 10ns。
        mdp_basename = use_mdp.name
        dst_mdp = voltage_dir / mdp_basename
        shutil.copy(str(use_mdp), str(dst_mdp))
        print(f"  复制 mdp -> {dst_mdp}")

        # grompp (-f 使用与源一致的文件名，保留 _50ns 后缀语义)
        run_command(
            [gmx, "grompp", "-f", mdp_basename, "-c", START_GRO,
             "-o", NVE_TPR, "-n", INDEX_NDX, "-maxwarn", "6"],
            cwd=str(voltage_dir),
            env=gmx_env,
        )

        # mdrun (首轮，无 cpt)
        mdrun_cmd = [
            gmx, "mdrun", "-deffnm", "nve",
            "-ntmpi", "1", "-ntomp", "32",
            "-tunepme", "no", "-v", "-pin", "on",
        ]
        if use_gpu:
            mdrun_cmd += ["-nb", "gpu", "-pme", "gpu", "-pmefft", "gpu"]
        print("  首轮模式 (无 cpt)")
        run_command(
            mdrun_cmd, cwd=str(voltage_dir), env=gmx_env,
            tolerate_artifacts=[voltage_dir / NVE_GRO, voltage_dir / NVE_XTC],
        )

    else:
        # ---- 续跑: convert-tpr + mdrun -cpi -append ----
        print(f"\n  --- 续跑 NVE (延长 {EXTEND_PS/1000}ns) ---")

        if not nve_tpr.is_file():
            fail(f"{voltage_dir} 缺少 {NVE_TPR}，无法续跑")
        if not nve_cpt.is_file():
            fail(f"{voltage_dir} 缺少 {NVE_CPT}，无法续跑")

        # 1. 重命名上一轮的 nve.tpr → nve_loop{loop-1}.tpr
        archived_tpr = voltage_dir / f"nve_loop{loop-1}.tpr"
        if nve_tpr.is_file():
            if archived_tpr.is_file():
                # 已存在同名归档，先删掉 (或加后缀)
                archived_tpr.unlink()
            shutil.move(str(nve_tpr), str(archived_tpr))
            print(f"  mv {NVE_TPR} -> {archived_tpr.name}")

        # 2. convert-tpr 延长模拟时间
        #    -extend <ps> 在原 tpr 末尾时间基础上延长
        print(f"  convert-tpr -extend {EXTEND_PS} (延长 {EXTEND_PS/1000}ns)")
        run_command(
            [gmx, "convert-tpr",
             "-s", str(archived_tpr),
             "-extend", str(EXTEND_PS),
             "-o", NVE_TPR],
            cwd=str(voltage_dir),
            env=gmx_env,
        )

        # 3. mdrun 续跑 (-cpi -append，轨迹追加到原 nve.xtc)
        mdrun_cmd = [
            gmx, "mdrun", "-deffnm", "nve",
            "-ntmpi", "1", "-ntomp", "32",
            "-tunepme", "no", "-v", "-pin", "on",
            "-cpi", NVE_CPT, "-append",
        ]
        if use_gpu:
            mdrun_cmd += ["-nb", "gpu", "-pme", "gpu", "-pmefft", "gpu"]
        print("  续跑模式 (-cpi nve.cpt -append)")
        run_command(
            mdrun_cmd, cwd=str(voltage_dir), env=gmx_env,
            tolerate_artifacts=[voltage_dir / NVE_GRO, voltage_dir / NVE_XTC],
        )

    # 5. 计算体相密度 (group 6)
    #    - 固定切片数 DENSITY_SL (不再按盒子高度动态计算)
    #    - 方向固定 Z (-d Z)
    #    - -s 明确指向当前轮的 tpr 文件 nve.tpr：续跑时旧 tpr 会先改名归档为
    #      nve_loop{N}.tpr，并用 convert-tpr 重新生成 nve.tpr，故当前轮恒为该名
    #    - -b 取整条累计轨迹的最后 DENSITY_WINDOW_NS ns：
    #        累计总时长 = FIRST_NVE_PS + (loop-1)*EXTEND_PS
    #        begin = 总时长 - 窗口 (首轮 50ns 时 begin=45000，每续跑10ns后移)
    sl = DENSITY_SL
    total_ps = FIRST_NVE_PS + (loop - 1) * EXTEND_PS

    # continue 模式下，优先读轨迹实测总时长作为密度窗口起点——
    # 它反映"这次模拟真正累计跑了多久"，比按轮次外推更贴近实际
    # (例如上一轮中途崩溃只跑了一部分时)。读取失败则回退到外推值。
    traj_total = None
    if args_mode == MODE_CONTINUE:
        traj_total = read_traj_total_time(voltage_dir, gmx, gmx_env)
        if traj_total is not None:
            print(f"  [continue] 轨迹实测总时长 = {traj_total:.1f} ps"
                  f" (按轮次外推 = {total_ps:.0f} ps)")
            total_ps = traj_total
        else:
            warn("读取轨迹实测总时长失败，回退到按轮次外推的时间口径")
    begin_ps = total_ps - DENSITY_WINDOW_PS

    nve_xtc = voltage_dir / NVE_XTC
    if not nve_xtc.is_file():
        fail(f"{voltage_dir} 缺少 {NVE_XTC}，无法计算密度")

    run_command(
        [gmx, "density", "-f", NVE_XTC, "-s", NVE_TPR, "-n", INDEX_NDX,
         "-sl", str(sl), "-d", "Z", "-b", str(begin_ps), "-o", DENSITY_XVG],
        cwd=str(voltage_dir),
        stdin_text=f"{DENSITY_GROUP}\n",
        env=gmx_env,
    )

    # 6. 算体相区间平均密度
    density = calc_average_density(
        voltage_dir / DENSITY_XVG,
        params["bulk_z_low"], params["bulk_z_high"]
    )
    print(f"  体相密度 (z={params['bulk_z_low']:.2f}~"
          f"{params['bulk_z_high']:.2f} nm) = {density:.4f}")

    # 7. 重新读电荷，计算 delta
    avg_charges = process_electrode_charge(
        charge_file, CHARGE_INTERVAL_STEPS
    )
    is_zero_voltage = voltage_dir.name == ZERO_V
    converged, delta = _charge_windows_converged(
        avg_charges, is_zero_voltage
    )

    # 8. 写日志 (当前轮) + 回填缺失的旧轮次
    append_density_log(log_file, loop, density)
    backfill_density_log(voltage_dir, loop - 1, gmx, gmx_env, params)

    # 9. 判收敛
    if delta is None:
        warn("续跑后电荷窗口仍 <2，无法判定收敛，下轮继续")
        return False
    if is_zero_voltage:
        print(f"  [续跑后] 0V 电荷窗口数={len(avg_charges)}, "
              f"|Δcharge|={abs(delta):.4f} e (阈<{CHARGE_ABS_CONV_TOL})")
    else:
        print(f"  [续跑后] 电荷窗口数={len(avg_charges)}, "
              f"delta={delta*100:.4f}%")
    if converged:
        write_equilibrium_log(
            voltage_dir, loop, avg_charges[-1], density
        )
        return True
    else:
        if is_zero_voltage:
            msg = f"|Δcharge|={abs(delta):.4f} e"
        else:
            msg = f"delta={delta*100:.4f}%"
        print(f"  {voltage_dir.name} 第 {loop} 轮未收敛 ({msg})，下轮继续")
        return False


# ============================================================
# 单电压点运行
# ============================================================

def run_single_voltage(system_root, voltage_name, gmx, params,
                       voltage_root, fine_gro, fine_ndx, fine_topol,
                       use_gpu=False, gmx_env=None):
    """只跑指定电压点，不进入其他电压点。"""
    vdir = system_root / voltage_name
    vdir.mkdir(parents=True, exist_ok=True)

    print(f"\n{'=' * 72}")
    print(f"单电压点模式：{voltage_name}")
    print(f"{'=' * 72}")

    # 检查是否已收敛
    if is_voltage_converged(vdir):
        print(f"{voltage_name} 已收敛 (存在 {EQUILIBRIUM_LOG})")
        print(f"如需重跑，请删除 {vdir / EQUILIBRIUM_LOG}")
        return

    # 准备输入结构
    print(f"\n--- 准备 {voltage_name} 输入结构 ---")
    if voltage_name == ZERO_V:
        # 0V: 从 fine 取
        prepare_input_structure(vdir, fine_gro, fine_ndx, fine_topol)
    else:
        # 1V-4V: 从 0V 的 nve.gro 取 (拓扑仍是 fine/topol.top)
        zero_v_gro = system_root / ZERO_V / NVE_GRO
        if not zero_v_gro.is_file():
            fail(
                f"0V 的 {NVE_GRO} 不存在：{zero_v_gro}\n"
                f"请先跑 0V：python {sys.argv[0]} {system_root} --voltage 0V --gmx {gmx}"
            )
        prepare_input_structure(vdir, zero_v_gro, fine_ndx, fine_topol)

    # 准备矩阵和控制文件
    print(f"\n--- 准备 {voltage_name} 矩阵和控制文件 ---")
    prepare_voltage_files(vdir, voltage_name, voltage_root)

    # 循环跑
    print(f"\n--- 运行 {voltage_name} 平衡循环 ---")
    loop_count = 0
    while loop_count < args_max_loops:
        loop_count += 1
        print(f"\n{'=' * 60}")
        print(f"{voltage_name} LOOP {loop_count} / {args_max_loops}")
        print(f"{'=' * 60}")

        converged = run_one_voltage(
            vdir, gmx, params, system_root,
            use_gpu=use_gpu, gmx_env=gmx_env,
        )
        if converged:
            break
    else:
        warn(
            f"{voltage_name} 在 {args_max_loops} 轮后仍未收敛，"
            f"请人工检查 {vdir}"
        )

    # 写单点汇总
    if is_voltage_converged(vdir):
        summary_path = system_root / EQUILIBRIUM_LOG
        elog = vdir / EQUILIBRIUM_LOG
        if elog.is_file():
            content = elog.read_text(encoding="utf-8")
            with open(summary_path, "w") as f:
                f.write(f"SINGLE VOLTAGE POINT: {voltage_name}\n\n")
                f.write(content)
            print(f"汇总结果已写入：{summary_path}")


# 全局变量 (由 main 设置，供 run_single_voltage 使用)
args_max_loops = DEFAULT_MAX_LOOPS
args_mode = MODE_AUTO


# ============================================================
# 主流程
# ============================================================

def main():
    global args_max_loops
    global args_mode

    parser = argparse.ArgumentParser(
        description="CPM 体系平衡检测循环 (0V 优先 + 多电压点 / 单电压点)"
    )
    parser.add_argument(
        "system_root",
        help="ACN 目录 (其下有 fine/, system_summary.json)",
    )
    parser.add_argument(
        "--voltage",
        choices=VOLTAGE_DIRS,
        default=None,
        help="只跑指定电压点 (0V/1V/2V/3V/4V)，不进入其他电压点",
    )
    parser.add_argument(
        "--gmx",
        default="gmx",
        help="gmx 可执行文件路径，默认 gmx",
    )
    parser.add_argument(
        "--max-loops",
        type=int,
        default=DEFAULT_MAX_LOOPS,
        help=f"每个电压点最大循环轮次，默认 {DEFAULT_MAX_LOOPS}",
    )
    parser.add_argument(
        "--gpu",
        action="store_true",
        default=False,
        help="启用 GPU 加速 mdrun：-nb gpu -pme gpu -pmefft gpu "
             "(默认关闭，纯 CPU 运行)",
    )
    parser.add_argument(
        "--mode",
        choices=[MODE_AUTO, MODE_RESTART, MODE_CONTINUE],
        default=MODE_AUTO,
        help="运行模式：auto=自动判定(默认)，restart=强制重新开始，"
             "continue=指定续跑(用 tpr 归档+density.log 双重判定，"
             "读取轨迹总时长并回填 density.log)",
    )
    args = parser.parse_args()

    args_max_loops = args.max_loops
    args_mode = args.mode

    # 一次性构造 gmx 子进程环境 (含 LD_LIBRARY_PATH 注入)
    gmx_env = build_gmx_env(args.gmx)

    system_root = Path(args.system_root).resolve()
    if not system_root.is_dir():
        fail(f"体系根目录不存在：{system_root}")

    # 自动检测 voltage_root
    # 优先级：
    #   1. candidate 下有 allMatrixA.bin (最可靠标记，不能是 ACN/ 工作目录)
    #   2. candidate 下 <电压值>/CPM_ControlFile.dat (例如 0V/CPM_ControlFile.dat)
    #      存在 (控制文件按 <voltage_root>/<V>/CPM_ControlFile.dat 组织)
    # 注意：不能用 candidate/0V 子目录是否存在判定，
    #       因为 ACN/0V 是脚本自动创建的工作目录，
    #       第一次重跑会误判 voltage_root=ACN。
    voltage_root = None
    candidates = [system_root, system_root.parent]
    # 先按 allMatrixA.bin 判定
    for candidate in candidates:
        if (candidate / MATRIX_FILE).is_file():
            voltage_root = candidate
            break
    # 再按控制文件路径判定 (allMatrixA.bin 没找到时)
    if voltage_root is None:
        for candidate in candidates:
            if (candidate / ZERO_V / CONTROL_FILE).is_file():
                voltage_root = candidate
                break

    if voltage_root is None:
        # 都找不到就默认 system_root.parent，避免后续报误导性错误
        voltage_root = system_root.parent

    print("=" * 72)
    print("CPM EQUILIBRIUM LOOP (NVE + convert-tpr 续跑)")
    print("=" * 72)
    print(f"System root  : {system_root}")
    print(f"Voltage root : {voltage_root}")
    print(f"GMX          : {args.gmx}")
    if gmx_env and gmx_env.get("LD_LIBRARY_PATH"):
        injected = gmx_env["LD_LIBRARY_PATH"].split(os.pathsep)[0]
        print(f"GMX lib path : {injected}  (LD_LIBRARY_PATH prepend)")
    else:
        print("GMX lib path : (未注入，使用进程继承环境)")
    print(f"Max loops    : {args.max_loops}")
    print(f"Run mode     : {args.mode}")
    if args.voltage:
        print(f"Mode         : 单电压点 ({args.voltage})")
    else:
        print(f"Mode         : 完整流程 (0V 优先 + 多电压点)")
        print(f"Voltage pts  : {VOLTAGE_DIRS}")

    # 1. 加载 system_summary，推导参数
    summary = load_system_summary(system_root)
    bulk_z_low, bulk_z_high = derive_density_region(summary)
    box_z_total = float(summary["box_dimensions_nm"]["z"])

    params = {
        "bulk_z_low": bulk_z_low,
        "bulk_z_high": bulk_z_high,
        "box_z_total": box_z_total,
    }

    print(
        f"\n体相密度区间 : z={bulk_z_low:.3f}~{bulk_z_high:.3f} nm"
    )
    print(f"盒子总长 z   : {box_z_total:.3f} nm")
    print(
        f"收敛判据     : |delta_charge| < "
        f"{CHARGE_CONVERGENCE_THRESHOLD*100:.1f}%"
    )
    print(f"续跑延长     : {EXTEND_PS/1000}ns/轮 (convert-tpr -extend)")

    # 2. 定位输入结构来源
    fine_dir = system_root / "fine"
    fine_gro = fine_dir / "first" / "nvt.gro"
    fine_ndx = fine_dir / INDEX_NDX
    fine_topol = fine_dir / TOPOL_TOP

    # 0V 模式或完整流程都需要 fine 结构。
    # topol.top 所有电压点共用，只要启动任务（包括单电压 1V-4V 模式）都要前置检查。
    if not fine_topol.is_file():
        fail(f"topol.top 源文件不存在：{fine_topol}\n请确认 fine 目录准备完整")

    if args.voltage == ZERO_V or args.voltage is None:
        if not fine_gro.is_file():
            fail(f"0V 输入结构不存在：{fine_gro}\n请确认 fine 收敛流程已完成")
        if not fine_ndx.is_file():
            fail(f"index.ndx 源文件不存在：{fine_ndx}")

    if args.voltage == ZERO_V:
        print(f"\n0V 输入结构来源 : {fine_gro}")
        print(f"index.ndx 来源  : {fine_ndx}")
        print(f"topol.top 来源  : {fine_topol}")

    # 3. 检查共用 mdp 文件
    mdp_10ns = system_root / GROMPP_MDP
    mdp_50ns = system_root / GROMPP_50NS_MDP
    if not mdp_10ns.is_file():
        fail(f"共用 mdp 不存在：{mdp_10ns}")
    if not mdp_50ns.is_file():
        warn(f"50ns mdp 不存在：{mdp_50ns}，首轮将使用 10ns mdp")

    # ============================================================
    # 单电压点模式
    # ============================================================
    if args.voltage:
        run_single_voltage(
            system_root, args.voltage, args.gmx, params,
            voltage_root, fine_gro, fine_ndx, fine_topol,
            use_gpu=args.gpu, gmx_env=gmx_env,
        )
        print(f"\n完成。")
        return

    # ============================================================
    # 完整流程 (0V 优先 + 多电压点)
    # ============================================================

    print(f"\n0V 输入结构来源 : {fine_gro}")
    print(f"index.ndx 来源  : {fine_ndx}")
    print(f"topol.top 来源  : {fine_topol}")

    # ---- Phase 0+1: 检测并运行 0V ----
    print(f"\n{'#' * 72}")
    print("# Phase 0: 检测 0V 是否已完成")
    print(f"{'#' * 72}")

    zero_v_dir = system_root / ZERO_V

    if is_voltage_converged(zero_v_dir):
        print(f"0V 已收敛 (存在 {EQUILIBRIUM_LOG})，跳过 0V，直接进入其他电压点")
    else:
        print(f"0V 未完成，开始准备并运行 0V")

        zero_v_dir.mkdir(parents=True, exist_ok=True)

        print(f"\n--- 准备 0V 输入结构 ---")
        prepare_input_structure(zero_v_dir, fine_gro, fine_ndx, fine_topol)

        print(f"\n--- 准备 0V 矩阵和控制文件 ---")
        prepare_voltage_files(zero_v_dir, ZERO_V, voltage_root)

        print(f"\n--- 运行 0V 平衡循环 ---")
        loop_count = 0
        while loop_count < args.max_loops:
            loop_count += 1
            print(f"\n{'=' * 60}")
            print(f"0V LOOP {loop_count} / {args.max_loops}")
            print(f"{'=' * 60}")

            converged = run_one_voltage(
                zero_v_dir, args.gmx, params, system_root,
                use_gpu=args.gpu, gmx_env=gmx_env,
            )
            if converged:
                break
        else:
            fail(
                f"0V 在 {args.max_loops} 轮后仍未收敛，"
                f"请人工检查 {zero_v_dir}"
            )

        print(f"\n0V 已收敛！")

    # ---- Phase 2: 准备 1V/2V/3V/4V ----
    print(f"\n{'#' * 72}")
    print("# Phase 2: 准备 1V/2V/3V/4V")
    print(f"{'#' * 72}")

    zero_v_gro = zero_v_dir / NVE_GRO
    if not zero_v_gro.is_file():
        fail(f"0V 的 {NVE_GRO} 不存在：{zero_v_gro}")

    other_voltages = [v for v in VOLTAGE_DIRS if v != ZERO_V]

    for vname in other_voltages:
        vdir = system_root / vname
        vdir.mkdir(parents=True, exist_ok=True)

        print(f"\n--- 准备 {vname} ---")
        prepare_input_structure(vdir, zero_v_gro, fine_ndx, fine_topol)
        prepare_voltage_files(vdir, vname, voltage_root)

    # ---- Phase 3: 循环跑 1V/2V/3V/4V ----
    print(f"\n{'#' * 72}")
    print("# Phase 3: 运行 1V/2V/3V/4V 平衡循环")
    print(f"{'#' * 72}")

    pending = []
    for vname in other_voltages:
        vdir = system_root / vname
        if is_voltage_converged(vdir):
            print(f"{vname} 已收敛，跳过")
        else:
            pending.append(vdir)

    loop_count = 0
    while pending and loop_count < args.max_loops:
        loop_count += 1
        print(f"\n{'#' * 72}")
        print(f"# GLOBAL LOOP {loop_count} / {args.max_loops}")
        print(f"# 待处理电压点: {[vd.name for vd in pending]}")
        print(f"{'#' * 72}")

        next_pending = []
        for vd in pending:
            converged = run_one_voltage(
                vd, args.gmx, params, system_root,
                use_gpu=args.gpu, gmx_env=gmx_env,
            )
            if not converged:
                next_pending.append(vd)

        pending = next_pending
        converged_count = len(other_voltages) - len(pending)
        print(
            f"\n第 {loop_count} 轮完成，"
            f"{converged_count}/{len(other_voltages)} "
            f"个电压点已收敛，{len(pending)} 个待续"
        )

    # ---- 汇总结果 ----
    print(f"\n{'=' * 72}")
    print("FINAL STATUS")
    print(f"{'=' * 72}")

    all_voltage_dirs = [system_root / v for v in VOLTAGE_DIRS]
    all_converged = all(
        is_voltage_converged(vd) for vd in all_voltage_dirs
    )

    summary_path = system_root / EQUILIBRIUM_LOG

    with open(summary_path, "w") as f:
        if all_converged:
            f.write("ALL VOLTAGE POINTS REACHED EQUILIBRIUM\n\n")
            f.write(
                f"{'Voltage':<10} {'Loops':<8} {'Avg_Charge':<14} "
                f"{'Bulk_Density':<14}\n"
            )
            for vd in all_voltage_dirs:
                elog = vd / EQUILIBRIUM_LOG
                if elog.is_file():
                    content = elog.read_text(encoding="utf-8")
                    loops = "?"
                    charge = "?"
                    density = "?"
                    for line in content.splitlines():
                        if line.startswith("loop="):
                            for part in line.split():
                                if part.startswith("loop="):
                                    loops = part.split("=")[1]
                                elif part.startswith("avg_charge="):
                                    charge = part.split("=")[1]
                                elif part.startswith("bulk_density="):
                                    density = part.split("=")[1]
                            break
                    f.write(
                        f"{vd.name:<10} {loops:<8} {charge:<14} "
                        f"{density:<14}\n"
                    )
            print(f"所有电压点均已收敛！")
            print(f"汇总结果已写入：{summary_path}")
        else:
            f.write("EQUILIBRIUM NOT REACHED FOR ALL POINTS\n\n")
            not_converged = [
                vd.name for vd in all_voltage_dirs
                if not is_voltage_converged(vd)
            ]
            f.write(f"未收敛电压点 ({len(not_converged)} 个):\n")
            for name in not_converged:
                f.write(f"  - {name}\n")
            print(f"仍有 {len(not_converged)} 个电压点未收敛：{not_converged}")
            print(f"汇总结果已写入：{summary_path}")
            print(
                "可手动继续跑："
                f"python {sys.argv[0]} {system_root} --gmx {args.gmx}"
            )

    print(f"\n完成。")


if __name__ == "__main__":
    enable_line_buffering()
    main()
