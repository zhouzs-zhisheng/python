#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
cpm_equilibrium_loop.py

恒电势 (CPM) 体系平衡检测循环脚本。

针对单个 MOF/溶剂体系下的多个电压点 (0V/1V/2V/3V/4V)，
循环运行 NVT 续跑 + 检测电极电荷变化率，直到所有电压点都达到平衡
(连续两个 5ns 窗口的平均电极电荷相对变化 < 0.5%)。

每个电压点独立平衡，无顺序依赖；并行 (轮内串行, 轮间循环)。

关键规则：
    - 单轮 NVT 默认 10 ns；首轮强制 50 ns (由 grompp.mdp 控制)
    - 最大循环 10 轮
    - 收敛判据：abs((q_last - q_prev) / q_prev) < 0.005
    - 电荷窗口：每 5000 步求一次平均 (假设 nstxout=1000 → 5 ns)
    - 体相密度：用 gmx density 选 index.ndx 的 group 6，
      z 区间从 system_summary.json 的 middle_vacuum 自动推导
    - 输出：
        * 每个电压点目录下 density.log (轮次记录)
        * 每个电压点目录下 cat.xvg (group 6 密度分布)
        * 体系根目录下 new_equilibrium_result.log (汇总平衡结果)

用法：
    python cpm_equilibrium_loop.py <system_root> [--gmx GMX] [--max-loops 10]

    system_root : 单个 MOF 的体系根目录
                   (其下应有 system_summary.json 和 0V/1V/2V/3V/4V 五个电压点目录)
    --gmx       : gmx 可执行文件路径，默认 gmx
    --max-loops : 最大循环轮次，默认 10
"""

import argparse
import json
import math
import os
import shutil
import subprocess
import sys
from pathlib import Path


# ============================================================
# 常量
# ============================================================

# 电压点子目录名 (相对 system_root)
VOLTAGE_DIRS = ["0V", "1V", "2V", "3V", "4V"]

# 电极电荷数据文件 (每个电压点目录下各一份)
CHARGE_FILE = "CPM_electrodeCharge.dat"

# 密度计算相关
DENSITY_GROUP = 6                 # index.ndx 中体相 density 的 group 编号
CHARGE_INTERVAL_STEPS = 5000       # 电荷求平均的窗口长度 (5 ns @ nstxout=1000)
CHARGE_CONVERGENCE_THRESHOLD = 0.005  # 电荷变化率 < 0.5% 视为收敛

# NVT 续跑
FIRST_LOOP_NS = 50                 # 首轮强制 50 ns
PER_LOOP_NS = 10                   # 后续每轮 10 ns
DEFAULT_MAX_LOOPS = 10             # 最大循环轮次

# 文件名
DENSITY_LOG = "density.log"
EQUILIBRIUM_LOG = "new_equilibrium_result.log"
START_GRO = "start.gro"
NVT_GRO = "nvt.gro"
NVT_TPR = "nvt.tpr"
NVT_CPT = "nvt.cpt"
NVT_XTC = "nvt.xtc"
GROMPP_MDP = "grompp.mdp"
INDEX_NDX = "index.ndx"
CAT_XVG = "cat.xvg"
SYSTEM_SUMMARY = "system_summary.json"


# ============================================================
# 工具函数
# ============================================================

def fail(message, code=1):
    """打印错误并退出。"""
    print(f"错误：{message}")
    sys.exit(code)


def warn(message):
    print(f"警告：{message}")


def run_command(args, cwd=None, stdin_text=None):
    """
    执行外部命令，捕获返回码。
    失败时直接 fail 退出 (CPM 流程中命令失败通常意味着体系出问题，
    不宜继续)。
    """
    print(f"[CMD] {' '.join(args)}")
    try:
        result = subprocess.run(
            args,
            cwd=cwd,
            input=stdin_text,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            check=True,
            text=True,
        )
        # 打印部分输出便于排查
        if result.stdout:
            lines = result.stdout.splitlines()
            tail = lines[-20:] if len(lines) > 20 else lines
            for line in tail:
                print(f"  {line}")
        return result
    except subprocess.CalledProcessError as exc:
        detail = exc.stdout or ""
        fail(
            f"命令执行失败 (返回码 {exc.returncode})：\n"
            f"  {' '.join(args)}\n{detail}"
        )
    except FileNotFoundError:
        fail(f"找不到外部命令：{args[0]}。请确认 gmx 路径正确。")


def load_system_summary(system_root):
    """加载 system_summary.json。"""
    summary_path = Path(system_root) / SYSTEM_SUMMARY
    if not summary_path.is_file():
        fail(f"找不到 {summary_path}")
    with open(summary_path, "r", encoding="utf-8") as f:
        return json.load(f)


def derive_density_region(summary):
    """
    从 system_summary.json 推导体相密度 z 区间。

    规则：取 middle_vacuum (type=vacuum 且非上下边界的那段)，
    取其中间 4 nm 作为体相区间 (与 adjust_density.py 的
    DENSITY_Z_MIN_NM=18.0, DENSITY_Z_MAX_NM=22.0 思路一致)。
    若 middle_vacuum 厚度 < 4 nm，则取整段中点 ± thickness/2。
    """
    regions = summary.get("z_structure_regions_nm", [])

    # 找 middle_vacuum：跳过最前和最后的 vacuum 段，取中间的 vacuum
    vacuum_regions = [
        r for r in regions if r.get("type") == "vacuum"
    ]
    if len(vacuum_regions) < 1:
        fail("system_summary.json 中找不到 vacuum 区域")

    # middle_vacuum 通常是 vacuum_regions 中间的那个
    # (结构：lower_vacuum, bottom_electrode, middle_vacuum,
    #        top_electrode, upper_vacuum)
    if len(vacuum_regions) >= 3:
        middle = vacuum_regions[len(vacuum_regions) // 2]
    else:
        middle = vacuum_regions[0]

    z_low = float(middle["z_low"])
    z_high = float(middle["z_high"])
    thickness = z_high - z_low

    # 取中间 4 nm 作为体相
    bulk_half = 2.0  # 4 nm / 2
    if thickness >= 4.0:
        center = (z_low + z_high) / 2.0
        return (center - bulk_half, center + bulk_half)
    else:
        # 太薄，用整段
        return (z_low, z_high)


def derive_electrode_regions(summary):
    """
    从 system_summary.json 推导电极区间 (用于密度记录，不用于收敛判定)。

    返回 [(low1, high1), (low2, high2)]，每端剔除 1 nm 边界
    (与参考脚本一致)。
    """
    regions = summary.get("z_structure_regions_nm", [])
    electrodes = [
        r for r in regions if r.get("type") == "electrode"
    ]
    if len(electrodes) < 2:
        fail("system_summary.json 中找不到 2 个 electrode 区域")

    result = []
    for e in sorted(electrodes, key=lambda x: float(x["z_low"])):
        z_low = float(e["z_low"]) + 1.0    # 剔除 1 nm 边界
        z_high = float(e["z_high"]) - 1.0
        result.append((z_low, z_high))
    return result


def process_electrode_charge(charge_file, interval):
    """
    读取电极电荷数据文件，每 interval 行求一次平均。

    文件格式：每行一个浮点数 (电极电荷)，# 开头的行跳过。
    返回 [avg_1, avg_2, ...]，每个元素是一个 interval 窗口的平均电荷。
    """
    data = []
    with open(charge_file, "r") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            try:
                # 取第一个数值列 (兼容多列格式)
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


def calc_average_density(xvg_path, z_low, z_high):
    """
    读取 gmx density 输出的 xvg，计算 z 在 [z_low, z_high] 区间内
    的平均密度。

    xvg 格式：前若干行 @/# 注释，数据行两列 (z, density)。
    """
    if not Path(xvg_path).is_file():
        fail(f"密度文件不存在：{xvg_path}")

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
        fail(f"密度文件无有效数据：{xvg_path}")

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
    """
    读取 density.log，返回 (next_loop_index, is_first_loop)。

    若文件不存在或仅有表头，返回 (1, True)。
    """
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
    """追加一行到 density.log。"""
    log_path = Path(log_file)
    # 首次写入时加表头
    if not log_path.is_file():
        with open(log_path, "w") as f:
            f.write("# Loop   bulk_density\n")
    with open(log_path, "a") as f:
        f.write(f"{loop:5d}    {density:10.4f}\n")


def write_equilibrium_log(voltage_dir, loop, avg_charge, density):
    """
    在电压点目录下写 new_equilibrium_result.log (单点收敛记录)。
    并在 system_root 下追加汇总 (所有点收敛时由 main 写最终行)。
    """
    log_path = Path(voltage_dir) / EQUILIBRIUM_LOG
    with open(log_path, "w") as f:
        f.write("reached equilibrium\n")
        f.write(
            f"loop={loop}  avg_charge={avg_charge:.4f}  "
            f"bulk_density={density:.4f}\n"
        )
    print(f"  -> {voltage_dir.name} 平衡记录已写入 {log_path}")


# ============================================================
# 核心流程
# ============================================================

def run_one_voltage(voltage_dir, gmx, params, shared_files_dir):
    """
    对单个电压点执行一轮平衡检测。

    返回 True 表示该电压点已收敛，False 表示未收敛需继续循环。
    """
    voltage_dir = Path(voltage_dir).resolve()
    print(f"\n{'=' * 60}")
    print(f"处理电压点：{voltage_dir.name} ({voltage_dir})")
    print(f"{'=' * 60}")

    # 1. 读取当前轮次
    log_file = voltage_dir / DENSITY_LOG
    loop, is_first = read_last_loop(log_file)

    # 2. 检查电荷文件
    charge_file = voltage_dir / CHARGE_FILE
    if not charge_file.is_file():
        fail(f"{voltage_dir} 缺少电荷文件：{charge_file}")

    # 3. 计算最近两个窗口的平均电荷
    avg_charges = process_electrode_charge(
        charge_file, CHARGE_INTERVAL_STEPS
    )

    if len(avg_charges) >= 2:
        delta = (
            (avg_charges[-1] - avg_charges[-2]) / avg_charges[-2]
        )
        print(
            f"  电荷窗口数={len(avg_charges)}, "
            f"最近两窗口: {avg_charges[-2]:.4f} -> {avg_charges[-1]:.4f}, "
            f"delta={delta*100:.4f}%"
        )
        if abs(delta) < CHARGE_CONVERGENCE_THRESHOLD:
            # 已收敛，不续跑，直接记录并返回
            # 算一次体相密度作为最终记录
            density = calc_average_density(
                voltage_dir / CAT_XVG,
                params["bulk_z_low"], params["bulk_z_high"]
            ) if (voltage_dir / CAT_XVG).is_file() else 0.0
            write_equilibrium_log(
                voltage_dir, loop, avg_charges[-1], density
            )
            return True
    else:
        print(f"  电荷窗口数={len(avg_charges)} (<2)，需继续跑")

    # 4. 未收敛，准备续跑 NVT
    #    首轮用已存在的 start.gro；后续用上一轮 nvt.gro -> start.gro
    start_gro = voltage_dir / START_GRO
    nvt_gro = voltage_dir / NVT_GRO

    if nvt_gro.is_file():
        # 用上一轮末态作为新起点
        shutil.copy(str(nvt_gro), str(start_gro))
        print(f"  cp {NVT_GRO} -> {START_GRO}")
    elif not start_gro.is_file():
        fail(
            f"{voltage_dir} 既无 {NVT_GRO} 也无 {START_GRO}，"
            f"无法启动 NVT"
        )

    # 5. 准备 mdp / index.ndx
    #    共用文件从 shared_files_dir 复制过来
    src_mdp = Path(shared_files_dir) / GROMPP_MDP
    src_ndx = Path(shared_files_dir) / INDEX_NDX
    if not src_mdp.is_file():
        fail(f"共用的 {GROMPP_MDP} 不存在：{src_mdp}")
    if not src_ndx.is_file():
        fail(f"共用的 {INDEX_NDX} 不存在：{src_ndx}")

    dst_mdp = voltage_dir / GROMPP_MDP
    dst_ndx = voltage_dir / INDEX_NDX
    shutil.copy(str(src_mdp), str(dst_mdp))
    shutil.copy(str(src_ndx), str(dst_ndx))

    # 6. grompp
    #    首轮 50ns (mdp 中已设置)，后续轮 10ns
    #    简化处理：要求用户在 mdp 中设 nsteps=50000000 (50 ns)
    #    后续轮通过 -nsteps 覆盖为 10 ns
    #    实际更稳妥的做法是修改 mdp 的 nsteps，但为避免改文件，
    #    这里用 grompp 不带 -nsteps，由 mdp 决定。
    #    若 mdp 已是 10 ns，首轮也是 10 ns (不满足 50ns 要求)
    #    -> 首轮需要单独的 50ns mdp。
    #
    #    决策：要求用户准备 grompp_50ns.mdp (首轮) 和 grompp.mdp (后续)
    #    首轮用 50ns 版，后续用 10ns 版。
    mdp_50ns = Path(shared_files_dir) / "grompp_50ns.mdp"
    mdp_10ns = Path(shared_files_dir) / GROMPP_MDP

    if is_first and mdp_50ns.is_file():
        # 首轮强制 50 ns
        use_mdp = mdp_50ns
        print(f"  首轮使用 50ns mdp：{use_mdp}")
    else:
        use_mdp = mdp_10ns
        print(f"  使用 10ns mdp：{use_mdp}")

    shutil.copy(str(use_mdp), str(dst_mdp))

    # grompp
    run_command(
        [gmx, "grompp", "-f", GROMPP_MDP, "-c", START_GRO,
         "-o", NVT_TPR, "-n", INDEX_NDX, "-maxwarn", "1"],
        cwd=str(voltage_dir),
    )

    # 7. mdrun
    #    续跑：若存在 nvt.cpt 则 -cpi -append，否则从头跑
    cpt = voltage_dir / NVT_CPT
    if cpt.is_file():
        mdrun_cmd = [
            gmx, "mdrun", "-deffnm", "nvt",
            "-ntmpi", "1", "-ntomp", "32",
            "-tunepme", "no", "-v", "-pin", "on", "-nstlist", "20",
            "-cpi", NVT_CPT, "-append",
        ]
        print("  续跑模式 (-cpi -append)")
    else:
        mdrun_cmd = [
            gmx, "mdrun", "-deffnm", "nvt",
            "-ntmpi", "1", "-ntomp", "32",
            "-tunepme", "no", "-v", "-pin", "on", "-nstlist", "20",
        ]
        print("  首轮模式 (无 cpt)")

    run_command(mdrun_cmd, cwd=str(voltage_dir))

    # 8. 计算体相密度 (group 6)
    box_z_total = params["box_z_total"]
    sl = math.ceil(box_z_total / 0.01)

    run_command(
        [gmx, "density", "-f", NVT_XTC, "-s", "nvt",
         "-sl", str(sl), "-o", CAT_XVG],
        cwd=str(voltage_dir),
        stdin_text=f"{DENSITY_GROUP}\n",
    )

    # 9. 算体相区间平均密度
    density = calc_average_density(
        voltage_dir / CAT_XVG,
        params["bulk_z_low"], params["bulk_z_high"]
    )
    print(f"  体相密度 (z={params['bulk_z_low']:.2f}~"
          f"{params['bulk_z_high']:.2f} nm) = {density:.4f}")

    # 10. 重新读电荷，计算 delta
    avg_charges = process_electrode_charge(
        charge_file, CHARGE_INTERVAL_STEPS
    )
    if len(avg_charges) >= 2:
        delta = (
            (avg_charges[-1] - avg_charges[-2]) / avg_charges[-2]
        )
        print(
            f"  续跑后电荷窗口数={len(avg_charges)}, "
            f"delta={delta*100:.4f}%"
        )
    else:
        delta = 1.0  # 强制未收敛
        warn("续跑后电荷窗口仍 <2，无法判定收敛")

    # 11. 写日志
    append_density_log(log_file, loop, density)

    # 12. 判收敛
    if abs(delta) < CHARGE_CONVERGENCE_THRESHOLD:
        write_equilibrium_log(
            voltage_dir, loop, avg_charges[-1], density
        )
        return True
    else:
        # 未收敛，nvt.gro 已生成，下轮 cp -> start.gro
        print(f"  {voltage_dir.name} 第 {loop} 轮未收敛 (delta="
              f"{delta*100:.4f}%)，下轮继续")
        return False


# ============================================================
# 主流程
# ============================================================

def main():
    parser = argparse.ArgumentParser(
        description="CPM 体系平衡检测循环 (多电压点并行)"
    )
    parser.add_argument(
        "system_root",
        help="单个 MOF 的体系根目录 "
             "(其下应有 system_summary.json 和 0V/1V/2V/3V/4V)",
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
        help=f"最大循环轮次，默认 {DEFAULT_MAX_LOOPS}",
    )
    args = parser.parse_args()

    system_root = Path(args.system_root).resolve()
    if not system_root.is_dir():
        fail(f"体系根目录不存在：{system_root}")

    print("=" * 72)
    print("CPM EQUILIBRIUM LOOP")
    print("=" * 72)
    print(f"System root : {system_root}")
    print(f"GMX         : {args.gmx}")
    print(f"Max loops   : {args.max_loops}")
    print(f"Voltage pts : {VOLTAGE_DIRS}")

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

    # 2. 确定电压点目录
    voltage_dirs = []
    for name in VOLTAGE_DIRS:
        vd = system_root / name
        if not vd.is_dir():
            warn(f"电压点目录不存在，跳过：{vd}")
            continue
        voltage_dirs.append(vd)

    if not voltage_dirs:
        fail(f"在 {system_root} 下找不到任何电压点目录 {VOLTAGE_DIRS}")

    print(f"\n将处理 {len(voltage_dirs)} 个电压点：")
    for vd in voltage_dirs:
        print(f"  - {vd}")

    # 3. 检查共用文件 (grompp.mdp / index.ndx)
    shared_files_dir = system_root
    for fname in [GROMPP_MDP, INDEX_NDX]:
        p = shared_files_dir / fname
        if not p.is_file():
            fail(f"共用文件不存在：{p}")

    # 4. 并行循环 (轮内串行, 轮间循环)
    pending = list(voltage_dirs)
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
                vd, args.gmx, params, shared_files_dir
            )
            if not converged:
                next_pending.append(vd)

        pending = next_pending
        print(
            f"\n第 {loop_count} 轮完成，"
            f"{len(voltage_dirs) - len(pending)}/{len(voltage_dirs)} "
            f"个电压点已收敛，{len(pending)} 个待续"
        )

    # 5. 汇总结果
    print(f"\n{'=' * 72}")
    print("FINAL STATUS")
    print(f"{'=' * 72}")

    all_converged = len(pending) == 0
    summary_path = system_root / EQUILIBRIUM_LOG

    with open(summary_path, "w") as f:
        if all_converged:
            f.write("ALL VOLTAGE POINTS REACHED EQUILIBRIUM\n\n")
            f.write(
                f"{'Voltage':<10} {'Loops':<8} {'Avg_Charge':<14} "
                f"{'Bulk_Density':<14}\n"
            )
            for vd in voltage_dirs:
                elog = vd / EQUILIBRIUM_LOG
                if elog.is_file():
                    content = elog.read_text(encoding="utf-8")
                    # 解析 loop / avg_charge / density
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
            f.write(f"未收敛电压点 ({len(pending)} 个):\n")
            for vd in pending:
                f.write(f"  - {vd.name}\n")
            print(f"仍有 {len(pending)} 个电压点未收敛，请检查：")
            for vd in pending:
                print(f"  - {vd.name}")
            print(f"汇总结果已写入：{summary_path}")
            print(
                "可手动继续跑："
                f"python {sys.argv[0]} {system_root} --gmx {args.gmx}"
            )

    print(f"\n完成。")


if __name__ == "__main__":
    main()
