#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
cpm_equilibrium_loop.py

恒电势 (CPM) 体系平衡检测循环脚本 (0V 优先 + 多电压点版本)。

流程：
    Phase 0: 检测 0V 是否已完成 (有 new_equilibrium_result.log)
    Phase 1: 若 0V 未完成，准备 0V (从 fine 取 nvt.gro + index.ndx)
             并循环跑 0V 直到电荷收敛
    Phase 2: 0V 完成后，用 0V 的 nvt.gro 作为 1V/2V/3V/4V 的输入结构，
             准备各电压点 (ln -s allMatrixA.bin + cp CPM_ControlFile.dat)
    Phase 3: 循环跑 1V/2V/3V/4V 直到全部电荷收敛

收敛判据：连续两个 5ns 窗口的平均电极电荷相对变化 < 0.5%。
每个电压点独立平衡，无顺序依赖 (但都必须等 0V 完成后才能开始)。

用法：
    python cpm_equilibrium_loop.py <system_root> [--gmx GMX] [--max-loops 10]

    system_root : ACN 目录 (其下有 fine/, system_summary.json)
    --gmx       : gmx 可执行文件路径，默认 gmx
    --max-loops : 每个电压点的最大循环轮次，默认 10

目录结构要求：
    qmof-xxx/
    ├── ACN/                      # = system_root
    │   ├── fine/
    │   │   ├── first/nvt.gro     # 0V 输入结构来源
    │   │   └── index.ndx
    │   ├── system_summary.json
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

VOLTAGE_DIRS = ["0V", "1V", "2V", "3V", "4V"]
ZERO_V = "0V"

CHARGE_FILE = "CPM_electrodeCharge.dat"
CONTROL_FILE = "CPM_ControlFile.dat"
MATRIX_FILE = "allMatrixA.bin"

DENSITY_GROUP = 6
CHARGE_INTERVAL_STEPS = 5000
CHARGE_CONVERGENCE_THRESHOLD = 0.005

DEFAULT_MAX_LOOPS = 10

DENSITY_LOG = "density.log"
EQUILIBRIUM_LOG = "new_equilibrium_result.log"
START_GRO = "start.gro"
NVT_GRO = "nvt.gro"
NVT_TPR = "nvt.tpr"
NVT_CPT = "nvt.cpt"
NVT_XTC = "nvt.xtc"
GROMPP_MDP = "grompp.mdp"
GROMPP_50NS_MDP = "grompp_50ns.mdp"
INDEX_NDX = "index.ndx"
CAT_XVG = "cat.xvg"
SYSTEM_SUMMARY = "system_summary.json"


# ============================================================
# 工具函数
# ============================================================

def fail(message, code=1):
    print(f"错误：{message}")
    sys.exit(code)


def warn(message):
    print(f"警告：{message}")


def run_command(args, cwd=None, stdin_text=None):
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
    文件格式：每行一个浮点数，# 开头跳过。
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
        # 用绝对路径创建符号链接
        os.symlink(str(matrix_src), str(link_path))
        print(f"  创建符号链接: {link_path} -> {matrix_src}")
    else:
        # 检查现有链接是否指向正确
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
    #    源在 voltage_root.parent/<电压值>/CPM_ControlFile.dat
    #    (对应原命令 ../../<电压值>/CPM_ControlFile.dat)
    ctrl_src = Path(voltage_root).resolve().parent / voltage_name / CONTROL_FILE
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


def prepare_input_structure(voltage_dir, source_gro, source_ndx):
    """
    准备电压点的输入结构 nvt.gro 和 index.ndx。

    - 0V: 从 fine/first/nvt.gro 复制
    - 1V/2V/3V/4V: 从 0V/nvt.gro 复制
    - index.ndx: 从 fine/index.ndx 复制 (所有电压点共用)
    """
    voltage_dir = Path(voltage_dir).resolve()

    # nvt.gro
    dst_gro = voltage_dir / NVT_GRO
    if not dst_gro.is_file():
        if not Path(source_gro).is_file():
            fail(f"输入结构源文件不存在：{source_gro}")
        shutil.copy(str(source_gro), str(dst_gro))
        print(f"  复制输入结构: {source_gro} -> {dst_gro}")
    else:
        print(f"  {NVT_GRO} 已存在，跳过复制")

    # index.ndx
    dst_ndx = voltage_dir / INDEX_NDX
    if not dst_ndx.is_file():
        if not Path(source_ndx).is_file():
            fail(f"index.ndx 源文件不存在：{source_ndx}")
        shutil.copy(str(source_ndx), str(dst_ndx))
        print(f"  复制 index.ndx: {source_ndx} -> {dst_ndx}")
    else:
        print(f"  {INDEX_NDX} 已存在，跳过复制")


# ============================================================
# 核心：单电压点单轮运行
# ============================================================

def run_one_voltage(voltage_dir, gmx, params, shared_files_dir):
    """
    对单个电压点执行一轮平衡检测。
    返回 True 表示已收敛，False 表示未收敛需继续循环。
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

    # 4. 准备续跑
    start_gro = voltage_dir / START_GRO
    nvt_gro = voltage_dir / NVT_GRO

    if nvt_gro.is_file():
        shutil.copy(str(nvt_gro), str(start_gro))
        print(f"  cp {NVT_GRO} -> {START_GRO}")
    elif not start_gro.is_file():
        fail(
            f"{voltage_dir} 既无 {NVT_GRO} 也无 {START_GRO}，"
            f"无法启动 NVT"
        )

    # 5. 准备 mdp (首轮 50ns，后续 10ns)
    mdp_50ns = Path(shared_files_dir) / GROMPP_50NS_MDP
    mdp_10ns = Path(shared_files_dir) / GROMPP_MDP

    if is_first and mdp_50ns.is_file():
        use_mdp = mdp_50ns
        print(f"  首轮使用 50ns mdp：{use_mdp}")
    elif mdp_10ns.is_file():
        use_mdp = mdp_10ns
        print(f"  使用 10ns mdp：{use_mdp}")
    else:
        fail(
            f"找不到 mdp 文件：{mdp_10ns} 或 {mdp_50ns}\n"
            f"请在 {shared_files_dir} 下准备 {GROMPP_MDP} 和 {GROMPP_50NS_MDP}"
        )

    dst_mdp = voltage_dir / GROMPP_MDP
    shutil.copy(str(use_mdp), str(dst_mdp))

    # 6. grompp
    run_command(
        [gmx, "grompp", "-f", GROMPP_MDP, "-c", START_GRO,
         "-o", NVT_TPR, "-n", INDEX_NDX, "-maxwarn", "1"],
        cwd=str(voltage_dir),
    )

    # 7. mdrun (续跑用 -cpi -append)
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
        delta = 1.0
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
        print(f"  {voltage_dir.name} 第 {loop} 轮未收敛 (delta="
              f"{delta*100:.4f}%)，下轮继续")
        return False


# ============================================================
# 主流程
# ============================================================

def main():
    parser = argparse.ArgumentParser(
        description="CPM 体系平衡检测循环 (0V 优先 + 多电压点)"
    )
    parser.add_argument(
        "system_root",
        help="ACN 目录 (其下有 fine/, system_summary.json)",
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
    args = parser.parse_args()

    system_root = Path(args.system_root).resolve()
    if not system_root.is_dir():
        fail(f"体系根目录不存在：{system_root}")

    # 自动检测 voltage_root：有 0V/ 和 allMatrixA.bin 的目录
    # 优先 system_root 本身，其次 system_root.parent
    voltage_root = None
    for candidate in [system_root, system_root.parent]:
        if (candidate / ZERO_V).is_dir() or (
            candidate / MATRIX_FILE
        ).exists():
            voltage_root = candidate
            break

    if voltage_root is None:
        # 都没有，默认用 system_root，让后续检查报具体错误
        voltage_root = system_root

    # 但电压点工作目录在 system_root 下
    # (ACN/0V, ACN/1V, ...)
    # allMatrixA.bin 和 CPM_ControlFile.dat 源在 voltage_root 下

    print("=" * 72)
    print("CPM EQUILIBRIUM LOOP (0V 优先 + 多电压点)")
    print("=" * 72)
    print(f"System root  : {system_root}")
    print(f"Voltage root : {voltage_root}")
    print(f"GMX          : {args.gmx}")
    print(f"Max loops    : {args.max_loops}")
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

    # 2. 定位输入结构来源
    fine_dir = system_root / "fine"
    fine_gro = fine_dir / "first" / NVT_GRO
    fine_ndx = fine_dir / INDEX_NDX

    if not fine_gro.is_file():
        fail(f"0V 输入结构不存在：{fine_gro}\n请确认 fine 收敛流程已完成")
    if not fine_ndx.is_file():
        fail(f"index.ndx 源文件不存在：{fine_ndx}")

    print(f"\n0V 输入结构来源 : {fine_gro}")
    print(f"index.ndx 来源  : {fine_ndx}")

    # 3. 检查共用 mdp 文件
    mdp_10ns = system_root / GROMPP_MDP
    mdp_50ns = system_root / GROMPP_50NS_MDP
    if not mdp_10ns.is_file():
        fail(f"共用 mdp 不存在：{mdp_10ns}")
    if not mdp_50ns.is_file():
        warn(f"50ns mdp 不存在：{mdp_50ns}，首轮将使用 10ns mdp")

    # ============================================================
    # Phase 0+1: 检测并运行 0V
    # ============================================================
    print(f"\n{'#' * 72}")
    print("# Phase 0: 检测 0V 是否已完成")
    print(f"{'#' * 72}")

    zero_v_dir = system_root / ZERO_V

    if is_voltage_converged(zero_v_dir):
        print(f"0V 已收敛 (存在 {EQUILIBRIUM_LOG})，跳过 0V，直接进入其他电压点")
    else:
        print(f"0V 未完成，开始准备并运行 0V")

        # 创建 0V 目录
        zero_v_dir.mkdir(parents=True, exist_ok=True)

        # 准备输入结构 (从 fine 取)
        print(f"\n--- 准备 0V 输入结构 ---")
        prepare_input_structure(zero_v_dir, fine_gro, fine_ndx)

        # 准备矩阵和控制文件
        print(f"\n--- 准备 0V 矩阵和控制文件 ---")
        prepare_voltage_files(zero_v_dir, ZERO_V, voltage_root)

        # 循环跑 0V
        print(f"\n--- 运行 0V 平衡循环 ---")
        loop_count = 0
        while loop_count < args.max_loops:
            loop_count += 1
            print(f"\n{'=' * 60}")
            print(f"0V LOOP {loop_count} / {args.max_loops}")
            print(f"{'=' * 60}")

            converged = run_one_voltage(
                zero_v_dir, args.gmx, params, system_root
            )
            if converged:
                break
        else:
            fail(
                f"0V 在 {args.max_loops} 轮后仍未收敛，"
                f"请人工检查 {zero_v_dir}"
            )

        print(f"\n0V 已收敛！")

    # ============================================================
    # Phase 2: 准备 1V/2V/3V/4V (用 0V 的 nvt.gro 作为输入)
    # ============================================================
    print(f"\n{'#' * 72}")
    print("# Phase 2: 准备 1V/2V/3V/4V")
    print(f"{'#' * 72}")

    # 0V 的 nvt.gro 作为其他电压点的输入
    zero_v_gro = zero_v_dir / NVT_GRO
    if not zero_v_gro.is_file():
        fail(f"0V 的 {NVT_GRO} 不存在：{zero_v_gro}")

    other_voltages = [v for v in VOLTAGE_DIRS if v != ZERO_V]

    for vname in other_voltages:
        vdir = system_root / vname
        vdir.mkdir(parents=True, exist_ok=True)

        print(f"\n--- 准备 {vname} ---")
        # 输入结构：从 0V 的 nvt.gro 复制
        prepare_input_structure(vdir, zero_v_gro, fine_ndx)
        # 矩阵和控制文件
        prepare_voltage_files(vdir, vname, voltage_root)

    # ============================================================
    # Phase 3: 循环跑 1V/2V/3V/4V
    # ============================================================
    print(f"\n{'#' * 72}")
    print("# Phase 3: 运行 1V/2V/3V/4V 平衡循环")
    print(f"{'#' * 72}")

    # 过滤掉已收敛的
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
                vd, args.gmx, params, system_root
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

    # ============================================================
    # 汇总结果
    # ============================================================
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
    main()
