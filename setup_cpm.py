#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import json
import shutil
import subprocess
import sys
from pathlib import Path

SUMMARY_NAME = "system_summary.json"
VOLTAGE_CONFIG = {
    "0v": (-0.0, 0.0),
    "1v": (-0.5, 0.5),
    "2v": (-1.0, 1.0),
    "3v": (-1.5, 1.5),
}
BULK = (19, 21)
OUTQ = 1000
NZERO = 0
MATRIX_DOUBLE = 1
MATRIX_3DC = 0

def fail(message, code=1):
    print(f"错误：{message}")
    sys.exit(code)

def load_summary(path):
    path = Path(path)
    if not path.is_file():
        fail(f"找不到 summary 文件：{path}")
    try:
        data = json.loads(path.read_text(encoding="utf-8"))
    except Exception as exc:
        fail(f"无法解析 summary：{path}\n{exc}")
    try:
        total_atoms = int(data["mof_electrode_statistics"]["total_mof_atoms"])
    except (KeyError, TypeError, ValueError):
        fail(f"{path} 中缺少有效的 mof_electrode_statistics.total_mof_atoms")
    if total_atoms <= 0:
        fail(f"{path} 中 total_mof_atoms 非法：{total_atoms}")
    return total_atoms

def check_command(command):
    if shutil.which(command) is None:
        fail(f"找不到命令：{command}，请确认环境已经正确加载。")

def run_command(args, cwd):
    print(f"\n[运行目录] {cwd}")
    print("[命令] " + " ".join(args))
    try:
        subprocess.run(args, cwd=str(cwd), check=True)
    except FileNotFoundError:
        fail(f"找不到命令：{args[0]}")
    except subprocess.CalledProcessError as exc:
        fail(f"命令执行失败，返回码 {exc.returncode}：\n" + " ".join(args))

def main():
    mof_dir = Path.cwd().resolve()

    acn_dir = mof_dir / "ACN"
    pc_dir = mof_dir / "PC"

    if not acn_dir.is_dir():
        fail(f"找不到 ACN 目录：{acn_dir}")
    if not pc_dir.is_dir():
        fail(f"找不到 PC 目录：{pc_dir}")

    acn_total_atoms = load_summary(acn_dir / SUMMARY_NAME)
    pc_total_atoms = load_summary(pc_dir / SUMMARY_NAME)

    if acn_total_atoms != pc_total_atoms:
        fail(
            "ACN 与 PC summary 中的 MOF 电极总原子数不一致：\n"
            f"  ACN = {acn_total_atoms}\n"
            f"  PC  = {pc_total_atoms}"
        )

    total_atoms = acn_total_atoms
    if total_atoms % 2 != 0:
        fail(
            f"MOF 电极总原子数 {total_atoms} 不是偶数，"
            "无法按上下两块等原子数电极平均分配。"
        )

    single_electrode_atoms = total_atoms // 2
    nvt_tpr = acn_dir / "first" / "nvt.tpr"

    if not nvt_tpr.is_file():
        fail(f"找不到 matrix 输入文件：{nvt_tpr}")

    check_command("itptools")

    print("=" * 72)
    print("CPM INITIALIZATION")
    print("=" * 72)
    print(f"MOF directory           : {mof_dir.name}")
    print(f"Matrix source TPR       : {nvt_tpr}")
    print(f"Total electrode atoms   : {total_atoms}")
    print(f"Atoms per electrode     : {single_electrode_atoms}")

    matrix_cmd = [
        "/home/itp_students/zhouzhisheng/software/itptools_20210723/itptools", "matrix",
        "-s", str(nvt_tpr.relative_to(mof_dir)),
        "-n", str(total_atoms),
        "-double", str(MATRIX_DOUBLE),
        "-3dc", str(MATRIX_3DC),
    ]
    run_command(matrix_cmd, cwd=mof_dir)

    for dirname, (left_pot, right_pot) in VOLTAGE_CONFIG.items():
        voltage_dir = mof_dir / dirname
        voltage_dir.mkdir(exist_ok=True)

        pot_arg = f"[{left_pot:.1f},{right_pot:.1f}]"
        npot_arg = f"[{single_electrode_atoms},{single_electrode_atoms}]"
        bulk_arg = f"[{BULK[0]},{BULK[1]}]"

        cmd = [
            "/home/itp_students/zhouzhisheng/software/itptools_20210723/itptools", "newCPM_File",
            "-pot", pot_arg,
            "-npot", npot_arg,
            "-bulk", bulk_arg,
            "-outQ", str(OUTQ),
            "-nzero", str(NZERO),
        ]
        run_command(cmd, cwd=voltage_dir)

    print("\nCPM 初始化完成。")

if __name__ == "__main__":
    main()
