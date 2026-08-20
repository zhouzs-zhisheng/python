#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
density_collect.py

用途
----
在指定的 NVT 模拟目录中运行：

    gmx density \
        -f nvt.xtc \
        -s nvt.tpr \
        -n <自动定位的 index.ndx> \
        -sl 1000 \
        -o density.xvg \
        -d Z \
        -b 7000

并自动选择 group 6。

与旧版不同，本脚本不会把 density.xvg 的全部 1000 个切片数据写入汇总文件，
而只统计 Z = 18.0 ~ 22.0 nm 区间内的密度平均值。

汇总文件不再位于单个 qmof-* 目录，而是位于所有 qmof-* 目录共同的上一级
MOF 总目录中：

    MOF/density_profiles.csv

每个体系在 CSV 中只占一行。

------------------------------------------------------------
使用示例
------------------------------------------------------------

当前目录：

    /home/.../MOF/

运行：

    python3 density_collect.py qmof-305c717/ACN/first/

也可以：

    python3 density_collect.py qmof-305c717/ACN/400/first/

或者使用绝对路径。

------------------------------------------------------------
目录识别
------------------------------------------------------------

例如输入：

    MOF/qmof-305c717/ACN/first/

自动识别：

    MOF name       = qmof-305c717
    Solvent        = ACN
    System         = 从最近的 system_summary.json 中读取 EMIM 数量
    Relative path  = qmof-305c717/ACN/first

例如输入：

    MOF/qmof-305c717/ACN/500/first/

则记录路径：

    qmof-305c717/ACN/500/first

------------------------------------------------------------
输出文件
------------------------------------------------------------

1. 当前模拟目录：

    density.xvg

    这是 GROMACS 原始输出，完整保留。

2. 所有 qmof-* 的共同父目录 MOF/：

    density_profiles.csv

CSV 每个体系一行，例如：

Directory,MOF,Solvent,System,Z_Min_nm,Z_Max_nm,Points,Density_Mean
qmof-305c717/ACN/first,qmof-305c717,ACN,400,18.0,22.0,201,xxx
qmof-305c717/PC/first,qmof-305c717,PC,400,18.0,22.0,201,xxx
qmof-7375a78/ACN/first,qmof-7375a78,ACN,400,18.0,22.0,201,xxx

------------------------------------------------------------
覆盖和排序规则
------------------------------------------------------------

唯一键：

    Directory

如果同一个目录再次运行，例如：

    qmof-305c717/ACN/first

已有记录会被新结果替换，不重复追加。

每次写入后，整个 CSV 按 Directory 字典序排序。

------------------------------------------------------------
默认参数
------------------------------------------------------------

density group = 6
slices        = 1000
direction     = Z
begin time    = 7000 ps

统计区间：
    18.0 <= Z <= 22.0 nm
"""

import argparse
import csv
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path


GROUP_ID = 6
N_SLICES = 1000
DIRECTION = "Z"
BEGIN_PS = 7000

Z_MIN_NM = 18.0
Z_MAX_NM = 22.0

XTC_NAME = "nvt.xtc"
TPR_NAME = "nvt.tpr"
XVG_NAME = "density.xvg"
CSV_NAME = "density_profiles.csv"


def fail(message, code=1):
    print(f"错误：{message}")
    sys.exit(code)


def warn(message):
    print(f"警告：{message}")


def find_solvent_and_mof(sim_dir):
    """
    从模拟目录向上寻找最近的 ACN 或 PC。

    返回：
        solvent_dir
        solvent_name
        mof_dir
    """
    current = sim_dir.resolve()

    while True:
        if current.name in ("ACN", "PC"):
            return current, current.name, current.parent

        if current.parent == current:
            break

        current = current.parent

    fail(
        f"无法从路径中识别 ACN/PC：{sim_dir}\n"
        "目录必须位于 .../<MOF>/ACN/... 或 .../<MOF>/PC/... 下。"
    )


def find_summary_for_system(sim_dir, solvent_dir):
    """
    从模拟目录开始向 solvent_dir 逐级向上寻找 system_summary.json。

    因此：
        ACN/first
        -> ACN/system_summary.json

        ACN/500/first
        -> 优先 ACN/500/system_summary.json
    """
    current = sim_dir.resolve()
    solvent_dir = solvent_dir.resolve()

    while True:
        candidate = current / "system_summary.json"

        if candidate.is_file():
            return candidate

        if current == solvent_dir or current.parent == current:
            break

        current = current.parent

    return None


def determine_system_name(sim_dir, solvent_dir, summary_path):
    """
    优先使用 summary 中 electrolyte_molecules.EMIM 作为体系编号。

    如果 summary 不存在，则：
        - 若 ACN/PC 下一级目录是数字，则用该数字；
        - 否则记作 base。
    """
    if summary_path is not None:
        try:
            data = json.loads(summary_path.read_text(encoding="utf-8"))
            emim = data.get("electrolyte_molecules", {}).get("EMIM")

            if isinstance(emim, int):
                return str(emim)

        except Exception as exc:
            warn(f"无法从 summary 读取 EMIM 数量：{exc}")

    try:
        relative = sim_dir.resolve().relative_to(solvent_dir.resolve())

        if relative.parts:
            first = relative.parts[0]

            if first.isdigit():
                return first

    except Exception:
        pass

    return "base"


def find_index_file(sim_dir, solvent_dir):
    """
    从 sim_dir 向 solvent_dir 搜索 index.ndx。

    支持：
        ACN/first -> ACN/index.ndx
        ACN/500/first -> ACN/500/index.ndx 或 ACN/index.ndx
    """
    current = sim_dir.resolve()
    solvent_dir = solvent_dir.resolve()

    while True:
        candidate = current / "index.ndx"

        if candidate.is_file():
            return candidate

        if current == solvent_dir or current.parent == current:
            break

        current = current.parent

    fail(
        f"从 {sim_dir} 到 {solvent_dir} 范围内找不到 index.ndx。"
    )


def determine_collection_root(mof_dir):
    """
    mof_dir 是：
        .../MOF/qmof-xxxx

    汇总文件要求放到：
        .../MOF/

    因而返回 mof_dir.parent。
    """
    return mof_dir.resolve().parent


def run_density(sim_dir, index_path):
    """
    执行 gmx density，stdin 自动输入 group 6。
    """
    sim_dir = sim_dir.resolve()

    xtc_path = sim_dir / XTC_NAME
    tpr_path = sim_dir / TPR_NAME
    xvg_path = sim_dir / XVG_NAME

    if not xtc_path.is_file():
        fail(f"找不到轨迹文件：{xtc_path}")

    if not tpr_path.is_file():
        fail(f"找不到 TPR 文件：{tpr_path}")

    if shutil.which("gmx") is None:
        fail("找不到 gmx 命令，请确认已经加载 GROMACS 环境。")

    index_rel = os.path.relpath(index_path, sim_dir)

    cmd = [
        "gmx",
        "density",
        "-f", XTC_NAME,
        "-s", TPR_NAME,
        "-n", index_rel,
        "-sl", str(N_SLICES),
        "-o", XVG_NAME,
        "-d", DIRECTION,
        "-b", str(BEGIN_PS),
    ]

    print("\n[运行目录]")
    print(sim_dir)

    print("\n[命令]")
    print(" ".join(cmd))

    print(f"\n[自动选择 group] {GROUP_ID}")

    try:
        subprocess.run(
            cmd,
            cwd=str(sim_dir),
            input=f"{GROUP_ID}\n",
            text=True,
            check=True,
        )

    except subprocess.CalledProcessError as exc:
        fail(f"gmx density 执行失败，返回码：{exc.returncode}")

    except FileNotFoundError:
        fail("找不到 gmx 命令。")

    if not xvg_path.is_file():
        fail(
            f"gmx density 已结束，但没有生成预期文件：{xvg_path}"
        )

    return xvg_path


def parse_xvg_region(xvg_path, z_min, z_max):
    """
    读取 density.xvg 数值区。

    跳过所有：
        #
        @

    仅保留：
        z_min <= Z <= z_max

    返回：
        selected_points
        mean_density
    """
    selected = []

    with open(
        xvg_path,
        "r",
        encoding="utf-8",
        errors="replace",
    ) as f:

        for lineno, line in enumerate(f, start=1):
            stripped = line.strip()

            if not stripped:
                continue

            if stripped.startswith("#") or stripped.startswith("@"):
                continue

            parts = stripped.split()

            if len(parts) < 2:
                continue

            try:
                z = float(parts[0])
                density = float(parts[1])

            except ValueError:
                warn(
                    f"{xvg_path}:{lineno} 无法解析数值，已忽略。"
                )
                continue

            if z_min <= z <= z_max:
                selected.append((z, density))

    if not selected:
        fail(
            f"{xvg_path} 中没有找到 Z={z_min}~{z_max} nm "
            "范围内的数据。\n"
            "请确认体系 Z 尺寸和 density.xvg 横坐标。"
        )

    mean_density = (
        sum(density for _, density in selected)
        / len(selected)
    )

    return selected, mean_density


def load_existing_csv(csv_path):
    """
    读取已有汇总 CSV。
    """
    if not csv_path.is_file():
        return []

    rows = []

    with open(
        csv_path,
        "r",
        newline="",
        encoding="utf-8",
    ) as f:

        reader = csv.DictReader(f)

        expected = {
            "Directory",
            "MOF",
            "Solvent",
            "System",
            "Z_Min_nm",
            "Z_Max_nm",
            "Points",
            "Density_Mean",
        }

        if (
            reader.fieldnames is None
            or not expected.issubset(set(reader.fieldnames))
        ):
            fail(
                f"已有 CSV 格式与当前脚本不兼容：{csv_path}\n"
                "如果这是旧版 density_collect.py 生成的完整 profile CSV，"
                "请先备份并删除/改名，然后重新运行。"
            )

        for row in reader:
            rows.append(row)

    return rows


def update_summary_csv(
    csv_path,
    directory_key,
    mof_name,
    solvent,
    system_name,
    selected_points,
    mean_density,
):
    """
    每个模拟目录仅保存一行。

    覆盖规则：
        Directory 相同 -> 删除旧记录，写入新记录。

    排序规则：
        按 Directory 字符串排序。
    """
    rows = load_existing_csv(csv_path)

    rows = [
        row
        for row in rows
        if row["Directory"] != directory_key
    ]

    rows.append({
        "Directory": directory_key,
        "MOF": mof_name,
        "Solvent": solvent,
        "System": system_name,
        "Z_Min_nm": f"{Z_MIN_NM:.3f}",
        "Z_Max_nm": f"{Z_MAX_NM:.3f}",
        "Points": str(len(selected_points)),
        "Density_Mean": f"{mean_density:.10g}",
    })

    rows.sort(
        key=lambda row: row["Directory"]
    )

    with open(
        csv_path,
        "w",
        newline="",
        encoding="utf-8",
    ) as f:

        fieldnames = [
            "Directory",
            "MOF",
            "Solvent",
            "System",
            "Z_Min_nm",
            "Z_Max_nm",
            "Points",
            "Density_Mean",
        ]

        writer = csv.DictWriter(
            f,
            fieldnames=fieldnames,
        )

        writer.writeheader()
        writer.writerows(rows)


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "运行 gmx density，仅统计 Z=18~22 nm 区间平均密度，"
            "并汇总到所有 MOF 目录的共同父目录 density_profiles.csv。"
        )
    )

    parser.add_argument(
        "simulation_dir",
        help=(
            "包含 nvt.xtc 和 nvt.tpr 的模拟目录，"
            "例如 qmof-305c717/ACN/first/"
        ),
    )

    return parser.parse_args()


def main():
    args = parse_args()

    sim_dir = Path(
        args.simulation_dir
    ).resolve()

    if not sim_dir.is_dir():
        fail(
            f"模拟目录不存在：{sim_dir}"
        )

    solvent_dir, solvent_name, mof_dir = (
        find_solvent_and_mof(sim_dir)
    )

    collection_root = determine_collection_root(
        mof_dir
    )

    summary_path = find_summary_for_system(
        sim_dir,
        solvent_dir,
    )

    system_name = determine_system_name(
        sim_dir,
        solvent_dir,
        summary_path,
    )

    index_path = find_index_file(
        sim_dir,
        solvent_dir,
    )

    try:
        directory_key = sim_dir.relative_to(
            collection_root
        ).as_posix()
    except ValueError:
        fail(
            f"模拟目录 {sim_dir} 不在汇总根目录 "
            f"{collection_root} 之下。"
        )

    print("=" * 72)
    print("DENSITY PROFILE COLLECTION")
    print("=" * 72)

    print(f"MOF              : {mof_dir.name}")
    print(f"Solvent          : {solvent_name}")
    print(f"System           : {system_name}")
    print(f"Simulation dir   : {sim_dir}")
    print(f"Directory key    : {directory_key}")
    print(f"Index file       : {index_path}")
    print(
        f"Summary          : "
        f"{summary_path if summary_path else '未找到'}"
    )
    print(
        f"Density range    : "
        f"{Z_MIN_NM:.3f} <= Z <= {Z_MAX_NM:.3f} nm"
    )

    xvg_path = run_density(
        sim_dir,
        index_path,
    )

    selected_points, mean_density = (
        parse_xvg_region(
            xvg_path,
            Z_MIN_NM,
            Z_MAX_NM,
        )
    )

    csv_path = (
        collection_root
        / CSV_NAME
    )

    update_summary_csv(
        csv_path=csv_path,
        directory_key=directory_key,
        mof_name=mof_dir.name,
        solvent=solvent_name,
        system_name=system_name,
        selected_points=selected_points,
        mean_density=mean_density,
    )

    print("\n" + "=" * 72)
    print("密度统计完成")
    print("=" * 72)

    print(
        f"原始 density.xvg   : {xvg_path}"
    )

    print(
        f"统计 Z 区间        : "
        f"{Z_MIN_NM:.3f} ~ {Z_MAX_NM:.3f} nm"
    )

    print(
        f"区间数据点数量     : "
        f"{len(selected_points)}"
    )

    print(
        f"区间平均密度       : "
        f"{mean_density:.10g}"
    )

    print(
        f"全局汇总 CSV       : "
        f"{csv_path}"
    )


if __name__ == "__main__":
    main()
