#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
adjust_density.py
=================

用途
----
针对单个 MOF 的单个溶剂体系（ACN 或 PC），在一个长时间 PBS 作业内连续调节
电解液密度，直到 Z=18~22 nm 区域平均密度达到目标精度。

脚本不提交新的 PBS 任务。它假定自己已经运行在一个获得计算资源的 PBS 作业中。

示例 PBS：
    #!/bin/bash
    #PBS -N density_adjust
    #PBS -l nodes=1:ppn=32
    #PBS -l walltime=720:00:00
    #PBS -q new
    #PBS -j oe

    cd $PBS_O_WORKDIR

    python3 scripts/adjust_density.py \
        qmof-305c717/ACN \
        --target-density 1186.5

核心调节规则
------------
1. ACN 体系每次新增严格满足：
       ΔEMIM : ΔBF4 : ΔACN = 1 : 1 : 5
       ΔPC = 0

2. PC 体系每次新增严格满足：
       ΔEMIM : ΔBF4 : ΔPC = 1 : 1 : 5
       ΔACN = 0

3. 只增加分子，不删除分子。

4. 粗调阈值默认 3%：
       |rho-target| / target <= 0.03

   第一次达到粗调阈值时，立刻保存：
       <solvent_root>/coarse_checkpoint.json

   保存完成之后，才切换到 fine 阶段。

5. 精调阈值默认 1%：
       |rho-target| / target <= 0.01

6. 精调发生过冲时，不从高密度体系删除分子，而是利用：
       lower bracket（低于目标）
       upper bracket（高于目标）
   从 lower bracket 对应的旧体系重新生成一个中间 EMIM 数量的体系。

目录假设
--------
MOF/
├── basic/
│   └── cpm_mdp/
│       ├── min.mdp
│       └── nvt.mdp
├── density_profiles.csv
├── scripts/
│   ├── adjust_density.py
│   ├── density_collect_region_summary.py
│   ├── incremental_add_molecules.py
│   ├── setup_cpm.py
│   ├── start_inp_top.py
│   ├── fix_top.py
│   ├── matrix.sh
│   ├── tijiao.sh
│   └── start_inp_top.sh
└── qmof-305c717/
    ├── ACN/
    │   ├── system.gro
    │   ├── topol.top
    │   ├── system_summary.json
    │   ├── index.ndx
    │   ├── first/
    │   │   ├── min.gro
    │   │   ├── nvt.gro
    │   │   ├── nvt.tpr
    │   │   └── nvt.xtc
    │   ├── 420/
    │   ├── 440/
    │   └── ...
    └── PC/
        └── ...

运行方式
--------
从 MOF 总目录运行：

    python3 scripts/adjust_density.py qmof-305c717/ACN --target-density 1186.5

或者：

    python3 scripts/adjust_density.py qmof-305c717/PC --target-density 1200.0

也可从某个已存在的数字体系开始：

    python3 scripts/adjust_density.py qmof-305c717/ACN/440 --target-density 1186.5

会生成/更新
----------
每个试验体系目录：
    index.ndx
    first/min.tpr
    first/min.gro
    first/nvt.tpr
    first/nvt.gro
    first/nvt.xtc
    first/density.xvg

ACN/ 或 PC/ 根目录：
    density_adjustment.json
    density_adjustment.csv
    coarse_checkpoint.json
    final_density_result.json
    coarse/                       # 第一次达到 coarse 标准的体系副本
    fine/                         # 最终达到 fine 标准的体系副本

MOF 总目录：
    density_profiles.csv

新增体系目录由 incremental_add_molecules.py 生成，例如：
    ACN/420/
    ACN/440/
    ACN/446/

密度相关文件的数据层级
--------------------
1. <system>/first/density.xvg
   GROMACS 原始 Z 方向 density profile。完整保留，不作为跨体系状态文件。

2. <solvent_root>/density_adjustment.json
   自动调密度主状态。保存目标密度、coarse/fine 阶段、所有已评价体系、
   当前 trial 与粗调 checkpoint。程序恢复时主要读取此文件。

3. <solvent_root>/density_adjustment.csv
   与 JSON 中 measurements 对应的人类可读历史表，每个被评价体系一行。

4. <solvent_root>/coarse_checkpoint.json
   第一次达到粗调容差后立即保存的状态快照。
   同时把该收敛体系复制到 <solvent_root>/coarse/。
   两者完成后才进入 fine。

5. <solvent_root>/final_density_result.json
   达到 fine 容差后的最终结果。
   同时把最终收敛体系复制到 <solvent_root>/fine/。

6. <MOF总目录>/density_profiles.csv
   与 density_collect_region_summary.py 共用的全局汇总格式。
   每个体系仅一行，保存 Z=18~22 nm 区域平均密度；Directory 相同则覆盖，
   最终按 Directory 排序。

逻辑顺序：
   NVT -> first/density.xvg -> 区间平均密度 -> density_adjustment.json/csv
       -> density_profiles.csv -> coarse_checkpoint 或 final_density_result

恢复机制
--------
即使长任务意外中止，已经生成的体系和日志仍保留。
重新运行本脚本时：
    - 如果数字体系目录只生成了一半，会先改名备份为 *.incomplete_backup_*，
      再从父体系重新生成该数字体系；
    - 已有 index.ndx 不重建；
    - 已完成 min.gro 不重跑 minimization；
    - 若 NVT 已有 nvt.gro，则视为完成；
    - 若只有 nvt.cpt，可从 checkpoint 继续 NVT；
    - 已记录过的 density 不重复计算。
"""

import argparse
import csv
import json
import math
import os
import shutil
import subprocess
import sys
import time
from pathlib import Path


COARSE_TOLERANCE = 0.03
FINE_TOLERANCE = 0.01

INITIAL_STEP_FRACTION = 0.05
COARSE_MAX_STEP_FRACTION = 0.10
FINE_MAX_STEP_FRACTION = 0.02

DENSITY_BEGIN_PS = 7000
DENSITY_Z_MIN_NM = 18.0
DENSITY_Z_MAX_NM = 22.0
DENSITY_SLICES = 1000
DENSITY_GROUP = 6

# 真空判据：
# 某个允许电解液存在的 Z 区域内，如果连续至少两个 density slice
# 低于目标密度的 1%，则认为存在局部真空/断层。
VACUUM_DENSITY_FRACTION = 0.02
VACUUM_MIN_CONSECUTIVE_POINTS = 5

# 根据真空切片比例估算补充分子数：
#   k ~= N_EMIM * f / (1-f)
# 为避免一次补得过猛，单轮最多增加当前 EMIM 的 20%。
VACUUM_MAX_STEP_FRACTION = 0.20

# 如果真空切片位于 MOF 电极区域，说明其可能与孔道/低可达体积有关，
# 对补充分子量的贡献只按普通电解液区域的 1/4 计。
ELECTRODE_VACUUM_WEIGHT = 0.25

NTMPI = 32
NTOMP = 1
MAX_ITERATIONS = 30

SUMMARY_NAME = "system_summary.json"
STATE_JSON = "density_adjustment.json"
STATE_CSV = "density_adjustment.csv"
COARSE_CHECKPOINT_JSON = "coarse_checkpoint.json"
FINAL_RESULT_JSON = "final_density_result.json"
GLOBAL_DENSITY_CSV = "density_profiles.csv"

# 收敛体系快照目录名：
#   <ACN或PC>/coarse/
#   <ACN或PC>/fine/
COARSE_SNAPSHOT_DIR = "coarse"
FINE_SNAPSHOT_DIR = "fine"

MOLECULES = ("EMIM", "BF4", "ACN", "PC")


def fail(message, code=1):
    print(f"\n错误：{message}")
    sys.exit(code)


def warn(message):
    print(f"警告：{message}")


def now_string():
    return time.strftime("%Y-%m-%d %H:%M:%S")


def run_command(args, cwd, stdin_text=None):
    print("\n[运行目录]")
    print(cwd)
    print("[命令]")
    print(" ".join(str(x) for x in args))

    try:
        subprocess.run(
            [str(x) for x in args],
            cwd=str(cwd),
            input=stdin_text,
            text=True if stdin_text is not None else False,
            check=True,
        )
    except FileNotFoundError:
        fail(f"找不到命令：{args[0]}")
    except subprocess.CalledProcessError as exc:
        fail(
            f"命令执行失败，返回码 {exc.returncode}：\n"
            + " ".join(str(x) for x in args)
        )


def check_command(name):
    if shutil.which(name) is None:
        fail(f"找不到命令：{name}。请确认 PBS 环境已经加载所需程序。")


def find_solvent_root(path):
    current = Path(path).resolve()

    while True:
        if current.name in ("ACN", "PC"):
            solvent_root = current
            mof_dir = current.parent
            collection_root = mof_dir.parent
            return solvent_root, current.name, mof_dir, collection_root

        if current.parent == current:
            break

        current = current.parent

    fail(
        f"无法从路径识别 ACN/PC：{path}\n"
        "输入应位于 .../<MOF>/ACN 或 .../<MOF>/PC 下。"
    )


def find_incremental_script(collection_root, explicit=None):
    if explicit:
        p = Path(explicit).resolve()
        if not p.is_file():
            fail(f"指定的增量脚本不存在：{p}")
        return p

    scripts_dir = collection_root / "scripts"

    candidates = [
        # 优先使用 Packmol 容错版：即使用户没有重命名，也可自动使用。
        scripts_dir / "incremental_add_molecules_packmol_safe.py",
        scripts_dir / "incremental_add_molecules.py",
        scripts_dir / "incremental_add_molecules_fixed_top.py",
    ]

    for p in candidates:
        if p.is_file():
            return p.resolve()

    fail(
        "在 MOF/scripts/ 中找不到可用的增量建模脚本。"
        "可用 --incremental-script 显式指定。"
    )


def find_mdp_files(collection_root):
    min_mdp = collection_root / "basic" / "cpm_mdp" / "min.mdp"
    nvt_mdp = collection_root / "basic" / "cpm_mdp" / "nvt.mdp"

    if not min_mdp.is_file():
        fail(f"找不到 min.mdp：{min_mdp}")

    if not nvt_mdp.is_file():
        fail(f"找不到 nvt.mdp：{nvt_mdp}")

    return min_mdp.resolve(), nvt_mdp.resolve()


def load_summary(system_dir):
    path = Path(system_dir) / SUMMARY_NAME

    if not path.is_file():
        fail(f"找不到体系 summary：{path}")

    try:
        data = json.loads(path.read_text(encoding="utf-8"))
    except Exception as exc:
        fail(f"summary 解析失败：{path}\n{exc}")

    raw = data.get("electrolyte_molecules")

    if not isinstance(raw, dict):
        fail(f"{path} 中缺少 electrolyte_molecules。")

    composition = {}

    for mol in MOLECULES:
        value = raw.get(mol, 0)

        if isinstance(value, bool) or not isinstance(value, int) or value < 0:
            fail(f"{path} 中 {mol} 数量非法：{value!r}")

        composition[mol] = value

    return data, composition


def increment_ratio(solvent):
    if solvent == "ACN":
        return {"EMIM": 1, "BF4": 1, "ACN": 5, "PC": 0}

    if solvent == "PC":
        return {"EMIM": 1, "BF4": 1, "ACN": 0, "PC": 5}

    fail(f"不支持的 solvent：{solvent}")


def build_target_composition(parent_comp, target_emim, solvent):
    current_emim = parent_comp["EMIM"]

    if target_emim <= current_emim:
        fail(
            f"目标 EMIM={target_emim} 必须大于父体系 EMIM={current_emim}。"
        )

    k = target_emim - current_emim
    ratio = increment_ratio(solvent)
    target = dict(parent_comp)

    for mol in MOLECULES:
        target[mol] += k * ratio[mol]

    return target, k


def new_state(solvent, mof_name, target_density, initial_system_dir):
    return {
        "version": 1,
        "created_at": now_string(),
        "updated_at": now_string(),
        "mof_name": mof_name,
        "solvent": solvent,
        "target_density": target_density,
        "coarse_tolerance": COARSE_TOLERANCE,
        "fine_tolerance": FINE_TOLERANCE,
        "density_region_nm": [DENSITY_Z_MIN_NM, DENSITY_Z_MAX_NM],
        "increment_rule": increment_ratio(solvent),
        "stage": "coarse",
        "status": "running",
        "initial_system_dir": str(Path(initial_system_dir).resolve()),
        "current_trial_dir": str(Path(initial_system_dir).resolve()),
        "coarse_checkpoint": None,
        "measurements": [],
    }


def state_path(solvent_root):
    return Path(solvent_root) / STATE_JSON


def save_state(solvent_root, state):
    state["updated_at"] = now_string()

    state_path(solvent_root).write_text(
        json.dumps(state, indent=4, ensure_ascii=False),
        encoding="utf-8",
    )

    write_state_csv(solvent_root, state)


def load_or_create_state(
    solvent_root,
    solvent,
    mof_name,
    target_density,
    initial_system_dir,
    reset=False,
):
    path = state_path(solvent_root)

    if reset and path.exists():
        backup = path.with_name(path.stem + f".backup_{int(time.time())}.json")
        shutil.copy2(path, backup)
        path.unlink()
        print(f"旧状态已备份：{backup}")

    if path.is_file():
        try:
            state = json.loads(path.read_text(encoding="utf-8"))
        except Exception as exc:
            fail(f"状态文件解析失败：{path}\n{exc}")

        old_target = float(state.get("target_density"))

        if abs(old_target - target_density) > 1e-12:
            fail(
                f"已有状态目标密度为 {old_target}，本次为 {target_density}。\n"
                "如要重新开始，请使用 --reset-state。"
            )

        return state

    state = new_state(
        solvent=solvent,
        mof_name=mof_name,
        target_density=target_density,
        initial_system_dir=initial_system_dir,
    )
    save_state(solvent_root, state)
    return state


def write_state_csv(solvent_root, state):
    path = Path(solvent_root) / STATE_CSV

    fieldnames = [
        "Iteration",
        "Directory",
        "EMIM",
        "BF4",
        "ACN",
        "PC",
        "Density",
        "Target",
        "RelativeError",
        "AbsRelativeError",
        "StageAtMeasurement",
        "VacuumDetected",
        "VacuumFraction",
        "ElectrolyteVacuumFraction",
        "ElectrodeVacuumFraction",
        "EffectiveVacuumFraction",
        "VacuumThreshold",
        "VacuumSegments",
        "MeasuredAt",
    ]

    rows = []

    for idx, m in enumerate(state["measurements"], start=1):
        comp = m["composition"]
        rows.append({
            "Iteration": idx,
            "Directory": m["directory"],
            "EMIM": comp["EMIM"],
            "BF4": comp["BF4"],
            "ACN": comp["ACN"],
            "PC": comp["PC"],
            "Density": m["density"],
            "Target": state["target_density"],
            "RelativeError": m["relative_error"],
            "AbsRelativeError": abs(m["relative_error"]),
            "StageAtMeasurement": m["stage"],
            "VacuumDetected": m.get("vacuum_detected", False),
            "VacuumFraction": m.get("vacuum_fraction", 0.0),
            "ElectrolyteVacuumFraction": m.get(
                "electrolyte_vacuum_fraction", 0.0
            ),
            "ElectrodeVacuumFraction": m.get(
                "electrode_vacuum_fraction", 0.0
            ),
            "EffectiveVacuumFraction": m.get(
                "effective_vacuum_fraction", 0.0
            ),
            "VacuumThreshold": m.get("vacuum_threshold", ""),
            "VacuumSegments": len(m.get("vacuum_segments", [])),
            "MeasuredAt": m["measured_at"],
        })

    with open(path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def find_measurement_by_dir(state, system_dir):
    target = str(Path(system_dir).resolve())

    for item in state["measurements"]:
        if str(Path(item["directory"]).resolve()) == target:
            return item

    return None


def ensure_index(system_dir):
    """
    创建 index.ndx。

    某些 GROMACS/批处理环境下，把 make_ndx 拆成两次调用并分别通过
    subprocess stdin 传入命令，可能触发：

        Fatal error:
        Error reading user input

    因此这里改成一次 make_ndx 调用，一次性输入：

        3 | 4 | 5
        name 6 IONS
        q

    等价于 shell 中已经验证可用的：

        printf "3 | 4 | 5\\nname 6 IONS\\nq\\n" | \
        gmx make_ndx -f system.gro -o index.ndx

    第一个命令会把 group 3、4、5 合并为新 group 6；
    第二个命令立即把 group 6 重命名为 IONS；
    最后 q 退出并写出 index.ndx。
    """
    system_dir = Path(system_dir)
    index_path = system_dir / "index.ndx"

    if index_path.is_file():
        print(f"index 已存在，跳过：{index_path}")
        return index_path

    system_gro = system_dir / "system.gro"

    if not system_gro.is_file():
        fail(f"找不到 system.gro：{system_gro}")

    make_ndx_input = (
        "3 | 4 | 5\n"
        "name 6 IONS\n"
        "q\n"
    )

    run_command(
        [
            "gmx",
            "make_ndx",
            "-f", "system.gro",
            "-o", "index.ndx",
        ],
        cwd=system_dir,
        stdin_text=make_ndx_input,
    )

    if not index_path.is_file():
        fail(f"make_ndx 后找不到：{index_path}")

    print(f"index 已成功生成：{index_path}")
    return index_path


def ensure_minimization(system_dir, min_mdp):
    system_dir = Path(system_dir)
    first = system_dir / "first"
    first.mkdir(exist_ok=True)

    min_gro = first / "min.gro"
    min_tpr = first / "min.tpr"

    if min_gro.is_file():
        print(f"min 已完成，跳过：{min_gro}")
        return min_gro

    if not min_tpr.is_file():
        run_command(
            [
                "gmx", "grompp",
                "-f", str(min_mdp),
                "-c", "system.gro",
                "-o", "first/min.tpr",
                "-maxwarn", "6",
            ],
            cwd=system_dir,
        )

    run_command(
        [
            "gmx", "mdrun",
            "-s", "first/min.tpr",
            "-deffnm", "first/min",
            "-v",
        ],
        cwd=system_dir,
    )

    if not min_gro.is_file():
        fail(f"minimization 结束后找不到：{min_gro}")

    return min_gro


def ensure_nvt(system_dir, nvt_mdp):
    """
    保证 NVT 完成，并支持意外中断后的 checkpoint 续跑。

    判定规则：
    - nvt.gro + nvt.tpr + nvt.xtc 都存在：视为已正常完成，直接跳过；
    - nvt.gro 不存在但 nvt.cpt 存在：视为上次 NVT 中途停止，使用
      -cpi first/nvt.cpt -append 从 checkpoint 继续；
    - 没有 checkpoint：从现有 nvt.tpr 正常开始。

    注意：仅有 nvt.xtc 并不代表完成，因为任务中断时 xtc 往往已经存在。
    最终 nvt.gro 是这里判断完整结束的关键文件。
    """
    system_dir = Path(system_dir)
    first = system_dir / "first"
    first.mkdir(exist_ok=True)

    nvt_gro = first / "nvt.gro"
    nvt_tpr = first / "nvt.tpr"
    nvt_cpt = first / "nvt.cpt"
    nvt_xtc = first / "nvt.xtc"

    if nvt_gro.is_file() and nvt_tpr.is_file() and nvt_xtc.is_file():
        print(f"NVT 已完成，跳过：{nvt_gro}")
        return nvt_gro

    if not nvt_tpr.is_file():
        run_command(
            [
                "gmx", "grompp",
                "-f", str(nvt_mdp),
                "-c", "first/min.gro",
                "-n", "index.ndx",
                "-o", "first/nvt.tpr",
                "-maxwarn", "6",
            ],
            cwd=system_dir,
        )

    cmd = [
        "gmx", "mdrun",
        "-s", "first/nvt.tpr",
        "-deffnm", "first/nvt",
        "-v",
        "-ntomp", str(NTOMP),
        "-ntmpi", str(NTMPI),
    ]

    if nvt_cpt.is_file():
        print(
            "发现 NVT checkpoint，但最终 nvt.gro 尚未生成；"
            "判定为上次 NVT 意外中断，将从 checkpoint 继续。"
        )
        cmd.extend(["-cpi", "first/nvt.cpt", "-append"])
    elif nvt_xtc.is_file() and not nvt_gro.is_file():
        warn(
            "检测到 first/nvt.xtc 已存在，但没有 nvt.gro/nvt.cpt。"
            "这可能是异常中断留下的轨迹。由于没有可恢复 checkpoint，"
            "本次将从 nvt.tpr 重新运行 NVT，并由 GROMACS 处理已有输出文件。"
        )

    run_command(cmd, cwd=system_dir)

    if not nvt_gro.is_file():
        fail(f"NVT 结束后找不到：{nvt_gro}")

    if not nvt_xtc.is_file():
        fail(f"NVT 结束后找不到：{nvt_xtc}")

    return nvt_gro


def ensure_simulation(system_dir, min_mdp, nvt_mdp):
    ensure_index(system_dir)
    ensure_minimization(system_dir, min_mdp)
    ensure_nvt(system_dir, nvt_mdp)


def run_density(system_dir):
    system_dir = Path(system_dir)
    first = system_dir / "first"

    xtc = first / "nvt.xtc"
    tpr = first / "nvt.tpr"
    xvg = first / "density.xvg"
    index = system_dir / "index.ndx"

    for p in (xtc, tpr, index):
        if not p.is_file():
            fail(f"density 输入文件不存在：{p}")

    run_command(
        [
            "gmx", "density",
            "-f", "nvt.xtc",
            "-s", "nvt.tpr",
            "-n", "../index.ndx",
            "-sl", str(DENSITY_SLICES),
            "-o", "density.xvg",
            "-d", "Z",
            "-b", str(DENSITY_BEGIN_PS),
        ],
        cwd=first,
        stdin_text=f"{DENSITY_GROUP}\n",
    )

    if not xvg.is_file():
        fail(f"density 命令后未生成：{xvg}")

    return xvg


def get_electrode_regions(summary):
    """
    从 system_summary.json 的 z_structure_regions_nm 中提取 MOF 电极 Z 区域。

    兼容此前 summary 的典型格式：
        {
            "region_id": 2,
            "name": "lower_electrode",
            "type": "electrode",
            "z_low": ...,
            "z_high": ...
        }

    判定优先使用 type == "electrode"；
    如果旧 summary 没有规范 type，则 name 中包含 "electrode" 也视为电极区。

    返回：
        [(z_low, z_high), ...]
    """
    regions = summary.get("z_structure_regions_nm")

    if not isinstance(regions, list) or not regions:
        fail(
            "system_summary.json 中缺少 z_structure_regions_nm，"
            "无法判断真空段位于电极区域还是电解液区域。"
        )

    electrode_regions = []

    for idx, region in enumerate(regions, start=1):
        if not isinstance(region, dict):
            continue

        region_type = str(region.get("type", "")).strip().lower()
        region_name = str(region.get("name", "")).strip().lower()

        is_electrode = (
            region_type == "electrode"
            or "electrode" in region_name
        )

        if not is_electrode:
            continue

        try:
            z_low = float(region["z_low"])
            z_high = float(region["z_high"])
        except (KeyError, TypeError, ValueError):
            fail(
                f"z_structure_regions_nm 第 {idx} 个电极区域缺少合法的 "
                "z_low/z_high。"
            )

        if z_high <= z_low:
            fail(
                f"z_structure_regions_nm 第 {idx} 个电极区域范围非法："
                f"[{z_low}, {z_high}]"
            )

        electrode_regions.append((z_low, z_high))

    if not electrode_regions:
        fail(
            "z_structure_regions_nm 中没有识别到任何 electrode 区域。"
        )

    return electrode_regions


def point_region_type(z, electrode_regions):
    """
    根据 Z 坐标判断一个 density slice 属于：
        electrode
        electrolyte

    五区模型中，除上下两个 electrode 区域外，其余 Z 区域都按
    electrolyte 区域处理。
    """
    for z_low, z_high in electrode_regions:
        if z_low <= z <= z_high:
            return "electrode"

    return "electrolyte"


def analyze_density_profile(xvg_path, summary, target_density):
    """
    同时完成：

    1. 体相平均密度
       仍然只取 18 <= Z <= 22 nm。

    2. 整盒真空检测
       直接扫描整个 density.xvg。
       连续至少 5 个点低于目标密度 2% 才认定为 vacuum segment。

    3. 真空位置分类
       利用 summary["z_structure_regions_nm"] 判断每一个真空切片是：
           - electrode
           - electrolyte

       注意这里不是用一个 segment 的中心点简单分类，而是逐点分类。
       因此如果真空段跨过电极/电解液边界，其两部分会分别统计。

    返回的真空比例：
        vacuum_fraction
            = 所有真空切片 / 全部有效 Z 切片

        electrolyte_vacuum_fraction
            = 电解液区域真空切片 / 全部有效 Z 切片

        electrode_vacuum_fraction
            = 电极区域真空切片 / 全部有效 Z 切片

        effective_vacuum_fraction
            = electrolyte_vacuum_fraction
              + ELECTRODE_VACUUM_WEIGHT * electrode_vacuum_fraction

    其中 ELECTRODE_VACUUM_WEIGHT 默认 0.25，
    即电极区域真空对补充分子量的贡献只有正常电解液区域的 1/4。
    """
    profile = []

    with open(xvg_path, "r", encoding="utf-8", errors="replace") as f:
        for line in f:
            s = line.strip()

            if not s or s.startswith("#") or s.startswith("@"):
                continue

            parts = s.split()

            if len(parts) < 2:
                continue

            try:
                z = float(parts[0])
                rho = float(parts[1])
            except ValueError:
                continue

            profile.append((z, rho))

    if not profile:
        fail(f"{xvg_path} 中没有解析到有效 density 数据。")

    # --------------------------------------------------------
    # 1. 18~22 nm 体相平均密度
    # --------------------------------------------------------
    bulk_values = [
        rho
        for z, rho in profile
        if DENSITY_Z_MIN_NM <= z <= DENSITY_Z_MAX_NM
    ]

    if not bulk_values:
        fail(
            f"{xvg_path} 中找不到 Z={DENSITY_Z_MIN_NM}~"
            f"{DENSITY_Z_MAX_NM} nm 数据。"
        )

    bulk_density = sum(bulk_values) / len(bulk_values)

    # --------------------------------------------------------
    # 2. 获取上下两个电极区域
    # --------------------------------------------------------
    electrode_regions = get_electrode_regions(summary)

    # --------------------------------------------------------
    # 3. 整盒连续低密度检测
    # --------------------------------------------------------
    vacuum_threshold = (
        float(target_density) * VACUUM_DENSITY_FRACTION
    )

    total_points = len(profile)
    vacuum_points = 0
    electrode_vacuum_points = 0
    electrolyte_vacuum_points = 0

    vacuum_segments = []
    current_run = []

    def flush_run():
        nonlocal vacuum_points
        nonlocal electrode_vacuum_points
        nonlocal electrolyte_vacuum_points
        nonlocal current_run

        if len(current_run) >= VACUUM_MIN_CONSECUTIVE_POINTS:
            segment_electrode_points = 0
            segment_electrolyte_points = 0

            classified_points = []

            for z, rho in current_run:
                location = point_region_type(
                    z,
                    electrode_regions,
                )

                classified_points.append((z, rho, location))

                if location == "electrode":
                    segment_electrode_points += 1
                else:
                    segment_electrolyte_points += 1

            segment_total = len(current_run)

            vacuum_points += segment_total
            electrode_vacuum_points += segment_electrode_points
            electrolyte_vacuum_points += segment_electrolyte_points

            # 给 segment 一个主要位置标签，便于日志阅读。
            if (
                segment_electrode_points > 0
                and segment_electrolyte_points == 0
            ):
                location = "electrode"
            elif (
                segment_electrolyte_points > 0
                and segment_electrode_points == 0
            ):
                location = "electrolyte"
            else:
                location = "mixed"

            vacuum_segments.append({
                "z_start": current_run[0][0],
                "z_end": current_run[-1][0],
                "points": segment_total,
                "location": location,
                "electrode_points": segment_electrode_points,
                "electrolyte_points": segment_electrolyte_points,
                "min_density": min(
                    rho for _, rho in current_run
                ),
                "mean_density": (
                    sum(rho for _, rho in current_run)
                    / segment_total
                ),
            })

        current_run = []

    # Python nonlocal 不能使用括号列表，因此上面的定义需要普通语法。

    for z, rho in profile:
        if rho < vacuum_threshold:
            current_run.append((z, rho))
        else:
            flush_run()

    flush_run()

    vacuum_fraction = (
        vacuum_points / total_points
        if total_points > 0
        else 0.0
    )

    electrode_vacuum_fraction = (
        electrode_vacuum_points / total_points
        if total_points > 0
        else 0.0
    )

    electrolyte_vacuum_fraction = (
        electrolyte_vacuum_points / total_points
        if total_points > 0
        else 0.0
    )

    effective_vacuum_fraction = (
        electrolyte_vacuum_fraction
        + ELECTRODE_VACUUM_WEIGHT
        * electrode_vacuum_fraction
    )

    vacuum_detected = len(vacuum_segments) > 0

    return {
        "bulk_density": bulk_density,
        "bulk_points": len(bulk_values),

        "vacuum_detected": vacuum_detected,
        "vacuum_fraction": vacuum_fraction,
        "vacuum_threshold": vacuum_threshold,
        "vacuum_points": vacuum_points,
        "accessible_points": total_points,
        "vacuum_segments": vacuum_segments,

        "electrode_regions_nm": [
            {
                "z_low": z_low,
                "z_high": z_high,
            }
            for z_low, z_high in electrode_regions
        ],

        "electrode_vacuum_points": electrode_vacuum_points,
        "electrolyte_vacuum_points": electrolyte_vacuum_points,

        "electrode_vacuum_fraction": electrode_vacuum_fraction,
        "electrolyte_vacuum_fraction": electrolyte_vacuum_fraction,

        "electrode_vacuum_weight": ELECTRODE_VACUUM_WEIGHT,
        "effective_vacuum_fraction": effective_vacuum_fraction,
    }


def update_global_density_csv(
    collection_root,
    mof_name,
    solvent,
    system_dir,
    system_name,
    density,
    points,
):
    path = Path(collection_root) / GLOBAL_DENSITY_CSV

    try:
        directory_key = (
            Path(system_dir).resolve()
            .relative_to(Path(collection_root).resolve())
            .as_posix()
            + "/first"
        )
    except ValueError:
        directory_key = str(Path(system_dir).resolve())

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

    rows = []

    if path.is_file():
        with open(path, "r", newline="", encoding="utf-8") as f:
            reader = csv.DictReader(f)

            if reader.fieldnames and set(fieldnames).issubset(reader.fieldnames):
                rows.extend(reader)
            else:
                warn(
                    f"已有 {path} 格式不是当前区域汇总格式，"
                    "本次不修改该文件。"
                )
                return

    rows = [r for r in rows if r["Directory"] != directory_key]

    rows.append({
        "Directory": directory_key,
        "MOF": mof_name,
        "Solvent": solvent,
        "System": str(system_name),
        "Z_Min_nm": f"{DENSITY_Z_MIN_NM:.3f}",
        "Z_Max_nm": f"{DENSITY_Z_MAX_NM:.3f}",
        "Points": str(points),
        "Density_Mean": f"{density:.10g}",
    })

    rows.sort(key=lambda r: r["Directory"])

    with open(path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def measure_system(
    system_dir,
    state,
    solvent_root,
    collection_root,
    mof_name,
    solvent,
):
    existing = find_measurement_by_dir(state, system_dir)

    if existing is not None:
        print(f"体系已有密度记录，复用：{existing['density']}")
        return existing

    summary, comp = load_summary(system_dir)

    xvg = run_density(system_dir)

    target = float(state["target_density"])

    density_analysis = analyze_density_profile(
        xvg_path=xvg,
        summary=summary,
        target_density=target,
    )

    density = density_analysis["bulk_density"]
    points = density_analysis["bulk_points"]
    error = (density - target) / target

    measurement = {
        "directory": str(Path(system_dir).resolve()),
        "composition": comp,
        "density": density,
        "relative_error": error,
        "stage": state["stage"],
        "measured_at": now_string(),
        "points": points,

        # 完整盒子 Z-profile 的局部真空检测信息
        "vacuum_detected": density_analysis["vacuum_detected"],
        "vacuum_fraction": density_analysis["vacuum_fraction"],
        "vacuum_threshold": density_analysis["vacuum_threshold"],
        "vacuum_points": density_analysis["vacuum_points"],
        "accessible_points": density_analysis["accessible_points"],
        "vacuum_segments": density_analysis["vacuum_segments"],

        # 真空区域位置分类
        "electrode_vacuum_points": density_analysis[
            "electrode_vacuum_points"
        ],
        "electrolyte_vacuum_points": density_analysis[
            "electrolyte_vacuum_points"
        ],
        "electrode_vacuum_fraction": density_analysis[
            "electrode_vacuum_fraction"
        ],
        "electrolyte_vacuum_fraction": density_analysis[
            "electrolyte_vacuum_fraction"
        ],
        "electrode_vacuum_weight": density_analysis[
            "electrode_vacuum_weight"
        ],
        "effective_vacuum_fraction": density_analysis[
            "effective_vacuum_fraction"
        ],
    }

    state["measurements"].append(measurement)

    update_global_density_csv(
        collection_root=collection_root,
        mof_name=mof_name,
        solvent=solvent,
        system_dir=system_dir,
        system_name=comp["EMIM"],
        density=density,
        points=points,
    )

    save_state(solvent_root, state)

    return measurement



def backup_existing_snapshot(snapshot_dir):
    """
    如果 coarse/ 或 fine/ 已存在，不静默覆盖。
    先改名为 *.previous_backup_<timestamp>，再生成新的快照。
    """
    snapshot_dir = Path(snapshot_dir)

    if not snapshot_dir.exists():
        return None

    stamp = int(time.time())
    backup = snapshot_dir.with_name(
        snapshot_dir.name + f".previous_backup_{stamp}"
    )

    counter = 1
    while backup.exists():
        backup = snapshot_dir.with_name(
            snapshot_dir.name
            + f".previous_backup_{stamp}_{counter}"
        )
        counter += 1

    snapshot_dir.rename(backup)

    print(
        f"已有收敛快照 {snapshot_dir}，"
        f"已备份为 {backup}"
    )

    return backup


def copy_system_snapshot(source_dir, solvent_root, label):
    """
    将达到 coarse/fine 标准的体系复制到：
        <solvent_root>/coarse
        <solvent_root>/fine

    两种情况分别处理：

    A. source_dir 是数字体系目录，例如 ACN/446：
       直接完整复制整个 446/ -> ACN/coarse/ 或 ACN/fine/。

    B. source_dir 就是 ACN/ 或 PC/ 根目录：
       不能把根目录复制到自己的子目录，否则会递归。
       此时只复制“当前物理体系本身”的核心文件和 first/ 目录，
       不复制数字试验目录、density_adjustment 状态、coarse/fine 等控制文件。

    返回快照目录绝对路径。
    """
    source_dir = Path(source_dir).resolve()
    solvent_root = Path(solvent_root).resolve()
    snapshot_dir = solvent_root / label

    if label not in (COARSE_SNAPSHOT_DIR, FINE_SNAPSHOT_DIR):
        fail(f"未知 snapshot label：{label}")

    backup_existing_snapshot(snapshot_dir)

    # --------------------------------------------------------
    # 数字体系或其他子体系：整目录复制。
    # --------------------------------------------------------
    if source_dir != solvent_root:
        shutil.copytree(source_dir, snapshot_dir)
        print(
            f"已复制收敛体系：\n"
            f"  source = {source_dir}\n"
            f"  target = {snapshot_dir}"
        )
        return snapshot_dir.resolve()

    # --------------------------------------------------------
    # 初始 base 体系就是 solvent_root 时：
    # 只复制体系核心内容，避免递归复制整个 ACN/PC。
    # --------------------------------------------------------
    snapshot_dir.mkdir(parents=False, exist_ok=False)

    core_files = [
        "system.gro",
        "topol.top",
        "system_summary.json",
        "packmol_input.inp",
        "packmol.log",
        "merged_system.pdb",
        "merged_system.gro",
        "electrode_for_packmol.pdb",
        "index.ndx",
        "target_ratio.json",
        "mdout.mdp",
    ]

    copied_any = False

    for name in core_files:
        src = source_dir / name

        if src.is_file():
            shutil.copy2(src, snapshot_dir / name)
            copied_any = True

    first_dir = source_dir / "first"

    if first_dir.is_dir():
        shutil.copytree(
            first_dir,
            snapshot_dir / "first",
        )
        copied_any = True

    if not copied_any:
        fail(
            f"无法从 base 体系 {source_dir} 复制任何核心体系文件。"
        )

    # 保存来源说明，便于之后确认这个 coarse/fine 是从 base 体系复制的。
    metadata = {
        "snapshot_label": label,
        "source_directory": str(source_dir),
        "created_at": now_string(),
        "note": (
            "Source was the ACN/PC root itself; only the physical-system "
            "files and first/ simulation directory were copied to avoid "
            "recursive copying of trial/control directories."
        ),
    }

    (snapshot_dir / "snapshot_info.json").write_text(
        json.dumps(metadata, indent=4, ensure_ascii=False),
        encoding="utf-8",
    )

    print(
        f"已复制 base 收敛体系核心文件：\n"
        f"  source = {source_dir}\n"
        f"  target = {snapshot_dir}"
    )

    return snapshot_dir.resolve()



def save_coarse_checkpoint(solvent_root, state, measurement):
    # 达到 coarse 标准后，先复制当前完整体系到 <solvent_root>/coarse/
    # 再写 coarse_checkpoint.json，并切换到 fine。
    snapshot_dir = copy_system_snapshot(
        source_dir=measurement["directory"],
        solvent_root=solvent_root,
        label=COARSE_SNAPSHOT_DIR,
    )

    checkpoint = {
        "saved_at": now_string(),
        "target_density": state["target_density"],
        "coarse_tolerance": state["coarse_tolerance"],
        "directory": measurement["directory"],
        "snapshot_directory": str(snapshot_dir),
        "composition": measurement["composition"],
        "density": measurement["density"],
        "relative_error": measurement["relative_error"],
        "vacuum_detected": measurement.get("vacuum_detected", False),
        "vacuum_fraction": measurement.get("vacuum_fraction", 0.0),
    }

    path = Path(solvent_root) / COARSE_CHECKPOINT_JSON

    path.write_text(
        json.dumps(checkpoint, indent=4, ensure_ascii=False),
        encoding="utf-8",
    )

    state["coarse_checkpoint"] = checkpoint
    state["stage"] = "fine"

    save_state(solvent_root, state)

    print("\n" + "=" * 72)
    print("COARSE CHECKPOINT SAVED")
    print("=" * 72)
    print(f"Directory : {measurement['directory']}")
    print(f"Snapshot  : {snapshot_dir}")
    print(f"Density   : {measurement['density']:.6f}")
    print(f"Error     : {measurement['relative_error'] * 100.0:.4f}%")


def save_final_result(solvent_root, state, measurement):
    # 达到 fine 标准后，将最终收敛体系复制到 <solvent_root>/fine/
    snapshot_dir = copy_system_snapshot(
        source_dir=measurement["directory"],
        solvent_root=solvent_root,
        label=FINE_SNAPSHOT_DIR,
    )

    result = {
        "status": "CONVERGED",
        "finished_at": now_string(),
        "snapshot_directory": str(snapshot_dir),
        "mof_name": state["mof_name"],
        "solvent": state["solvent"],
        "target_density": state["target_density"],
        "fine_tolerance": state["fine_tolerance"],
        "directory": measurement["directory"],
        "composition": measurement["composition"],
        "density": measurement["density"],
        "relative_error": measurement["relative_error"],
        "vacuum_detected": measurement.get("vacuum_detected", False),
        "vacuum_fraction": measurement.get("vacuum_fraction", 0.0),
        "coarse_checkpoint": state.get("coarse_checkpoint"),
    }

    path = Path(solvent_root) / FINAL_RESULT_JSON

    path.write_text(
        json.dumps(result, indent=4, ensure_ascii=False),
        encoding="utf-8",
    )

    state["status"] = "converged"
    state["current_trial_dir"] = measurement["directory"]
    save_state(solvent_root, state)


def distinct_measurements_by_emim(state):
    mapping = {}

    for m in state["measurements"]:
        mapping[int(m["composition"]["EMIM"])] = m

    return [mapping[n] for n in sorted(mapping)]


def bracket_measurements(state):
    target = float(state["target_density"])
    ms = distinct_measurements_by_emim(state)

    lower_candidates = [m for m in ms if float(m["density"]) < target]
    upper_candidates = [m for m in ms if float(m["density"]) > target]

    lower = None
    upper = None

    if lower_candidates:
        lower = max(
            lower_candidates,
            key=lambda m: int(m["composition"]["EMIM"])
        )

    if upper_candidates:
        upper = min(
            upper_candidates,
            key=lambda m: int(m["composition"]["EMIM"])
        )

    if (
        lower is not None
        and upper is not None
        and int(lower["composition"]["EMIM"])
        >= int(upper["composition"]["EMIM"])
    ):
        upper = None

    return lower, upper


def linear_prediction(n1, rho1, n2, rho2, target):
    denom = rho2 - rho1

    if abs(denom) < 1e-12:
        return None

    pred = n2 + (target - rho2) * (n2 - n1) / denom

    if not math.isfinite(pred):
        return None

    return pred


def choose_vacuum_fill_trial(current_measurement):
    """
    真空优先补充分子，并区分真空位于电解液区域还是电极区域。

    普通电解液区域真空：
        权重 = 1.0

    电极区域真空：
        权重 = ELECTRODE_VACUUM_WEIGHT = 0.25

    因此先计算等效真空比例：

        f_eff
        = f_electrolyte
          + 0.25 * f_electrode

    再按原来的体积缺失估计：

        k_raw = N_EMIM * f_eff / (1 - f_eff)

    这样：
    - 真空全部位于电解液区域 -> 使用正常增加量；
    - 真空全部位于电极区域 -> 贡献约为正常增加量的 1/4；
    - 同时存在两类真空 -> 分别加权后合并。

    最终仍受 VACUUM_MAX_STEP_FRACTION 限制，
    且真正添加的组分继续由 1:1:5 规则控制。
    """
    current_n = int(
        current_measurement["composition"]["EMIM"]
    )

    total_f = float(
        current_measurement.get("vacuum_fraction", 0.0)
    )

    electrolyte_f = float(
        current_measurement.get(
            "electrolyte_vacuum_fraction",
            0.0,
        )
    )

    electrode_f = float(
        current_measurement.get(
            "electrode_vacuum_fraction",
            0.0,
        )
    )

    effective_f = float(
        current_measurement.get(
            "effective_vacuum_fraction",
            (
                electrolyte_f
                + ELECTRODE_VACUUM_WEIGHT
                * electrode_f
            ),
        )
    )

    if total_f <= 0.0:
        fail(
            "choose_vacuum_fill_trial 被调用，"
            "但 vacuum_fraction <= 0。"
        )

    if effective_f <= 0.0:
        # 理论上不会发生，但为了防止旧状态文件缺字段。
        effective_f = (
            total_f * ELECTRODE_VACUUM_WEIGHT
        )

    safe_f = min(effective_f, 0.95)

    raw_step = int(
        math.ceil(
            current_n
            * safe_f
            / (1.0 - safe_f)
        )
    )

    raw_step = max(1, raw_step)

    max_step = max(
        1,
        int(
            math.ceil(
                current_n
                * VACUUM_MAX_STEP_FRACTION
            )
        ),
    )

    step = min(raw_step, max_step)

    reason = (
        "vacuum weighted fill: "
        f"total={total_f:.6f}, "
        f"electrolyte={electrolyte_f:.6f}*1.0, "
        f"electrode={electrode_f:.6f}"
        f"*{ELECTRODE_VACUUM_WEIGHT:.2f}, "
        f"effective={effective_f:.6f}, "
        f"raw_step={raw_step}, "
        f"capped_step={step}"
    )

    return (
        current_measurement,
        current_n + step,
        reason,
    )


def choose_next_trial(state, current_measurement):
    target_density = float(state["target_density"])
    stage = state["stage"]

    # 真空检测优先级高于平均密度插值。
    # 只要当前体系存在局部真空，就先补分子，不使用当前点做普通 bracket 决策。
    if current_measurement.get("vacuum_detected", False):
        return choose_vacuum_fill_trial(current_measurement)

    lower, upper = bracket_measurements(state)

    # 已有 bracket：优先夹逼
    if lower is not None and upper is not None:
        n_low = int(lower["composition"]["EMIM"])
        n_high = int(upper["composition"]["EMIM"])
        rho_low = float(lower["density"])
        rho_high = float(upper["density"])

        if n_high - n_low <= 1:
            fail(
                "精调区间只剩相邻两个整数 EMIM 数量，且均未达到 fine tolerance。\n"
                f"lower: N={n_low}, rho={rho_low}\n"
                f"upper: N={n_high}, rho={rho_high}"
            )

        pred = linear_prediction(
            n_low, rho_low,
            n_high, rho_high,
            target_density,
        )

        if pred is None:
            candidate = (n_low + n_high) // 2
            reason = "bracket midpoint fallback"
        else:
            candidate = int(round(pred))
            reason = "lower-upper interpolation"

        candidate = max(n_low + 1, candidate)
        candidate = min(n_high - 1, candidate)

        return lower, candidate, reason

    # 没有 upper，只能继续从低侧向上
    current_n = int(current_measurement["composition"]["EMIM"])
    current_rho = float(current_measurement["density"])

    if current_rho >= target_density:
        fail(
            "当前密度已经高于目标，但没有找到 lower bracket。"
            "建议从更低密度的初始体系重新开始。"
        )

    all_ms = distinct_measurements_by_emim(state)
    lower_ms = [m for m in all_ms if float(m["density"]) < target_density]

    max_fraction = (
        FINE_MAX_STEP_FRACTION
        if stage == "fine"
        else COARSE_MAX_STEP_FRACTION
    )

    max_step = max(1, int(math.ceil(current_n * max_fraction)))

    if len(lower_ms) < 2:
        step = max(1, int(math.ceil(current_n * INITIAL_STEP_FRACTION)))
        step = min(step, max_step)

        return (
            current_measurement,
            current_n + step,
            "initial proportional step",
        )

    m1 = lower_ms[-2]
    m2 = lower_ms[-1]

    n1 = int(m1["composition"]["EMIM"])
    n2 = int(m2["composition"]["EMIM"])
    rho1 = float(m1["density"])
    rho2 = float(m2["density"])

    pred = linear_prediction(
        n1, rho1,
        n2, rho2,
        target_density,
    )

    if pred is None:
        candidate = current_n + max(
            1,
            int(math.ceil(current_n * INITIAL_STEP_FRACTION))
        )
        reason = "secant denominator fallback"
    else:
        candidate = int(math.ceil(pred))
        reason = "secant extrapolation"

    if candidate <= current_n:
        candidate = current_n + 1

    if candidate - current_n > max_step:
        candidate = current_n + max_step
        reason += " + step cap"

    return current_measurement, candidate, reason


def write_target_ratio(solvent_root, target_emim, target_comp):
    path = Path(solvent_root) / f".density_target_{target_emim}.json"

    path.write_text(
        json.dumps(target_comp, indent=4, ensure_ascii=False),
        encoding="utf-8",
    )

    return path


def ensure_target_system(
    solvent_root,
    solvent,
    parent_measurement,
    target_emim,
    incremental_script,
    collection_root,
):
    target_dir = Path(solvent_root) / str(target_emim)

    required = [
        target_dir / "system.gro",
        target_dir / "topol.top",
        target_dir / SUMMARY_NAME,
    ]

    if target_dir.is_dir() and all(p.is_file() for p in required):
        print(f"目标体系已存在，直接复用：{target_dir}")
        return target_dir.resolve()

    if target_dir.exists():
        # 上一次任务可能在增量建模/Packmol/topology/summary 阶段意外中断。
        # 此时数字目录已经存在，但还不是一个可用于 MD 的完整体系。
        # 不直接报错，也不静默删除，而是先整体备份后重新生成。
        backup_dir = target_dir.with_name(
            target_dir.name
            + f".incomplete_backup_{int(time.time())}"
        )

        counter = 1
        while backup_dir.exists():
            backup_dir = target_dir.with_name(
                target_dir.name
                + f".incomplete_backup_{int(time.time())}_{counter}"
            )
            counter += 1

        print(
            "检测到上一次中断留下的不完整目标体系目录："
            f"{target_dir}"
        )
        print(
            "将其备份为："
            f"{backup_dir}"
        )
        target_dir.rename(backup_dir)
        print("随后将从父体系重新执行本轮增量建模。")

    parent_dir = Path(parent_measurement["directory"]).resolve()
    parent_gro = parent_dir / "system.gro"

    if not parent_gro.is_file():
        fail(f"父体系缺少 system.gro：{parent_gro}")

    _, parent_comp = load_summary(parent_dir)

    target_comp, k = build_target_composition(
        parent_comp=parent_comp,
        target_emim=target_emim,
        solvent=solvent,
    )

    print("\n" + "=" * 72)
    print("GENERATE NEW SYSTEM")
    print("=" * 72)
    print(f"Parent          : {parent_dir}")
    print(f"Parent EMIM     : {parent_comp['EMIM']}")
    print(f"Target EMIM     : {target_emim}")
    print(f"Increment k     : {k}")

    if solvent == "ACN":
        print(f"Add             : EMIM +{k}, BF4 +{k}, ACN +{5*k}")
    else:
        print(f"Add             : EMIM +{k}, BF4 +{k}, PC +{5*k}")

    print(f"Target comp     : {target_comp}")

    ratio_file = write_target_ratio(
        solvent_root,
        target_emim,
        target_comp,
    )

    run_command(
        [
            sys.executable,
            str(incremental_script),
            str(parent_gro),
            "--ratio", str(ratio_file),
        ],
        cwd=collection_root,
    )

    if not target_dir.is_dir():
        fail(f"增量脚本结束后没有找到预期目录：{target_dir}")

    if not all(p.is_file() for p in required):
        fail(f"增量脚本生成的体系不完整：{target_dir}")

    shutil.copy2(ratio_file, target_dir / "target_ratio.json")

    return target_dir.resolve()


def process_measurement(solvent_root, state, measurement):
    error = float(measurement["relative_error"])
    abs_error = abs(error)

    print("\n" + "=" * 72)
    print("DENSITY RESULT")
    print("=" * 72)
    print(f"Directory        : {measurement['directory']}")
    print(f"Density          : {measurement['density']:.6f}")
    print(f"Target           : {state['target_density']:.6f}")
    print(f"Relative error   : {error * 100.0:.4f}%")
    print(f"Current stage    : {state['stage']}")

    vacuum_detected = bool(
        measurement.get("vacuum_detected", False)
    )

    vacuum_fraction = float(
        measurement.get("vacuum_fraction", 0.0)
    )

    print(
        f"Vacuum detected  : "
        f"{'YES' if vacuum_detected else 'NO'}"
    )
    print(
        f"Vacuum fraction  : "
        f"{vacuum_fraction * 100.0:.4f}%"
    )
    print(
        "  electrolyte    : "
        f"{float(measurement.get('electrolyte_vacuum_fraction', 0.0))*100.0:.4f}%"
    )
    print(
        "  electrode      : "
        f"{float(measurement.get('electrode_vacuum_fraction', 0.0))*100.0:.4f}%"
        f" (weight={ELECTRODE_VACUUM_WEIGHT:.2f})"
    )
    print(
        "  effective      : "
        f"{float(measurement.get('effective_vacuum_fraction', 0.0))*100.0:.4f}%"
    )
    print(
        f"Vacuum threshold : "
        f"{float(measurement.get('vacuum_threshold', 0.0)):.6f}"
    )

    if vacuum_detected:
        print(
            "Vacuum segments  : "
            f"{len(measurement.get('vacuum_segments', []))}"
        )

        for idx, seg in enumerate(
            measurement.get("vacuum_segments", []),
            start=1,
        ):
            print(
                f"  - segment {idx}: "
                f"Z={seg['z_start']:.4f}~{seg['z_end']:.4f} nm, "
                f"points={seg['points']}, "
                f"min_rho={seg['min_density']:.6f}, "
                f"mean_rho={seg['mean_density']:.6f}"
            )

        print(
            "局部真空存在：本轮禁止进入 coarse/fine 收敛判定，"
            "下一步将优先按真空比例补充分子。"
        )

        return "continue"

    # 必须先保存 coarse checkpoint，再允许进入 fine 判定
    if (
        state["stage"] == "coarse"
        and abs_error <= float(state["coarse_tolerance"])
    ):
        save_coarse_checkpoint(
            solvent_root,
            state,
            measurement,
        )

    if (
        state["coarse_checkpoint"] is not None
        and abs_error <= float(state["fine_tolerance"])
    ):
        save_final_result(
            solvent_root,
            state,
            measurement,
        )
        return "converged"

    return "continue"


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "在长时间 PBS 作业内，对单个 ACN/PC 体系按固定 1:1:5 "
            "增量自动进行密度粗调和精调。"
        )
    )

    parser.add_argument(
        "system_dir",
        help=(
            "初始体系目录，例如 qmof-305c717/ACN、"
            "qmof-305c717/PC 或 qmof-305c717/ACN/440。"
        ),
    )

    parser.add_argument(
        "--target-density",
        type=float,
        required=True,
        help="目标密度，单位与 gmx density 输出一致，通常为 kg/m^3。",
    )

    parser.add_argument(
        "--incremental-script",
        default=None,
        help=(
            "显式指定 incremental_add_molecules.py 路径。"
            "默认从 MOF/scripts/ 中寻找。"
        ),
    )

    parser.add_argument(
        "--reset-state",
        action="store_true",
        help="备份并重置已有 density_adjustment.json。",
    )

    parser.add_argument(
        "--max-iterations",
        type=int,
        default=MAX_ITERATIONS,
        help=f"最多评价体系数，默认 {MAX_ITERATIONS}。",
    )

    return parser.parse_args()


def main():
    args = parse_args()

    initial_system_dir = Path(args.system_dir).resolve()

    if not initial_system_dir.is_dir():
        fail(f"输入体系目录不存在：{initial_system_dir}")

    (
        solvent_root,
        solvent,
        mof_dir,
        collection_root,
    ) = find_solvent_root(initial_system_dir)

    if args.target_density <= 0:
        fail("--target-density 必须大于 0。")

    incremental_script = find_incremental_script(
        collection_root,
        args.incremental_script,
    )

    min_mdp, nvt_mdp = find_mdp_files(collection_root)

    check_command("gmx")

    # 初始体系合法性检查
    load_summary(initial_system_dir)

    state = load_or_create_state(
        solvent_root=solvent_root,
        solvent=solvent,
        mof_name=mof_dir.name,
        target_density=args.target_density,
        initial_system_dir=initial_system_dir,
        reset=args.reset_state,
    )

    print("=" * 72)
    print("AUTOMATIC DENSITY ADJUSTMENT")
    print("=" * 72)
    print(f"MOF                : {mof_dir.name}")
    print(f"Solvent            : {solvent}")
    print(f"Target density     : {args.target_density}")
    print(
        "Increment rule     : "
        + (
            "EMIM:BF4:ACN = 1:1:5"
            if solvent == "ACN"
            else "EMIM:BF4:PC = 1:1:5"
        )
    )
    print(
        f"Tolerance          : "
        f"coarse={COARSE_TOLERANCE*100:.1f}% "
        f"fine={FINE_TOLERANCE*100:.1f}%"
    )
    print(
        f"Density region     : "
        f"Z={DENSITY_Z_MIN_NM}~{DENSITY_Z_MAX_NM} nm"
    )
    print(f"Increment script   : {incremental_script}")
    print(f"State file         : {state_path(solvent_root)}")

    if state.get("status") == "converged":
        print("\n已有状态显示体系已经收敛。")
        print(f"结果文件：{Path(solvent_root) / FINAL_RESULT_JSON}")
        return

    current_trial = Path(
        state.get("current_trial_dir", str(initial_system_dir))
    ).resolve()

    for loop_index in range(1, args.max_iterations + 1):
        print("\n" + "#" * 72)
        print(f"# CONTROL LOOP {loop_index}")
        print("#" * 72)
        print(f"Current trial : {current_trial}")

        # 1. 当前体系完成 index/min/NVT
        ensure_simulation(
            system_dir=current_trial,
            min_mdp=min_mdp,
            nvt_mdp=nvt_mdp,
        )

        # 2. 测量密度
        measurement = measure_system(
            system_dir=current_trial,
            state=state,
            solvent_root=solvent_root,
            collection_root=collection_root,
            mof_name=mof_dir.name,
            solvent=solvent,
        )

        # 3. coarse checkpoint / fine convergence
        result = process_measurement(
            solvent_root,
            state,
            measurement,
        )

        if result == "converged":
            print("\n" + "=" * 72)
            print("DENSITY ADJUSTMENT CONVERGED")
            print("=" * 72)
            print(f"Final directory  : {measurement['directory']}")
            print(f"Final density    : {measurement['density']:.6f}")
            print(
                f"Final error      : "
                f"{measurement['relative_error']*100.0:.4f}%"
            )
            print(
                f"Final result file: "
                f"{Path(solvent_root) / FINAL_RESULT_JSON}"
            )
            return

        # 4. 预测下一体系
        parent_measurement, target_emim, reason = choose_next_trial(
            state,
            measurement,
        )

        print("\n[NEXT TRIAL]")
        print(
            f"Parent EMIM     : "
            f"{parent_measurement['composition']['EMIM']}"
        )
        print(f"Target EMIM     : {target_emim}")
        print(f"Prediction      : {reason}")

        # 5. 从父体系按严格 1:1:5 生成下一体系
        next_dir = ensure_target_system(
            solvent_root=solvent_root,
            solvent=solvent,
            parent_measurement=parent_measurement,
            target_emim=target_emim,
            incremental_script=incremental_script,
            collection_root=collection_root,
        )

        state["current_trial_dir"] = str(next_dir)
        save_state(solvent_root, state)

        current_trial = next_dir

    state["status"] = "max_iterations_reached"
    save_state(solvent_root, state)

    fail(
        f"已达到最大体系评价次数 {args.max_iterations}，"
        "但尚未达到 fine tolerance。"
    )


if __name__ == "__main__":
    main()
