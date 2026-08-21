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

fine 收敛后的 5ns 续跑校验（POST-CONVERGENCE 5NS CHECK）
-----------------------------------------------------
fine 阶段达到容差后，并不立即结束，而是继续对最终收敛体系做稳定性校验：

    1. 用 gmx convert-tpr -extend 5000 + mdrun -cpi 从 first/nvt.cpt 续算 5ns，
       生成 first/nvt_check_r{N}.tpr / .xtc / .gro / .cpt。
    2. 只对“新增的 5ns 段”跑 gmx density（-b 跳过原 NVT 段），
       生成 first/density_check_r{N}.xvg，取 Z=18~22 nm 平均密度 ρ_check。
    3. 判据：
           |ρ_check - ρ_fine| / ρ_fine > 1%   -> 重跑（再续 5ns）
           否则                              -> check_passed，记录
       其中 ρ_fine 是 fine 收敛时记录的密度。
    4. 最多重跑 CHECK_MAX_RETRIES（默认 3）轮，仍未稳定则告警退出。
    5. 第 N+1 轮基于第 N 轮的 nvt_check_rN.cpt 续算，
       即累计延伸 (N+1)*5ns，每轮独立统计该轮新增 5ns 段的密度。

可重复性（重要）
----------------
重新运行本脚本时，对已经 fine 收敛的体系不会重复调密度、不会重写
final_density_result.json，而是按 state["check_status"] 处理：

    - check_status == "passed"              直接报告完成
    - check_status == "max_retries_exceeded" 告警退出（提示人工检查）
    - check_status in (pending, running)    从 state["check_round"]+1 继续

每轮续算/密度统计的中间文件按 round 索引独立命名
（nvt_check_r1.xtc / density_check_r1.xvg ...），已存在的文件直接复用，
不会因任务中断而需要从头重跑整个 5ns check 流程。

调度器无关性
------------
本脚本本身只依赖 gmx / packmol / itptools 在 PATH 中可执行，与具体
批处理调度器（PBS / Slurm / LSF / 本地）无关。PBS 字样仅出现在示例
文档中，可替换为对应的 sbatch / bsub / 直接 python 运行。
唯一约束：调密度主循环需要在单个长任务内连续运行，依赖
density_adjustment.json 维护迭代状态。
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

# ============================================================
# 达到 fine 收敛后的“5ns 续跑 + 重新校验密度”机制
#
# 流程：
#   1. 在 fine 收敛体系上，用 gmx convert-tpr -extend + mdrun -cpi
#      从 first/nvt.cpt 续算 5ns，生成 first/nvt_check_rN.{tpr,xtc,gro,cpt}。
#   2. 只对“新增的 5ns 段”跑 gmx density（-b 跳过原 NVT 段）。
#   3. 比较 ρ_check 与 ρ_fine：
#        |ρ_check - ρ_fine| / ρ_fine > CHECK_TOLERANCE  -> 重跑（再续 5ns）
#        否则记录为 check_passed。
#   4. 最多重跑 CHECK_MAX_RETRIES 次；仍未稳定则告警退出。
#
# 中断恢复：
#   - 每轮独立的 nvt_check_rN.* 与 density_check_rN.xvg 文件，
#     重新运行脚本时已存在的文件不重跑。
#   - state["check_status"] / state["check_round"] 持久化进度。
# ============================================================
CHECK_EXTRA_NS = 5.0          # 每轮续跑的额外模拟时长（ns）
CHECK_TOLERANCE = 0.01        # |ρ_check - ρ_fine| / ρ_fine 阈值（1%）
CHECK_MAX_RETRIES = 3         # 最多重跑次数

# check 阶段产物命名前缀（按 round 索引，后缀由 check_round_files 拼接 r{N}）
CHECK_TPR_PREFIX = "nvt_check_"
CHECK_XTC_PREFIX = "nvt_check_"
CHECK_GRO_PREFIX = "nvt_check_"
CHECK_CPT_PREFIX = "nvt_check_"
CHECK_XVG_PREFIX = "density_check_"

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
        # === fine 收敛后的 5ns 续跑校验状态 ===
        # check_status:
        #   "pending"           已达 fine，尚未开始 5ns check
        #   "running"           正在执行某一轮 5ns 续算/密度统计
        #   "passed"            某轮 check 密度变化 <= 1%，已确认稳定
        #   "max_retries_exceeded" 已达最大重跑次数仍未稳定
        "check_status": "pending",
        "check_round": 0,            # 已开始/已完成的 round 数
        "fine_density": None,        # fine 收敛时记录的 ρ_fine
        "fine_directory": None,      # fine 收敛体系目录
        "check_tolerance": CHECK_TOLERANCE,
        "check_max_retries": CHECK_MAX_RETRIES,
        "check_extra_ns": CHECK_EXTRA_NS,
        "check_measurements": [],    # 每轮 check 的密度记录
    }


def ensure_check_fields(state):
    """兼容旧 state：补齐 5ns check 机制引入的字段。"""
    defaults = {
        "check_status": "pending",
        "check_round": 0,
        "fine_density": None,
        "fine_directory": None,
        "check_tolerance": CHECK_TOLERANCE,
        "check_max_retries": CHECK_MAX_RETRIES,
        "check_extra_ns": CHECK_EXTRA_NS,
        "check_measurements": [],
    }
    changed = False
    for k, v in defaults.items():
        if k not in state:
            state[k] = v
            changed = True
    return changed


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

        # 兼容旧 state：补齐 5ns check 机制引入的字段
        if ensure_check_fields(state):
            save_state(solvent_root, state)

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
        # === 5ns 续跑校验初始状态 ===
        # fine 收敛刚达成时，check 阶段尚未开始，标记为 pending。
        # 后续 run_post_convergence_check 会把这些字段逐步更新。
        "check_status": "pending",
        "check_round": 0,
        "fine_density": float(measurement["density"]),
        "fine_directory": str(snapshot_dir),
        "check_tolerance": state.get("check_tolerance", CHECK_TOLERANCE),
        "check_max_retries": state.get("check_max_retries", CHECK_MAX_RETRIES),
        "check_extra_ns": state.get("check_extra_ns", CHECK_EXTRA_NS),
        "check_measurements": [],
    }

    path = Path(solvent_root) / FINAL_RESULT_JSON

    path.write_text(
        json.dumps(result, indent=4, ensure_ascii=False),
        encoding="utf-8",
    )

    state["status"] = "converged"
    state["current_trial_dir"] = measurement["directory"]
    # 同步 state 中的 fine_density / fine_directory / check_status
    # fine_directory 指向 fine/ 快照目录（续算在此进行），
    # 不是原始 trial 目录（measurement["directory"]）。
    state["fine_density"] = float(measurement["density"])
    state["fine_directory"] = str(snapshot_dir)
    state["check_status"] = "pending"
    state["check_round"] = 0
    state["check_measurements"] = []
    save_state(solvent_root, state)


# ============================================================
# 5ns 续跑 + 重新校验密度：辅助函数
# ============================================================

def parse_mdp_nsteps_dt(mdp_path):
    """
    从 GROMACS mdp 中读取 nsteps 与 dt，返回 (nsteps, dt_ps)。

    用于推算原 NVT 的总时长 T0 = nsteps * dt（ps），
    作为续算后只统计新增 5ns 段密度的 -b 起算时间。

    mdp 中 nsteps 通常为整数；dt 可能是 0.001 / 0.002 / 0.005 等。
    """
    mdp_path = Path(mdp_path)
    if not mdp_path.is_file():
        fail(f"找不到 mdp 文件：{mdp_path}")

    nsteps = None
    dt = None

    with open(mdp_path, "r", encoding="utf-8", errors="replace") as f:
        for raw in f:
            line = raw.split(";", 1)[0].strip()
            if not line or "=" not in line:
                continue
            key, _, value = line.partition("=")
            key = key.strip().lower()
            value = value.strip()

            if key == "nsteps":
                try:
                    nsteps = int(value)
                except ValueError:
                    fail(f"mdp 中 nsteps 非法：{value!r}")
            elif key == "dt":
                try:
                    dt = float(value)
                except ValueError:
                    fail(f"mdp 中 dt 非法：{value!r}")

            if nsteps is not None and dt is not None:
                break

    if nsteps is None:
        fail(f"mdp 中未找到 nsteps：{mdp_path}")
    if dt is None:
        fail(f"mdp 中未找到 dt：{mdp_path}")
    if nsteps <= 0:
        fail(f"mdp 中 nsteps 必须为正：{nsteps}")
    if dt <= 0:
        fail(f"mdp 中 dt 必须为正：{dt}")

    return nsteps, dt


def original_nvt_end_ps(nvt_mdp):
    """
    计算原 NVT tpr 的总时长（ps）。

    若 mdp 中 nsteps = -1（无限制），则取一个保守的默认值
    DENSITY_BEGIN_PS（7000 ps）作为兜底，并打印警告。
    """
    nsteps, dt = parse_mdp_nsteps_dt(nvt_mdp)
    if nsteps < 0:
        warn(
            f"nvt.mdp 中 nsteps={nsteps}（无限制），"
            f"使用默认 {DENSITY_BEGIN_PS} ps 作为原 NVT 时长。"
        )
        return float(DENSITY_BEGIN_PS)
    return nsteps * dt


def check_round_files(system_dir, round_idx):
    """
    返回第 round_idx 轮 check 的所有产物路径。

    每轮独立命名，便于中断后恢复：
        first/nvt_check_r1.tpr / .xtc / .gro / .cpt
        first/density_check_r1.xvg
    """
    first = Path(system_dir) / "first"
    suffix = f"r{round_idx}"
    return {
        "tpr":  first / f"{CHECK_TPR_PREFIX}{suffix}.tpr",
        "xtc":  first / f"{CHECK_XTC_PREFIX}{suffix}.xtc",
        "gro":  first / f"{CHECK_GRO_PREFIX}{suffix}.gro",
        "cpt":  first / f"{CHECK_CPT_PREFIX}{suffix}.cpt",
        "xvg":  first / f"{CHECK_XVG_PREFIX}{suffix}.xvg",
    }


def extend_nvt_for_check(
    system_dir,
    nvt_mdp,
    round_idx,
    extra_ns,
    base_cpt=None,
):
    """
    在 fine 收敛体系的 NVT 基础上续算 extra_ns 纳秒。

    使用 gmx convert-tpr -extend 把 tpr 的总步数延长 extra_ns/dt 步，
    再用 mdrun -cpi 从上一轮的 cpt（或原 nvt.cpt）续算。

    参数：
        system_dir : fine 收敛体系目录
        nvt_mdp    : 原 nvt.mdp（用于解析 dt）
        round_idx  : 当前轮次（1 起）
        extra_ns   : 本轮续算时长（ns）
        base_cpt   : 续算起点的 cpt；None 时按 round_idx 自动选择
                     round 1 -> first/nvt.cpt
                     round N -> first/nvt_check_r(N-1).cpt

    可重复性：
        若本轮 tpr+xtc+gro+cpt 都已存在，直接跳过续算；
        若只 tpr 缺失但其它都存在，视为异常，fail。
    """
    system_dir = Path(system_dir)
    first = system_dir / "first"
    first.mkdir(exist_ok=True)

    paths = check_round_files(system_dir, round_idx)

    # 已完整生成则跳过
    if (
        paths["tpr"].is_file()
        and paths["xtc"].is_file()
        and paths["gro"].is_file()
        and paths["cpt"].is_file()
    ):
        print(
            f"第 {round_idx} 轮 NVT 续算已完成，跳过："
            f"{paths['xtc']}"
        )
        return paths

    # 解析 dt
    _, dt = parse_mdp_nsteps_dt(nvt_mdp)
    extend_ps = float(extra_ns) * 1000.0
    extend_steps = int(round(extend_ps / dt))

    # 自动选择 base cpt
    if base_cpt is None:
        if round_idx == 1:
            base_cpt = first / "nvt.cpt"
        else:
            prev = check_round_files(system_dir, round_idx - 1)
            base_cpt = prev["cpt"]
    else:
        base_cpt = Path(base_cpt)

    if not base_cpt.is_file():
        fail(
            f"找不到续算起点 cpt：{base_cpt}\n"
            f"请确认原 NVT 已完成（first/nvt.cpt 存在）。"
        )

    # gmx convert-tpr -s first/nvt.tpr -extend <ps> -o first/nvt_check_rN.tpr
    # 注意：convert-tpr 的 -extend 单位是 ps。
    if not paths["tpr"].is_file():
        print(
            f"\n[convert-tpr] 续算 {extend_ps:.1f} ps "
            f"({extra_ns:.1f} ns / dt={dt}) -> {paths['tpr']}"
        )
        run_command(
            [
                "gmx", "convert-tpr",
                "-s", "first/nvt.tpr",
                "-extend", f"{extend_ps:.1f}",
                "-o", str(paths["tpr"].relative_to(system_dir)),
            ],
            cwd=system_dir,
        )
    else:
        print(f"第 {round_idx} 轮 tpr 已存在，跳过 convert-tpr：{paths['tpr']}")

    # gmx mdrun -s first/nvt_check_rN.tpr -cpi <base_cpt>
    #            -deffnm first/nvt_check_rN
    # -cpi 必须用相对于 cwd(system_dir) 的路径，否则会找不到 checkpoint
    # 而从 tpr 起始 step 重新开始模拟（即从头跑完整时长，而非续算）。
    try:
        cpt_arg = str(base_cpt.relative_to(system_dir))
    except ValueError:
        # base_cpt 不在 system_dir 下，用绝对路径
        cpt_arg = str(base_cpt)

    cmd = [
        "gmx", "mdrun",
        "-s", f"first/{CHECK_TPR_PREFIX}r{round_idx}.tpr",
        "-cpi", cpt_arg,
        "-deffnm", f"first/{CHECK_TPR_PREFIX}r{round_idx}",
        "-v",
        "-ntomp", str(NTOMP),
        "-ntmpi", str(NTMPI),
    ]

    print(
        f"[mdrun] 续算 {extend_ps:.1f} ps，"
        f"基于 cpt {base_cpt} -> {paths['gro']}"
    )
    run_command(cmd, cwd=system_dir)

    for key in ("gro", "xtc", "cpt"):
        if not paths[key].is_file():
            fail(
                f"第 {round_idx} 轮 NVT 续算结束后未生成：{paths[key]}"
            )

    return paths


def run_check_density(system_dir, round_idx, begin_ps):
    """
    对第 round_idx 轮续算得到的 nvt_check_rN.xtc 跑 gmx density，
    只统计 t >= begin_ps 之后的密度。

    自动选择 group 6（与原 measure_system 保持一致）。

    可重复性：
        若 density_check_rN.xvg 已存在，直接复用。
    """
    system_dir = Path(system_dir)
    first = system_dir / "first"
    paths = check_round_files(system_dir, round_idx)

    xtc = paths["xtc"]
    tpr = paths["tpr"]
    xvg = paths["xvg"]
    index = system_dir / "index.ndx"

    if xvg.is_file():
        print(f"第 {round_idx} 轮 check 密度已存在，复用：{xvg}")
        return xvg

    for p in (xtc, tpr, index):
        if not p.is_file():
            fail(f"check density 输入文件不存在：{p}")

    run_command(
        [
            "gmx", "density",
            "-f", f"first/{CHECK_XTC_PREFIX}r{round_idx}.xtc",
            "-s", f"first/{CHECK_TPR_PREFIX}r{round_idx}.tpr",
            "-n", "../index.ndx",
            "-sl", str(DENSITY_SLICES),
            "-o", f"first/{CHECK_XVG_PREFIX}r{round_idx}.xvg",
            "-d", "Z",
            "-b", f"{begin_ps:.1f}",
        ],
        cwd=first,
        stdin_text=f"{DENSITY_GROUP}\n",
    )

    if not xvg.is_file():
        fail(f"gmx density 后未生成：{xvg}")

    return xvg


def evaluate_check(fine_density, check_density, tolerance):
    """
    判断 5ns 续跑后密度是否相对 fine 收敛密度变化超过 tolerance。

    返回 (passed, relative_change)。
    """
    if fine_density is None or fine_density == 0:
        fail("fine_density 缺失或为 0，无法评估 check 密度变化。")

    relative_change = abs(
        float(check_density) - float(fine_density)
    ) / float(fine_density)

    return relative_change <= float(tolerance), relative_change


def save_check_result(
    solvent_root,
    state,
    round_idx,
    check_density,
    begin_ps,
    relative_change,
    passed,
    final=False,
):
    """
    把第 round_idx 轮 check 结果写回 final_density_result.json 与 state。

    final=True 表示该轮通过（passed=True）或已超过最大重试次数，
    会同步更新 state["check_status"]。
    """
    measurement = {
        "round": round_idx,
        "measured_at": now_string(),
        "directory": state["fine_directory"],
        "xtc": f"first/{CHECK_XTC_PREFIX}r{round_idx}.xtc",
        "xvg": f"first/{CHECK_XVG_PREFIX}r{round_idx}.xvg",
        "begin_ps": round(begin_ps, 3),
        "density": float(check_density),
        "fine_density": float(state["fine_density"]),
        "relative_change": float(relative_change),
        "tolerance": float(state["check_tolerance"]),
        "passed": bool(passed),
    }

    state["check_measurements"].append(measurement)
    state["check_round"] = round_idx

    if passed:
        state["check_status"] = "passed"
    elif final:
        state["check_status"] = "max_retries_exceeded"
    else:
        state["check_status"] = "running"

    save_state(solvent_root, state)

    # 同步写回 final_density_result.json
    result_path = Path(solvent_root) / FINAL_RESULT_JSON
    if result_path.is_file():
        try:
            result = json.loads(
                result_path.read_text(encoding="utf-8")
            )
        except Exception:
            result = {}
    else:
        result = {}

    result["check_status"] = state["check_status"]
    result["check_round"] = state["check_round"]
    result["check_tolerance"] = state["check_tolerance"]
    result["check_max_retries"] = state["check_max_retries"]
    result["check_extra_ns"] = state["check_extra_ns"]
    result["fine_density"] = state["fine_density"]
    result["check_measurements"] = state["check_measurements"]
    result["last_check"] = measurement

    result_path.write_text(
        json.dumps(result, indent=4, ensure_ascii=False),
        encoding="utf-8",
    )

    print("\n" + "=" * 72)
    print(f"CHECK ROUND {round_idx} {'PASSED' if passed else 'RETRY'}")
    print("=" * 72)
    print(f"Fine density    : {state['fine_density']:.6f}")
    print(f"Check density   : {check_density:.6f}")
    print(f"Relative change : {relative_change * 100.0:.4f}%")
    print(f"Tolerance       : {state['check_tolerance'] * 100.0:.2f}%")
    print(f"Check status    : {state['check_status']}")

    return measurement


def run_post_convergence_check(
    solvent_root,
    state,
    nvt_mdp,
    collection_root,
    mof_name,
    solvent,
):
    """
    fine 收敛后的 5ns 续跑 + 重新校验密度主循环。

    状态机：
        check_status == "passed"           -> 直接返回，不重跑
        check_status == "max_retries_exceeded" -> 告警退出
        check_status in (pending, running)  -> 从 check_round 开始继续

    每轮：
        1. 调 extend_nvt_for_check 续算 5ns
        2. 调 run_check_density 统计新增 5ns 段密度
        3. |ρ_check - ρ_fine| / ρ_fine > 1% -> 下一轮（基于本轮 cpt 续算）
        4. 否则 -> check_passed，更新 density_profiles.csv 与 final_density_result.json
        5. 最多 CHECK_MAX_RETRIES 轮
    """
    check_status = state.get("check_status", "pending")
    max_retries = int(state.get("check_max_retries", CHECK_MAX_RETRIES))

    if check_status == "passed":
        print("\n" + "=" * 72)
        print("POST-CONVERGENCE CHECK ALREADY PASSED")
        print("=" * 72)
        print(f"Check rounds    : {state.get('check_round', 0)}")
        print(f"Fine density    : {state.get('fine_density')}")
        if state.get("check_measurements"):
            last = state["check_measurements"][-1]
            print(f"Last check rho  : {last['density']:.6f}")
            print(
                f"Last change     : "
                f"{last['relative_change'] * 100.0:.4f}%"
            )
        print(f"Result file     : {Path(solvent_root) / FINAL_RESULT_JSON}")
        return

    if check_status == "max_retries_exceeded":
        warn(
            "此前 5ns check 已达到最大重试次数仍未稳定。\n"
            f"如需重新开始 check，请删除或编辑 "
            f"{Path(solvent_root) / FINAL_RESULT_JSON} 与 "
            f"{Path(solvent_root) / STATE_JSON} 中的 check_* 字段。"
        )
        return

    # 兼容旧 state（改造前已收敛的体系）：fine_density / fine_directory 缺失。
    # 优先从 final_density_result.json 回填，其次从 measurements 最后一条回填。
    # fine_directory 必须指向 <solvent_root>/fine/ 快照目录（续算在此进行），
    # 不是原始 trial 目录。
    if state.get("fine_density") is None or state.get("fine_directory") is None:
        backfilled = False

        # 1. 从 final_density_result.json 回填
        result_path = Path(solvent_root) / FINAL_RESULT_JSON
        if result_path.is_file():
            try:
                result = json.loads(
                    result_path.read_text(encoding="utf-8")
                )
            except Exception:
                result = {}

            if state.get("fine_density") is None:
                rho = result.get("density")
                if rho is not None:
                    try:
                        state["fine_density"] = float(rho)
                        backfilled = True
                    except (TypeError, ValueError):
                        pass

            # fine_directory 优先用 snapshot_directory（指向 fine/ 快照），
            # 而不是 directory（指向原始 trial 目录）。
            if state.get("fine_directory") is None:
                d = result.get("snapshot_directory")
                if d:
                    state["fine_directory"] = str(d)
                    backfilled = True

        # 2. 从 measurements 最后一条回填 fine_density
        if state.get("fine_density") is None and state.get("measurements"):
            last = state["measurements"][-1]
            rho = last.get("density")
            if rho is not None:
                try:
                    state["fine_density"] = float(rho)
                    backfilled = True
                except (TypeError, ValueError):
                    pass

        # 3. fine_directory 兜底：从 solvent_root 拼 <solvent_root>/fine/
        if state.get("fine_directory") is None:
            candidate = Path(solvent_root) / FINE_SNAPSHOT_DIR
            if candidate.is_dir():
                state["fine_directory"] = str(candidate.resolve())
                backfilled = True

        if backfilled:
            save_state(solvent_root, state)
            print(
                "\n[backfill] 检测到旧 state 缺 fine_density/fine_directory，"
                "已从 final_density_result.json 或 solvent_root/fine/ 回填："
            )
            print(f"  fine_density    = {state.get('fine_density')}")
            print(f"  fine_directory  = {state.get('fine_directory')}")

    if state.get("fine_density") is None:
        fail(
            "state.fine_density 缺失，且无法从 final_density_result.json "
            "或 measurements 回填。请检查 "
            f"{Path(solvent_root) / FINAL_RESULT_JSON} 与 "
            f"{Path(solvent_root) / STATE_JSON}。"
        )

    fine_dir = Path(state["fine_directory"]) if state.get("fine_directory") else None
    if fine_dir is None or not fine_dir.is_dir():
        fail(
            f"fine 收敛体系目录不存在或未指定：{fine_dir}\n"
            "请检查 final_density_result.json 中的 directory 字段。"
        )

    # 原始 NVT 总时长（ps），作为第 1 轮 check 的 -b 起算时间
    t0 = original_nvt_end_ps(nvt_mdp)
    extra_ns = float(state.get("check_extra_ns", CHECK_EXTRA_NS))
    extra_ps = extra_ns * 1000.0

    # 校验：T0 应不小于原 measure_system 使用的 DENSITY_BEGIN_PS，
    # 否则说明原 NVT 比预期短，统计区间与原 fine 不一致。
    if t0 < DENSITY_BEGIN_PS:
        warn(
            f"原 NVT 总时长 {t0:.1f} ps 小于原密度统计起始 "
            f"{DENSITY_BEGIN_PS} ps。这可能导致 check 段密度统计区间"
            "与原 fine 收敛时不一致，请检查 nvt.mdp 的 nsteps。"
        )

    print("\n" + "=" * 72)
    print("POST-CONVERGENCE 5NS CHECK")
    print("=" * 72)
    print(f"Fine directory  : {fine_dir}")
    print(f"Fine density     : {state['fine_density']:.6f}")
    print(f"Original NVT end : {t0:.1f} ps")
    print(f"Extra per round  : {extra_ns:.1f} ns ({extra_ps:.1f} ps)")
    print(f"Tolerance        : {state['check_tolerance'] * 100.0:.2f}%")
    print(f"Max retries      : {max_retries}")

    # 从 check_round 继续；若 check_round=0 表示尚未开始，从 1 开始
    start_round = int(state.get("check_round", 0)) + 1

    for round_idx in range(start_round, max_retries + 1):
        print("\n" + "#" * 72)
        print(f"# CHECK ROUND {round_idx} / {max_retries}")
        print("#" * 72)

        # 1. 续算 5ns
        paths = extend_nvt_for_check(
            system_dir=fine_dir,
            nvt_mdp=nvt_mdp,
            round_idx=round_idx,
            extra_ns=extra_ns,
        )

        # 2. 跑密度：begin = t0 + (round_idx - 1) * extra_ps
        #    第 1 轮只统计原 NVT 之后的 5ns 段；
        #    第 N 轮只统计第 N 轮新增的 5ns 段（基于上轮 cpt 续算）。
        begin_ps = t0 + (round_idx - 1) * extra_ps

        xvg = run_check_density(
            system_dir=fine_dir,
            round_idx=round_idx,
            begin_ps=begin_ps,
        )

        # 3. 解析密度（只取 18~22 nm 区间，与 measure_system 一致）
        summary, _ = load_summary(fine_dir)
        density_analysis = analyze_density_profile(
            xvg_path=xvg,
            summary=summary,
            target_density=float(state["target_density"]),
        )
        check_density = density_analysis["bulk_density"]

        # 4. 判断
        passed, relative_change = evaluate_check(
            fine_density=float(state["fine_density"]),
            check_density=check_density,
            tolerance=state["check_tolerance"],
        )

        is_final = passed or (round_idx >= max_retries)
        save_check_result(
            solvent_root=solvent_root,
            state=state,
            round_idx=round_idx,
            check_density=check_density,
            begin_ps=begin_ps,
            relative_change=relative_change,
            passed=passed,
            final=is_final,
        )

        # 同步更新全局 density_profiles.csv，便于跨体系查看
        update_global_density_csv(
            collection_root=collection_root,
            mof_name=mof_name,
            solvent=solvent,
            system_dir=fine_dir,
            system_name=str(state.get("check_round", 0)),
            density=check_density,
            points=density_analysis["bulk_points"],
        )

        if passed:
            print("\n" + "=" * 72)
            print("POST-CONVERGENCE CHECK PASSED")
            print("=" * 72)
            print(
                f"5ns 续跑后密度变化 {relative_change * 100.0:.4f}% "
                f"<= {state['check_tolerance'] * 100.0:.2f}%，"
                "体系已确认稳定。"
            )
            print(f"Result file     : {Path(solvent_root) / FINAL_RESULT_JSON}")
            return

    # 已达最大重试仍未通过
    print("\n" + "=" * 72)
    print("POST-CONVERGENCE CHECK MAX RETRIES EXCEEDED")
    print("=" * 72)
    warn(
        f"已续跑 {max_retries} 轮 5ns（共 {max_retries * extra_ns:.1f} ns），"
        "密度仍未相对 fine 收敛值稳定在 1% 以内。\n"
        "结果已记录为 max_retries_exceeded，请人工检查体系。"
    )
    print(f"Result file     : {Path(solvent_root) / FINAL_RESULT_JSON}")


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
        # 已收敛的体系，重新运行时进入 5ns check 流程：
        #   - check 已通过：直接报告完成
        #   - check 已达最大重试：告警退出
        #   - check 未完成：从中断点继续
        print("\n已有状态显示体系已经收敛，进入 5ns 续跑校验流程。")
        print(f"结果文件：{Path(solvent_root) / FINAL_RESULT_JSON}")
        run_post_convergence_check(
            solvent_root=solvent_root,
            state=state,
            nvt_mdp=nvt_mdp,
            collection_root=collection_root,
            mof_name=mof_dir.name,
            solvent=solvent,
        )
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

            # fine 收敛后，进入 5ns 续跑 + 重新校验密度流程
            run_post_convergence_check(
                solvent_root=solvent_root,
                state=state,
                nvt_mdp=nvt_mdp,
                collection_root=collection_root,
                mof_name=mof_dir.name,
                solvent=solvent,
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
