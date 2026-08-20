#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
incremental_add_molecules.py

以任意名称的已有 GRO 文件作为母体系，保留其中全部原子坐标，
仅增加目标组成相对于当前组成多出来的 EMIM/BF4/ACN/PC 分子。

核心原则：
1. 输入 GRO 文件名不固定，由命令行显式指定。
2. 当前总分子数只读取同目录 summary 的 electrolyte_molecules。
3. 新 summary 只保存最终体系汇总状态，不保存会影响下一轮判断的增量历史。
4. 父 GRO 转为 PDB 后作为 Packmol fixed structure；Packmol 只添加 delta molecules。
5. 新 topology 保留父 [ molecules ] 顺序，并在其后追加本轮新增 molecule 条目，
   以匹配 Packmol 输出坐标顺序。
6. 父 topol.top 中的相对 #include 路径会重写为相对于新子目录的路径。
7. 新目录默认命名为目标 EMIM 数量。

典型用法：
python incremental_add_molecules.py ./any_name.gro --ratio new_mix_ratio.json

可选显式指定伴随文件：
python incremental_add_molecules.py ./any_name.gro \
    --ratio new_mix_ratio.json \
    --summary system_summary.json \
    --inp packmol_input.inp \
    --topol topol.top
"""

import argparse
import copy
import json
import math
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path


MOLECULE_ORDER = ["EMIM", "BF4", "ACN", "PC"]
DEFAULT_SUMMARY = "system_summary.json"
DEFAULT_INP = "packmol_input.inp"
DEFAULT_TOPOL = "topol.top"
DEFAULT_OUTPUT_GRO = "system.gro"
PACKMOL_TIMEOUT = 3600


def fail(message, code=1):
    print(f"错误：{message}")
    sys.exit(code)


def warn(message):
    print(f"警告：{message}")


def norm_posix(path):
    return os.path.normpath(str(path)).replace("\\", "/")


def run_command(args, cwd=None, stdin_path=None, stdout_path=None, timeout=None):
    """不使用 shell=True 执行外部命令。"""
    stdin_handle = None
    stdout_handle = None
    try:
        if stdin_path is not None:
            stdin_handle = open(stdin_path, "rb")
        if stdout_path is not None:
            stdout_handle = open(stdout_path, "wb")

        return subprocess.run(
            args,
            cwd=cwd,
            stdin=stdin_handle,
            stdout=stdout_handle if stdout_handle else subprocess.PIPE,
            stderr=subprocess.STDOUT if stdout_handle else subprocess.PIPE,
            check=True,
            timeout=timeout,
        )
    except FileNotFoundError:
        fail(f"找不到外部命令：{args[0]}。请确认已经安装并加入 PATH。")
    except subprocess.TimeoutExpired:
        fail(f"命令运行超时：{' '.join(args)}")
    except subprocess.CalledProcessError as exc:
        if stdout_handle:
            detail = f"请检查日志文件：{stdout_path}"
        else:
            out = (exc.stdout or b"").decode(errors="replace")
            err = (exc.stderr or b"").decode(errors="replace")
            detail = (out + "\n" + err).strip()
        fail(f"命令执行失败：{' '.join(args)}\n{detail}")
    finally:
        if stdin_handle:
            stdin_handle.close()
        if stdout_handle:
            stdout_handle.close()


def choose_companion_file(parent_dir, explicit, default_name, suffix, description):
    """优先显式指定；否则默认名；最后尝试同目录唯一候选。"""
    parent_dir = Path(parent_dir)

    if explicit:
        p = Path(explicit)
        if not p.is_absolute():
            p = parent_dir / p
        p = p.resolve()
        if not p.is_file():
            fail(f"指定的{description}不存在：{p}")
        return p

    default_path = parent_dir / default_name
    if default_path.is_file():
        return default_path.resolve()

    candidates = [
        p for p in parent_dir.iterdir()
        if p.is_file() and p.name.lower().endswith(suffix.lower())
    ]

    if len(candidates) == 1:
        print(f"未找到默认 {default_name}，自动使用：{candidates[0].name}")
        return candidates[0].resolve()

    if not candidates:
        fail(f"在 {parent_dir} 中找不到{description}，请使用命令行参数显式指定。")

    fail(
        f"在 {parent_dir} 中找到多个可能的{description}："
        + ", ".join(p.name for p in candidates)
        + "。请使用命令行参数显式指定。"
    )


def load_target_ratio(path):
    """读取目标最终分子总数；支持对象或只含一个对象的列表。"""
    path = Path(path).resolve()
    if not path.is_file():
        fail(f"目标配比文件不存在：{path}")

    try:
        data = json.loads(path.read_text(encoding="utf-8"))
    except Exception as exc:
        fail(f"目标配比 JSON 解析失败：{exc}")

    if isinstance(data, list):
        if len(data) != 1 or not isinstance(data[0], dict):
            fail("ratio 若使用列表格式，必须只包含一个 JSON 对象。")
        data = data[0]

    if not isinstance(data, dict):
        fail("目标配比文件必须是 JSON 对象。")

    missing = [mol for mol in MOLECULE_ORDER if mol not in data]
    if missing:
        fail("目标配比缺少键：" + ", ".join(missing))

    result = {}
    for mol in MOLECULE_ORDER:
        value = data[mol]
        if isinstance(value, bool) or not isinstance(value, int) or value < 0:
            fail(f"{mol} 的目标数量必须为非负整数，当前值：{value!r}")
        result[mol] = value

    return result


def load_parent_summary(path):
    path = Path(path)
    try:
        data = json.loads(path.read_text(encoding="utf-8"))
    except Exception as exc:
        fail(f"summary 解析失败：{path}\n{exc}")

    if "electrolyte_molecules" not in data:
        fail(f"summary 缺少 electrolyte_molecules：{path}")

    current = data["electrolyte_molecules"]
    if not isinstance(current, dict):
        fail("summary 中 electrolyte_molecules 必须是对象。")

    normalized = {}
    for mol in MOLECULE_ORDER:
        value = current.get(mol, 0)
        if isinstance(value, bool) or not isinstance(value, int) or value < 0:
            fail(f"summary 中 {mol} 数量非法：{value!r}")
        normalized[mol] = value

    if "box_dimensions_nm" not in data:
        fail("summary 缺少 box_dimensions_nm。")

    if "packmol_insertion_regions_nm" not in data:
        fail("summary 缺少 packmol_insertion_regions_nm。")

    regions = data["packmol_insertion_regions_nm"]
    if not isinstance(regions, list) or len(regions) != 3:
        fail("summary 中 packmol_insertion_regions_nm 必须包含三个区域。")

    for idx, region in enumerate(regions, start=1):
        for key in ("z_low", "z_high", "thickness"):
            if key not in region:
                fail(f"Packmol 区域 {idx} 缺少 {key}。")

    return data, normalized


def compute_additions(current, target):
    """计算目标总数减当前总数；只允许增加。"""
    added = {}
    for mol in MOLECULE_ORDER:
        delta = target[mol] - current[mol]
        if delta < 0:
            fail(
                f"{mol} 当前数量为 {current[mol]}，目标数量为 {target[mol]}。"
                "本工具只支持增添分子，不支持删除已有分子。"
            )
        added[mol] = delta

    if sum(added.values()) == 0:
        fail("目标组成与当前体系完全相同，没有需要新增的分子。")

    if target["EMIM"] != target["BF4"]:
        warn(
            f"目标 EMIM={target['EMIM']}，BF4={target['BF4']}。"
            "如果二者分别为 +1/-1 离子，请确认是否有意构建非电中性体系。"
        )

    return added


_STRUCTURE_RE = re.compile(r"^\s*structure\s+(.+?)\s*$", re.IGNORECASE)
_NUMBER_RE = re.compile(r"^\s*number\s+(\d+)\s*$", re.IGNORECASE)
_INSIDE_RE = re.compile(
    r"^\s*inside\s+box\s+"
    r"([-+0-9.eE]+)\s+([-+0-9.eE]+)\s+([-+0-9.eE]+)\s+"
    r"([-+0-9.eE]+)\s+([-+0-9.eE]+)\s+([-+0-9.eE]+)\s*$",
    re.IGNORECASE,
)
_TOL_RE = re.compile(r"^\s*tolerance\s+([-+0-9.eE]+)\s*$", re.IGNORECASE)


def parse_packmol_input(inp_path):
    """
    旧 inp 只提供可继承的 Packmol 参数和 PDB 路径线索。
    不把旧 inp 中的 number 当成当前总数量。
    """
    inp_path = Path(inp_path)
    lines = inp_path.read_text(encoding="utf-8", errors="replace").splitlines()

    tolerance = None
    blocks = []
    current = None

    for raw in lines:
        m = _TOL_RE.match(raw)
        if m:
            tolerance = float(m.group(1))
            continue

        m = _STRUCTURE_RE.match(raw)
        if m:
            current = {
                "path": m.group(1).strip().strip('"').strip("'"),
                "number": None,
                "inside_box": None,
            }
            blocks.append(current)
            continue

        if current is None:
            continue

        m = _NUMBER_RE.match(raw)
        if m:
            current["number"] = int(m.group(1))
            continue

        m = _INSIDE_RE.match(raw)
        if m:
            current["inside_box"] = tuple(float(m.group(i)) for i in range(1, 7))
            continue

        if raw.strip().lower() == "end structure":
            current = None

    if tolerance is None:
        warn("旧 Packmol inp 中未找到 tolerance，将使用 2.0 Å。")
        tolerance = 2.0

    inside_blocks = [b for b in blocks if b["inside_box"] is not None]
    if not inside_blocks:
        fail("旧 Packmol inp 中没有找到任何 inside box，无法继承 X/Y 插入范围。")

    first_box = inside_blocks[0]["inside_box"]
    x_low, y_low, _, x_high, y_high, _ = first_box

    molecule_paths = {}
    for block in blocks:
        stem = Path(block["path"]).stem.upper()
        if stem in MOLECULE_ORDER:
            molecule_paths[stem] = block["path"]

    return {
        "tolerance_a": tolerance,
        "x_low_a": x_low,
        "x_high_a": x_high,
        "y_low_a": y_low,
        "y_high_a": y_high,
        "molecule_paths": molecule_paths,
    }


def resolve_old_relative_path(parent_dir, path_text):
    p = Path(path_text)
    if p.is_absolute():
        return p.resolve()
    return (Path(parent_dir) / p).resolve()


def find_molecule_pdb(parent_dir, mol, packmol_info):
    """优先沿用旧 inp 路径；否则由同级 PDB 或 basic/str 推断。"""
    parent_dir = Path(parent_dir).resolve()

    old = packmol_info["molecule_paths"].get(mol)
    if old:
        p = resolve_old_relative_path(parent_dir, old)
        if p.is_file():
            return p

    for old_path in packmol_info["molecule_paths"].values():
        known = resolve_old_relative_path(parent_dir, old_path)
        candidate = known.parent / f"{mol}.pdb"
        if candidate.is_file():
            return candidate.resolve()

    cur = parent_dir
    for _ in range(8):
        candidate = cur / "basic" / "str" / f"{mol}.pdb"
        if candidate.is_file():
            return candidate.resolve()
        if cur.parent == cur:
            break
        cur = cur.parent

    fail(
        f"无法找到 {mol}.pdb。旧 inp 没有可用路径，"
        "且父目录链中也没有找到 basic/str。"
    )


def largest_remainder_allocate(total, weights):
    if total == 0:
        return [0] * len(weights)

    total_weight = sum(weights)
    if total_weight <= 0:
        fail("三个 Packmol 插入区域总厚度不大于 0。")

    raw = [total * w / total_weight for w in weights]
    base = [math.floor(x) for x in raw]
    remainder = total - sum(base)

    order = sorted(
        range(len(weights)),
        key=lambda i: raw[i] - base[i],
        reverse=True,
    )

    for i in range(remainder):
        base[order[i % len(order)]] += 1

    return base


def allocate_additions(added, summary):
    """仅对本轮新增数量按三个插入区有效厚度比例分配。"""
    regions = summary["packmol_insertion_regions_nm"]
    weights = [float(r["thickness"]) for r in regions]
    return {
        mol: largest_remainder_allocate(added[mol], weights)
        for mol in MOLECULE_ORDER
    }


def determine_series_root(parent_dir):
    """数字父目录的下一代输出为其兄弟目录；否则输出到当前父目录下。"""
    parent_dir = Path(parent_dir).resolve()
    if parent_dir.name.isdigit():
        return parent_dir.parent
    return parent_dir


def prepare_output_dir(series_root, target_emim, overwrite=False):
    out_dir = Path(series_root) / str(target_emim)

    if out_dir.exists():
        if not overwrite:
            fail(
                f"目标目录已存在：{out_dir}\n"
                "为避免覆盖已有体系，请修改目标 EMIM 数量，或使用 --overwrite。"
            )
        shutil.rmtree(out_dir)

    out_dir.mkdir(parents=True, exist_ok=False)
    return out_dir.resolve()


def convert_parent_gro_to_pdb(input_gro, out_dir):
    output_pdb = Path(out_dir) / "parent_system.pdb"

    run_command([
        "gmx", "editconf",
        "-f", str(Path(input_gro).resolve()),
        "-o", str(output_pdb),
    ])

    if not output_pdb.is_file():
        fail(f"gmx editconf 未生成：{output_pdb}")

    return output_pdb


def write_incremental_packmol_input(
    out_dir,
    parent_pdb,
    summary,
    packmol_info,
    assignment,
    molecule_pdbs,
):
    out_dir = Path(out_dir)
    regions = summary["packmol_insertion_regions_nm"]

    lines = [
        f"tolerance {packmol_info['tolerance_a']:.6f}",
        "output merged_system.pdb",
        "filetype pdb",
        "",
        f"structure {Path(parent_pdb).name}",
        "  number 1",
        "  fixed 0.0 0.0 0.0 0.0 0.0 0.0",
        "end structure",
        "",
    ]

    for mol in MOLECULE_ORDER:
        mol_counts = assignment[mol]
        if sum(mol_counts) == 0:
            continue

        pdb_rel = norm_posix(os.path.relpath(molecule_pdbs[mol], out_dir))

        for region_idx, count in enumerate(mol_counts):
            if count <= 0:
                continue

            region = regions[region_idx]
            z_low_a = float(region["z_low"]) * 10.0
            z_high_a = float(region["z_high"]) * 10.0

            lines.extend([
                f"structure {pdb_rel}",
                f"  number {count}",
                (
                    "  inside box "
                    f"{packmol_info['x_low_a']:.3f} "
                    f"{packmol_info['y_low_a']:.3f} "
                    f"{z_low_a:.3f}  "
                    f"{packmol_info['x_high_a']:.3f} "
                    f"{packmol_info['y_high_a']:.3f} "
                    f"{z_high_a:.3f}"
                ),
                "end structure",
                "",
            ])

    inp_path = out_dir / "packmol_input.inp"
    inp_path.write_text("\n".join(lines), encoding="utf-8")
    return inp_path


def count_pdb_atoms(pdb_path):
    """统计 PDB 中 ATOM/HETATM 记录数。"""
    count = 0
    with open(pdb_path, "r", encoding="utf-8", errors="replace") as f:
        for line in f:
            record = line[:6].strip().upper()
            if record in ("ATOM", "HETATM"):
                count += 1

    if count <= 0:
        fail(f"PDB 中没有读取到任何 ATOM/HETATM：{pdb_path}")

    return count


def calculate_expected_packmol_atoms(parent_pdb, added, molecule_pdbs):
    """
    计算本轮 Packmol 输出理论原子总数：

        parent_system.pdb 原子数
        + Σ(新增分子数 × 对应 molecule PDB 原子数)

    这里直接使用实际参与 Packmol 的 PDB 模板计数，避免 ITP/PDB 原子数
    定义不一致导致误判。
    """
    parent_atoms = count_pdb_atoms(parent_pdb)
    expected = parent_atoms
    details = {
        "parent_atoms": parent_atoms,
        "molecules": {},
    }

    for mol in MOLECULE_ORDER:
        count = int(added.get(mol, 0))
        if count <= 0:
            continue

        pdb_path = molecule_pdbs[mol]
        atoms_per_molecule = count_pdb_atoms(pdb_path)
        added_atoms = count * atoms_per_molecule
        expected += added_atoms

        details["molecules"][mol] = {
            "molecule_count": count,
            "atoms_per_molecule": atoms_per_molecule,
            "added_atoms": added_atoms,
        }

    details["expected_total_atoms"] = expected
    return expected, details


def run_packmol(out_dir, inp_path, expected_atoms):
    """
    运行 Packmol，并对“达到循环上限但已经输出完整结构”的情况做容错。

    关键原则：
    1. Packmol return code == 0：仍然必须检查输出原子数；
    2. Packmol return code != 0：不立即退出；
    3. 若 merged_system.pdb 不存在，则真正失败；
    4. 若文件存在，统计 ATOM/HETATM 数：
           actual == expected -> 输出完整，允许继续；
           actual != expected -> 输出不完整，停止。

    因而不会因为 Packmol 仅仅达到内部 loop limit 就丢弃一个已经完整生成的结构，
    同时也不会错误接受缺少部分分子的半成品。
    """
    out_dir = Path(out_dir).resolve()
    inp_path = Path(inp_path).resolve()
    log_path = out_dir / "packmol.log"
    merged = out_dir / "merged_system.pdb"

    stdin_handle = None
    stdout_handle = None

    try:
        stdin_handle = open(inp_path, "rb")
        stdout_handle = open(log_path, "wb")

        try:
            result = subprocess.run(
                ["packmol"],
                cwd=str(out_dir),
                stdin=stdin_handle,
                stdout=stdout_handle,
                stderr=subprocess.STDOUT,
                check=False,
                timeout=PACKMOL_TIMEOUT,
            )
        except FileNotFoundError:
            fail("找不到外部命令：packmol。请确认已经安装并加入 PATH。")
        except subprocess.TimeoutExpired:
            fail(f"Packmol 运行超时，请检查日志：{log_path}")

    finally:
        if stdin_handle:
            stdin_handle.close()
        if stdout_handle:
            stdout_handle.close()

    return_code = int(result.returncode)

    if not merged.is_file():
        fail(
            f"Packmol 返回码为 {return_code}，且未生成 {merged}。\n"
            f"请检查日志：{log_path}"
        )

    actual_atoms = count_pdb_atoms(merged)

    print("\n[ Packmol output validation ]")
    print(f"  Packmol return code : {return_code}")
    print(f"  Expected atoms      : {expected_atoms}")
    print(f"  Actual PDB atoms    : {actual_atoms}")

    if actual_atoms != expected_atoms:
        fail(
            "Packmol 虽然生成了 merged_system.pdb，但输出结构不完整。\n"
            f"  Expected atoms : {expected_atoms}\n"
            f"  Actual atoms   : {actual_atoms}\n"
            f"  Difference     : {actual_atoms - expected_atoms}\n"
            f"请检查 Packmol 日志：{log_path}\n"
            "为避免坐标与 topol.top 分子数不匹配，程序不会继续。"
        )

    if return_code != 0:
        warn(
            "Packmol 返回非零退出码（例如达到 loop limit），"
            "但 merged_system.pdb 的实际原子数与理论原子数完全一致。\n"
            "因此判定结构完整，继续执行 GRO/topology/summary 后处理。\n"
            f"Packmol 日志：{log_path}"
        )
    else:
        print("  Packmol 状态       : 正常结束且结构完整")

    return merged


def convert_merged_to_gro_and_fix_box(out_dir, summary, output_name=DEFAULT_OUTPUT_GRO):
    out_dir = Path(out_dir)
    merged_pdb = out_dir / "merged_system.pdb"
    merged_gro = out_dir / "merged_system.gro"

    run_command([
        "gmx", "editconf",
        "-f", str(merged_pdb),
        "-o", str(merged_gro),
    ])

    if not merged_gro.is_file():
        fail(f"未生成中间 GRO：{merged_gro}")

    box = summary["box_dimensions_nm"]
    try:
        x = float(box["x"])
        y = float(box["y"])
        z = float(box["z"])
    except Exception:
        fail("summary 中 box_dimensions_nm 必须含合法 x/y/z。")

    lines = merged_gro.read_text(encoding="utf-8", errors="replace").splitlines(True)
    if not lines:
        fail(f"GRO 文件为空：{merged_gro}")

    lines[-1] = f"{x:10.5f}{y:10.5f}{z:10.5f}\n"

    final_gro = out_dir / output_name
    final_gro.write_text("".join(lines), encoding="utf-8")
    return final_gro


_QUOTED_INCLUDE_RE = re.compile(r'^(\s*#include\s+")([^"]+)(".*)$')


def resolve_topology_includes(parent_topol, child_dir):
    """把父 topology quoted include 重写为相对于新目录的路径。"""
    parent_topol = Path(parent_topol).resolve()
    parent_dir = parent_topol.parent
    child_dir = Path(child_dir).resolve()

    lines = parent_topol.read_text(
        encoding="utf-8", errors="replace"
    ).splitlines(True)

    rewritten = []

    for line in lines:
        raw = line.rstrip("\n")
        m = _QUOTED_INCLUDE_RE.match(raw)

        if not m:
            rewritten.append(line)
            continue

        old_path = m.group(2)
        p = Path(old_path)
        absolute = p if p.is_absolute() else (parent_dir / p).resolve()

        if not absolute.exists():
            warn(
                f"父 topology include 目标当前不存在：{old_path}\n"
                f"解析后：{absolute}"
            )

        new_rel = norm_posix(os.path.relpath(absolute, child_dir))
        newline = "\n" if line.endswith("\n") else ""
        rewritten.append(f'{m.group(1)}{new_rel}{m.group(3)}{newline}')

    return rewritten


def append_added_molecules_to_topology(parent_topol, child_dir, added):
    """
    保留父 [ molecules ] 顺序，并按 Packmol 新增结构顺序追加本轮新增条目。
    不把父数量和新增数量聚合，以避免破坏坐标-拓扑的分子顺序对应。
    """
    lines = resolve_topology_includes(parent_topol, child_dir)

    mol_start = None
    next_section = None

    for i, line in enumerate(lines):
        if line.strip().lower().startswith("[ molecules ]"):
            mol_start = i
            break

    if mol_start is None:
        fail(f"父 topology 中找不到 [ molecules ]：{parent_topol}")

    for i in range(mol_start + 1, len(lines)):
        stripped = lines[i].strip()
        if stripped.startswith("[") and stripped.endswith("]"):
            next_section = i
            break

    insert_at = next_section if next_section is not None else len(lines)

    additions = [
        f"{mol:8s}  {added[mol]}\n"
        for mol in MOLECULE_ORDER
        if added[mol] > 0
    ]

    block = [
        "\n",
        "; incremental molecules appended for this coordinate file\n",
        *additions,
    ]

    new_lines = lines[:insert_at] + block + lines[insert_at:]

    out_top = Path(child_dir) / "topol.top"
    out_top.write_text("".join(new_lines), encoding="utf-8")
    return out_top


def parse_itp_molar_mass(itp_path):
    """从 ITP [ atoms ] 段累加最后一列质量。"""
    lines = Path(itp_path).read_text(
        encoding="utf-8", errors="replace"
    ).splitlines()

    in_atoms = False
    atom_count = 0
    mass = 0.0

    for line in lines:
        stripped = line.strip()
        if not stripped or stripped.startswith(";"):
            continue

        if stripped.lower().startswith("[ atoms"):
            in_atoms = True
            continue

        if in_atoms and stripped.startswith("["):
            break

        if not in_atoms:
            continue

        data = stripped.split(";", 1)[0].split()
        if len(data) < 8:
            continue

        try:
            atom_mass = float(data[-1])
        except ValueError:
            continue

        atom_count += 1
        mass += atom_mass

    if atom_count == 0:
        fail(f"ITP 的 [ atoms ] 中没有读取到带质量原子：{itp_path}")

    return {
        "atom_count_per_molecule": atom_count,
        "molar_mass_g_mol": mass,
    }


def parse_topology_include_targets(topol_path):
    topol_path = Path(topol_path).resolve()
    result = []

    for line in topol_path.read_text(
        encoding="utf-8", errors="replace"
    ).splitlines():
        m = _QUOTED_INCLUDE_RE.match(line)
        if not m:
            continue

        p = Path(m.group(2))
        if not p.is_absolute():
            p = (topol_path.parent / p).resolve()
        result.append(p)

    return result


def find_molecule_itp(parent_topol, mol):
    targets = parse_topology_include_targets(parent_topol)

    for p in targets:
        if p.name.upper() == f"{mol}.ITP" and p.is_file():
            return p

    cur = Path(parent_topol).parent.resolve()
    for _ in range(8):
        candidate = cur / "basic" / "FF" / f"{mol}.itp"
        if candidate.is_file():
            return candidate.resolve()
        if cur.parent == cur:
            break
        cur = cur.parent

    return None


def collect_molecule_properties(parent_topol, parent_summary, target):
    old_stats = parent_summary.get("electrolyte_statistics", {})
    old_species = old_stats.get("species", {}) if isinstance(old_stats, dict) else {}

    props = {}

    for mol in MOLECULE_ORDER:
        if target[mol] <= 0:
            continue

        itp = find_molecule_itp(parent_topol, mol)

        if itp is not None:
            info = parse_itp_molar_mass(itp)
            info["itp_file"] = str(itp)
            props[mol] = info
            continue

        fallback = old_species.get(mol, {}) if isinstance(old_species, dict) else {}
        mass = fallback.get("molar_mass_g_mol")
        atoms = fallback.get("atom_count_per_molecule")

        if mass is None:
            fail(
                f"无法找到 {mol}.itp，且父 summary 中也没有可复用的 {mol} 摩尔质量。"
            )

        warn(f"未定位到 {mol}.itp，沿用父 summary 中的 {mol} 物性。")
        props[mol] = {
            "atom_count_per_molecule": atoms,
            "molar_mass_g_mol": float(mass),
            "itp_file": None,
        }

    return props


def gcd_many(values):
    values = [abs(int(v)) for v in values if int(v) != 0]
    if not values:
        return 1

    g = values[0]
    for value in values[1:]:
        g = math.gcd(g, value)
    return g


def build_electrolyte_statistics(target, properties):
    active = [mol for mol in MOLECULE_ORDER if target[mol] > 0]
    total_particles = sum(target[mol] for mol in active)

    if total_particles <= 0:
        return {
            "present_species": [],
            "number_of_species": 0,
            "total_electrolyte_particles": 0,
            "simplified_count_ratio": {},
            "species": {},
            "average_molar_mass_g_mol": 0.0,
            "average_molar_mass_definition": "sum(N_i*M_i)/sum(N_i)",
        }

    g = gcd_many([target[mol] for mol in active])

    weighted_mass_sum = sum(
        target[mol] * float(properties[mol]["molar_mass_g_mol"])
        for mol in active
    )

    species = {}

    for mol in active:
        count = target[mol]
        molar_mass = float(properties[mol]["molar_mass_g_mol"])
        mole_fraction = count / total_particles
        count_x_mass = count * molar_mass
        mass_fraction = count_x_mass / weighted_mass_sum if weighted_mass_sum > 0 else 0.0

        species[mol] = {
            "molecule_count": count,
            "atom_count_per_molecule": properties[mol].get("atom_count_per_molecule"),
            "molar_mass_g_mol": round(molar_mass, 6),
            "mole_fraction": round(mole_fraction, 10),
            "percentage": round(mole_fraction * 100.0, 6),
            "count_x_molar_mass": round(count_x_mass, 6),
            "mass_fraction": round(mass_fraction, 10),
        }

    average = weighted_mass_sum / total_particles

    return {
        "present_species": active,
        "number_of_species": len(active),
        "total_electrolyte_particles": total_particles,
        "simplified_count_ratio": {mol: target[mol] // g for mol in active},
        "species": species,
        "average_molar_mass_g_mol": round(average, 6),
        "average_molar_mass_definition": "sum(N_i*M_i)/sum(N_i)",
    }


_HISTORY_KEYS = {
    "generation_info",
    "molecule_changes",
    "added_molecules",
    "incremental_addition",
    "molecule_region_distribution",
}


def build_final_summary(parent_summary, out_dir, output_gro_name, target, electrolyte_stats):
    """
    只保存最终体系快照。
    主动移除可能存在的历史/增量字段，避免下一轮误用。
    """
    summary = copy.deepcopy(parent_summary)

    for key in list(summary.keys()):
        if key in _HISTORY_KEYS:
            summary.pop(key, None)

    summary["system_name"] = Path(out_dir).name
    summary["final_structure_file"] = output_gro_name
    summary["electrolyte_molecules"] = dict(target)
    summary["electrolyte_statistics"] = electrolyte_stats

    return summary


def write_summary_files(out_dir, summary):
    out_dir = Path(out_dir)

    json_path = out_dir / "system_summary.json"
    json_path.write_text(
        json.dumps(summary, indent=4, ensure_ascii=False),
        encoding="utf-8"
    )

    txt_path = out_dir / "system_summary.txt"
    lines = []
    lines.append("=" * 66)
    lines.append(f" SYSTEM SUMMARY REPORT: {summary.get('system_name', out_dir.name)}")
    lines.append("=" * 66)
    lines.append(f"MOF Name: {summary.get('mof_name', 'N/A')}")
    lines.append(f"Final Structure: {summary.get('final_structure_file', 'N/A')}")

    box = summary.get("box_dimensions_nm", {})
    if box:
        lines.append(
            "Box Size (nm): "
            f"X={float(box.get('x', 0.0)):.3f}, "
            f"Y={float(box.get('y', 0.0)):.3f}, "
            f"Z={float(box.get('z', 0.0)):.3f}"
        )

    regions = summary.get("z_structure_regions_nm", [])
    if regions:
        lines.append("")
        lines.append("[ Z-direction Five-Part Structure (nm) ]")
        for r in regions:
            lines.append(
                f"  - Region {r.get('region_id', '?')} "
                f"{r.get('name', ''):18s} "
                f"({r.get('type', ''):9s}) : "
                f"Z = [{float(r.get('z_low', 0.0)):.3f}, "
                f"{float(r.get('z_high', 0.0)):.3f}], "
                f"Thickness = {float(r.get('thickness', 0.0)):.3f}"
            )

    insertion = summary.get("packmol_insertion_regions_nm", [])
    if insertion:
        lines.append("")
        lines.append("[ Packmol Electrolyte Insertion Regions (nm) ]")
        for i, r in enumerate(insertion, start=1):
            rid = r.get("vacuum_region_id", i)
            lines.append(
                f"  - Vacuum Region {rid}: "
                f"Z = [{float(r['z_low']):.3f}, {float(r['z_high']):.3f}], "
                f"Usable Thickness = {float(r['thickness']):.3f}"
            )

    mof_info = summary.get("mof_electrode_statistics", {})
    if mof_info:
        lines.append("")
        lines.append("[ MOF Electrode Information ]")
        for label, key in [
            ("Total MOF Molecules in System", "total_mof_molecules"),
            ("Atoms per Single Molecule", "single_molecule_atoms"),
            ("Total MOF Atoms", "total_mof_atoms"),
            ("Total MOF Mass (g/mol or amu)", "total_mof_mass"),
        ]:
            if key in mof_info:
                lines.append(f"  - {label:34s}: {mof_info[key]}")

    composition = summary.get("electrolyte_molecules", {})
    lines.append("")
    lines.append("[ Electrolyte Composition - Final Total Counts ]")
    for mol in MOLECULE_ORDER:
        lines.append(f"  - {mol:6s}: {int(composition.get(mol, 0))}")

    estats = summary.get("electrolyte_statistics", {})
    if estats:
        lines.append("")
        lines.append("[ Electrolyte Statistics ]")
        lines.append("  Present species: " + ", ".join(estats.get("present_species", [])))
        lines.append(f"  Number of species       : {estats.get('number_of_species', 0)}")
        lines.append(f"  Total particles         : {estats.get('total_electrolyte_particles', 0)}")
        lines.append(
            f"  Average molar mass      : "
            f"{float(estats.get('average_molar_mass_g_mol', 0.0)):.6f} g/mol"
        )

        ratio = estats.get("simplified_count_ratio", {})
        if ratio:
            ratio_text = " : ".join(
                f"{mol}={ratio[mol]}"
                for mol in MOLECULE_ORDER
                if mol in ratio
            )
            lines.append(f"  Simplified count ratio  : {ratio_text}")

        species = estats.get("species", {})
        if species:
            lines.append("")
            lines.append(
                "  Species      Count     MolarMass(g/mol)    "
                "MoleFraction      Percent(%)      MassFraction"
            )
            for mol in MOLECULE_ORDER:
                if mol not in species:
                    continue
                s = species[mol]
                lines.append(
                    f"  {mol:8s} "
                    f"{int(s['molecule_count']):8d} "
                    f"{float(s['molar_mass_g_mol']):18.6f} "
                    f"{float(s['mole_fraction']):16.8f} "
                    f"{float(s['percentage']):15.6f} "
                    f"{float(s['mass_fraction']):16.8f}"
                )

    lines.append("=" * 66)
    txt_path.write_text("\n".join(lines) + "\n", encoding="utf-8")

    return json_path, txt_path


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "在任意已有 GRO 母体系基础上仅增加分子，"
            "并生成新的 system.gro/topol.top/summary。"
        )
    )

    parser.add_argument(
        "gro",
        help="母体系 GRO 文件；文件名不限，例如 system.gro、npt.gro、80.gro 等。"
    )
    parser.add_argument(
        "--ratio",
        required=True,
        help="目标最终分子总数量 JSON，例如 new_mix_ratio.json。"
    )
    parser.add_argument(
        "--summary",
        default=None,
        help="母体系 summary；默认同目录 system_summary.json，否则尝试唯一 *_summary.json。"
    )
    parser.add_argument(
        "--inp",
        default=None,
        help="母体系 Packmol inp；默认同目录 packmol_input.inp，否则尝试唯一 .inp。"
    )
    parser.add_argument(
        "--topol",
        default=None,
        help="母体系 topology；默认同目录 topol.top，否则尝试唯一 .top。"
    )
    parser.add_argument(
        "--output-gro-name",
        default=DEFAULT_OUTPUT_GRO,
        help=f"新体系最终 GRO 文件名，默认 {DEFAULT_OUTPUT_GRO}。"
    )
    parser.add_argument(
        "--output-root",
        default=None,
        help=(
            "新数字目录放置位置。默认：若母目录名为数字，则放到其父目录；"
            "否则放到母目录。"
        )
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="允许删除并重建已存在的目标 EMIM 数量目录。"
    )

    return parser.parse_args()


def main():
    args = parse_args()

    input_gro = Path(args.gro).resolve()
    if not input_gro.is_file():
        fail(f"输入 GRO 不存在：{input_gro}")

    parent_dir = input_gro.parent

    print("=" * 72)
    print("增量构建 MOF-电解液体系")
    print("=" * 72)
    print(f"母体系 GRO      : {input_gro}")
    print(f"母体系目录      : {parent_dir}")

    summary_path = choose_companion_file(
        parent_dir, args.summary, DEFAULT_SUMMARY, "_summary.json", "summary JSON"
    )
    inp_path = choose_companion_file(
        parent_dir, args.inp, DEFAULT_INP, ".inp", "Packmol inp"
    )
    topol_path = choose_companion_file(
        parent_dir, args.topol, DEFAULT_TOPOL, ".top", "topology"
    )

    print(f"Summary         : {summary_path}")
    print(f"Parent inp      : {inp_path}")
    print(f"Parent topology : {topol_path}")

    parent_summary, current = load_parent_summary(summary_path)
    target = load_target_ratio(args.ratio)
    added = compute_additions(current, target)

    print("\n[ 当前 -> 目标 -> 本轮需要增加 ]")
    for mol in MOLECULE_ORDER:
        print(
            f"  {mol:6s}: "
            f"{current[mol]:8d} -> "
            f"{target[mol]:8d} -> "
            f"+{added[mol]:8d}"
        )

    packmol_info = parse_packmol_input(inp_path)
    assignment = allocate_additions(added, parent_summary)

    print("\n[ 本轮新增分子的三区分配 ]")
    for mol in MOLECULE_ORDER:
        if added[mol] > 0:
            a = assignment[mol]
            print(
                f"  {mol:6s}: lower={a[0]}, middle={a[1]}, upper={a[2]}, total={sum(a)}"
            )

    if args.output_root:
        series_root = Path(args.output_root).resolve()
    else:
        series_root = determine_series_root(parent_dir)

    series_root.mkdir(parents=True, exist_ok=True)
    out_dir = prepare_output_dir(series_root, target["EMIM"], overwrite=args.overwrite)
    print(f"\n新体系目录      : {out_dir}")

    molecule_pdbs = {}
    for mol in MOLECULE_ORDER:
        if added[mol] > 0:
            molecule_pdbs[mol] = find_molecule_pdb(parent_dir, mol, packmol_info)
            print(f"{mol:6s} PDB      : {molecule_pdbs[mol]}")

    # 1. 任意名称母 GRO -> parent_system.pdb
    parent_pdb = convert_parent_gro_to_pdb(input_gro, out_dir)

    # Packmol 完整性校验所需的理论总原子数。
    # 直接根据实际参与 Packmol 的 PDB 文件计算。
    expected_atoms, atom_validation = calculate_expected_packmol_atoms(
        parent_pdb=parent_pdb,
        added=added,
        molecule_pdbs=molecule_pdbs,
    )

    print("\n[ Packmol expected atom count ]")
    print(f"  Parent atoms       : {atom_validation['parent_atoms']}")
    for mol, info in atom_validation["molecules"].items():
        print(
            f"  {mol:6s}: "
            f"{info['molecule_count']} molecules × "
            f"{info['atoms_per_molecule']} atoms = "
            f"{info['added_atoms']} atoms"
        )
    print(f"  Expected total     : {expected_atoms}")

    # 2. 生成只包含 delta molecules 的 Packmol 输入
    new_inp = write_incremental_packmol_input(
        out_dir=out_dir,
        parent_pdb=parent_pdb,
        summary=parent_summary,
        packmol_info=packmol_info,
        assignment=assignment,
        molecule_pdbs=molecule_pdbs,
    )

    # 3. 固定父体系，仅增加新分子
    run_packmol(out_dir, new_inp, expected_atoms=expected_atoms)

    # 4. 输出新 GRO 并恢复父体系 box
    final_gro = convert_merged_to_gro_and_fix_box(
        out_dir,
        parent_summary,
        output_name=args.output_gro_name,
    )

    # 5. topology：父分子顺序 + 本轮追加顺序
    new_topol = append_added_molecules_to_topology(
        parent_topol=topol_path,
        child_dir=out_dir,
        added=added,
    )

    # 6. 重新按最终总组成统计物性
    properties = collect_molecule_properties(
        parent_topol=topol_path,
        parent_summary=parent_summary,
        target=target,
    )
    electrolyte_stats = build_electrolyte_statistics(target, properties)

    # 7. summary 只保留最终快照，不保留影响下一轮的增量历史
    final_summary = build_final_summary(
        parent_summary=parent_summary,
        out_dir=out_dir,
        output_gro_name=Path(final_gro).name,
        target=target,
        electrolyte_stats=electrolyte_stats,
    )

    summary_json, summary_txt = write_summary_files(out_dir, final_summary)

    print("\n" + "=" * 72)
    print("新体系构建完成")
    print("=" * 72)
    print(f"Final GRO       : {final_gro}")
    print(f"Topology        : {new_topol}")
    print(f"Summary JSON    : {summary_json}")
    print(f"Summary TXT     : {summary_txt}")
    print(f"Packmol input   : {new_inp}")
    print(f"Packmol log     : {out_dir / 'packmol.log'}")
    print("")
    print(
        "新 summary 中 electrolyte_molecules 保存的是当前最终总数量；"
        "下一次可直接指定本次目录中的任意 GRO 文件作为新的母体系。"
    )


if __name__ == "__main__":
    main()
