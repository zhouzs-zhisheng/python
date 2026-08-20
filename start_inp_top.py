#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
build_electrolyte_system.py (Dynamic Relative Path & Symlink Friendly Version)

功能：
1. 智能动态探测 basic 文件夹的相对路径（支持软链接和多级跳转，摒弃绝对路径）；
2. 智能兼容 mix_ratio.json（支持双组参数）；
3. 严格使用 PDB 格式输入给 Packmol，输出 PDB 并通过 gmx editconf 转为 GRO；
4. 将盒子沿 Z 方向严格划分为下真空/下电极/中真空/上电极/上真空五部分，Packmol 仅向三个真空区插入电解液；
5. 精准修正盒子的最后一行尺寸（nm），最终结构命名为 system.gro；
6. 按规范格式生成 topol.top：每个 #include 独立成行，并重建 [ system ]/[ molecules ]；
7. 解析 MOF 的 itp 文件，计算总原子数与总质量，并导出结构化 JSON 及 TXT 统计报告。
"""

import os
import sys
import json
import subprocess
from math import floor, gcd

# ========== 用户可调参数（内部换算为 Å 传入 Packmol） ==========
MARGIN_Z_NM = 0.2          # nm，Z方向缓冲，防止分子与电极原子重叠
MARGIN_XY_NM = 0.2         # nm，XY方向缓冲，防止分子贴壁
TOLERANCE_NM = 0.2         # nm，Packmol 最小原子间距
PACKMOL_TIMEOUT = 3600     # 秒，Packmol 运行超时时间
MOLECULE_TYPES = ("EMIM", "BF4", "ACN", "PC")
# ==========================================================

def get_mof_name():
    """返回当前目录名，作为 MOF 名称"""
    return os.path.basename(os.getcwd())

def find_basic_directory():
    """
    智能动态查找 basic 文件夹（支持软链接和相对路径），
    避免使用绝对路径。依次检查当前目录、上级目录、上上级目录。
    """
    possible_paths = [
        "basic",
        os.path.join("..", "basic"),
        os.path.join("..", "..", "basic"),
        os.path.join("..", "..", "..", "basic")
    ]
    
    for path in possible_paths:
        if os.path.isdir(path):
            # 返回规范化的相对路径
            rel_path = os.path.normpath(path)
            print(f"成功找到 basic 目录路径: {rel_path}")
            return rel_path
            
    print("错误：在当前目录及其父目录中均找不到 basic 文件夹（或软链接），请检查路径！")
    sys.exit(1)

def get_basic_to_target_prefix(basic_dir):
    """
    根据 basic_dir 的相对路径层级，动态计算在生成的 topol.top 中
    访问 basic 目录所需的相对路径前缀（例如 "../" 或 "../../../"）
    """
    # 计算 basic_dir 包含多少个 ".."
    depth = basic_dir.count("..")
    if depth == 0:
        return "basic"
    elif depth == 1:
        return "../basic"
    elif depth == 2:
        return "../../basic"
    else:
        return "../../../basic"

def get_basic_path_from_workdir(basic_dir):
    """
    计算从 ACN/PC 工作子目录访问 basic 目录的相对路径。

    basic_dir 是以“脚本启动时的 MOF 当前目录”为基准得到的相对路径；
    Packmol 和生成后的 topol.top 位于其下一层 work_dir 中，因此需要
    先从 work_dir 返回 MOF 当前目录一级，再沿 basic_dir 访问 basic。

    例如：
        basic_dir = "basic"       -> "../basic"
        basic_dir = "../basic"    -> "../../basic"
        basic_dir = "../../basic" -> "../../../basic"
    """
    return os.path.normpath(os.path.join("..", basic_dir)).replace("\\", "/")


def check_prerequisites(basic_dir, mof):
    """检查所有必需的模板文件是否存在"""
    required_files = [
        os.path.join(basic_dir, "electrode_lhr", f"{mof}.gro"),
        os.path.join(basic_dir, "topol_lhr", f"{mof}.top"),
        os.path.join(basic_dir, "itp_lhr", f"{mof}.itp"),
        os.path.join(basic_dir, "FF", "ffatomtype_IL_ACN.itp"),
        os.path.join(basic_dir, "FF", "ffatomtype_MOF.itp"),
    ]
    for mol in MOLECULE_TYPES:
        required_files.append(os.path.join(basic_dir, "str", f"{mol}.pdb"))
        required_files.append(os.path.join(basic_dir, "FF", f"{mol}.itp"))

    missing = [f for f in required_files if not os.path.isfile(f)]
    if missing:
        print("错误：以下必需文件不存在：")
        for f in missing:
            print(f"  {f}")
        sys.exit(1)
    print("所有模板文件检查通过。")

def convert_gro_to_pdb(basic_dir, mof, work_dir):
    """使用 gmx editconf 将电极 gro 转换为 pdb，并放入指定的子文件夹中"""
    input_gro = os.path.join(basic_dir, "electrode_lhr", f"{mof}.gro")
    output_pdb = os.path.join(work_dir, "electrode_for_packmol.pdb")
    cmd = f"gmx editconf -f {input_gro} -o {output_pdb}"
    try:
        subprocess.run(cmd, shell=True, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        print(f"电极文件已转换为 {output_pdb}")
    except subprocess.CalledProcessError as e:
        print(f"gmx editconf 执行失败：{e.stderr.decode()}")
        sys.exit(1)
    return "electrode_for_packmol.pdb"

def parse_electrode_gro(basic_dir, mof):
    """
    解析电极 GRO 文件，提取盒子尺寸和上下两块 MOF 电极的完整 Z 边界。

    返回：
        x_len, y_len, z_len,
        z_bottom_min, z_bottom_max,
        z_top_min, z_top_max,
        mol_max

    所有坐标和盒长单位均为 nm。

    仍沿用原脚本的残基编号假设：
        residue <= mol_max 为下部电极；
        residue >  mol_max 为上部电极。
    """
    gro_path = os.path.join(basic_dir, "electrode_lhr", f"{mof}.gro")
    with open(gro_path, 'r') as f:
        lines = f.readlines()

    atom_lines = lines[2:-1]
    box_line = lines[-1].strip().split()
    x_len, y_len, z_len = map(float, box_line[:3])

    res_ids = []
    z_coords = []
    for line in atom_lines:
        res_no = int(line[0:5].strip())
        z = float(line[36:44].strip())
        res_ids.append(res_no)
        z_coords.append(z)

    if not res_ids:
        print(f"错误：电极 GRO 文件中没有读取到原子：{gro_path}")
        sys.exit(1)

    max_res = max(res_ids)
    mol_max = max_res // 2

    bottom_z = [z for rid, z in zip(res_ids, z_coords) if rid <= mol_max]
    top_z = [z for rid, z in zip(res_ids, z_coords) if rid > mol_max]

    if not bottom_z or not top_z:
        print("错误：未能根据残基编号正确分离底部和顶部电极。")
        sys.exit(1)

    z_bottom_min = min(bottom_z)
    z_bottom_max = max(bottom_z)
    z_top_min = min(top_z)
    z_top_max = max(top_z)

    if not (0.0 <= z_bottom_min <= z_bottom_max <= z_top_min <= z_top_max <= z_len):
        print("错误：解析得到的电极 Z 边界与盒子范围/上下电极顺序不一致。")
        print(f"  box Z       : [0.000, {z_len:.3f}] nm")
        print(f"  bottom MOF  : [{z_bottom_min:.3f}, {z_bottom_max:.3f}] nm")
        print(f"  top MOF     : [{z_top_min:.3f}, {z_top_max:.3f}] nm")
        sys.exit(1)

    print(f"盒子尺寸: x={x_len:.3f}, y={y_len:.3f}, z={z_len:.3f} nm")
    print(f"底部电极 Z 范围 = [{z_bottom_min:.3f}, {z_bottom_max:.3f}] nm")
    print(f"顶部电极 Z 范围 = [{z_top_min:.3f}, {z_top_max:.3f}] nm")
    print(f"MOF 单层分子数 = {mol_max}")

    return (x_len, y_len, z_len,
            z_bottom_min, z_bottom_max,
            z_top_min, z_top_max,
            mol_max)

def parse_mof_itp(basic_dir, mof, total_mof_molecules):
    """
    解析 MOF 的 itp 文件，提取每个原子的质量，计算总原子数与总质量。
    """
    itp_path = os.path.join(basic_dir, "itp_lhr", f"{mof}.itp")
    if not os.path.isfile(itp_path):
        print(f"错误：找不到 MOF itp 文件：{itp_path}")
        sys.exit(1)

    with open(itp_path, 'r') as f:
        lines = f.readlines()

    in_atoms_section = False
    single_molecule_atom_count = 0
    single_molecule_mass = 0.0

    for line in lines:
        stripped = line.strip()
        if not stripped or stripped.startswith(";"):
            continue
        if "[ atoms ]" in stripped:
            in_atoms_section = True
            continue
        if in_atoms_section:
            if stripped.startswith("["):
                break 
            parts = stripped.split()
            if len(parts) >= 8:
                try:
                    mass = float(parts[-1])
                    single_molecule_mass += mass
                    single_molecule_atom_count += 1
                except ValueError:
                    continue

    total_mof_atoms = single_molecule_atom_count * total_mof_molecules
    total_mof_mass = single_molecule_mass * total_mof_molecules

    mof_info = {
        "single_molecule_atoms": single_molecule_atom_count,
        "single_molecule_mass": round(single_molecule_mass, 4),
        "total_mof_molecules": total_mof_molecules,
        "total_mof_atoms": total_mof_atoms,
        "total_mof_mass": round(total_mof_mass, 4)
    }
    return mof_info


def parse_molecule_itp_mass(basic_dir, molecule_name):
    """
    从 basic/FF/<molecule_name>.itp 的 [ atoms ] 区域读取单个分子的摩尔质量。

    默认 GROMACS [ atoms ] 行最后一列为 mass：
        nr type resnr residue atom cgnr charge mass

    返回：
        {
            "molecule": 分子名称,
            "atom_count": 单分子原子数,
            "molar_mass_g_mol": 单分子摩尔质量
        }
    """
    itp_path = os.path.join(basic_dir, "FF", f"{molecule_name}.itp")

    if not os.path.isfile(itp_path):
        print(f"错误：找不到 {molecule_name} 的 ITP 文件：{itp_path}")
        sys.exit(1)

    with open(itp_path, "r") as f:
        lines = f.readlines()

    in_atoms_section = False
    atom_count = 0
    molar_mass = 0.0

    for line in lines:
        stripped = line.strip()

        if not stripped or stripped.startswith(";"):
            continue

        if stripped.lower() == "[ atoms ]":
            in_atoms_section = True
            continue

        if in_atoms_section:
            if stripped.startswith("["):
                break

            # 去掉行尾注释，避免注释内容干扰 split
            data_part = stripped.split(";", 1)[0].strip()
            if not data_part:
                continue

            parts = data_part.split()
            if len(parts) < 8:
                print(
                    f"错误：{itp_path} 的 [ atoms ] 中存在少于 8 列的原子行：\n"
                    f"  {line.rstrip()}"
                )
                sys.exit(1)

            try:
                mass = float(parts[7])
            except ValueError:
                print(
                    f"错误：无法解析 {itp_path} 的原子质量：\n"
                    f"  {line.rstrip()}"
                )
                sys.exit(1)

            atom_count += 1
            molar_mass += mass

    if not in_atoms_section:
        print(f"错误：{itp_path} 中找不到 [ atoms ] 区域。")
        sys.exit(1)

    if atom_count == 0:
        print(f"错误：{itp_path} 的 [ atoms ] 区域没有解析到有效原子。")
        sys.exit(1)

    return {
        "molecule": molecule_name,
        "atom_count": atom_count,
        "molar_mass_g_mol": round(molar_mass, 6),
    }


def compute_electrolyte_statistics(basic_dir, molecule_totals):
    """
    统计当前体系实际插入的电解液组分。

    仅 count > 0 的组分定义为“当前体系插入的分子种类”。

    对每种组分计算：
      - molecule_count：插入数量
      - molar_mass_g_mol：单分子摩尔质量
      - mole_fraction：按插入粒子数计算的摩尔分数
      - percentage：mole_fraction * 100
      - mass_weight：count * molar_mass，作为组成加权量

    电解液平均摩尔质量定义为粒子数/摩尔分数加权平均：
        M_avg = sum(N_i * M_i) / sum(N_i)

    注意：
    对 EMIM、BF4 这类离子，这里的“molecule”实际表示 Packmol/GROMACS
    中独立插入的物种粒子；平均摩尔质量因此是“所有插入物种粒子的平均摩尔质量”，
    不是把阳离子+阴离子预先合并成一个离子对后的摩尔质量。
    """
    inserted_species = [
        mol for mol in MOLECULE_TYPES
        if molecule_totals.get(mol, 0) > 0
    ]

    total_count = sum(molecule_totals[mol] for mol in inserted_species)

    if total_count <= 0:
        print("错误：当前体系没有任何待插入的电解液分子。")
        sys.exit(1)

    species_stats = {}
    weighted_mass_sum = 0.0

    # 计算最简整数计数比，例如 100:100:500 -> 1:1:5
    common_divisor = 0
    for mol in inserted_species:
        common_divisor = gcd(common_divisor, molecule_totals[mol])
    simplified_ratio = {
        mol: molecule_totals[mol] // common_divisor
        for mol in inserted_species
    }

    for mol in inserted_species:
        count = molecule_totals[mol]
        mass_info = parse_molecule_itp_mass(basic_dir, mol)
        molar_mass = mass_info["molar_mass_g_mol"]
        mole_fraction = count / total_count
        mass_weight = count * molar_mass

        weighted_mass_sum += mass_weight

        species_stats[mol] = {
            "molecule_count": count,
            "atom_count_per_molecule": mass_info["atom_count"],
            "molar_mass_g_mol": molar_mass,
            "mole_fraction": round(mole_fraction, 8),
            "percentage": round(mole_fraction * 100.0, 4),
            "count_x_molar_mass": round(mass_weight, 6),
        }

    average_molar_mass = weighted_mass_sum / total_count

    # 基于 count*M 的质量分数，便于同时查看“数量比例”和“质量贡献比例”
    for mol in inserted_species:
        species_stats[mol]["mass_fraction"] = round(
            species_stats[mol]["count_x_molar_mass"] / weighted_mass_sum,
            8
        )

    return {
        "inserted_species": inserted_species,
        "number_of_species": len(inserted_species),
        "total_inserted_particles": total_count,
        "simplified_count_ratio": simplified_ratio,
        "species": species_stats,
        "average_molar_mass_g_mol": round(average_molar_mass, 6),
        "average_molar_mass_definition": "sum(N_i*M_i)/sum(N_i)",
    }


def compute_system_regions(z_len, z_bottom_min, z_bottom_max, z_top_min, z_top_max):
    """
    根据盒子和两块电极的完整 Z 边界构造五个物理区域，并计算三个
    真空区域中真正允许 Packmol 插入分子的有效区域。

    五个物理区域（nm）：
        1. lower_vacuum     : 0              -> z_bottom_min
        2. bottom_electrode : z_bottom_min   -> z_bottom_max
        3. middle_vacuum    : z_bottom_max   -> z_top_min
        4. top_electrode    : z_top_min      -> z_top_max
        5. upper_vacuum     : z_top_max      -> z_len

    Packmol 只使用三个 vacuum 区域，并在真空区两侧各保留 MARGIN_Z_NM。
    返回：
        physical_regions_nm : 五个物理区域，单位 nm
        packmol_regions_a   : 三个可插入区域，单位 Å，格式 (low, high, thickness)
    """
    physical_regions_nm = [
        {"region_id": 1, "name": "lower_vacuum", "type": "vacuum",
         "z_low": 0.0, "z_high": z_bottom_min,
         "thickness": max(0.0, z_bottom_min)},
        {"region_id": 2, "name": "bottom_electrode", "type": "electrode",
         "z_low": z_bottom_min, "z_high": z_bottom_max,
         "thickness": max(0.0, z_bottom_max - z_bottom_min)},
        {"region_id": 3, "name": "middle_vacuum", "type": "vacuum",
         "z_low": z_bottom_max, "z_high": z_top_min,
         "thickness": max(0.0, z_top_min - z_bottom_max)},
        {"region_id": 4, "name": "top_electrode", "type": "electrode",
         "z_low": z_top_min, "z_high": z_top_max,
         "thickness": max(0.0, z_top_max - z_top_min)},
        {"region_id": 5, "name": "upper_vacuum", "type": "vacuum",
         "z_low": z_top_max, "z_high": z_len,
         "thickness": max(0.0, z_len - z_top_max)},
    ]

    vacuum_bounds_nm = [
        (0.0, z_bottom_min),
        (z_bottom_max, z_top_min),
        (z_top_max, z_len),
    ]

    packmol_regions_a = []
    for region_index, (physical_low, physical_high) in enumerate(vacuum_bounds_nm, start=1):
        low_nm = physical_low + MARGIN_Z_NM
        high_nm = physical_high - MARGIN_Z_NM
        thickness_nm = high_nm - low_nm

        if thickness_nm <= 0.0:
            print(
                f"错误：第 {region_index} 个真空区域扣除两侧 Z 缓冲 {MARGIN_Z_NM:.3f} nm 后没有可用空间。\n"
                f"  物理范围 = [{physical_low:.3f}, {physical_high:.3f}] nm\n"
                f"  可用范围 = [{low_nm:.3f}, {high_nm:.3f}] nm"
            )
            sys.exit(1)

        packmol_regions_a.append((low_nm * 10.0, high_nm * 10.0, thickness_nm * 10.0))

    return physical_regions_nm, packmol_regions_a

def load_molecule_ratios(json_path="mix_ratio.json"):
    """智能读取分子比例文件，兼容标准 JSON 列表或多行对象格式"""
    if not os.path.isfile(json_path):
        print(f"错误：找不到分子比例文件 {json_path}")
        sys.exit(1)
        
    with open(json_path, 'r') as f:
        content = f.read().strip()
        
    try:
        data = json.loads(content)
    except json.JSONDecodeError:
        try:
            fixed_content = "[" + content.replace("}{", "},{") + "]"
            data = json.loads(fixed_content)
        except Exception as e:
            print(f"错误：mix_ratio.json 格式不正确，解析失败。详细错误：{e}")
            sys.exit(1)
            
    if not isinstance(data, list) or len(data) < 2:
        print("错误：mix_ratio.json 必须包含至少两组参数配置")
        sys.exit(1)
        
    required_keys = {"EMIM", "BF4", "ACN", "PC"}
    for idx, item in enumerate(data[:2]):
        if not required_keys.issubset(item.keys()):
            print(f"错误：mix_ratio.json 中第 {idx+1} 行必须包含 EMIM, BF4, ACN, PC 四个键")
            sys.exit(1)
        for key in required_keys:
            if not isinstance(item[key], int) or item[key] < 0:
                print(f"错误：第 {idx+1} 行的 {key} 数量必须为非负整数")
                sys.exit(1)
                
    return data[0], data[1]

def assign_molecules_by_volume(molecule_totals, region_thicknesses):
    """按区域厚度比例分配每种分子到三个区域"""
    total_thick = sum(region_thicknesses)
    if total_thick <= 0:
        print("错误：所有区域厚度均为零，无法放置分子")
        sys.exit(1)

    assignment = {}
    for mol, total in molecule_totals.items():
        if total == 0:
            assignment[mol] = [0, 0, 0]
            continue
        raw_counts = [total * thick / total_thick for thick in region_thicknesses]
        int_counts = [floor(c) for c in raw_counts]
        frac_parts = [c - floor(c) for c in raw_counts]
        remainder = total - sum(int_counts)
        if remainder > 0:
            indices = sorted(range(3), key=lambda i: frac_parts[i], reverse=True)
            for i in range(remainder):
                int_counts[indices[i % 3]] += 1
        assignment[mol] = int_counts
    return assignment

def write_packmol_input(basic_dir, mof, x_len, y_len, regions, assignment, work_dir):
    """生成 packmol_input.inp 文件，分子pdb路径采用相对 basic_dir 的纯相对路径"""
    lines = []
    lines.append(f"tolerance {TOLERANCE_NM * 10.0}")
    lines.append("output merged_system.pdb")
    lines.append("filetype pdb")
    lines.append("")
    
    lines.append("structure electrode_for_packmol.pdb")
    lines.append("  number 1")
    lines.append("  fixed 0.0 0.0 0.0 0.0 0.0 0.0")
    lines.append("end structure")
    lines.append("")

    x_len_a = x_len * 10.0
    y_len_a = y_len * 10.0
    m_xy = MARGIN_XY_NM * 10.0

    # Packmol 在 work_dir 中执行，因此 PDB 路径必须以 work_dir 为基准。
    # basic_dir 本身则是以脚本启动时的 MOF 当前目录为基准。
    # work_dir 比 MOF 当前目录深一级，所以统一先退一级，再拼接 basic_dir。
    basic_from_workdir = get_basic_path_from_workdir(basic_dir)

    for mol in MOLECULE_TYPES:
        counts = assignment[mol]
        for region_idx, count in enumerate(counts):
            if count == 0:
                continue
            low, high, thick = regions[region_idx]
            if low >= high:
                continue
            
            pdb_path = os.path.join(basic_from_workdir, "str", f"{mol}.pdb")
            pdb_path = os.path.normpath(pdb_path).replace("\\", "/")

            # 写入 Packmol 输入前，再从 work_dir 视角验证一次实际路径。
            resolved_pdb = os.path.normpath(os.path.join(work_dir, pdb_path))
            if not os.path.isfile(resolved_pdb):
                print(
                    f"错误：Packmol 分子模板相对路径无效：{pdb_path}\n"
                    f"  从工作目录 {work_dir} 解析后为：{resolved_pdb}"
                )
                sys.exit(1)

            lines.append(f"structure {pdb_path}")
            lines.append(f"  number {count}")
            x_low = m_xy
            x_high = x_len_a - m_xy
            y_low = m_xy
            y_high = y_len_a - m_xy
            lines.append(f"  inside box {x_low:.3f} {y_low:.3f} {low:.3f}  {x_high:.3f} {y_high:.3f} {high:.3f}")
            lines.append("end structure")
            lines.append("")

    inp_path = os.path.join(work_dir, "packmol_input.inp")
    with open(inp_path, 'w') as f:
        f.write("\n".join(lines))

def run_packmol(work_dir):
    """在对应的子文件夹内执行 Packmol"""
    cwd = os.getcwd()
    os.chdir(work_dir)
    cmd = "packmol < packmol_input.inp > packmol.log 2>&1"
    try:
        subprocess.run(cmd, shell=True, check=True, timeout=PACKMOL_TIMEOUT)
    except Exception as e:
        print(f"Packmol 运行出错于 {work_dir}: {e}")
        os.chdir(cwd)
        sys.exit(1)
    os.chdir(cwd)

def convert_merged_pdb_to_gro(work_dir):
    """使用 gmx editconf 将 Packmol 输出的 merged_system.pdb 转换为 merged_system.gro"""
    cwd = os.getcwd()
    os.chdir(work_dir)
    cmd = "gmx editconf -f merged_system.pdb -o merged_system.gro"
    try:
        subprocess.run(cmd, shell=True, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    except Exception as e:
        print(f"gmx editconf 转换失败于 {work_dir}: {e}")
        os.chdir(cwd)
        sys.exit(1)
    os.chdir(cwd)

def fix_gro_box_and_write_system(work_dir, x_len, y_len, z_len):
    """
    读取 gmx editconf 生成的 merged_system.gro，使用原电极 GRO 的准确盒子尺寸
    替换最后一行，并将修正后的最终结构写为 system.gro。

    merged_system.gro 保留为中间文件；system.gro 是后续模拟应使用的最终结构。
    """
    input_gro = os.path.join(work_dir, "merged_system.gro")
    output_gro = os.path.join(work_dir, "system.gro")

    with open(input_gro, 'r') as f:
        lines = f.readlines()

    if not lines:
        print(f"错误：GRO 文件为空：{input_gro}")
        sys.exit(1)

    lines[-1] = f"{x_len:10.5f}{y_len:10.5f}{z_len:10.5f}\n"

    with open(output_gro, 'w') as f:
        f.writelines(lines)

    print(f"已修正盒子尺寸并生成最终结构文件：{output_gro}")

def create_system_top(basic_dir, mof, molecule_totals, work_dir):
    """
    生成格式规范的 GROMACS topol.top。

    输出格式严格重建为：
        1. 每个 #include 独占一行；
        2. 空行；
        3. [ system ]；
        4. MOF 名称；
        5. 空行；
        6. [ molecules ] 及其内容。

    这样生成出的 topol.top 已经包含 fix_top.py 的核心修正规则，
    不需要再额外执行 fix_top.py。

    同时：
        - #include 路径统一以 work_dir 为基准；
        - 保留原始 MOF topology 中 [ molecules ] 的非电解液条目；
        - 删除原有 EMIM/BF4/ACN/PC 条目；
        - 根据当前 molecule_totals 重新写入电解液数量。
    """
    top_orig = os.path.join(basic_dir, "topol_lhr", f"{mof}.top")

    if not os.path.isfile(top_orig):
        print(f"错误：找不到原始 MOF topology：{top_orig}")
        sys.exit(1)

    with open(top_orig, "r") as f:
        orig_lines = f.readlines()

    # ----------------------------------------------------------
    # 1. 找到原 topology 的 [ molecules ]。
    #    与 fix_top.py 一样，只保留 [ molecules ] 及其后真正需要的内容，
    #    [ system ] 之前的旧内容不直接复制，而是重新规范生成。
    # ----------------------------------------------------------
    mol_section_index = None

    for i, line in enumerate(orig_lines):
        if line.strip().lower().startswith("[ molecules ]"):
            mol_section_index = i
            break

    if mol_section_index is None:
        print(f"错误：原始 topology 中找不到 [ molecules ]：{top_orig}")
        sys.exit(1)

    molecules_rest = orig_lines[mol_section_index:]

    # ----------------------------------------------------------
    # 2. 从 work_dir 访问 basic 的统一相对路径。
    #    Packmol 与 topol.top 使用同一层级基准。
    # ----------------------------------------------------------
    prefix = get_basic_path_from_workdir(basic_dir)

    include_specs = [
        ("FF", "ffatomtype_IL_ACN.itp"),
        ("FF", "ffatomtype_MOF.itp"),
        ("itp_lhr", f"{mof}.itp"),
        ("FF", "EMIM.itp"),
        ("FF", "BF4.itp"),
        ("FF", "ACN.itp"),
        ("FF", "PC.itp"),
    ]

    include_lines = []
    include_targets = []

    for subdir, filename in include_specs:
        rel_target = os.path.normpath(
            os.path.join(prefix, subdir, filename)
        ).replace("\\", "/")

        include_targets.append(rel_target)
        include_lines.append(f'#include "{rel_target}"\n')

    # ----------------------------------------------------------
    # 3. 在真正写 topol.top 前逐个验证 include 路径。
    # ----------------------------------------------------------
    for rel_target in include_targets:
        resolved_target = os.path.normpath(
            os.path.join(work_dir, rel_target)
        )

        if not os.path.isfile(resolved_target):
            print(
                "错误：topol.top 的 include 相对路径无效："
                f"{rel_target}\n"
                f"  从工作目录 {work_dir} 解析后为：{resolved_target}"
            )
            sys.exit(1)

    # ----------------------------------------------------------
    # 4. 整理 [ molecules ]。
    #
    #    只在 [ molecules ] section 内删除旧的
    #    EMIM/BF4/ACN/PC 行，避免误删其他 section 中的文本。
    # ----------------------------------------------------------
    cleaned_molecules = []
    in_molecules = False
    next_section_index = None

    for line in molecules_rest:
        stripped = line.strip()

        if stripped.lower().startswith("[ molecules ]"):
            in_molecules = True
            cleaned_molecules.append("[ molecules ]\n")
            continue

        # 如果 [ molecules ] 后还有其他 section，则不再按 molecule 行处理。
        if in_molecules and stripped.startswith("[") and stripped.endswith("]"):
            in_molecules = False

        if in_molecules:
            # 空行/注释原样保留
            if not stripped or stripped.startswith(";"):
                cleaned_molecules.append(line)
                continue

            first_token = stripped.split()[0]

            if first_token in MOLECULE_TYPES:
                # 删除原体系中的旧电解液数量，稍后统一重写。
                continue

        cleaned_molecules.append(line)

    # ----------------------------------------------------------
    # 5. 确定新电解液条目应该插到哪里：
    #    放在 [ molecules ] 中所有原有非电解液条目之后，
    #    若后面还有其他 section，则插在下一个 section 之前。
    # ----------------------------------------------------------
    insert_pos = len(cleaned_molecules)

    for i in range(1, len(cleaned_molecules)):
        stripped = cleaned_molecules[i].strip()

        if stripped.startswith("[") and stripped.endswith("]"):
            insert_pos = i
            break

    electrolyte_lines = []

    for mol in MOLECULE_TYPES:
        total = molecule_totals.get(mol, 0)

        if total > 0:
            electrolyte_lines.append(
                f"{mol:8s}  {total}\n"
            )

    # 保证 molecule section 尾部格式干净。
    molecule_prefix = cleaned_molecules[:insert_pos]
    molecule_suffix = cleaned_molecules[insert_pos:]

    while (
        molecule_prefix
        and molecule_prefix[-1].strip() == ""
    ):
        molecule_prefix.pop()

    if electrolyte_lines:
        molecule_prefix.append("\n")
        molecule_prefix.extend(electrolyte_lines)

    if molecule_suffix:
        molecule_prefix.append("\n")

    final_molecules = molecule_prefix + molecule_suffix

    # ----------------------------------------------------------
    # 6. 按 fix_top.py 的最终目标格式重新构建整个 topology：
    #
    #    #include ...
    #    #include ...
    #
    #    [ system ]
    #    MOF_NAME
    #
    #    [ molecules ]
    #    ...
    # ----------------------------------------------------------
    new_lines = []

    new_lines.extend(include_lines)
    new_lines.append("\n")

    new_lines.append("[ system ]\n")
    new_lines.append(f"{mof}\n")
    new_lines.append("\n")

    new_lines.extend(final_molecules)

    out_top = os.path.join(work_dir, "topol.top")

    with open(out_top, "w") as f:
        f.writelines(new_lines)

    print(f"已生成规范化 topology：{out_top}")
    print(f"  [ system ] 名称：{mof}")
    print(f"  #include 数量：{len(include_lines)}")

def export_system_statistics(work_dir, mof, x_len, y_len, z_len,
                             physical_regions_nm, packmol_regions_a,
                             molecule_totals, mof_info, electrolyte_info):
    """导出系统统计信息的 JSON 和 TXT 报告，包括空间结构、MOF 与电解液组成统计。"""
    packmol_vacuum_regions_nm = [
        {
            "vacuum_region_id": i + 1,
            "z_low": r[0] / 10.0,
            "z_high": r[1] / 10.0,
            "thickness": r[2] / 10.0,
        }
        for i, r in enumerate(packmol_regions_a)
    ]

    electrode_regions = [r for r in physical_regions_nm if r["type"] == "electrode"]

    stats = {
        "system_name": work_dir,
        "mof_name": mof,
        "final_structure_file": "system.gro",
        "box_dimensions_nm": {
            "x": x_len,
            "y": y_len,
            "z": z_len
        },
        "z_structure_regions_nm": physical_regions_nm,
        "electrode_z_boundaries_nm": {
            "bottom_electrode": {
                "z_min": electrode_regions[0]["z_low"],
                "z_max": electrode_regions[0]["z_high"]
            },
            "top_electrode": {
                "z_min": electrode_regions[1]["z_low"],
                "z_max": electrode_regions[1]["z_high"]
            }
        },
        "packmol_insertion_regions_nm": packmol_vacuum_regions_nm,
        "packmol_z_margin_nm": MARGIN_Z_NM,
        "mof_electrode_statistics": mof_info,
        "electrolyte_molecules": molecule_totals,
        "electrolyte_statistics": electrolyte_info
    }

    json_path = os.path.join(work_dir, "system_summary.json")
    with open(json_path, 'w', encoding='utf-8') as f:
        json.dump(stats, f, indent=4, ensure_ascii=False)

    txt_path = os.path.join(work_dir, "system_summary.txt")
    with open(txt_path, 'w', encoding='utf-8') as f:
        f.write("==================================================\n")
        f.write(f" SYSTEM SUMMARY REPORT: {work_dir}\n")
        f.write("==================================================\n")
        f.write(f"MOF Name: {mof}\n")
        f.write(f"Final Structure: system.gro\n")
        f.write(f"Box Size (nm): X={x_len:.3f}, Y={y_len:.3f}, Z={z_len:.3f}\n\n")

        f.write("[ Z-direction Five-Part Structure (nm) ]\n")
        for r in physical_regions_nm:
            f.write(
                f"  - Region {r['region_id']} {r['name']:18s} "
                f"({r['type']:9s}) : Z = [{r['z_low']:.3f}, {r['z_high']:.3f}], "
                f"Thickness = {r['thickness']:.3f}\n"
            )

        f.write("\n[ Packmol Electrolyte Insertion Regions (nm) ]\n")
        f.write(f"  Z margin on each vacuum boundary: {MARGIN_Z_NM:.3f} nm\n")
        for r in packmol_vacuum_regions_nm:
            f.write(
                f"  - Vacuum Region {r['vacuum_region_id']} : "
                f"Z = [{r['z_low']:.3f}, {r['z_high']:.3f}], "
                f"Usable Thickness = {r['thickness']:.3f}\n"
            )

        f.write("\n[ MOF Electrode Information ]\n")
        f.write(f"  - Total MOF Molecules in System : {mof_info['total_mof_molecules']}\n")
        f.write(f"  - Atoms per Single Molecule     : {mof_info['single_molecule_atoms']}\n")
        f.write(f"  - Total MOF Atoms               : {mof_info['total_mof_atoms']}\n")
        f.write(f"  - Total MOF Mass (g/mol or amu) : {mof_info['total_mof_mass']:.4f}\n\n")

        f.write("[ Electrolyte Composition ]\n")
        f.write(
            "  Inserted species: "
            + ", ".join(electrolyte_info["inserted_species"])
            + "\n"
        )
        f.write(
            f"  Number of inserted species : "
            f"{electrolyte_info['number_of_species']}\n"
        )
        f.write(
            f"  Total inserted particles   : "
            f"{electrolyte_info['total_inserted_particles']}\n"
        )
        ratio_text = ":".join(
            str(electrolyte_info["simplified_count_ratio"][mol])
            for mol in electrolyte_info["inserted_species"]
        )
        ratio_names = ":".join(electrolyte_info["inserted_species"])
        f.write(
            f"  Simplified count ratio      : "
            f"{ratio_names} = {ratio_text}\n"
        )
        f.write(
            f"  Average molar mass          : "
            f"{electrolyte_info['average_molar_mass_g_mol']:.6f} g/mol\n"
        )
        f.write(
            "  Average definition          : "
            f"{electrolyte_info['average_molar_mass_definition']}\n\n"
        )

        f.write(
            "  Species     Count    MolarMass(g/mol)    "
            "MoleFraction    Percentage(%)\n"
        )
        for mol in electrolyte_info["inserted_species"]:
            info = electrolyte_info["species"][mol]
            f.write(
                f"  {mol:8s} "
                f"{info['molecule_count']:8d} "
                f"{info['molar_mass_g_mol']:18.6f} "
                f"{info['mole_fraction']:15.8f} "
                f"{info['percentage']:15.4f}\n"
            )

        f.write("\n[ Raw Molecule Counts ]\n")
        for mol in MOLECULE_TYPES:
            f.write(f"  - {mol:6s} : {molecule_totals.get(mol, 0)}\n")

        f.write("==================================================\n")

def process_system_pipeline(basic_dir, mof, molecule_totals, work_dir,
                            x_len, y_len, z_len,
                            z_bottom_min, z_bottom_max,
                            z_top_min, z_top_max,
                            mol_max):
    """针对单组参数的完整执行流水线。"""
    print(f"\n================ 开始处理子文件夹: {work_dir} ================")
    os.makedirs(work_dir, exist_ok=True)

    convert_gro_to_pdb(basic_dir, mof, work_dir)

    physical_regions_nm, packmol_regions_a = compute_system_regions(
        z_len,
        z_bottom_min, z_bottom_max,
        z_top_min, z_top_max
    )
    region_thicknesses = [r[2] for r in packmol_regions_a]

    assignment = assign_molecules_by_volume(molecule_totals, region_thicknesses)
    write_packmol_input(basic_dir, mof, x_len, y_len, packmol_regions_a, assignment, work_dir)
    run_packmol(work_dir)

    convert_merged_pdb_to_gro(work_dir)
    fix_gro_box_and_write_system(work_dir, x_len, y_len, z_len)

    create_system_top(basic_dir, mof, molecule_totals, work_dir)

    total_mof_molecules = mol_max * 2
    mof_info = parse_mof_itp(basic_dir, mof, total_mof_molecules)

    # 读取当前实际插入物种的 ITP，计算摩尔质量、摩尔分数和平均摩尔质量
    electrolyte_info = compute_electrolyte_statistics(
        basic_dir,
        molecule_totals
    )

    export_system_statistics(
        work_dir, mof, x_len, y_len, z_len,
        physical_regions_nm, packmol_regions_a,
        molecule_totals, mof_info, electrolyte_info
    )

    print(f"================ {work_dir} 目录处理完成 ================\n")

def main():
    print("===== 构建多体系电解液-MOF 混合系统（纯相对路径与软链接自适应版） =====")
    mof = get_mof_name()
    print(f"MOF 名称: {mof}")
    
    # 动态探测 basic 目录位置，完全抛弃绝对路径
    basic_dir = find_basic_directory()

    check_prerequisites(basic_dir, mof)
    (x_len, y_len, z_len,
     z_bottom_min, z_bottom_max,
     z_top_min, z_top_max,
     mol_max) = parse_electrode_gro(basic_dir, mof)

    acn_totals, pc_totals = load_molecule_ratios("mix_ratio.json")

    # 分别处理 ACN 和 PC 体系
    process_system_pipeline(
        basic_dir, mof, acn_totals, "ACN",
        x_len, y_len, z_len,
        z_bottom_min, z_bottom_max,
        z_top_min, z_top_max,
        mol_max
    )
    process_system_pipeline(
        basic_dir, mof, pc_totals, "PC",
        x_len, y_len, z_len,
        z_bottom_min, z_bottom_max,
        z_top_min, z_top_max,
        mol_max
    )

    print("所有体系构建步骤及信息统计全部完成！")

if __name__ == "__main__":
    main()
