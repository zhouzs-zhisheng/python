#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
analyze_capacitance.py

从 CPM 恒电势模拟结果计算 MOF 电极的质量电容（积分电容 + 微分电容）。

数据来源（每个电压点目录下，目录名形如 -2V / -1.5V / 0V / 0.5V / 1V / 2V，
 由脚本自动发现并按电压数值从小到大排序）：
  - CPM_electrodeCharge.dat : 每行 2 列 = (负极电荷, 正极电荷)，单位 |e|，无时间列。
                              负极电荷为负值是正常的，保留负号输出。
  - CPM_potential.dat       : 每行 3 列 = (Bulk, 负极侧电极电位, 正极侧电极电位)，
                              单位 V。Bulk 列是字符串 (如 not calc)，忽略。
  - system_summary.json     : mof_electrode_statistics.total_mof_mass 为电极质量
                              (Da，数值上等于 g/mol)，只算 MOF 电极。

时间换算（电荷/电势文件无时间列，靠头注释的输出频率换算）：
  头注释形如 "calculate frequency: 10, output frequency: 200"
  每行时间间隔 row_interval_fs = dt_fs * output_frequency
  电荷与电势两个文件各自读自己的 output_frequency（可能不同），各自取末窗口。
  dt_fs (计算间隔, fs) 由用户给出。

积分电容（每个非 0V 电压点）：
  C_int [F/g] = 4 * Q_abs_avg * F / (|V_rel| * M)
  Q_abs_avg   = (|Q_neg| + |Q_pos|) / 2   （正负电极等效电荷，末窗口均值，|e|）
  V_rel       = 该点实测电极间电位差 − 0V 实测电极间电位差 (V)，取绝对值
  M           = 电极质量 (total_mof_mass，Da 数值直接用)
  F           = 法拉第常数 96485 C/mol
  公式来源（用户定义）：C_int = 4·Q∞/(M·V)，Q∞ 为该电压下平衡电荷 (|e|)，
  V 取 |V_rel| 使负电压侧也得到正值；乘 4 为电极对几何换算因子。
  推导：q_C = Q_e*F/N_A，m_g = M/N_A  →  C_g = Q_e*F/(V*M)（N_A 抵消）

微分电容（电压点按数值排序后取相邻点对，正极/负极各自计算）：
  排序：目录名自动解析电压数值（-2V, -1.5V, -1V, -0.5V, 0V, 0.5V, 1V, 1.5V, 2V ...），
       从小到大排列，相邻两个电压点构成一个区间，每个电压点对应一份电荷量。
  ΔV          = V_rel,i+1 − V_rel,i
  ΔQ_pos      = Q_pos,i+1 − Q_pos,i
  ΔQ_neg      = |Q_neg,i+1| − |Q_neg,i|    （负极电荷为负，取幅值差分）
  C_diff,pos  = ΔQ_pos * F / (ΔV * M)
  C_diff,neg  = ΔQ_neg * F / (ΔV * M)
  注：对称扫描时负电压侧 |Q_neg| 随 |V| 减小而减小，ΔQ_neg/C_diff,neg 在负侧
      可能为负（负极是正极的镜像），属正常；如需恒为正，可改用电池总电荷
      ΔQ_cell = (Q_pos − Q_neg) 差分（两侧均为正）。
  同时输出电荷量随电压值的变化 (Q_neg/Q_pos/Q_abs_avg 随 V，见 CSV)。

输出（--output-dir 下，前缀 --prefix）：
  <prefix>.json : 结构化结果
  <prefix>.csv  : 一份 CSV，kind 列区分积分表(voltage)与微分表(pair)

用法：
  python analyze_capacitance.py --system-dir ACN --dt-fs 2.0
  python analyze_capacitance.py --system-dir ACN --dt-fs 2.0 --window-ns 5 \
        --output-dir ./results --prefix cap
"""

import argparse
import json
import math
import re
import sys
from pathlib import Path

# ============================================================
# 常量
# ============================================================

# 电压点目录命名规则：带符号十进制数 + 可选 'V' 后缀，如 -2V / -1.5V / 0V / 0.5V / 1V / 2V。
# 脚本自动发现并排序，不依赖固定列表（0 也可写作 "0"）。
VOLTAGE_DIR_RE = re.compile(r"^([+-]?\d+(?:\.\d+)?)\s*V?$", re.IGNORECASE)

CHARGE_FILE = "CPM_electrodeCharge.dat"
POTENTIAL_FILE = "CPM_potential.dat"
SYSTEM_SUMMARY = "system_summary.json"

FARADAY = 96485.0          # C/mol
FREQ_RE = re.compile(
    r"output\s*frequency\s*:\s*(\d+)", re.IGNORECASE
)


# ============================================================
# 通用工具
# ============================================================

def fail(msg):
    print(f"[FATAL] {msg}", file=sys.stderr)
    sys.exit(1)


def warn(msg):
    print(f"[WARN] {msg}")


def load_mof_mass(system_dir, override=None):
    """
    读取电极质量统计并输出明细 (便于自我检查)，返回 (total_mof_mass, stats)。
    只算 MOF 电极，总质量取自 mof_electrode_statistics.total_mof_mass (Da=g/mol)。
    交叉验证：single_molecule_mass × total_mof_molecules ≈ total_mof_mass。
    """
    if override is not None:
        stats = {
            "source": "命令行覆盖 --mof-mass",
            "total_mof_mass": override,
        }
        print(f"电极质量: 使用命令行覆盖值 {override} (Da)")
        return override, stats

    sf = Path(system_dir) / SYSTEM_SUMMARY
    if not sf.is_file():
        fail(f"{sf} 不存在，无法读取电极质量")
    with open(sf, "r") as f:
        summary = json.load(f)
    try:
        stats = summary["mof_electrode_statistics"]
        mass = float(stats["total_mof_mass"])
    except (KeyError, TypeError):
        fail(f"{sf} 中缺少 mof_electrode_statistics.total_mof_mass")

    # 交叉验证：单分子质量 × 分子数 是否等于 总质量
    single = stats.get("single_molecule_mass")
    nmol = stats.get("total_mof_molecules")
    if single is not None and nmol is not None:
        check = float(single) * float(nmol)
        if abs(check - mass) > 1e-4 * max(1.0, abs(mass)):
            warn(f"电极质量交叉验证不符：single_molecule_mass({single}) × "
                 f"total_mof_molecules({nmol}) = {check:.4f} "
                 f"≠ total_mof_mass({mass})")

    print("\n电极质量统计 (mof_electrode_statistics):")
    for k in ("single_molecule_atoms", "single_molecule_mass",
              "total_mof_molecules", "total_mof_atoms"):
        if k in stats:
            print(f"  {k:24s}: {stats[k]}")
    print(f"  {'total_mof_mass (用于电容计算)':24s}: {mass} Da (=g/mol)")
    return mass, stats


def read_output_frequency(path):
    """
    从头注释行提取 output frequency。
    注释区在数据行开始前结束；返回 int，找不到返回 None。
    """
    with open(path, "r") as f:
        for line in f:
            if line.startswith("#"):
                m = FREQ_RE.search(line)
                if m:
                    return int(m.group(1))
            else:
                break
    return None


def resolve_row_interval(path, dt_fs, sample_fs_fallback):
    """返回 (row_interval_fs, output_frequency)。头注释解析失败时回退 sample_fs。"""
    freq = read_output_frequency(path)
    if freq is not None:
        return dt_fs * freq, freq
    if sample_fs_fallback is not None:
        warn(f"{Path(path).name} 未找到 output frequency，回退 --sample-fs "
             f"{sample_fs_fallback} fs")
        return sample_fs_fallback, None
    fail(f"{Path(path).name} 未找到 output frequency，且未提供 --sample-fs 兜底")


def read_charge_file(path):
    """每行 2 列 = (负极电荷, 正极电荷)，单位 |e|。跳过 # 与空行。"""
    negs, poss = [], []
    with open(path, "r") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) < 2:
                continue
            try:
                negs.append(float(parts[0]))
                poss.append(float(parts[1]))
            except ValueError:
                continue
    return negs, poss


def read_potential_file(path):
    """
    每行形如 '  not calc  -0.251866   0.748134'。
    Bulk 列是字符串 (可能占 1~2 个 token，如 'not calc')，电极电位是行内数字。
    稳健做法：收集行内所有可转 float 的 token，取最后两个 = (负极侧电位, 正极侧电位)。
    """
    v_negs, v_poss = [], []
    with open(path, "r") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            nums = []
            for tok in line.split():
                try:
                    nums.append(float(tok))
                except ValueError:
                    continue
            if len(nums) < 2:
                continue
            v_negs.append(nums[-2])
            v_poss.append(nums[-1])
    return v_negs, v_poss


def last_window_mean(series, row_interval_fs, window_fs):
    """取序列最后 window_fs 对应行数的均值；行数不足则取全部。"""
    if not series:
        return None
    rows = math.ceil(window_fs / row_interval_fs)
    rows = min(rows, len(series))
    return sum(series[-rows:]) / rows


# ============================================================
# 逐电压点收集
# ============================================================

def discover_voltage_dirs(system_dir):
    """
    扫描体系目录下所有电压点目录（如 -2V, -1.5V, 0V, 0.5V, 1V, 2V），
    按电压数值从小到大排序，返回目录名列表；找不到任何电压点时 fail。
    """
    found = []
    for child in Path(system_dir).iterdir():
        if not child.is_dir():
            continue
        m = VOLTAGE_DIR_RE.match(child.name)
        if m:
            found.append((float(m.group(1)), child.name))
    found.sort(key=lambda t: t[0])
    if not found:
        fail(f"{system_dir} 下未发现任何电压点目录（命名如 -2V / 0V / 1V）")
    return [name for _, name in found]


def collect_voltage_points(system_dir, dt_fs, window_fs, sample_fs_fallback):
    """
    自动发现电压点目录（-2V, -1.5V, 0V, 0.5V, 1V ...），按电压数值排序后
    逐一解析电荷与电势文件，取末窗口均值。
    返回点列表，每点 dict（含 name 目录名、voltage 数值）；
    缺目录/文件/数据的点带 error 标记。
    """
    vnames = discover_voltage_dirs(system_dir)
    points = []
    for vname in vnames:
        vdir = Path(system_dir) / vname
        m = VOLTAGE_DIR_RE.match(vname)
        voltage = float(m.group(1)) if m else None
        rec = {"name": vname, "voltage": voltage, "error": None}
        if not vdir.is_dir():
            rec["error"] = "目录缺失"
            points.append(rec)
            warn(f"跳过 {vname}：目录不存在 {vdir}")
            continue

        charge_file = vdir / CHARGE_FILE
        pot_file = vdir / POTENTIAL_FILE
        if not charge_file.is_file():
            rec["error"] = f"缺少 {CHARGE_FILE}"
            points.append(rec)
            warn(f"跳过 {vname}：{charge_file} 不存在")
            continue
        if not pot_file.is_file():
            rec["error"] = f"缺少 {POTENTIAL_FILE}"
            points.append(rec)
            warn(f"跳过 {vname}：{pot_file} 不存在")
            continue

        # 电荷文件：行间隔 + 末窗口均值
        negs, poss = read_charge_file(charge_file)
        row_c, freq_c = resolve_row_interval(charge_file, dt_fs, sample_fs_fallback)
        rec["row_interval_fs_charge"] = row_c
        rec["output_frequency_charge"] = freq_c
        rec["rows_last_charge"] = math.ceil(window_fs / row_c)
        rec["Q_neg"] = last_window_mean(negs, row_c, window_fs)
        rec["Q_pos"] = last_window_mean(poss, row_c, window_fs)

        # 电势文件：行间隔 + 末窗口均值
        v_negs, v_poss = read_potential_file(pot_file)
        row_p, freq_p = resolve_row_interval(pot_file, dt_fs, sample_fs_fallback)
        rec["row_interval_fs_pot"] = row_p
        rec["output_frequency_pot"] = freq_p
        rec["rows_last_pot"] = math.ceil(window_fs / row_p)
        rec["V_neg"] = last_window_mean(v_negs, row_p, window_fs)
        rec["V_pos"] = last_window_mean(v_poss, row_p, window_fs)

        # 派生量
        if rec["Q_neg"] is None or rec["Q_pos"] is None:
            rec["error"] = "电荷数据为空"
            warn(f"跳过 {vname}：电荷数据为空")
        elif rec["V_neg"] is None or rec["V_pos"] is None:
            rec["error"] = "电势数据为空"
            warn(f"跳过 {vname}：电势数据为空")
        else:
            rec["Q_abs_neg"] = abs(rec["Q_neg"])
            rec["Q_abs_avg"] = (rec["Q_abs_neg"] + abs(rec["Q_pos"])) / 2.0
            rec["V_diff"] = rec["V_pos"] - rec["V_neg"]
            rec["V_rel"] = None  # 待 0V 零点校正后填入
        points.append(rec)
    return points


def assign_v_rel(points):
    """V_rel = V_diff − V_diff(0V)。0V 缺失时退回 V_diff 并 warn。"""
    zero = next((p for p in points
                 if not p["error"] and p.get("voltage") == 0.0
                 and p["V_diff"] is not None), None)
    if zero is None:
        warn("未找到有效的 0V 点，V_rel 退回 V_diff（无 0V 零点校正）")
        for p in points:
            if not p["error"] and p["V_diff"] is not None:
                p["V_rel"] = p["V_diff"]
        return
    z_diff = zero["V_diff"]
    for p in points:
        if not p["error"] and p["V_diff"] is not None:
            p["V_rel"] = p["V_diff"] - z_diff


# ============================================================
# 电容计算
# ============================================================

def compute_integral(points, mof_mass):
    """积分电容：非 0V 点，C_int = 4*Q_abs_avg*F/(|V_rel|*M)（用户定义 4Q∞/(M·V)）。"""
    for p in points:
        p["C_int_F_per_g"] = None
        if p["error"]:
            continue
        if p["Q_abs_avg"] is None or p["V_rel"] is None or p["V_rel"] == 0:
            p["C_int_F_per_g"] = float("nan")
            continue
        p["C_int_F_per_g"] = (4.0 * p["Q_abs_avg"] * FARADAY
                              / (mof_mass * abs(p["V_rel"])))


def compute_differential(points, mof_mass):
    """微分电容：按电压数值排序后，相邻有效点对，正/负极各自计算，ΔQ 取幅值。"""
    ok = [p for p in points if not p["error"]]
    ok.sort(key=lambda p: p["voltage"] if p.get("voltage") is not None
            else float("-inf"))
    diffs = []
    for i in range(len(ok) - 1):
        p0, p1 = ok[i], ok[i + 1]
        d = {
            "pair": f"{p0['name']}→{p1['name']}",
            "dV": None,
            "dQ_pos": None,
            "dQ_neg": None,
            "C_diff_pos_F_per_g": None,
            "C_diff_neg_F_per_g": None,
        }
        if p0["V_rel"] is not None and p1["V_rel"] is not None:
            d["dV"] = p1["V_rel"] - p0["V_rel"]
        if p0["Q_pos"] is not None and p1["Q_pos"] is not None:
            d["dQ_pos"] = p1["Q_pos"] - p0["Q_pos"]
        if p0["Q_neg"] is not None and p1["Q_neg"] is not None:
            d["dQ_neg"] = abs(p1["Q_neg"]) - abs(p0["Q_neg"])  # 取幅值差分
        if (d["dV"] is not None and d["dV"] != 0
                and d["dQ_pos"] is not None):
            d["C_diff_pos_F_per_g"] = d["dQ_pos"] * FARADAY / (d["dV"] * mof_mass)
        else:
            d["C_diff_pos_F_per_g"] = float("nan")
        if (d["dV"] is not None and d["dV"] != 0
                and d["dQ_neg"] is not None):
            d["C_diff_neg_F_per_g"] = d["dQ_neg"] * FARADAY / (d["dV"] * mof_mass)
        else:
            d["C_diff_neg_F_per_g"] = float("nan")
        diffs.append(d)
    return diffs


def safe_avg(values):
    vals = [v for v in values if v is not None and not math.isnan(v)]
    return sum(vals) / len(vals) if vals else float("nan")


def build_summary(points, diffs):
    c_ints = [p["C_int_F_per_g"] for p in points
              if p.get("C_int_F_per_g") is not None]
    return {
        "n_voltages_total": len(points),
        "n_voltages_ok": sum(1 for p in points if not p["error"]),
        "n_pairs": len(diffs),
        "C_int_avg_F_per_g": safe_avg(c_ints),
        "C_int_max_F_per_g": max(c_ints) if c_ints else float("nan"),
        "C_diff_pos_avg_F_per_g": safe_avg(
            [d["C_diff_pos_F_per_g"] for d in diffs]),
        "C_diff_neg_avg_F_per_g": safe_avg(
            [d["C_diff_neg_F_per_g"] for d in diffs]),
    }


# ============================================================
# 输出
# ============================================================

def fmt(v, width, nd=4):
    if v is None or (isinstance(v, float) and math.isnan(v)):
        return f"{'-':>{width}}"
    return f"{v:>{width}.{nd}f}"


def print_integral_table(points):
    print("\n==================== 积分电容 ====================")
    print(f"{'点':<4} {'V_rel(V)':>10} {'Q_neg(e)':>12} {'Q_pos(e)':>12} "
          f"{'Q_abs_avg(e)':>13} {'C_int(F/g)':>12}")
    for p in points:
        if p["error"]:
            print(f"{p['name']:<4}  (跳过: {p['error']})")
            continue
        print(f"{p['name']:<4} {fmt(p['V_rel'], 10)} "
              f"{fmt(p['Q_neg'], 12)} {fmt(p['Q_pos'], 12)} "
              f"{fmt(p['Q_abs_avg'], 13)} {fmt(p['C_int_F_per_g'], 12)}")


def print_differential_table(diffs):
    print("\n==================== 微分电容 ====================")
    print(f"{'区间':<14} {'dV(V)':>8} {'dQ_pos(e)':>12} {'dQ_neg(e)':>12} "
          f"{'C_diff_pos(F/g)':>16} {'C_diff_neg(F/g)':>16}")
    for d in diffs:
        print(f"{d['pair']:<14} {fmt(d['dV'], 8)} {fmt(d['dQ_pos'], 12)} "
              f"{fmt(d['dQ_neg'], 12)} {fmt(d['C_diff_pos_F_per_g'], 16)} "
              f"{fmt(d['C_diff_neg_F_per_g'], 16)}")


def write_json(result, out_path):
    with open(out_path, "w") as f:
        json.dump(result, f, ensure_ascii=False, indent=2, allow_nan=True)
    print(f"JSON 结果写入: {out_path}")


CSV_COLUMNS = [
    "kind", "label", "V_rel", "dV",
    "Q_neg", "Q_abs_neg", "Q_pos", "Q_abs_avg",
    "dQ_pos", "dQ_neg",
    "C_int_F_per_g", "C_diff_pos_F_per_g", "C_diff_neg_F_per_g",
]


def write_csv(points, diffs, out_path):
    import csv as csv_mod
    rows = []
    for p in points:
        if p["error"]:
            continue
        rows.append({
            "kind": "voltage", "label": p["name"],
            "V_rel": p["V_rel"], "dV": None,
            "Q_neg": p["Q_neg"], "Q_abs_neg": p["Q_abs_neg"],
            "Q_pos": p["Q_pos"], "Q_abs_avg": p["Q_abs_avg"],
            "dQ_pos": None, "dQ_neg": None,
            "C_int_F_per_g": p["C_int_F_per_g"],
            "C_diff_pos_F_per_g": None, "C_diff_neg_F_per_g": None,
        })
    for d in diffs:
        rows.append({
            "kind": "pair", "label": d["pair"],
            "V_rel": None, "dV": d["dV"],
            "Q_neg": None, "Q_abs_neg": None,
            "Q_pos": None, "Q_abs_avg": None,
            "dQ_pos": d["dQ_pos"], "dQ_neg": d["dQ_neg"],
            "C_int_F_per_g": None,
            "C_diff_pos_F_per_g": d["C_diff_pos_F_per_g"],
            "C_diff_neg_F_per_g": d["C_diff_neg_F_per_g"],
        })
    with open(out_path, "w", newline="") as f:
        writer = csv_mod.DictWriter(f, fieldnames=CSV_COLUMNS)
        writer.writeheader()
        writer.writerows(rows)
    print(f"CSV 结果写入: {out_path}")


# ============================================================
# 主流程
# ============================================================

def main():
    parser = argparse.ArgumentParser(
        description="CPM 恒电势结果 -> MOF 电极质量电容分析 (积分+微分)"
    )
    parser.add_argument("--system-dir", required=True,
                        help="体系目录 (其下有 0V..4V 与 system_summary.json)")
    parser.add_argument("--dt-fs", type=float, required=True,
                        help="MD 计算间隔 (fs)，用于换算每行时间，如 2.0")
    parser.add_argument("--window-ns", type=float, default=5.0,
                        help="末窗口时长 (ns)，默认 5")
    parser.add_argument("--output-dir", default=None,
                        help="输出目录，默认 = 体系目录")
    parser.add_argument("--prefix", default="capacitance",
                        help="输出文件名前缀，默认 capacitance")
    parser.add_argument("--sample-fs", type=float, default=None,
                        help="兜底：每行时间间隔 (fs)，头注释解析失败时使用")
    parser.add_argument("--mof-mass", type=float, default=None,
                        help="覆盖电极质量 (Da)，默认读 system_summary.json")
    args = parser.parse_args()

    system_dir = Path(args.system_dir).resolve()
    if not system_dir.is_dir():
        fail(f"体系目录不存在：{system_dir}")
    # 若 --system-dir 指向的是电压点目录 (如 ACN/1V)，自动回退到其父目录
    # (system_summary.json 与 0V/1V.. 同级，位于 ACN 根)。
    if VOLTAGE_DIR_RE.match(system_dir.name):
        warn(f"{system_dir.name} 是电压点目录，自动使用父目录作为体系根："
             f"{system_dir.parent}")
        system_dir = system_dir.parent
    if args.dt_fs <= 0:
        fail(f"--dt-fs 必须为正数：{args.dt_fs}")

    mof_mass, electrode_stats = load_mof_mass(system_dir, args.mof_mass)
    window_fs = args.window_ns * 1e6

    print(f"\n体系目录 : {system_dir}")
    print(f"dt(计算间隔) : {args.dt_fs} fs | 末窗口 : {args.window_ns} ns")
    print(f"电极质量(用于电容计算) : {mof_mass} Da (=g/mol) | F = {FARADAY} C/mol")

    # 1. 逐电压点收集
    points = collect_voltage_points(
        system_dir, args.dt_fs, window_fs, args.sample_fs
    )

    # 2. 0V 零点校正 V_rel
    assign_v_rel(points)

    # 3. 积分电容
    compute_integral(points, mof_mass)

    # 4. 微分电容
    diffs = compute_differential(points, mof_mass)

    # 5. 汇总
    summary = build_summary(points, diffs)

    # 6. 终端输出
    print_integral_table(points)
    print_differential_table(diffs)
    print(f"\n==================== 汇总 ====================")
    for k, v in summary.items():
        if isinstance(v, float):
            print(f"  {k:28s}: {v:.4f}")
        else:
            print(f"  {k:28s}: {v}")

    # 7. 文件输出
    out_dir = Path(args.output_dir).resolve() if args.output_dir \
        else system_dir
    out_dir.mkdir(parents=True, exist_ok=True)
    json_path = out_dir / f"{args.prefix}.json"
    csv_path = out_dir / f"{args.prefix}.csv"

    result = {
        "system_name": system_dir.name,
        "dt_fs": args.dt_fs,
        "window_ns": args.window_ns,
        "mof_mass_g_mol": mof_mass,
        "mof_electrode_statistics": electrode_stats,
        "faraday_C_mol": FARADAY,
        "metadata": {
            p["name"]: {
                "output_frequency_charge": p["output_frequency_charge"],
                "row_interval_fs_charge": p["row_interval_fs_charge"],
                "rows_last_charge": p["rows_last_charge"],
                "output_frequency_pot": p["output_frequency_pot"],
                "row_interval_fs_pot": p["row_interval_fs_pot"],
                "rows_last_pot": p["rows_last_pot"],
            }
            for p in points if not p["error"] and "row_interval_fs_charge" in p
        },
        "voltages": [
            {k: p[k] for k in
             ("name", "voltage", "V_neg", "V_pos", "V_diff", "V_rel",
              "Q_neg", "Q_abs_neg", "Q_pos", "Q_abs_avg", "C_int_F_per_g")}
            if not p["error"] else {"name": p["name"], "voltage": p.get("voltage"),
                                    "error": p["error"]}
            for p in points
        ],
        "differential": diffs,
        "summary": summary,
    }
    write_json(result, json_path)
    write_csv(points, diffs, csv_path)

    print("\n完成。")


if __name__ == "__main__":
    main()
