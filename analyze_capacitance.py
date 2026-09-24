#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
analyze_capacitance.py

从 CPM 恒电势模拟结果计算 MOF 电极的质量电容（积分电容 + 微分电容）。

数据来源（每个电压点目录 0V/1V/2V/3V/4V 下）：
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

积分电容（每个非 0V 电压点，用户定义 C_int = 4Q∞/(M·V)）：
  C_int [F/g] = 4 * Q_abs_avg * F / (V_rel * M)
  Q_abs_avg   = (|Q_neg| + |Q_pos|) / 2   （正负电极等效电荷，末窗口均值，|e|）
  V_rel       = 该点实测电极间电位差 − 0V 实测电极间电位差 (V)
  M           = 电极质量 (total_mof_mass，Da 数值直接用)
  F           = 法拉第常数 96485 C/mol
  推导：q_C = Q_e*F/N_A，m_g = M/N_A  →  C_g = Q_e*F/(V*M)（N_A 抵消）；
  乘 4 为电极对几何换算因子（Q∞ 即该电压下平衡电荷）。

微分电容（9 点电极电位曲线，每点一个电荷量）：
  每个模拟目录 (0V..4V) 的 CPM_potential.dat / CPM_electrodeCharge.dat 给出两个电极：
    (V_neg, Q_neg) 与 (V_pos, Q_pos)   （V 为实测电极电位，Q 为末窗口均值电荷）
  对应关系（9 个电极电位点来自 0-4V 的 5 个模拟）：
    4V → -2V    3V → -1.5V  2V → -1V  1V → -0.5V
    0V → 0   （0V 两电极电位≈0，合并为一个点）
    1V → +0.5V  2V → +1V    3V → +1.5V  4V → +2V
  按电极电位排序后取相邻点差分：
    dV          = V_i+1 − V_i
    dQ          = Q_i+1 − Q_i              （带符号；曲线单调递增 → dQ 为正值）
    C_diff      = factor * dQ * F / (dV * M)，factor 默认 1（可用 --diff-factor 改，
                 如与积分 4Q∞/(M·V) 一致可设 4）
  同时输出电荷量随电压值的变化（电极电位曲线表，见 CSV kind=point）。

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

VOLTAGE_DIRS = ["0V", "1V", "2V", "3V", "4V"]

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

def collect_voltage_points(system_dir, dt_fs, window_fs, sample_fs_fallback):
    """
    遍历 0V..4V，解析电荷与电势文件，取末窗口均值。
    返回点列表，每点 dict；缺目录/文件/数据的点带 error 标记。
    """
    points = []
    for vname in VOLTAGE_DIRS:
        vdir = Path(system_dir) / vname
        rec = {"name": vname, "error": None}
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
    zero = next((p for p in points if p["name"] == "0V" and not p["error"]), None)
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
    """积分电容：非 0V 点，C_int = 4*Q_abs_avg*F/(M*V_rel)（用户定义 4Q∞/(M·V)）。"""
    for p in points:
        p["C_int_F_per_g"] = None
        if p["error"] or p["name"] == "0V":
            continue
        if p["Q_abs_avg"] is None or p["V_rel"] is None or p["V_rel"] <= 0:
            p["C_int_F_per_g"] = float("nan")
            continue
        p["C_int_F_per_g"] = (4.0 * p["Q_abs_avg"] * FARADAY
                              / (mof_mass * p["V_rel"]))


def safe_avg(values):
    vals = [v for v in values if v is not None and not math.isnan(v)]
    return sum(vals) / len(vals) if vals else float("nan")


def build_summary(points, diffs, curve):
    c_ints = [p["C_int_F_per_g"] for p in points
              if p.get("C_int_F_per_g") is not None]
    c_diffs = [d["C_diff_F_per_g"] for d in diffs
               if d.get("C_diff_F_per_g") is not None]
    return {
        "n_voltages_total": len(points),
        "n_voltages_ok": sum(1 for p in points if not p["error"]),
        "n_curve_points": len(curve),
        "n_pairs": len(diffs),
        "C_int_avg_F_per_g": safe_avg(c_ints),
        "C_int_max_F_per_g": max(c_ints) if c_ints else float("nan"),
        "C_diff_avg_F_per_g": safe_avg(c_diffs),
    }


# ============================================================
# 微分电容：9 点电极曲线
# ============================================================

def build_electrode_curve(points, zero_charge_mode="abs_avg"):
    """
    由各模拟目录 (0V..4V) 构建 9 点电极电位曲线。

    每个模拟目录的 CPM_potential.dat / CPM_electrodeCharge.dat 给出两个电极：
      (V_neg, Q_neg) 与 (V_pos, Q_pos)   （V 为实测电极电位，Q 为末窗口均值电荷）
    对应关系（9 个电极电位点来自 0-4V 的 5 个模拟）：
      4V -> -2V   3V -> -1.5V  2V -> -1V  1V -> -0.5V
      0V -> 0
      1V -> +0.5V 2V -> +1V    3V -> +1.5V 4V -> +2V
    其中 0V 目录的两个电极电位都≈0，合并为一个 "0" 点。

    返回按电极电位升序排列的点列表，每点 dict:
      potential, charge, sim(来源目录), electrode('neg'/'pos'/'mid')
    """
    curve = []
    zero = None
    for p in points:
        if p["error"]:
            continue
        if p["name"] == "0V":
            # 0V：两个电极电位≈0，合并为一个中点
            pot = (p["V_neg"] + p["V_pos"]) / 2.0
            if zero_charge_mode == "abs_avg":
                q = p["Q_abs_avg"]
            elif zero_charge_mode == "zero":
                q = 0.0
            elif zero_charge_mode == "neg":
                q = p["Q_neg"]
            else:  # default abs_avg
                q = p["Q_abs_avg"]
            zero = {"potential": pot, "charge": q,
                    "sim": p["name"], "electrode": "mid"}
            continue
        curve.append({"potential": p["V_neg"], "charge": p["Q_neg"],
                      "sim": p["name"], "electrode": "neg"})
        curve.append({"potential": p["V_pos"], "charge": p["Q_pos"],
                      "sim": p["name"], "electrode": "pos"})
    if zero is not None:
        curve.append(zero)
    curve.sort(key=lambda c: c["potential"])
    return curve


def compute_differential_curve(curve, mof_mass, factor=1.0):
    """
    微分电容：9 点电极曲线中相邻两点的 dQ/dV。
    返回 8 个区间，每区间 dict:
      pair, side(neg/pos/neg-pos), V_from, V_to, dV, dQ, C_diff_F_per_g
    """
    diffs = []
    for i in range(len(curve) - 1):
        c0, c1 = curve[i], curve[i + 1]
        side = c0["electrode"]
        if c0["electrode"] != c1["electrode"]:
            side = f"{c0['electrode']}-{c1['electrode']}"
        d = {
            "pair": f"{c0['potential']:.3f}V→{c1['potential']:.3f}V",
            "side": side,
            "V_from": c0["potential"],
            "V_to": c1["potential"],
            "dV": c1["potential"] - c0["potential"],
            "dQ": c1["charge"] - c0["charge"],
            "C_diff_F_per_g": None,
        }
        if (d["dV"] is not None and d["dV"] != 0
                and d["dQ"] is not None):
            d["C_diff_F_per_g"] = (factor * d["dQ"] * FARADAY
                                   / (d["dV"] * mof_mass))
        else:
            d["C_diff_F_per_g"] = float("nan")
        diffs.append(d)
    return diffs


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


def print_curve_table(curve):
    print("\n==== 电极电位曲线 (9 点，电荷量随电压值的变化) ====")
    print(f"{'电位(V)':>10} {'电荷(e)':>12} {'来源':>6} {'电极':>6}")
    for c in curve:
        print(f"{fmt(c['potential'], 10, 3)} {fmt(c['charge'], 12)} "
              f"{c['sim']:>6} {c['electrode']:>6}")


def print_differential_table(diffs):
    print("\n==================== 微分电容 ====================")
    print(f"{'区间':<16} {'侧':>10} {'dV(V)':>8} {'dQ(e)':>10} "
          f"{'C_diff(F/g)':>14}")
    for d in diffs:
        print(f"{d['pair']:<16} {d['side']:>10} {fmt(d['dV'], 8)} "
              f"{fmt(d['dQ'], 10)} {fmt(d['C_diff_F_per_g'], 14)}")


def write_json(result, out_path):
    with open(out_path, "w") as f:
        json.dump(result, f, ensure_ascii=False, indent=2, allow_nan=True)
    print(f"JSON 结果写入: {out_path}")


CSV_COLUMNS = [
    "kind", "label", "V_rel", "dV",
    "Q_neg", "Q_abs_neg", "Q_pos", "Q_abs_avg",
    "potential", "charge", "sim", "electrode",
    "dQ",
    "C_int_F_per_g", "C_diff_F_per_g",
]


def write_csv(points, diffs, curve, out_path):
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
            "potential": None, "charge": None, "sim": None, "electrode": None,
            "dQ": None,
            "C_int_F_per_g": p["C_int_F_per_g"],
            "C_diff_F_per_g": None,
        })
    for c in curve:
        rows.append({
            "kind": "point", "label": f"{c['potential']:.3f}V",
            "V_rel": None, "dV": None,
            "Q_neg": None, "Q_abs_neg": None,
            "Q_pos": None, "Q_abs_avg": None,
            "potential": c["potential"], "charge": c["charge"],
            "sim": c["sim"], "electrode": c["electrode"],
            "dQ": None,
            "C_int_F_per_g": None, "C_diff_F_per_g": None,
        })
    for d in diffs:
        rows.append({
            "kind": "pair", "label": d["pair"],
            "V_rel": None, "dV": d["dV"],
            "Q_neg": None, "Q_abs_neg": None,
            "Q_pos": None, "Q_abs_avg": None,
            "potential": None, "charge": None,
            "sim": None, "electrode": d["side"],
            "dQ": d["dQ"],
            "C_int_F_per_g": None,
            "C_diff_F_per_g": d["C_diff_F_per_g"],
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
    parser.add_argument("--zero-charge", choices=["abs_avg", "zero", "neg"],
                        default="abs_avg",
                        help="微分曲线 '0' 点电荷取值：abs_avg=正负电极幅值均值"
                             "(默认)，zero=0，neg=负电极实测电荷")
    parser.add_argument("--diff-factor", type=float, default=1.0,
                        help="微分电容乘数，默认 1（即 ΔQ*F/(ΔV*M)）；"
                             "若需与积分 4Q∞/(M·V) 一致可设 4")
    args = parser.parse_args()

    system_dir = Path(args.system_dir).resolve()
    if not system_dir.is_dir():
        fail(f"体系目录不存在：{system_dir}")
    # 若 --system-dir 指向的是电压点目录 (如 ACN/1V)，自动回退到其父目录
    # (system_summary.json 与 0V/1V.. 同级，位于 ACN 根)。
    if system_dir.name in VOLTAGE_DIRS:
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

    # 1. 逐电压点收集 (0V..4V)
    points = collect_voltage_points(
        system_dir, args.dt_fs, window_fs, args.sample_fs
    )

    # 2. 0V 零点校正 V_rel
    assign_v_rel(points)

    # 3. 积分电容 (4Q∞/(M·V))
    compute_integral(points, mof_mass)

    # 4. 微分电容：9 点电极曲线 -> 相邻点差分
    curve = build_electrode_curve(points, args.zero_charge)
    diffs = compute_differential_curve(curve, mof_mass, args.diff_factor)

    # 5. 汇总
    summary = build_summary(points, diffs, curve)

    # 6. 终端输出
    print_integral_table(points)
    print_curve_table(curve)
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
        "integral_formula": "C_int = 4 * Q_abs_avg * F / (M * V_rel)",
        "differential_formula": f"C_diff = {args.diff_factor} * dQ * F / (dV * M)",
        "zero_charge_mode": args.zero_charge,
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
             ("name", "V_neg", "V_pos", "V_diff", "V_rel",
              "Q_neg", "Q_abs_neg", "Q_pos", "Q_abs_avg", "C_int_F_per_g")}
            if not p["error"] else {"name": p["name"], "error": p["error"]}
            for p in points
        ],
        "electrode_curve": curve,
        "differential": diffs,
        "summary": summary,
    }
    write_json(result, json_path)
    write_csv(points, diffs, curve, csv_path)

    print("\n完成。")


if __name__ == "__main__":
    main()
