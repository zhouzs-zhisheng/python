# 成分调整脚本（0V 密度判定）设计方案

> 目标：通过调整体系分子数（只能在现有基础上增量增加），让电解液体系（如 ACN 或 PC）的**体相密度**达到目标设定值。
> 核心不同点：**密度判定不用 NVT 轨迹，而用 0V 平衡模拟（NVE）轨迹**。
> 本文档为**设计思路存档**，尚未实现代码，仅保存当前确定的需求与设计，供下次继续开发。

---

## 一、总体思路

调度"成分（分子数）"的外循环，套一层"每个候选成分先做普通驰豫（min/NVT），再跑一段 0V 平衡（NVE），最后在 0V 轨迹上测密度"的流程。

两套旧代码各取其一：
- 增量和续跑（NVT 阶段） ← 参考 `adjust_density.py`
- 0V 模拟 + 0V 密度判定 ← 参考 `cpm_equilibrium_loop.py`

---

## 二、已确认的关键决策点

1. **MOL1 / MOL2 命名与 `topol.top` 的 `[moleculetype]` 完全一致**。
2. 每次生长时 `topol.top` 的 `[molecules]` 更新**仍由 `incremental_add_molecules.py` 负责**，脚本只负责拼接参数并调用它。
3. 过密度微冲时的裁量：**允许极小正值容差**（见收敛判据一节）。

---

## 三、从 adjust_density 摘取的重要思路

### 1. 增量生长（从无到有、逐步加分子）
- 体系不能删分子，只能**增加**。每次生长在上一体系基础上按比例插入指定分子。
- 每个候选成分 = 一个数字目录 `<N>`（N 为 MOL1 数）。首轮需要一个**种子体系**（最低成分），之后每次从"父体系"增量长出更高成分候选。
- 生长来源（本次修订）：**父体系的 `0V/nve.gro`** → `gmx editconf` 转 `parent_system.pdb` → 交给增量脚本按比例插入新分子、重写 `topol.top` 与 `system_summary.json`，产出 `system.gro / topol.top / system_summary.json`。
- 中途中断容错：目标数字目录已存在但不完整 → 整体备份为 `*.incomplete_backup_*` → 从父体系重新生成。

### 2. 幂等 + 体系续跑（checkpoint 恢复）
- `index.ndx` 已存在 → 跳过 `make_ndx`。
- `first/min.gro` 已存在 → 跳过 minimization。
- NVT：`nvt.gro + nvt.tpr + nvt.xtc` 齐全 = 完成，跳过；只有 `nvt.cpt` → `mdrun -cpi -append` 续跑；否则从 `nvt.tpr` 重新起跑。
- 最终判定以 `nvt.gro` 存在为准（仅 `nvt.xtc` 不算完成，因为中断常会留下轨迹）。

### 3. 状态持久化（重启不丢进度）
- 根目录维护 `composition_adjustment.json` + `.csv`，记录目标密度、coarse/fine 阶段、全部已测 `measurements`（成分+密度）、当前 trial、失败记录。
- 已测过的密度不重复计算；同时写一份全局 `density_profiles.csv` 汇总。

### 4. 收敛分两段 + 上下夹逼
- coarse 区间、fine 区间均为**单侧**（见收敛判据一节）。
- 精调避免过渡：用 **bracket（下界 < 目标 < 上界）** 线性内插估算中间 MOL 数。

### 5. 失败二分回退
- 模拟失败（分子过多等原因）→ 备份失败目录 → 记录到 `failed_trials` → 用 `new_target = parent + step//2` 从父体系重试。

> 不摘取：`adjust_density` 的"5ns 续跑稳定性校验"（`run_post_convergence_check`）——它被"0V 密度判定"取代。

---

## 四、从 cpm_equilibrium_loop 摘取的重要思路（0V 判定）

### 1. 0V 平衡（NVE）—— 单段固定时长，无续跑
- 0V 平衡**只跑一段**：用 `grompp_20ns.mdp` 生成 `nve.tpr`，`mdrun` 一次跑到底（20 ns）。
- **删除** 一切 `convert-tpr -extend` / `mdrun -cpi -append` 续跑逻辑（本次修订）。
- 因此只需两个可执行动作：`grompp`（`-maxwarn 2`）+ `mdrun`。产物齐全（`nve.gro / nve.xtc / nve.tpr`）即视为完成。
- **0V 不做电荷收敛**：完全不读 `CPM_electrodeCharge.dat`，不做任何电荷/时间收敛判定，固定时长跑完即 `gmx density`。
- 0V 只需要一个 mdp：**`grompp_20ns.mdp`（20ns）**；不再需要续跑用的 `grompp.mdp`（10ns）。

### 2. 0V 密度测量
- `gmx density -f nve.xtc -s nve.tpr -n index.ndx -d Z -b <取轨迹后段窗口> -sl <切片数> -o density.xvg`。
- 用 `system_summary.json` 的体相 Z 区间对 profile 求平均 = `ρ_0V`，作为该候选体系的**官方密度**。

---

## 五、收敛判据（单侧区间）

- 记 `error = (ρ_0V − target) / target`。
- **coarse 区间**：`-0.03 ≤ error ≤ 0`；首次落到该区间 → 记 coarse checkpoint。
- **fine 区间**：`-0.01 ≤ error ≤ 0`；落到该区间 → 精调收敛，结束。
- **不允许过密度（`error > 0` 原则上不算达标）**；但因"只增不减且密度随分子数大致单调"，最后一步可能微过冲。
- **过密度微冲裁量**（已确认）：**允许极小正值容差**通过，即存在一个很小的正容差阈值（如 `+OVER_TOL`），当 `0 < error ≤ OVER_TOL` 时视为"最佳贴近"，允许算作完成；超出该容差的正过冲则需二分缩小增量、从 lower 体系重生长逼近。

---

## 六、主流程

1. 读取文件开头的全局配置（目标密度、MOL1/MOL2 与路径、RATIO、coarse/fine 区间、0V 20ns、Z 区间、两个 gmx 路径、maxwarn 等）；加载或创建状态。
2. 确定当前成分（从状态或最低种子开始）。
3. 若种子（或当前候选）尚无 `0V/nve.gro` → 单段 20ns 0V，得到模板。
4. **生长候选**：父 `0V/nve.gro` → pdb → 按 `RATIO / DELTA_MAX` 插入 MOL1/MOL2（及溶剂）→ `<N>/system.gro, topol.top, system_summary.json`。
5. **驰豫**：index → min → NVT（用 `--gmx-minnvt`，幂等 + checkpoint 续跑）。
6. **0V 平衡**：`<N>/0V/`（用 `--gmx-0v`）单段 20ns，固定时长即停，无续跑、无电荷判定。
7. **0V 测密度**：`gmx density` → 体相 Z 区间平均 `ρ_0V`。
8. **单侧判定**：见收敛判据一节。
9. 未收敛 → 从 lower 体系按比例缩小 step 重生长 → 回 4；迭代上限未果则告警。
10. 每步写状态 json/csv；输出汇总与最终结果 json。

---

## 七、全局参数（建议集中放在文件开头，改造即改这里）

```python
# ---- 分子/插入规则 ----
MOL1 = "EMIM"              # 与 topol.top 的 [moleculetype] 完全一致
MOL2 = "BF4"
MOL1_DIR = "/path/to/EMIM.pdb"   # mol1 的单分子结构 pdb 路径
MOL2_DIR = "/path/to/BF4.pdb"    # mol2 的单分子结构 pdb 路径
RATIO = (1, 1)                    # mol1 : mol2

LIGAND = "ACN"             # 溶剂分子（可选，0 则不插溶剂）
LIGAND_DIR = "/path/to/ACN.pdb"
LIGAND_RATIO = (1, 1, 5)   # 每组 mol1/mol2 附带 5 个溶剂

DELTA_MAX = (1, 1, 5)      # 单步最大增量批次

# ---- 目标与判据 ----
TARGET_DENSITY = 1.20      # g/cm3
COARSE_LOW,  COARSE_HIGH  = -0.03, 0.0
FINE_LOW,    FINE_HIGH    = -0.01, 0.0
OVER_TOL = 0.001            # 过密度微冲允许的小正容差（已确认）

# ---- 0V 模拟 ----
ZERO_V_MDP = "grompp_20ns.mdp"   # 单段 20ns，无续跑
ZERO_V_GROUPS = 6
DENSITY_GROUP = 6
DENSITY_START_NS = 10.0          # 取轨迹后段窗口起点
DENSITY_SLICES = 4141
Z_STRUCTURE_REGIONS = (18.0, 22.0)   # 体相 Z 区间（取自 system_summary.json）

# ---- 环境 / 命令 ----
GMX_MIN_NVT = "gmx"          # min/NVT 用的 GROMACS
GMX_ZERO_V  = "gmx"          # 0V 用的 GROMACS（两个版本不同路径）
MAXWARN = 2
NTMPI = 1
NTOMP = 32
GPU_ARGS = ("-nb gpu -pme gpu -pmefft gpu -tunepme no")

# ---- 流程 ----
MAX_LOOPS = 30               # 迭代上限
```

---

## 八、需要用户提供的文件与目录路径

```
MOF总目录/
├── basic/cpm_mdp/
│   ├── min.mdp                    # ① min 参数
│   └── nvt.mdp                    # ① NVT 参数
├── allMatrixA.bin                 # ⑥ 0V 恒电势矩阵
├── 0V/CPM_ControlFile.dat         # ⑥ 0V 控制文件
└── qmof-305c717/
    ├── ACN/                       # solvent_root
    │   ├── system.gro             # ② 种子体系坐标（最低成分）
    │   ├── topol.top              # ③ 种子拓扑
    │   ├── system_summary.json    # ④ 体相 Z 区间/盒子
    │   ├── index.ndx              # ⑤ 组定义（或脚本生成）
    │   ├── grompp_20ns.mdp        # ⑥ 0V NVE 20ns（唯一 0V mdp，无续跑）
    │   ├── fine/                  # 已有可复用，否则种子生长
    │   ├── <N>/
    │   │   ├── system.gro / topol.top / system_summary.json  （脚本生成）
    │   │   ├── first/min.*, nvt.*                            （脚本生成）
    │   │   └── 0V/index.ndx, nve.*, density.xvg              （脚本生成）
    │   └── composition_adjustment.json / .csv + density_profiles.csv（脚本生成）
    └── scripts/
        └── incremental_add_molecules.py   # ⑦ 增量加分子
```

**必填清单：**
1. `min.mdp`、`nvt.mdp`（驰豫参数）
2. 种子体系 `system.gro`（最低成分坐标）
3. `topol.top`（含 MOL1/MOL2（及溶剂）定义与力场，`[moleculetype]` 名与 MOL1/MOL2 一致）
4. `system_summary.json`（体相 Z 区间、盒子尺寸）
5. `index.ndx`（或给出密度组定义让脚本 `make_ndx`）
6. 0V mdp：`grompp_20ns.mdp`（20ns，唯一）+ 0V 配套（`CPM_ControlFile.dat`、`allMatrixA.bin`）
7. `incremental_add_molecules.py`
8. **两个 GROMACS 路径**：`--gmx-minnvt`（min/NVT）、`--gmx-0v`（0V）
9. **插入分子模板 pdb**：`MOL1_DIR`（如 EMIM.pdb）、`MOL2_DIR`（如 TFO.pdb）；溶剂则给 `LIGAND_DIR`（如 ACN.pdb）

---

## 九、其他约定与开放事项

- 5ns 稳定性校验不摘取；由 0V 密度判定取代。
- 0V 恒电势配套（`CPM_ControlFile.dat` + `allMatrixA.bin`）由脚本软链/复制到各候选 `0V/` 目录（沿用 cpm 脚本思路）。
- 仅一处未完全敲定：溶液是 ACN 溶剂恒定，还是也参与增量（由 `LIGAND` 与 `LIGAND_RATIO` 是否=0 决定）。
- 剩余未实现：代码正文。本文档为设计存档，供后续开发实现。