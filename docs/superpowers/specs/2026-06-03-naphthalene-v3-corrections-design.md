# 萘晶体 v3 纠错计算 — 设计规范

> 基于 naphthalene-transport-v2 的错误分析文档
> 目标: 修正4个关键公式/流程错误，重新计算完整的转移积分→迁移率管线
> 方案A: 仅修复Python后处理代码，复用v2全部VASP原始数据

---

## 1. 错误总览

基于 `docs/映射模型哈密顿量计算错误分析.md` 的深度审计，识别出4个关键错误：

| # | 位置 | 描述 | 影响 | 修复难度 |
|---|------|------|------|----------|
| **E1** | Stage 6 | 时间积分缺少因子2 (偶函数) | μ偏小2× | 简单 |
| **E2** | Stage 5 | A₀/B求和遗漏反向近邻路径 | μ偏小2× | 简单 |
| **E3** | Stage 1 | H_AB(k)未使用P2₁/a对称性 | TB拟合失真, t值偏差 | 中等 |
| **E4** | Stage 3 | 有限差分中的符号简并风险 | ∂t/∂Q可能灾难性错误 | 简单 |

E1+E2叠加使迁移率系统性偏小4倍。E3影响TB参数的所有下游结果。E4可能导致电声耦合导数完全错误。

---

## 2. 项目架构

### 2.1 目录结构

```
naphthalene-transport-v3/
├── archive_v2/                                  # v2项目完整备份
│   └── naphthalene-transport-v2/                # 原v2所有文件
├── original_docs/                               # 原始文献 (符号链接)
│   └── model_hamiltonian/ -> ../../naphthalene-transport-v2/original_docs/model_hamiltonian/
├── code/                                        # v3修改后的代码
│   ├── config.yaml                              # 全局配置
│   ├── validate_pipeline.py                     # 端到端验证脚本 (更新)
│   ├── fix_poscar.py                            # POSCAR修复 (不变,从v2复制)
│   ├── naph_transport/                          # Python包 v3.0.0
│   │   ├── __init__.py                          # 版本号更新
│   │   ├── constants.py                         # 物理常数 (不变)
│   │   ├── crystal_utils.py                     # 晶体结构 (不变,已有正确H_AB)
│   │   ├── io_utils.py                          # I/O解析 (不变)
│   │   ├── validator.py                         # 验证工具 (更新checkpoints)
│   │   ├── stage1_transfer.py                   # [修改-E3] 使用build_tb_H_AB
│   │   ├── stage2_phonon.py                     # 声子后处理 (不变)
│   │   ├── stage3_eph.py                        # [修改-E4] 连续性校验符号对齐
│   │   ├── stage4_assemble.py                   # 参数组装 (不变)
│   │   ├── stage5_mobility_assemble.py          # [修改-E2] A₀/B ×2因子
│   │   └── stage6_mobility_integrate.py         # [修改-E1] I₀/I_λ ×2因子
│   ├── slurm/                                   # Slurm脚本 (预留给后续VASP重算)
│   │   ├── 01_phonopy_dfpt.slurm                # DFPT (从v2复制)
│   │   └── 03_eph_array_v3.slurm                # e-ph array (从v2复制)
│   └── vasp_inputs/                             # VASP输入模板 (从v2复制)
│       ├── INCAR_dfpt
│       ├── INCAR_eph
│       └── KPOINTS_444
├── data/                                        # 计算数据
│   ├── crystal/ -> ../archive_v2/naphthalene-transport-v2/data/crystal/
│   ├── phonon/  -> ../archive_v2/naphthalene-transport-v2/data/phonon/
│   ├── eph/     -> ../archive_v2/naphthalene-transport-v2/data/eph/
│   └── results/                                 # v3新结果输出
├── docs/
│   ├── v3_error_fixes.md                        # 纠错说明 (中文)
│   └── superpowers/
│       ├── specs/2026-06-03-naphthalene-v3-corrections-design.md  # 本文档
│       └── plans/2026-06-03-naphthalene-v3-corrections-plan.md    # 实现计划
└── PIPELINE_COMPLETION_REPORT_V3.md             # v3完成报告
```

### 2.2 从v2复用的数据

| 数据 | 路径 | 用途 | 方式 |
|------|------|------|------|
| 优化结构 | `data/crystal/CONTCAR` | 所有VASP计算的结构输入 | 符号链接 |
| DFT能带 | `data/crystal/band/EIGENVAL` | Stage 1 TB拟合参考数据 | 符号链接 |
| WAVECAR | `data/crystal/WAVECAR` | SCF加速 (如需) | 符号链接 |
| CHGCAR | `data/crystal/CHGCAR` | SCF加速 (如需) | 符号链接 |
| DFPT力常数 | `data/phonon/FORCE_CONSTANTS` | Stage 2声子分析 | 符号链接 |
| 176×EIGENVAL | `data/eph/eph_mode_*/` | Stage 3 TB重拟合 | 符号链接 |
| 176×OUTCAR | `data/eph/eph_mode_*/` | Stage 3 力提取 | 符号链接 |

### 2.3 不再需要的v2计算结果

v2的全部Python计算结果（stage1-6的.npz文件）将被v3新结果替换，但保留在archive_v2中作为对比参考。

---

## 3. 错误修复详细设计

### 3.1 E1: 时间积分偶函数修正 (Stage 6)

**物理原理**: 根据Kubo公式 (Wang 2007 Eq. 2-3)，函数 $f(T)$ 和 $h_q(T)$ 的时间积分域为 $(-\infty, \infty)$。被积函数 $C(t)\cos[S(t)]$ 和 hopping项的被积函数均为 $t$ 的**偶函数**，因此 $\int_{-\infty}^{\infty} = 2\int_{0}^{\infty}$。

**代码位置**: `stage6_mobility_integrate.py:compute_mobility_at_T()`

**修复**:
```python
# 修改前
I0 = simpson(C * np.cos(S), x=t_array)
I_lam = simpson(integrand_hop, x=t_array)

# 修改后
I0 = 2.0 * simpson(C * np.cos(S), x=t_array)
I_lam = 2.0 * simpson(integrand_hop, x=t_array)
```

**验证**: 修复后纯论文输入验证 (t_paper + paper_G) 应得到 μ_c'(300K) 放大2倍。

---

### 3.2 E2: 近邻求和因子 (Stage 5)

**物理原理**: 论文中 $A_\alpha$ 和 $B_{\alpha q}$ 的求和遍历**所有近邻分子** ($\sum_{n \neq m}$)。萘晶体中提取了6个独立dimer方向，但每个方向有正负2个等效跳跃路径 ($+R_n$ 和 $-R_n$)。因此对所有近邻求和 = 对独立对求和的 **2倍**。

**代码位置**: `stage5_mobility_assemble.py:compute_A0_alpha()` 和 `compute_B_alpha_lambda()`

**修复**:
```python
# 修改前
A0_alpha = (1/(2*HBAR_SI**2)) * np.sum(R_proj**2 * t_values**2)
B_alpha = (omega_rads**2 / 4) * np.sum(R_proj**2 * g_nonlocal**2, axis=1)

# 修改后
A0_alpha = 2.0 * (1/(2*HBAR_SI**2)) * np.sum(R_proj**2 * t_values**2)
B_alpha = 2.0 * (omega_rads**2 / 4) * np.sum(R_proj**2 * g_nonlocal**2, axis=1)
```

**注意**: 这里乘以2是因为每个独立dimer对 (如 a方向的 +R_a 和 -R_a) 都对A₀/B有相同的贡献。此因子独立于E3的TB修正。

**验证**: 纯论文t值 + paper_G测试中，A₀和B应比v2大2倍。

---

### 3.3 E3: H_AB(k) P2₁/a对称性修正 (Stage 1)

**物理原理**: P2₁/a空间群的2₁螺旋轴意味着交叉子晶格(A→B)的跳跃不仅对应一个矢量R，还包含与子晶格相对位移τ_AB相关的对称等效路径。

**当前问题**: `stage1_transfer.py:fit_tb_parameters()` 中内联了简化的H_AB:
```python
H_AB = tac * exp(ik·R_ac) + tab * exp(ik·R_ab) + tabc * exp(ik·R_abc)
```
忽略了τ_AB相位因子。

**正确形式**: `crystal_utils.py:build_tb_H_AB()` 已经实现了正确的P2₁/a对称性形式:
```python
H_AB(k) = Σ_n t_n · [exp(ik·(R_n+τ_AB)) + exp(ik·(R_n-τ_AB))]
        = Σ_n 2·t_n · cos(k·τ_AB) · exp(ik·R_n)
```

**对|H_AB|的影响**:
- 修改前: |H_AB(k)| = |Σ t_n·exp(ik·R_n)| — 不存在k→-k对称性破缺（已验证E(k)=E(-k)），但缺少τ_AB调制
- 修改后: |H_AB(k)| = |Σ 2·t_n·cos(k·τ_AB)·exp(ik·R_n)| — 包含cos(k·τ_AB)的k-dependent调制

这改变了|H_AB(k)|的k依赖函数形式，因此拟合出的t_mn值会发生变化。

**代码位置**: `stage1_transfer.py:fit_tb_parameters()` -> `residuals()`

**修复**: 将内联的H_AB替换为调用`crystal_utils.build_tb_H_AB()`:
```python
# 修改前 (内联简化版)
H_AB = (tac * np.exp(1j * np.dot(k_cart, R_cross[0])) +
        tab * np.exp(1j * np.dot(k_cart, R_cross[1])) +
        tabc * np.exp(1j * np.dot(k_cart, R_cross[2])))

# 修改后 (使用已有正确函数)
from .crystal_utils import build_tb_H_AB
t_cross_arr = np.array([tac, tab, tabc])
H_AB = build_tb_H_AB(k_cart, R_cross, t_cross_arr, tau_AB)
```

**备选路径 (注释到代码中)**:
若±τ_AB形式不足以完全捕获P2₁/a对称性，需进一步扩展为4符号变体：
```python
# t_ab连接 (1/2, ±1/2, 0) (±τ_AB)
# H_AB += t_ab * [exp(ik·(R_ab+τ_AB)) + exp(ik·(R_ab-τ_AB)) 
#                  + exp(ik·(-R_ab+τ_AB)) + exp(ik·(-R_ab-τ_AB))]
# 注意: 备选路径仅在方案A修复后RMSE未改善时启用
```

**验证检查点**: 
- CP1: 修正后RMSE是否 ≤ 修正前的25.8 meV
- CP1b: κ(condition number)是否保持良好 (< 100)
- CP1c: TB拟合能带是否满足 E(k) = E(-k) (验证对称性)

---

### 3.4 E4: 符号连续性校验 (Stage 3)

**物理原理**: 2×2矩阵的本征值 $E_\pm(k) = \varepsilon_0 + H_{AA} \pm |H_{AB}|$ 仅依赖$|H_{AB}|$的绝对值。对于小位移 ΔQ=0.01 √amu·Å，t值的变化仅有几个百分点。符号翻转不可能是物理的——它只是最小二乘拟合中进入错误简并分支的数学artifact。

**风险分析**:
- 单个t的符号翻转 → |H_AB(k)|的k依赖改变 → 代价函数惩罚 → 通常不会发生
- 交叉子晶格t的**整体**符号翻转 → |H_AB(k)|不变 → 完全简并 → 这是真正的风险
- 同一个mode的 +ΔQ 和 -ΔQ 可能进入相反符号的分支 → 有限差商错误

**代码位置**: `stage3_eph.py:process_one_displacement()`

**修复**: 以v3平衡态拟合结果为参考进行连续性校验:
```python
def enforce_sign_continuity(t_fitted, t_equilibrium, mode_idx, verbose=True):
    """Ensure sign continuity for cross-sublattice parameters.

    For small displacements, t values should not flip sign. If sign
    changes but magnitude change is small (< 50%), it's a fitting
    artifact — correct the sign. If magnitude change is large, warn.
    """
    warnings = []
    cross_keys = ['ac', 'ab', 'abc']  # cross-sublattice only
    for key in cross_keys:
        t_eq = t_equilibrium[key]
        t_fit = t_fitted[key]

        if np.sign(t_fit) != np.sign(t_eq):
            rel_change = abs(t_fit - t_eq) / max(abs(t_eq), 1e-6)
            if abs(t_fit) < 1e-8:
                # Near-zero: sign is meaningless, keep as-is
                continue
            elif rel_change < 0.5:
                # Small change, sign flip is artifact → correct
                t_fitted[key] = np.sign(t_eq) * abs(t_fit)
                if verbose:
                    print(f"  [Mode {mode_idx}] {key}: sign flip corrected "
                          f"({t_fit*1000:.2f} → {t_fitted[key]*1000:.2f} meV, "
                          f"Δrel={rel_change*100:.0f}%)")
            else:
                # Large change — could be real physics or fitting failure
                warnings.append(
                    f"Mode {mode_idx}, {key}: sign flip with "
                    f"large Δ|t| ({rel_change*100:.0f}%), t_eq={t_eq*1000:.1f}, "
                    f"t_fit={t_fit*1000:.1f} meV")
    return t_fitted, warnings
```

**与简单符号强制的区别**: 
- ❌ 简单强制: `t = sign(t_init) * abs(t)` — 使用论文值作为符号参考，可能强制错误约定
- ✅ 连续性校验: 以**自己拟合的平衡态**为基准，仅修正小位移下的artifact，大变化时发出警告

**验证检查点**:
- CP3: 修正后检查所有88个mode的符号翻转警告数量
- CP3b: 对于发出警告的mode，检查 ∂t/∂Q 是否物理合理

---

## 4. 计算流程图 (修改后)

```
Stage 0 (复用v2) → Stage 1 [修改-E3] → Stage 2 (复用v2) → Stage 3 [修改-E4]
        ↓                                                                    ↓
    CONTCAR + EIGENVAL                                            176×EIGENVAL (v2复用)
        ↓                                                                    ↓
    build_tb_H_AB()  ← 使用P2₁/a对称性                     连续性校验 + central FD
        ↓                                                                    ↓
    t_mn(v3), RMSE(v3)                                          ∂ε₀/∂Q(v3), ∂t/∂Q(v3)
                                                                               ↓
Stage 6 [修改-E1] ← Stage 5 [修改-E2] ← Stage 4 (复用)
       ↓                    ↓                  ↓
    I₀×2, I_λ×2         A₀×2, B×2          G_λ(v3)
```

---

## 5. 代码修改清单

| 文件 | 修改内容 | 对应的错误 |
|------|----------|-----------|
| `naph_transport/stage1_transfer.py` | `fit_tb_parameters()`: H_AB改用`build_tb_H_AB()` | E3 |
| `naph_transport/stage3_eph.py` | `process_one_displacement()`: 添加符号连续性校验 | E4 |
| `naph_transport/stage5_mobility_assemble.py` | `compute_A0_alpha()`, `compute_B_alpha_lambda()`: 添加×2因子 | E2 |
| `naph_transport/stage6_mobility_integrate.py` | `compute_mobility_at_T()`: I₀和I_λ添加×2因子 | E1 |
| `naph_transport/__init__.py` | 版本号更新为 3.0.0 | — |
| `validate_pipeline.py` | 更新checkpoint引用路径和验证阈值 | — |
| `config.yaml` | 更新版本号和注释 | — |

---

## 6. 验证策略

### 6.1 修复后的验证管线

执行顺序: Stage 1 → Stage 3 Post → Stage 4 → Stage 5 → Stage 6

每个阶段完成后运行对应的检查点:

| 检查点 | 验证内容 | 通过标准 |
|--------|----------|----------|
| CP1 | TB拟合质量 | RMSE ≤ v2 (25.8 meV), κ < 100, E(k)=E(-k) |
| CP3 | 电声耦合导数 | 符号翻转警告 < 5个mode, ∂t/∂Q分布合理 |
| CP4 | 参数组装 | G_λ > 0, 单位正确 |
| CP5 | 迁移率预因子 | A₀/B与论文量级一致 (考虑×4因子) |
| CP6 | 迁移率积分 | 收敛性检查 (dt-halving Δ<5%), 量级与论文可比 |
| CP7 | 纯论文验证 | 论文t+G → μ_c'(300K) ~ 0.2-0.4 cm²/Vs (考虑×2因子) |

### 6.2 最终验证

1. **纯论文输入测试**: t_paper + paper_G (3 modes) → 验证代码正确性
2. **v3完整输入**: v3 t + v3 G (88 modes) → 获取最终μ_α(T)
3. **与v2对比**: 量化E1-E4修复对结果的影响
4. **与论文对比**: 量级、温度趋势、传输机制

---

## 7. 已知限制与未来工作

### 7.1 E3备选路径

若使用 `build_tb_H_AB()` (±τ_AB形式) 后TB拟合RMSE未改善或恶化，则需要:
1. 进行完整的P2₁/a晶体学对称性分析
2. 对每个交叉子晶格dimer对枚举所有4个符号变体
3. 必要时考虑是否需重算VASP能带 (检查DFT数据质量)

### 7.2 G值偏小问题

v2中已知v2 DFPT的G值 (尤其是低频模式) 偏小于论文经验力场值。E1-E4的修复不解决此问题。若修复后μ仍偏离论文，需考虑:
- DFPT力常数 vs 经验力场的系统性差异
- 是否需要更精确的交换关联泛函 (hybrid functional)

### 7.3 GPU使用

当前方案A全部为Python后处理，不需要GPU。若后续因E3备选路径或其他原因需要重新运行VASP:
- VASP为CPU密集型代码，GPU加速有限
- 如需提交，使用Slurm CPU分区 (32-64核)
- A100/V100仅在需要大规模并行或GPU加速的后处理时使用

---

## 8. 参考文档

| 文档 | 路径 |
|------|------|
| v2完成报告 | `archive_v2/naphthalene-transport-v2/PIPELINE_COMPLETION_REPORT.md` |
| 错误分析 | `archive_v2/naphthalene-transport-v2/docs/映射模型哈密顿量计算错误分析.md` |
| 计算流程 | `archive_v2/naphthalene-transport-v2/docs/CALCULATION_WORKFLOW.md` |
| v2设计规范 | `archive_v2/naphthalene-transport-v2/docs/superpowers/specs/2026-06-02-naphthalene-transfer-integral-coupling-design.md` |
| v2实现计划 | `archive_v2/naphthalene-transport-v2/docs/superpowers/plans/2026-06-02-naphthalene-transfer-integral-coupling-plan.md` |
| 原始论文 | `original_docs/model_hamiltonian/044506_1_online.pdf` |
| 文献阅读笔记 | `original_docs/model_hamiltonian/映射模型哈密顿量-文献阅读.md` |
