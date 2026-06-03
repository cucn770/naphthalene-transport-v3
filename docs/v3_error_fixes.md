# Naphthalene Transport v3 — 纠错说明

> 基于 v2 的深度物理审计 (`docs/映射模型哈密顿量计算错误分析.md`)
> v3.0.0 | 2026-06-03

## 修复清单

| # | 位置 | 修复 | 影响 |
|---|------|------|------|
| E1 | Stage 6 时间积分 | 添加 ×2 因子 (偶函数) | μ 放大 2× |
| E2 | Stage 5 A₀/B | 添加 ×2 因子 (近邻求和) | μ 放大 2× |
| E3 | Stage 1 H_AB | 使用 build_tb_H_AB() P2₁/a 对称性 | t_mn 值变化 |
| E4 | Stage 3 符号 | 连续性校验 (以平衡态为参考) | ∂t/∂Q 稳定性 |

## E1+E2 净效应

迁移率系统性放大 4 倍（相对于 v2 错误结果）。

## E3 对称性改进

H_AB(k) 从 `t·exp(ik·R)` 改为 `2t·cos(k·τ)·exp(ik·R)`，
包含 P2₁/a 螺旋轴的 ±τ 等效路径。

备选路径 (若 RMSE 未改善):
  对 ab、abc 方向使用完整的 4 符号变体展开。

## 数据依赖

所有 VASP 原始数据 (EIGENVAL, FORCE_CONSTANTS, OUTCAR) 复用 v2。
仅修改 Python 后处理管线。

## 参考

- 原始论文: J. Chem. Phys. 127, 044506 (2007)
- v2 完成报告: archive_v2/PIPELINE_COMPLETION_REPORT.md
- v2 错误分析: archive_v2/docs/映射模型哈密顿量计算错误分析.md
- v3 设计规范: docs/superpowers/specs/2026-06-03-naphthalene-v3-corrections-design.md
