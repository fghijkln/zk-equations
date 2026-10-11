# 实验分支封存说明：exp/complex-v0.4.3

**封存日期**：2026-10-11
**状态**：研究完成，分支冻结，不合并，不发版。

## 分支内容

本分支在 v0.4.2 基础上，探索了**微分方程的零知识证明**，完成了 2×2 矩阵：

| | ODE | PDE |
|---|---|---|
| 线性 | 精确证明 | 精确证明 |
| 非线性 | 定点残差证明 | 定点残差证明 |

### 提交记录（2026-10-11）

- `b66d29c` 非线性 ODE 原型：parser 接受 `y^2`，电路乘法门，实数域 Newton 求解器
- `5f22795` 非线性 ODE：明确标注 prove 未实现，solver 返回近似解
- `7af176b` 2D Poisson 求解器：Chebyshev tau，零 Dirichlet，精确分数
- `6203967` 2D Poisson 完整 prove/verify：`pde2d` 电路 + API
- `ae062d0` 非线性 ODE 定点残差证明：`cheb_fp` 电路 + API
- `0669daf` CI：云端测试 workflow（pytest + pyflakes + bandit）
- `3771607` 非线性 PDE 定点证明：`pde2d_fp` 电路 + 求解器 + API
- `c13f263` 文档：非线性 ODE/PDE 定点证明设计

## 关键结论

### 1. 精确非线性证明不可能

`x² - 2 = 0` 有实数解但 mod 3 无解。没有任何定理把 ODE 适定性转移到固定素数 `p` 的有限域解。约 37% 的二次系统在给定 `p` 下无解。

**试过的死路**：有限域 Newton（无收敛理论）、Hensel（方向反了）、复数同伦（无 `ℂ→𝔽_p` 映射）、Gröbner（算得出但解可能不存在）。

### 2. 改证残差有界

Statement 改为：*我知道整数 `z`，使得 `|F_i(z/S)| ≤ ε`*。

- `S = 10⁶`（定点 scale），`ε = 10⁻⁵`（容差）
- 电路：整数缩放残差 `R_i = D·S²·F_i(z/S)`，bit 分解 range proof
- Prover：float Newton → 量化 → 验界 → 证明

### 3. 线性保持精确

线性系统的有理数解 reduce mod p 必定是解，精确证明继续用。

## API 速查

```python
# 非线性 ODE：y' - y + y² = 0
api.prove_equation('cheb(y\' - y + y^2 = 0, y(0) = 0.5, deg = 3)', '', '')

# 非线性 PDE：u_xx + u_yy + u² = 2
api.prove_equation('pde2d(u_xx + u_yy + u^2 = 2, deg = 3)', '', '')

# 线性 PDE：u_xx + u_yy = 2
api.prove_equation('pde2d(u_xx + u_yy = 2, deg = 3)', '', '')
```

## 测试状态

- 本地：`python3 -m tests.test_zk` → **116/116 通过**（2026-10-11）
- 云端：GitHub Actions `ZK Core Tests`（test + lint）→ **success**（run 38126729466）

## 未做事项

- **App 集成**：无 Java bridge，无 UI。最新 APK（witness-ode2.apk）仅含 Taylor ODE。
- **合并**：分支未合进 main，未打 tag，未发版。
- **F-Droid**：MR !51202 仍在 v0.4.2/versionCode 7 线，未动。

## 以后捡起来

1. 看本文档 + `docs/nonlinear-ode-pde.md`
2. 跑测试：`python3 -m tests.test_zk`
3. 关键文件：
   - `core/circuit.py`：`_build_cheb_fp`、`_build_pde_fp`
   - `core/api.py`：`cheb_fp_prove`、`pde2d_nl_fp_prove`
   - `core/api_pde_nl.py`：2D Newton 求解器
4. 调研原文：`~/workspace/research_notes/finite-field-polynomial-roots-20261011-0840/report.md`

## 限制（原型）

- ODE 非线性仅二次；PDE 非线性仅 `u²` 项
- PDE 仅奇数 `deg`，区域 `[-1,1]²`，零 Dirichlet
- ODE 边值点仅限 `-1, 0, 1`
