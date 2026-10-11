# 非线性 ODE 与 PDE 的零知识证明

## 概述

Witness 的微分方程证明支持 2×2 矩阵：

| | ODE | PDE |
|---|---|---|
| **线性** | 精确证明 | 精确证明 |
| **非线性** | 定点残差证明 | 定点残差证明 |

- **线性**：Chebyshev tau 离散化得到有理系数线性系统，精确求解后 reduce mod p，证明残差**恰好为零**。
- **非线性**：精确有限域解**可能不存在**（见下），改为证明**残差有界**：`‖F(z/S)‖ ≤ ε`。

## 为什么非线性不能精确证明

### 根本障碍：解的存在性

`x² - 2 = 0` 有实数解 `±√2`，但在 mod 3 下无解。没有任何定理把实数解（或 ODE 适定性）转移到固定素数 `p` 下的有限域解。

启发式估计：约 37% 的二次多项式系统在给定素数 `p` 下无解（Poisson(1) 模型）。如果命题本身为假，诚实 prover 也证不出来。

### 试过的死路

- **有限域 Newton**：有限域没有拓扑/度量，收敛理论不成立。
- **Hensel 提升**：方向反了——需要 mod p 解作为输入往上提，不能从实数近似得到 mod p 解。
- **复数同伦延拓**：没有 `ℂ → 𝔽_p` 的映射，帮不上忙。
- **Gröbner 基**：算得出来（F4/F5 对 4–16 个未知数是秒级），但算出来之前不知道解是否存在。

### 标准 ZK 范式的启示

标准 ZK（circom 的 `<--` hint）允许 out-of-circuit 任意计算，in-circuit 只验证。我们的 witness 生成（float Newton → 量化）符合这个范式，区别在于：标准 ZK 的 witness 按构造一定存在，我们的可能不存在。所以 statement 必须改弱。

## 定点残差证明

### Statement

> 我知道整数向量 `z`，使得对所有 tau 残差 `F_i`，有 `|F_i(z/S)| ≤ ε`。

其中 `S` 是定点 scale（如 `10⁶`），`ε` 是容差（如 `10⁻⁵`）。

### 电路构造

缩放整数残差：
```
R_i(z) = D·S²·F_i(z/S)
```
其中 `D` 清除有理系数的分母。则 `|F_i| ≤ ε ⟺ |R_i| ≤ D·S²·ε =: B`。

Range proof：`R_i + B = Σ_{t=0}^{m-1} 2ᵗ·bit_{i,t}`，其中 `2ᵐ > 2B`，每个 bit 用 Hadamard 约束 `b·(1-b) = 0`。

### 复杂度

- ODE（`deg=N`）：`(N+1)` 个系数，`O(N²)` 个乘积项，`(N+1)·m` 个 bit。
- PDE（`deg=N`）：`(N+1)²` 个系数，`O(N⁴)` 个乘积项（`u²` 的 2D 卷积），`M·m` 个 bit（`M` 为约束数）。

`N=3` 时：ODE 约 60 个乘法门，PDE 约 200 个乘法门。Bulletproofs 轻松处理。

## API

### 非线性 ODE

```python
from core import api

# y' - y + y² = 0, y(0) = 0.5, deg = 3
proof = api.prove_equation('cheb(y\' - y + y^2 = 0, y(0) = 0.5, deg = 3)', '', '')
ok = api.verify_equation('cheb(y\' - y + y^2 = 0, y(0) = 0.5, deg = 3)', proof, '')
# ok == True
```

证明的 canonical 形式包含参数：
```
cheb_fp(-*y+y*y+y'=0,y(0)=1/2,deg=3,S=1000000,eps=1/100000)
```

自定义参数：
```python
proof = api.cheb_fp_prove(eq_str, S=10**7, eps=Fraction(1, 10**6))
```

### 非线性 PDE

```python
# u_xx + u_yy + u² = 2, 零 Dirichlet, deg = 3
proof = api.prove_equation('pde2d(u_xx + u_yy + u^2 = 2, deg = 3)', '', '')
ok = api.verify_equation('pde2d(u_xx + u_yy + u^2 = 2, deg = 3)', proof, '')
```

### 线性（精确证明，不变）

```python
# 线性 ODE（Taylor）
api.prove_equation("ode(y' - y = 0, y(0) = 1, deg = 5)", '', '')

# 线性 Chebyshev ODE
api.prove_equation('cheb(y\' = y, y(0) = 1, deg = 7)', '', '')

# 线性 PDE（2D Poisson）
api.prove_equation('pde2d(u_xx + u_yy = 2, deg = 3)', '', '')
```

## 限制

- ODE 非线性：仅支持**二次**（`y²`, `y·y'` 等）；三次及以上拒绝。
- PDE 非线性：仅支持 `u²` 项；`u_xx + u_yy + u² = f` 形式。
- PDE 仅支持**奇数** `deg`（偶数因边界行奇异）。
- PDE 区域固定 `[-1,1]²`，零 Dirichlet 边界。
- ODE 区域固定 `[-1,1]`，边值点仅限 `-1, 0, 1`。

## 文件

- `core/circuit.py`：
  - `Circuit._build_cheb(..., _frac_mode)`：Fraction 约束生成
  - `Circuit._build_cheb_fp(spec, S, eps)`：ODE 定点电路
  - `Circuit._build_pde_fp(spec, S, eps)`：PDE 定点电路
- `core/api.py`：
  - `cheb_fp_prove(eq_str, S, eps)` / verify 集成
  - `pde2d_nl_fp_prove(eq_str, S, eps)` / `_verify_pde_nl_fp`
- `core/api_pde_nl.py`：2D 非线性 float Newton 求解器

## 参考

- 调研报告：`~/workspace/research_notes/finite-field-polynomial-roots-20261011-0840/report.md`
- Gröbner 基：Faugère F4 (1999), F5 (2002)；开源实现 msolve (ISSAC 2021)
- 存在性反例：`x²-2=0` mod 3；随机系统期望根数 `E[#roots]=1`
- ZKML 量化思路：arXiv:2502.18535
