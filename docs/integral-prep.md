# Witness 积分证明：开工前准备

日期 2026-10-08。用户已定方向（`x = int(x^2, 0, 1)` 形式，证书 `equation`
字段只回显符号形式），未开工；两个候选起手点（Simpson 精度测算 /
int 解析器草稿）都已备好，任选其一可立即开干。另有一个测算中发现的
诚实问题（§3），需用户拍板。

本文件未提交、未推送（repo 是公开的，设计定稿前保持 untracked）。

## 1. Simpson tolerance 测算（已跑完）

脚本 `docs/simpson_measure.py`（untracked），结论：

**编译期 Simpson 没有定点舍入误差。** 关键澄清：tolerance 模式的电路
在系数上做的是精确整数运算（`A_e = round(10^70·c_e/10^{4e})`，域内无回绕，
`bound < N//4` 已 guard），唯一的舍入是编译期系数取整与 witness 量化，
两者都远小于 1/k。Simpson 求和本身在编译期用 `Fraction` 精确算，
所以要测的只有**截断误差** `e_n = |S_n − I_true|`。

**预算分配**（电路要求 `|x − S_n| < 1/k`，严格）：
- Simpson 截断误差：`< 1/(4k)`（取 1/4 预算）
- witness 与真值距离：`< 1/(k·10^4)`（沿用 roots 的严格定义）
- witness 量化 `|X/10^4 − w| ≤ 0.5·10^{−4}`：k=10000 时这是最大头，
  占预算一半——1/(4k) 的 Simpson 规则正是为了给它留出余量
- 系数取整：`0.5/10^70` 量级，可忽略

三项相加 `< 1/(4k) + 1/(k·10^4) + 0.5·10^{−4} < 1/k`，对 k ≤ 10000 成立。

**实测**（8 个被积函数，n 为子区间数，误差精确值）：

| 被积函数 | k=10 | k=100 | k=1000 | k=10000 |
|---|---|---|---|---|
| x^2（Simpson 对 ≤3 次精确） | n=2 | n=2 | n=2 | n=2 |
| x^3+2x+5 | n=2 | n=2 | n=2 | n=2 |
| x^10 | n=4 | n=8 | n=16 | n=32 |
| sin(x)（13 次 Chebyshev） | n=4 | n=8 | n=16 | n=16 |
| exp(x)（11 次） | n=2 | n=2 | n=4 | n=4 |
| cos(x)（14 次） | n=2 | n=2 | n=4 | n=8 |
| ln(x)（9 次） | n=2 | n=2 | n=2 | n=4 |
| exp(x^2)（复合，22 次） | n=2 | n=4 | n=8 | n=16 |

规律干净：k 每 ×10，n 约 ×1.8（Simpson `O(h^4)` 的理论值），
最坏情况 k=10000 也只要 n=32。编译期耗时：n=256 的 22 次多项式
0.02 秒——手机端编译毫无压力。

**推荐规则**：不用查表，用**自适应 Richardson**——从 n=4 起每次翻倍，
`|S_{2n} − S_n|/15 < 1/(4k)` 即停（上限 n=1024，超限报"被积函数对该 k
太振荡"）。对任意用户手写的被积函数都自洽，且 prover/verifier 从同一
字符串出发得到完全相同的 S（确定性，见 §2）。

## 2. `int(f, a, b)` 解析器草稿

### 2.1 改动点（core/circuit.py）

1. **tokenize**：加 `COMMA` token（`,` 目前报 "unexpected character"）。
2. **parse_primary**：`IDENT == "int"` 时走特殊分支（在查 `_APPROX`
   之前；"int" 现在是 unknown function，无冲突）：
   ```
   int '(' expr ',' expr ',' expr ')'
   ```
   - 第一个 expr：被积函数（x 的多项式，超越函数走现有 Chebyshev 替换）
   - 第二、三个 expr：下限/上限，必须是**常数**（含 x 的项 → 报错
     "integration bounds must be numbers"；小数按 Fraction 精确处理）
   - 返回**常数多项式** `{0: S}`，S = 自适应 Simpson（§1 规则，
     `Fraction` 精确算）
3. **形式限制**：normalize（LHS−RHS）后必须是严格的 `±(x − S)`，
   即只允许 `x = int(...)` / `int(...) = x` 两种写法；其他形式
   （如 `x + int(...) = 1`）直接报错。理由：这个功能的证明对象就是
   "我知道这个积分的值"，别的写法没有保密语义，只会把数值泄进
   canonical（见 §3）。
4. **canonical（红线所在）**：tolerance 模式现在是
   `canonical = _canon_poly_frac(poly) + ";k=%d"`，对积分方程必须
   **绕过**，改走符号形式：
   ```
   x=int(<canon 被积函数>,<canon 下限>,<canon 上限>);k=100
   ```
   被积函数 canonical 复用 `_canon_poly_frac`，上下限 canonical 化小数
   写法。prover/verifier 从同一字符串得到同一 canonical，
   `verify_equation` 的 `proof["equation"] == circ.canonical` 检查照常
   通过，且证书里**永远不出现数值**。`proof["equation"]` 赋值处
   （api.py）无需改动——它本来就取 `circ.canonical`。
   实现上：parse 阶段多返回一个 `integral_spec`（被积/canonical 三件套），
   `Circuit.__init__` 收到非 None 就走符号 canonical 分支。
5. **mode**：S 通常是分数 → 自动进 tolerance 模式（现有判定逻辑不用动）；
   万一积分为整数（如 `int(x,0,2)=2`）会进 exact 模式精确证明，无害。
6. **定义域**：`|S| ≤ 5`（witness 范围 `|X| ≤ 5·10^4` 的直接推论），超限
   编译期报错；超越函数的区间诚实注记（"证明对象是 Chebyshev 多项式
   本人"）原样沿用。

### 2.2 prover 侧（core/api.py，草稿）

- `I_true`：被积多项式的**解析原函数**，`Fraction` 精确值——
  比 roots 那边还干净，不需要数值 solver。
- "查看积分值"按钮：`I_true` 格式化到 12 位小数，用户点选即 witness。
- witness 纪律（沿用严格定义）：`|w − I_true| < 1/(k·10^4)`，
  `Fraction` 精确比较，不满足直接 refuse。诚实 prover 必过：
  `|w − S_n| ≤ |w − I| + |I − S_n| < 1/(k·10^4) + 1/(4k)`，
  再加量化 `0.5·10^{−4}`，仍 `< 1/k`（§1 预算）。
- `solve.py` 不用动——积分不需要求解器。

### 2.3 电路规模

`x = S` 是 1 次多项式：`n_real = n_pow + m + m2` 很小，
证明/验证开销与现有 tolerance 方程同量级。

## 3. 一个诚实发现：Simpson 版的数值保密边界

§1–§2 的设计（Simpson 和在**编译期**算成常数 `S_n` 进电路）有一个
必须摆到台面上的性质：

**电路是公开的，`S_n` 可从电路精确恢复。**
tolerance 约束行是 `Σ a_j·w(X^j) − Σ 2^i·b_i = −(a_0 + B)`，
其中 `a_0 = round(−10^70·S_n)` 就躺在公开的 `c` 向量里；
`verify_equation` 本来就要从方程字符串**重新编译**电路，
任何一个 verifier 写一行 Python 就能算出 `S_n = −a_0/10^70`，
精确到约 70 位有效数字。

也就是说：
- ✅ 用户定的红线（证书 `equation` 字段只回显符号形式）守得住；
- ✅ 走 App 流程（UI 从不展示电路内部）实际看不到数值；
- ❌ 但"积分数值为秘密"在电路层不成立——对抗性的 verifier
  拿着源码一定能恢复数值。保密只在**输出层**成立。

用户原话是"数值只活在 witness 与电路约束里"——Simpson 版字面符合
这句话（数值确实只在 witness 和电路约束里，不在证书里）。
如果这个输出层保密模型就是想要的，§1–§2 可直接开工。

如果想要**电路层也保密**（verifier 从公开物完全推不出数值），
有另一条路线（FTC/原函数版），草稿如下：

- witness 改为原函数的**系数** `D_e = round(10^4·d_e)`（秘密 wire）
  加 `W`（秘密 wire）；公开的只有被积函数系数。
- 约束全是**线性**的：`(e+1)·D_{e+1} − C_e ≈ 0`（导数关系，
  小容差 `B_e = ceil((e+2)/2)` 覆盖取整误差）、
  `W − Σ D_e·(b^e−a^e) ≈ 0`。verifier 只能看到线性关系，
  解不出 `W`——信息论意义上的隐藏。
- 副作用全是正面的：**没有 Simpson、没有截断误差、没有 n/k 权衡**，
  原函数是精确的；且零乘法门（只有 range proof 的位约束需要乘法门）。
- 代价：① 电路 builder 要支持**多 witness 变量**（现在只支持单个 x，
  是中等改动）；② 精度固定在约 `10^{−4}` 量级（k 参数失去意义，
  对用户反而更简单）；③ 系数取整误差分析要重做（Chebyshev 高次小系数
  在 `10^4` 尺度下会直接归零，需证明自洽）。

## 4. 待用户拍板

1. **路线**：(a) Simpson 版——改动小（解析器+编译期求和+符号 canonical，
   约 100 行），保密止于输出层；还是 (b) FTC 版——数值真隐藏，
   精确无误差，但要写多变量线性电路 builder。
2. 若选 (a)：n 用自适应 Richardson（§1 推荐）还是固定查表
   `{10:4, 100:8, 1000:16, 10000:32}`；k 默认值沿用现有交互即可。
