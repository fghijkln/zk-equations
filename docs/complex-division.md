# 复除法：设计与安全红线

> 基于 2026-10-10 的论文调研
> （`research_notes/complex-division-zk-circuits-20261010-1152/report.md`）。
> v0.4.3 未实现复除法；本文档为未来工作存档安全要求。

## 1. 标准做法：hint + 验证

全行业一致（gnark、arkworks、circom、0xPolygonZero）：

1. Prover 在电路外算好商，作为 witness（hint）。
2. 电路只验证 `商 × 除数 = 被除数`。

复数除法 `(a+bi)/(c+di) = (e+fi)`：
- hint `(e, f)`；
- 约束 `a = e·c − f·d`，`b = e·d + f·c`；
- **4 个乘法约束**（R1CS）/ 4 个乘法门。

显式公式 `[(ac+bd)+(bc−ad)i]/(c²+d²)` **只用于算 hint**，
不进电路（进电路要 8 个门，比直接法贵一倍）。

## 2. 安全红线：零因子（secp256k1 必读）

secp256k1 标量域 `n ≡ 1 (mod 4)`，故 `Fr[x]/(x²+1)` **不是域**，
有零因子：存在非零 `(c,d)` 使 `c²+d² = 0`（如 `(s,1)`，`s²=−1`）。

后果：
- **"除数 ≠ (0,0)" 的检查是不够的**。恶意 prover 可用零因子做除数，
  此时商不唯一（约束欠定），可伪造任意商。
- **正确做法：检查范数** `c²+d² ≠ 0`：
  hint `invN = 1/(c²+d²)`，约束 `(c²+d²)·invN = 1`（3 个门）。
- 带检查的复除法总共 **7 个门**（4 + 3）。

替代方案：换不可约定义多项式（如 `x²−β`，β 为二次非剩余），
使扩展为真正的域——这是所有配对库（gnark 等）的做法。
但那就不是 `i²=−1` 编码了。

**本项目当前状态**：v0.4.3 的多项式求值只用加乘，不用除法，
零因子问题目前不影响。**未来任何引入复除法的工作必须做范数检查**，
绝不能只查 `(c,d) ≠ (0,0)`。

## 3. 现有实现

| 库 | 复除法 | 零因子处理 |
|---|---|---|
| gnark `Ext2.DivUnchecked` | 有（hint + `x == div*y`） | 无（文档明示不管零除数） |
| gnark `Ext2.Inverse` | 有（hint + `inv*x == 1`） | 隐式（`inv*x==1` 在 x=0 时无解） |
| arkworks r1cs-std | 无 | — |
| circomlib | 无 | — |
| halo2-lib | 无 | — |

注意：gnark 的 E2 全用不可约塔多项式（`u²+5` 等），是真正的域，
天然没有零因子——它们的"不检查"在我们这里**不适用**。

## 4. 定点 rescale 说明

复数定点数除以**实数** scale 时，拆成两个独立实除法即可
（`(a+bi)/s = (a/s) + (b/s)i`），不需要复除法 machinery。
只有除以**真正的复数**除数时，才需要 §1–§2 的 7 门版本。

## 5. 参考

- gnark `e2.go`（DivUnchecked/Inverse，实现已验证）：
  https://github.com/Consensys-Incorporated/gnark/blob/master/std/algebra/emulated/fields_bls12381/e2.go
- El Housni "Pairings in Rank-1 Constraint Systems" (ACNS 2023)：
  https://eprint.iacr.org/2022/1162
- arkworks r1cs-std PR #70（`mul_by_inverse` 可靠性修复）：
  https://github.com/arkworks-rs/r1cs-std/blob/HEAD/CHANGELOG.md
