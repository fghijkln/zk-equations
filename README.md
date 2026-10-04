# zk-equations

Zero-knowledge proofs for arithmetic equations, shipped as a standalone Android APK.

**The idea.** You type an equation like `x^3 + 2*x + 5 = 38` and prove *"I know an x
that satisfies this"* — the verifier learns the equation holds, but never learns x.

**How.** Bulletproofs arithmetic-circuit proofs (Bünz et al.), implemented in pure
Python (standard library only): no trusted setup, no bilinear pairings — just group
operations on secp256k1, so it runs fine on a phone. The equation is compiled into
an arithmetic circuit (multiplication gates + linear constraints), then proved.

**Status.** 🚧 Early development. Core cryptography is being implemented and reviewed.

**Honesty note.** This is a learning-grade implementation: the math is real, but the
code has not been audited. Do not use it to protect real assets.

## Layout

- `core/` — pure-Python ZK core: curve ops, Pedersen commitments, inner-product
  argument, equation compiler, prover/verifier, public API.
- `tests/` — correctness + tamper-resistance tests (`python3 -m tests.test_zk`).
- `app/` — the **Witness** Android app (Gradle project root): Chaquopy Python
  shell around `core/` (copied in by CI), WebView UI, release signed via CI.
  Built by GitHub Actions on every push to `main`.

## Equation language (v2)

Single variable `x`, integer/decimal constants, `+ - * ^` (non-negative
integer exponent), parentheses, `=`, implicit multiplication (`2x`),
`**` as `^` alias, and transcendental functions `sin`, `cos`, `exp`,
`ln`. Examples: `x^3 + 2*x + 5 = 38`, `x = 0.5`, `sin(x) = 0.5`.

Two proving modes, chosen automatically:

- **Exact.** Integer polynomial equations are proved exactly (mod the
  secp256k1 group order).
- **Tolerance.** Equations with decimals or transcendental functions
  need a precision denominator `k` (final precision `1/k`). Values use
  4 fixed decimal places; the circuit proves `|T(x) − y| < 1/k`, where
  `T` is the equation's polynomial — for transcendental functions,
  their fixed **degree-12 Taylor polynomial** (sin/cos/exp around 0,
  ln around x=1). The proved statement is exactly about that Taylor
  polynomial; for small |x| it genuinely approximates the true
  function. Accuracy guide: sin/cos for |x| ≤ 2π, exp for |x| ≤ 2,
  ln for 0.5 ≤ x ≤ 1.5 (the degree-12 ln series is already ~0.04 off
  at x=2 — the prover honestly refuses a witness that misses the
  requested precision).

Proofs are bound to the **canonical** equation (polynomial normal
form, with `k` for tolerance mode), also hashed into the Fiat–Shamir
transcript. Term-moved or side-swapped spellings of the same equation
(`x^3+2*x+5=38` vs `x^3+2*x=33` vs `38=x^3+2*x+5`) verify each other's
proofs; genuinely different equations do not.

## Quick start (core)

```bash
python3 -m core.demo
python3 -m tests.test_zk   # 60 tests: correctness + tamper resistance
```

```python
from core import api
proof = api.prove_equation("x^3 + 2*x + 5 = 38", 3)
api.verify_equation("x^3 + 2*x = 33", proof)          # True (moved term)
proof = api.prove_equation("sin(x) = 0.5", "0.5236", "1000")
api.verify_equation("sin(x) = 0.5", proof, "1000")    # True
```
