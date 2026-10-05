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

**Known limitations.**
- *No side-channel resistance.* Blinding factors use `secrets` (CSPRNG),
  but the curve arithmetic itself is plain Python big-int code — CPython's
  integer ops, branches and allocations are not constant-time, so timing
  attacks are not defended against. This is an architectural property of
  "pure standard-library Python", not a bug to be fixed here.
- *Transcendental domain.* Chebyshev interpolation is only faithful inside
  its interval (e.g. sin/cos on [−π, π]); the solver searches |x| ≤ 5, so
  the band [π, 5] is a gray zone where the polynomial has "run away" and
  grows no roots there. Empirically (9 equations, 13 roots) nothing escapes,
  but that is measurement, not a proof — treat out-of-interval witnesses
  with suspicion.

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
  need a precision denominator `k`. Strict definition: **the witness
  must be within `1/(k·10⁴)` of an exact root** — the prover
  numerically solves the equation and refuses (no proof) if the
  witness is farther away. Witnesses support up to 8 decimal places.
  The circuit itself works in 4 fixed decimal places and proves
  `|T(x) − y| < 1/k`, where `T` is the equation's polynomial — for
  transcendental functions,
  their fixed **Chebyshev-interpolation polynomial** (near-minimax; per
  Trefethen, strictly better than Taylor at the same degree, and the
  consensus choice in the ZK literature, e.g. Kurik & Laud IACR
  2024/859). Measured max errors: sin 2.3e-9 on [−π, π] (deg 13),
  cos 2.9e-10 on [−π, π] (deg 14), exp 5.3e-9 on [−2, 2] (deg 11),
  ln 6.1e-9 on [1, 2] (deg 9). No in-circuit range reduction is done
  (it needs division/comparisons); outside the interval the proof stays
  valid but is about the polynomial, not the true function. ln is only
  meaningful on [1, 2] (singular at 0 — hard restriction).

- **Solver.** The app has a "find roots" button: pure-Python numerical
  solver (Aberth–Ehrlich simultaneous iteration + Brent bracketing +
  Newton polish, zero dependencies), shows up to 5 real roots in the
  meaningful interval, 8 decimals by default (1–12 selectable). Tap a
  root to use it as the witness. `api.solve_equation(eq, decimals=8)`.

Proofs are bound to the **canonical** equation (polynomial normal
form, with `k` for tolerance mode), also hashed into the Fiat–Shamir
transcript. Term-moved or side-swapped spellings of the same equation
(`x^3+2*x+5=38` vs `x^3+2*x=33` vs `38=x^3+2*x+5`) verify each other's
proofs; genuinely different equations do not.

## Quick start (core)

```bash
python3 -m core.demo
python3 -m tests.test_zk   # 90 tests: correctness + tamper resistance + certificate
```

```python
from core import api
proof = api.prove_equation("x^3 + 2*x + 5 = 38", 3)
api.verify_equation("x^3 + 2*x = 33", proof)          # True (moved term)
proof = api.prove_equation("sin(x) = 0.5", "0.5236", "1000")
api.verify_equation("sin(x) = 0.5", proof, "1000")    # True
```
