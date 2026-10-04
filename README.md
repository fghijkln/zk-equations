# zk-equations

Zero-knowledge proofs for arithmetic equations, shipped as a standalone Android APK.

**The idea.** You type an equation like `x^3 + 2*x + 5 = 42` and prove *"I know an x
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

## Quick start (core)

```bash
python3 -m core.demo
```
