# Committed inputs (debug/experimental)

Branch `exp/committed-inputs`. Extends the circuit protocol with one
publicly committed input (m = 1): the equation's variable x (wire 0,
`aL[0]`) is committed as `V = G_0^v · H^γ`, published before proving.
The prover shows: "I know (v, γ) with V = G_0^v·H^γ and f(v) = 0" —
without revealing v or γ.

## Protocol delta vs the base protocol

Notation: `v = aL[0]`, `G_0 = gvec[0]`, `H` the blinding generator.

**Prover** (`bulletproof.prove_committed`, deltas marked `DELTA`):

1. Transcript uses `"zkeq-circuit-v1-committed"` (domain separation:
   base proofs can never verify as committed and vice versa).
2. `V = v·G_0 + γ·H` is appended to the transcript **before** `A_I`
   (hence before the y, z, x challenges). Binding via Fiat-Shamir:
   the prover cannot tune v after seeing challenges.
3. `A_I = α·H + Σ_{i≥1} aL[i]·G_i + <aR, hvec>` — wire 0 excluded.
   Its `v·G_0` term is supplied by the public V instead.
4. `l(X), r(X), t(X)` are computed from the **full** witness
   (unchanged); check (a) is unaffected (it involves no A_I/μ/V).
5. `μ = α·x + β·x² + ρ·x³ + γ·x` — absorbs the extra `x·γ·H` term
   the verifier introduces via `x·V`.
6. The prover's `Pip` gets the same `+ x·V` term as the verifier's,
   so both sides prove `Pip = <lx,gvec> + <rx,hp> + t̂·U`.
7. Proof dict gains `"V"` (compressed point hex).

**Verifier** (`bulletproof.verify_detail_committed`):

1. Recomputes the transcript with V (domain-separated).
2. `Pvef' = Pvef + x·V`, where Pvef is built from the proof's A_I
   exactly as in the base protocol.
3. Checks (a) and (b) run unchanged; (b)'s IPA input is
   `Pip' = Pvef' − μ·h`.

## Why it is sound

- `x·A_I' + x·V = x·A_I(base) + x·γ·H`; the prover's `μ` absorbs the
  `x·γ·H`, so an honest prover's IPA input is byte-identical to the
  base protocol's. Completeness follows from base completeness.
- V is fixed before y, z, x. In check (b) the IPA binds `lx[0]`; its
  `G_0`-coefficient must equal `x·(v* + y⁻¹·wR[0]) + x²·aO[0] +
  x³·sL[0]` where `v*` is V's committed value. Since A_I'/A_O/S/V are
  fixed before x, coefficient comparison forces `aL[0] = v*` (DLOG
  binding; guessing x is infeasible). Check (a) then forces the wires
  to satisfy the circuit (base Schwartz–Zippel argument), so `v*`
  satisfies the equation.
- Knowledge of γ is enforced: the H-coefficient of Pvef' carries
  `x·γ*`; only someone knowing γ* can supply the matching μ'
  (DLOG in H). A prover that doesn't know the opening of V cannot
  complete check (b).
- Cross-protocol substitution is blocked by the transcript domain
  separation string.

## What is NOT covered

- This is experimental: the delta above is derived by the author,
  not independently reviewed. It is isolated on this branch precisely
  so the stable protocol on `main` is untouched.
- m = 1 only (single committed input = the equation variable).
- No range proof on v: the commitment hides v, but v is not shown
  to lie in any interval beyond what the circuit constrains.
- Side-channel note from the base protocol applies unchanged.

## API

- `api.commit_input(eq, witness, precision="")` → `{"V", "gamma"}`
- `api.prove_committed_equation(eq, witness, gamma_hex, precision="")`
- `api.verify_committed_equation(eq, V_hex, proof, precision="")`
