"""Human-readable proof certificates, styled like traditional math proofs.

A Bulletproofs proof is a transcript of a 3-act protocol; this module
renders it as 命题 / 证明 / 证毕 (Proposition / Proof / QED), with each
step labeled in plain language. Readable does not mean understandable:
the reader can follow the logical structure without knowing
elliptic-curve math.

The certificate is derived from the proof dict alone (the circuit is
recompiled from proof["equation"]); the witness never appears.
"""

from . import circuit as circuit_mod
from . import curve
from . import transcript as transcript_mod


def _trunc(h, head=8):
    """Truncate a hex string for display: head chars + … + last 4."""
    h = str(h)
    if len(h) <= head + 4:
        return h
    return h[:head] + "…" + h[-4:]


def _replay_challenges(n, q, eq_canonical, proof):
    """Replay the Fiat-Shamir transcript to recover y, z, x."""
    tr = transcript_mod.Transcript()
    tr.append_str("proto", "zkeq-circuit-v1")
    tr.append_int("n", n)
    tr.append_int("Q", q)
    tr.append_str("eq", eq_canonical)
    tr.append_point("AI", curve.decompress(bytes.fromhex(proof["AI"])))
    tr.append_point("AO", curve.decompress(bytes.fromhex(proof["AO"])))
    tr.append_point("S", curve.decompress(bytes.fromhex(proof["S"])))
    y = tr.challenge("y")
    z = tr.challenge("z")
    for i in (1, 3, 4, 5, 6):
        tr.append_point("T%d" % i,
                        curve.decompress(bytes.fromhex(proof["T"][str(i)])))
    x = tr.challenge("x")
    return y, z, x


def _short_eq(eq, limit=80):
    if len(eq) <= limit:
        return eq
    return eq[:limit] + "…"


def certificate(proof):
    """Build a bilingual proof certificate from a proof dict.

    Returns a JSON-serializable dict; every prose field is
    {"zh": ..., "en": ...}, math notation and values are shared.
    """
    import json
    eq = proof["equation"]
    k = proof.get("precision")
    ipa = proof["ipa"]
    rounds = len(ipa["Ls"])
    # n/q: stored in new proofs; legacy proofs recompile (may fail on
    # rational-coefficient canonical forms) or derive n from the IPA.
    n = proof.get("n")
    q = proof.get("q")
    if n is None:
        try:
            circ = circuit_mod.compile(eq, k if k is not None else "")
            n, q = circ.n, circ.q
        except ValueError:
            n, q = 2 ** rounds, None
    challenges = None
    if q is not None:
        try:
            challenges = _replay_challenges(n, q, eq, proof)
        except (ValueError, KeyError):
            challenges = None
    size = len(json.dumps(proof).encode("utf-8"))
    eq_short = _short_eq(eq)

    if k is None:
        prec_zh = "精确方程：witness 为整数"
        prec_en = "exact equation: integer witness"
    else:
        prec_zh = "精度 k=%d：witness 与精确根相差 < 1/%d" % (k, k * 10000)
        prec_en = ("precision k=%d: witness within 1/%d of the exact root"
                   % (k, k * 10000))

    T = lambda zh, en: {"zh": zh, "en": en}

    steps = [
        {
            "num": 1,
            "title": T("承诺", "Commit"),
            "math": "A_I = αH + ⟨a_L, G⟩ + ⟨a_R, H⟩,  "
                    "A_O = βH + ⟨a_O, G⟩,  "
                    "S = ρH + ⟨s_L, G⟩ + ⟨s_R, H⟩",
            "values": [
                ("A_I", _trunc(proof["AI"])),
                ("A_O", _trunc(proof["AO"])),
                ("S", _trunc(proof["S"])),
            ],
            "note": T(
                "证明者把秘密线值 a_L, a_R, a_O 封入三个 Pedersen 承诺。"
                "验证者无法从承诺反推 x。",
                "The prover seals the secret wires a_L, a_R, a_O into three "
                "Pedersen commitments. The verifier cannot recover x from them."),
        },
        {
            "num": 2,
            "title": T("挑战", "Challenge"),
            "math": "y, z ← transcript",
            "values": ([
                ("y", _trunc("%064x" % challenges[0])),
                ("z", _trunc("%064x" % challenges[1])),
            ] if challenges else []),
            "note": T(
                "双方由 Fiat-Shamir transcript 共同算出挑战值，无需交互。",
                "Both sides derive the challenges from the Fiat-Shamir "
                "transcript — no interaction needed."),
        },
        {
            "num": 3,
            "title": T("多项式承诺", "Polynomial commitments"),
            "math": "t(X) = ⟨l(X), r(X)⟩,  T_i = t_i·G + τ_i·H",
            "values": [("T_%d" % i, _trunc(proof["T"][str(i)]))
                       for i in (1, 3, 4, 5, 6)],
            "note": T(
                "证明者承诺多项式 t(X) 的系数（t₂ = 0，故无 T₂）。"
                "t(X) 编码了电路的全部乘法门约束。",
                "The prover commits to the coefficients of t(X) "
                "(t₂ = 0, so no T₂). t(X) encodes all multiplication-gate "
                "constraints of the circuit."),
        },
        {
            "num": 4,
            "title": T("求值", "Evaluate"),
            "math": "x ← transcript,  t̂ = ⟨l(x), r(x)⟩",
            "values": (([("x", _trunc("%064x" % challenges[2]))] if challenges
                        else []) + [
                ("t̂", _trunc(proof["t_hat"])),
                ("τ_x", _trunc(proof["tau_x"])),
                ("μ", _trunc(proof["mu"])),
            ]),
            "note": T(
                "在挑战点 x 处求值并公开盲化组合。τ_x, μ 隐藏了随机盲化元，"
                "t̂ 是双方后面要核对的目标值。",
                "Evaluate at challenge x and reveal the blinded combinations. "
                "t̂ is the target value both sides will check against."),
        },
        {
            "num": 5,
            "title": T("内积论证（%d 轮折半）" % rounds,
                       "Inner-product argument (%d rounds)" % rounds),
            "math": "⟨l(x), r(x)⟩ ≟ t̂",
            "values": [("L_%d" % (j + 1), _trunc(ipa["Ls"][j]))
                       for j in range(rounds)] +
                      [("R_%d" % (j + 1), _trunc(ipa["Rs"][j]))
                       for j in range(rounds)] +
                      [("a", _trunc(ipa["a"])),
                       ("b", _trunc(ipa["b"]))],
            "note": T(
                "每轮把向量长度折半，最终归约为两个数 a, b，"
                "验证者只需检查 a·b = t̂。",
                "Each round halves the vector length, reducing to scalars "
                "a, b; the verifier only checks a·b = t̂."),
        },
    ]

    checks = [
        T("g^t̂·h^τ_x = g^{δ·x²}·∏T_i^{xⁱ}（多项式恒等式）",
          "g^t̂·h^τ_x = g^{δ·x²}·∏T_i^{xⁱ} (polynomial identity)"),
        T("内积论证接受：承诺 P 打开为 (a, b)，且 a·b = t̂",
          "Inner-product argument accepts: commitment P opens to (a, b) "
          "with a·b = t̂"),
    ]

    q_txt_zh = str(q) if q is not None else "?"
    q_txt_en = str(q) if q is not None else "unknown"
    return {
        "proposition": T(
            "证明者知道 x，使得 %s 成立（%s），"
            "且 x 的值不被透露。" % (eq_short, prec_zh),
            "The prover knows x satisfying %s (%s), "
            "without revealing x." % (eq_short, prec_en)),
        "stats": {"n": n, "q": q, "rounds": rounds, "bytes": size},
        "stats_text": T(
            "电路：%d 条线，%s 个乘法门；证明 %d 字节" % (n, q_txt_zh, size),
            "Circuit: %d wires, %s multiplication gates; proof %d bytes"
            % (n, q_txt_en, size)),
        "steps": steps,
        "checks": checks,
        "qed_note": T(
            "以上各式皆成立，故验证者确信证明者知道满足方程的 x，"
            "而 x 本身在证明中从未出现。",
            "All checks hold; the verifier is convinced the prover knows "
            "an x satisfying the equation, while x itself never appears "
            "in the proof."),
    }
