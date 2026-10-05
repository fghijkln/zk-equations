"""Improved inner-product argument (Bulletproofs paper, Section 3).

Protocol 2 proves the paper's relation (2):
    P = g^a * h^b   and   <a, b> = c
i.e. knowledge of vectors with a given inner product, where P carries
no u^c term -- exactly as in the paper. (Relation (1),
P = g^a*h^b*u^{<a,b>}, is what the inner Protocol 1 proves.)

Reduction (the paper's Protocol 2, implemented unchanged): verifier
sends w <- Z_p^* (after P and the scalar c are fixed in the
transcript); both compute P' = P * u^{w*c}, u' = u^w; then
Protocol 1 proves relation (1):
    P' = g^a * h^b * u'^{<a,b>}.
Soundness: P' = g^a*h^b*u^{w*c} matches g^a*h^b*(u^w)^{<a,b>} iff
w*c = w*<a,b> iff c = <a,b> (w != 0). Since P and c are fixed before
w is drawn, a cheating prover cannot tune them to a lucky w
(rewinding argument: two distinct challenges force the extracted
(a, b) to coincide by DLOG binding, hence <a,b> = c and P = g^a*h^b).

Protocol 1 then proves relation (2):  P = g^a * h^b * u^<a,b>
with 2*log2(n) group elements + 2 scalars of communication.

Non-interactive via the Fiat-Shamir transcript (paper Section 4.4).
n must be a power of two (the circuit compiler pads to this).
"""

from . import curve


def _dot(a, b):
    return sum(x * y for x, y in zip(a, b)) % curve.N


def _fold_gens(gens, x, x_inv, n2):
    """g'_i = x^{-1}*g_i + x*g_{i+n'}   (paper Protocol 1, step 23)"""
    return [curve.add(curve.mul(x_inv, gens[i]), curve.mul(x, gens[n2 + i]))
            for i in range(n2)]


def _fold_gens_h(gens, x, x_inv, n2):
    """h'_i = x*h_i + x^{-1}*h_{i+n'}   (paper Protocol 1, step 25)"""
    return [curve.add(curve.mul(x, gens[i]), curve.mul(x_inv, gens[n2 + i]))
            for i in range(n2)]


def _vec_commit(scalars, gens):
    return curve.msm(scalars, gens)


def prove_relation2(g, h, u, P, a, b, tr):
    """Protocol 1 (prover): P = g^a h^b u^<a,b>. Returns proof dict."""
    n = len(a)
    assert n == len(b) == len(g) == len(h) and n >= 1 and (n & (n - 1)) == 0
    Ls, Rs = [], []
    while n > 1:
        n2 = n // 2
        a_lo, a_hi = a[:n2], a[n2:]
        b_lo, b_hi = b[:n2], b[n2:]
        cL = _dot(a_lo, b_hi)
        cR = _dot(a_hi, b_lo)
        # L = g_hi^a_lo * h_lo^b_hi * u^cL   (paper step 18)
        L = curve.add(curve.add(_vec_commit(a_lo, g[n2:]),
                                _vec_commit(b_hi, h[:n2])),
                      curve.mul(cL, u))
        # R = g_lo^a_hi * h_hi^b_lo * u^cR   (paper step 19)
        R = curve.add(curve.add(_vec_commit(a_hi, g[:n2]),
                                _vec_commit(b_lo, h[n2:])),
                      curve.mul(cR, u))
        tr.append_point("ipa-L", L)
        tr.append_point("ipa-R", R)
        Ls.append(L)
        Rs.append(R)
        x = tr.challenge("ipa-x")
        x_inv = pow(x, curve.N - 2, curve.N)
        x2 = (x * x) % curve.N
        x2_inv = (x_inv * x_inv) % curve.N
        g = _fold_gens(g, x, x_inv, n2)
        h = _fold_gens_h(h, x, x_inv, n2)
        # P' = L^{x^2} * P * R^{x^{-2}}   (paper step 26)
        P = curve.add(curve.add(curve.mul(x2, L), P), curve.mul(x2_inv, R))
        # a' = a_lo*x + a_hi*x^{-1}; b' = b_lo*x^{-1} + b_hi*x  (steps 28-29)
        a = [(a_lo[i] * x + a_hi[i] * x_inv) % curve.N for i in range(n2)]
        b = [(b_lo[i] * x_inv + b_hi[i] * x) % curve.N for i in range(n2)]
        n = n2
    # n == 1: reveal a, b (paper step 7)
    return {"Ls": Ls, "Rs": Rs, "a": a[0] % curve.N, "b": b[0] % curve.N}


def verify_relation2(g, h, u, P, proof, tr):
    """Protocol 1 (verifier)."""
    n = len(g)
    assert n == len(h) and n >= 1 and (n & (n - 1)) == 0
    Ls, Rs = proof["Ls"], proof["Rs"]
    rounds = len(Ls)
    assert rounds == len(Rs) and n == 2 ** rounds
    for L, R in zip(Ls, Rs):
        tr.append_point("ipa-L", L)
        tr.append_point("ipa-R", R)
        x = tr.challenge("ipa-x")
        x_inv = pow(x, curve.N - 2, curve.N)
        x2 = (x * x) % curve.N
        x2_inv = (x_inv * x_inv) % curve.N
        n2 = n // 2
        g = _fold_gens(g, x, x_inv, n2)
        h = _fold_gens_h(h, x, x_inv, n2)
        P = curve.add(curve.add(curve.mul(x2, L), P), curve.mul(x2_inv, R))
        n = n2
    a, b = proof["a"] % curve.N, proof["b"] % curve.N
    # P =? g^a * h^b * u^{a*b}  (paper step 9)
    expect = curve.add(curve.add(curve.mul(a, g[0]), curve.mul(b, h[0])),
                       curve.mul((a * b) % curve.N, u))
    return P == expect


def prove(g, h, u, P, c, a, b, tr):
    """Protocol 2 (prover): prove  P = g^a * h^b  and  <a, b> = c.

    Reduction: draw w, set P' = P * u^{w*c} and u' = u^w, then run
    Protocol 1 to prove  P' = g^a * h^b * u'^{<a,b>}.
    (See the module docstring for why this reduction is sound.)
    """
    tr.append_point("ip1-P", P)
    tr.append_scalar("ip1-c", c)
    w = tr.challenge("ip1-w")
    # P' = P * u^{w*c}; run Protocol 1 with blinding generator u' = u^w.
    # This is the paper's Protocol 2 reduction, implemented unchanged.
    P2 = curve.add(P, curve.mul((w * c) % curve.N, u))
    return prove_relation2(g, h, curve.mul(w, u), P2, a, b, tr)


def verify(g, h, u, P, c, proof, tr):
    """Protocol 2 (verifier)."""
    tr.append_point("ip1-P", P)
    tr.append_scalar("ip1-c", c)
    w = tr.challenge("ip1-w")
    P2 = curve.add(P, curve.mul((w * c) % curve.N, u))
    return verify_relation2(g, h, curve.mul(w, u), P2, proof, tr)


# --- (de)serialization: SEC-compressed hex for points, hex for scalars ---

def _p2hex(Pt):
    return curve.compress(Pt).hex()


def _hex2p(s):
    return curve.decompress(bytes.fromhex(s))


def _s2hex(s):
    return "%064x" % (s % curve.N)


def _hex2s(s):
    v = int(s, 16)
    if not 0 <= v < curve.N:
        raise ValueError("scalar out of range")
    return v


def to_jsonable(proof):
    return {
        "Ls": [_p2hex(L) for L in proof["Ls"]],
        "Rs": [_p2hex(R) for R in proof["Rs"]],
        "a": _s2hex(proof["a"]),
        "b": _s2hex(proof["b"]),
    }


def from_jsonable(d):
    return {
        "Ls": [_hex2p(s) for s in d["Ls"]],
        "Rs": [_hex2p(s) for s in d["Rs"]],
        "a": _hex2s(d["a"]),
        "b": _hex2s(d["b"]),
    }
