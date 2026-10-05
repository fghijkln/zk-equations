"""Zero-knowledge proof for arithmetic-circuit satisfiability.

Implements paper Section 5 (Protocol 3), with the Section 5.2
logarithmic compression: the linear-size (l, r) transfer is replaced by
the inner-product argument from Section 3.

Relation proved (paper eq. 66, with m = 0 committed inputs, i.e. no V):
    a_L o a_R = a_O                                    (Hadamard product)
    W_L*a_L + W_R*a_R + W_O*a_O = c                    (Q linear constraints)

Protocol (non-interactive via Fiat-Shamir, cf. paper Section 4.4):

  P -> V: A_I = h^a * gvec^a_L * hvec^a_R
          A_O = h^b * gvec^a_O
          S   = h^r * gvec^s_L * hvec^s_R
  V -> P: y, z
  both:   w_L = sum_q z^q W_L[q]  (likewise w_R, w_O)
          k(y,z)  = <y^{-n} o w_R, w_L>
          d(y,z)  = <z_{1:Q}, c> + k(y,z)
  P defines (paper Protocol 3, Part 1):
          l(X) = (a_L + y^{-n} o w_R)*X + a_O*X^2 + s_L*X^3
          r(X) = y^n o (a_R*X - 1^n + s_R*X^3) + w_L*X + w_O
          t(X) = <l(X), r(X)> = sum_{i=1..6} t_i X^i
  P -> V: T_i = g^{t_i} h^{tau_i}  for i in {1,3,4,5,6}   (no T_2)
  V -> P: x
  P -> V: tau_x = sum_{i!=2} tau_i x^i,  mu = a*x + b*x^2 + r*x^3,
          t_hat = <l(x), r(x)>
  both:   h'_i = hvec_i^{y^{-i}}
          Pvef = A_I^x * A_O^{x^2} * S^{x^3}
                 * gvec^{x*(y^{-n} o w_R)} * hvec^{-1^n} * h'^{x*w_L + w_O}
  V checks (paper eqs. 81-84, with 75 replaced per Section 5.2):
    (a) g^{t_hat} h^{tau_x} =?= g^{d(y,z)*x^2} * prod_{i!=2} T_i^{x^i}
    (b) inner-product argument on (gvec, h', u) proving
        Pvef * h^{-mu} = gvec^l * h'^r * u^{t_hat}, i.e. <l, r> = t_hat.

Why it works: (b) binds (l, r) to the committed wires, so
t_hat = t(x); (a) then forces t_2 = d(y,z); and
t_2 - d(y,z) = <a_L o a_R - a_O, y^n> + <z_{1:Q}, w - c>, which is the
zero polynomial in (y, z) iff the witness satisfies every constraint
(Schwartz-Zippel; the y^n factor stops a prover hiding a Hadamard
error e != 0 with <e, 1^n> = 0).
"""

import secrets

from . import curve
from . import ipa
from .transcript import Transcript

N = curve.N


def _rand():
    return secrets.randbelow(N)


def _rand_vec(n):
    return [_rand() for _ in range(n)]


def _dot(a, b):
    return sum(x * y for x, y in zip(a, b)) % N


def _vec_commit(scalars, gens):
    return curve.msm(scalars, gens)


def _add_vec(a, b):
    return [(x + y) % N for x, y in zip(a, b)]


def _hadamard(a, b):
    return [(x * y) % N for x, y in zip(a, b)]


def _scale(s, a):
    return [(s * x) % N for x in a]


def prove(circuit, aL, aR, aO):
    """Prove knowledge of wires satisfying the circuit. Returns a
    JSON-serializable proof dict."""
    n, Q = circuit.n, circuit.q
    assert len(aL) == len(aR) == len(aO) == n
    gvec = [curve.Gvec(i) for i in range(n)]
    hvec = [curve.Hvec(i) for i in range(n)]
    h = curve.H()
    g = curve.G

    tr = Transcript()
    tr.append_str("proto", "zkeq-circuit-v1")
    tr.append_int("n", n)
    tr.append_int("Q", Q)
    tr.append_str("eq", circuit.canonical)

    # ---- Part 1: commit to the wires ----
    alpha, beta, rho = _rand(), _rand(), _rand()
    sL, sR = _rand_vec(n), _rand_vec(n)
    AI = curve.add(curve.add(curve.mul(alpha, h), _vec_commit(aL, gvec)),
                             _vec_commit(aR, hvec))
    AO = curve.add(curve.mul(beta, h), _vec_commit(aO, gvec))
    S = curve.add(curve.add(curve.mul(rho, h), _vec_commit(sL, gvec)),
                  _vec_commit(sR, hvec))
    tr.append_point("AI", AI)
    tr.append_point("AO", AO)
    tr.append_point("S", S)
    y = tr.challenge("y")
    z = tr.challenge("z")

    # ---- shared values from y, z ----
    ypow = [1] * (n + 1)      # ypow[i] = y^i
    for i in range(1, n + 1):
        ypow[i] = (ypow[i - 1] * y) % N
    y_n = ypow[1:n + 1]                       # (y^1..y^n)
    y_inv = pow(y, N - 2, N)
    yinvpow = [1] * (n + 1)
    for i in range(1, n + 1):
        yinvpow[i] = (yinvpow[i - 1] * y_inv) % N
    y_neg_n = yinvpow[1:n + 1]                # (y^{-1}..y^{-n})
    zpow = [1] * (Q + 1)
    for i in range(1, Q + 1):
        zpow[i] = (zpow[i - 1] * z) % N
    z_1Q = zpow[1:Q + 1]                      # (z^1..z^Q)

    def wsum(W):
        w = [0] * n
        for q in range(Q):
            zq = z_1Q[q]
            row = W[q]
            for i in range(n):
                w[i] = (w[i] + zq * row[i]) % N
        return w

    wL, wR, wO = wsum(circuit.WL), wsum(circuit.WR), wsum(circuit.WO)
    k_yz = _dot(_hadamard(y_neg_n, wR), wL)
    delta_yz = (_dot(z_1Q, circuit.c) + k_yz) % N

    # ---- polynomials l(X), r(X), t(X) ----
    # l(X) = (aL + y^{-n} o wR) X + aO X^2 + sL X^3
    l1 = _add_vec(aL, _hadamard(y_neg_n, wR))
    l2 = list(aO)
    l3 = list(sL)
    # r(X) = y^n o (aR X - 1 + sR X^3) + wL X + wO
    r0 = _add_vec(_scale(N - 1, y_n), wO)
    r1 = _add_vec(_hadamard(y_n, aR), wL)
    r3 = _hadamard(y_n, sR)
    # t(X) = <l(X), r(X)> = sum_{i=1..6} t_i X^i
    lcoeffs = [[0] * n, l1, l2, l3]
    rcoeffs = [r0, r1, [0] * n, r3]
    t = [0] * 7
    for i in range(4):
        for j in range(4):
            if i + j >= 1:
                t[i + j] = (t[i + j] + _dot(lcoeffs[i], rcoeffs[j])) % N

    taus = {}
    Ts = {}
    for i in (1, 3, 4, 5, 6):
        tau = _rand()
        taus[i] = tau
        Ts[i] = curve.add(curve.mul(t[i], g), curve.mul(tau, h))
        tr.append_point("T%d" % i, Ts[i])
    x = tr.challenge("x")
    xpow = [1] * 7
    for i in range(1, 7):
        xpow[i] = (xpow[i - 1] * x) % N

    # ---- evaluate at x ----
    lx = [0] * n
    rx = [0] * n
    for i in range(n):
        lx[i] = (l1[i] * x + l2[i] * xpow[2] + l3[i] * xpow[3]) % N
        rx[i] = (r0[i] + r1[i] * x + r3[i] * xpow[3]) % N
    t_hat = _dot(lx, rx)
    tau_x = sum(taus[i] * xpow[i] for i in (1, 3, 4, 5, 6)) % N
    mu = (alpha * x + beta * xpow[2] + rho * xpow[3]) % N
    tr.append_scalar("tau_x", tau_x)
    tr.append_scalar("mu", mu)
    tr.append_scalar("t_hat", t_hat)

    # ---- transformed generators + verifier-computable P ----
    hp = [curve.mul(y_neg_n[i], hvec[i]) for i in range(n)]
    Pvef = curve.mul(x, AI)
    Pvef = curve.add(Pvef, curve.mul(xpow[2], AO))
    Pvef = curve.add(Pvef, curve.mul(xpow[3], S))
    Pvef = curve.add(Pvef, _vec_commit(_scale(x, _hadamard(y_neg_n, wR)), gvec))
    Pvef = curve.add(Pvef, _vec_commit([(N - 1)] * n, hvec))      # hvec^{-1^n}
    Pvef = curve.add(Pvef, _vec_commit(_add_vec(_scale(x, wL), wO), hp))

    # ---- inner-product argument: <lx, rx> = t_hat ----
    Pip = curve.add(Pvef, curve.mul((N - mu) % N, h))  # Pvef * h^{-mu}
    ipa_proof = ipa.prove(gvec, hp, curve.U(), Pip, t_hat, lx, rx, tr)

    return {
        "AI": curve.compress(AI).hex(),
        "AO": curve.compress(AO).hex(),
        "S": curve.compress(S).hex(),
        "T": {str(i): curve.compress(Ts[i]).hex() for i in (1, 3, 4, 5, 6)},
        "tau_x": "%064x" % tau_x,
        "mu": "%064x" % mu,
        "t_hat": "%064x" % t_hat,
        "ipa": ipa.to_jsonable(ipa_proof),
    }


def prove_committed(circuit, aL, aR, aO, gamma):
    """Prove with wire 0 (aL[0]) as a publicly committed input.

    Public: V = G_0^v * H^gamma (v = aL[0]). The prover knows (v, gamma)
    and the rest of the witness. Deltas vs prove() are marked DELTA.
    See docs/committed-input.md for the protocol and soundness note.
    Returns a JSON-serializable proof dict (includes "V").
    """
    n, Q = circuit.n, circuit.q
    assert len(aL) == len(aR) == len(aO) == n
    gvec = [curve.Gvec(i) for i in range(n)]
    hvec = [curve.Hvec(i) for i in range(n)]
    h = curve.H()
    g = curve.G
    v = aL[0] % N
    gamma = gamma % N

    tr = Transcript()
    tr.append_str("proto", "zkeq-circuit-v1-committed")  # DELTA: domain sep
    tr.append_int("n", n)
    tr.append_int("Q", Q)
    tr.append_str("eq", circuit.canonical)

    # DELTA: V committed BEFORE any challenge (binding via Fiat-Shamir).
    V = curve.add(curve.mul(v, gvec[0]), curve.mul(gamma, h))
    tr.append_point("V", V)

    # ---- Part 1: commit to the wires ----
    alpha, beta, rho = _rand(), _rand(), _rand()
    sL, sR = _rand_vec(n), _rand_vec(n)
    # DELTA: AI excludes wire 0 (v*G_0 comes from the public V instead).
    AI = curve.add(curve.mul(alpha, h), _vec_commit(aR, hvec))
    AI = curve.add(AI, _vec_commit(aL[1:], gvec[1:]))
    AO = curve.add(curve.mul(beta, h), _vec_commit(aO, gvec))
    S = curve.add(curve.add(curve.mul(rho, h), _vec_commit(sL, gvec)),
                  _vec_commit(sR, hvec))
    tr.append_point("AI", AI)
    tr.append_point("AO", AO)
    tr.append_point("S", S)
    y = tr.challenge("y")
    z = tr.challenge("z")

    # ---- shared values from y, z ---- (identical to prove())
    ypow = [1] * (n + 1)
    for i in range(1, n + 1):
        ypow[i] = (ypow[i - 1] * y) % N
    y_n = ypow[1:n + 1]
    y_inv = pow(y, N - 2, N)
    yinvpow = [1] * (n + 1)
    for i in range(1, n + 1):
        yinvpow[i] = (yinvpow[i - 1] * y_inv) % N
    y_neg_n = yinvpow[1:n + 1]
    zpow = [1] * (Q + 1)
    for i in range(1, Q + 1):
        zpow[i] = (zpow[i - 1] * z) % N
    z_1Q = zpow[1:Q + 1]

    def wsum(W):
        w = [0] * n
        for q in range(Q):
            zq = z_1Q[q]
            row = W[q]
            for i in range(n):
                w[i] = (w[i] + zq * row[i]) % N
        return w

    wL, wR, wO = wsum(circuit.WL), wsum(circuit.WR), wsum(circuit.WO)
    k_yz = _dot(_hadamard(y_neg_n, wR), wL)
    delta_yz = (_dot(z_1Q, circuit.c) + k_yz) % N

    # ---- polynomials l(X), r(X), t(X) ---- (identical; full aL used)
    l1 = _add_vec(aL, _hadamard(y_neg_n, wR))
    l2 = list(aO)
    l3 = list(sL)
    r0 = _add_vec(_scale(N - 1, y_n), wO)
    r1 = _add_vec(_hadamard(y_n, aR), wL)
    r3 = _hadamard(y_n, sR)
    lcoeffs = [[0] * n, l1, l2, l3]
    rcoeffs = [r0, r1, [0] * n, r3]
    t = [0] * 7
    for i in range(4):
        for j in range(4):
            if i + j >= 1:
                t[i + j] = (t[i + j] + _dot(lcoeffs[i], rcoeffs[j])) % N

    taus = {}
    Ts = {}
    for i in (1, 3, 4, 5, 6):
        tau = _rand()
        taus[i] = tau
        Ts[i] = curve.add(curve.mul(t[i], g), curve.mul(tau, h))
        tr.append_point("T%d" % i, Ts[i])
    x = tr.challenge("x")
    xpow = [1] * 7
    for i in range(1, 7):
        xpow[i] = (xpow[i - 1] * x) % N

    # ---- evaluate at x ---- (identical)
    lx = [0] * n
    rx = [0] * n
    for i in range(n):
        lx[i] = (l1[i] * x + l2[i] * xpow[2] + l3[i] * xpow[3]) % N
        rx[i] = (r0[i] + r1[i] * x + r3[i] * xpow[3]) % N
    t_hat = _dot(lx, rx)
    tau_x = sum(taus[i] * xpow[i] for i in (1, 3, 4, 5, 6)) % N
    # DELTA: mu absorbs the x*gamma*H term the verifier adds via x*V.
    mu = (alpha * x + beta * xpow[2] + rho * xpow[3] + gamma * x) % N
    tr.append_scalar("tau_x", tau_x)
    tr.append_scalar("mu", mu)
    tr.append_scalar("t_hat", t_hat)

    # ---- transformed generators + verifier-computable P ---- (identical)
    hp = [curve.mul(y_neg_n[i], hvec[i]) for i in range(n)]
    Pvef = curve.mul(x, AI)
    Pvef = curve.add(Pvef, curve.mul(xpow[2], AO))
    Pvef = curve.add(Pvef, curve.mul(xpow[3], S))
    Pvef = curve.add(Pvef, _vec_commit(_scale(x, _hadamard(y_neg_n, wR)), gvec))
    Pvef = curve.add(Pvef, _vec_commit([(N - 1)] * n, hvec))
    Pvef = curve.add(Pvef, _vec_commit(_add_vec(_scale(x, wL), wO), hp))

    # ---- inner-product argument: <lx, rx> = t_hat ---- (identical)
    Pip = curve.add(Pvef, curve.mul((N - mu) % N, h))
    # DELTA: the prover adds the same x*V term as the verifier, so both
    # sides' Pip equals <lx,gvec> + <rx,hp> + t_hat*U.
    Pip = curve.add(Pip, curve.mul(x, V))
    ipa_proof = ipa.prove(gvec, hp, curve.U(), Pip, t_hat, lx, rx, tr)

    return {
        "V": curve.compress(V).hex(),  # DELTA
        "AI": curve.compress(AI).hex(),
        "AO": curve.compress(AO).hex(),
        "S": curve.compress(S).hex(),
        "T": {str(i): curve.compress(Ts[i]).hex() for i in (1, 3, 4, 5, 6)},
        "tau_x": "%064x" % tau_x,
        "mu": "%064x" % mu,
        "t_hat": "%064x" % t_hat,
        "ipa": ipa.to_jsonable(ipa_proof),
    }
def verify_detail_committed(circuit, proof):
    """Verify a committed-input proof dict against the circuit.

    Deltas vs verify_detail() are marked DELTA. The proof must carry
    "V"; the caller must additionally check proof["V"] equals the
    commitment the verifier was given (binding the proof to V).
    Returns (valid, checks) with checks "poly"/"ipa" as usual.
    """
    checks = []
    try:
        n, Q = circuit.n, circuit.q
        gvec = [curve.Gvec(i) for i in range(n)]
        hvec = [curve.Hvec(i) for i in range(n)]
        h = curve.H()
        g = curve.G

        V = curve.decompress(bytes.fromhex(proof["V"]))  # DELTA
        AI = curve.decompress(bytes.fromhex(proof["AI"]))
        AO = curve.decompress(bytes.fromhex(proof["AO"]))
        S = curve.decompress(bytes.fromhex(proof["S"]))
        Ts = {i: curve.decompress(bytes.fromhex(proof["T"][str(i)]))
              for i in (1, 3, 4, 5, 6)}
        tau_x = int(proof["tau_x"], 16)
        mu = int(proof["mu"], 16)
        t_hat = int(proof["t_hat"], 16)
        for v in (tau_x, mu, t_hat):
            if not 0 <= v < N:
                return False, checks

        tr = Transcript()
        tr.append_str("proto", "zkeq-circuit-v1-committed")  # DELTA
        tr.append_int("n", n)
        tr.append_int("Q", Q)
        tr.append_str("eq", circuit.canonical)
        tr.append_point("V", V)  # DELTA: V before any challenge
        tr.append_point("AI", AI)
        tr.append_point("AO", AO)
        tr.append_point("S", S)
        y = tr.challenge("y")
        z = tr.challenge("z")

        y_inv = pow(y, N - 2, N)
        ypow = [1] * (n + 1)
        for i in range(1, n + 1):
            ypow[i] = (ypow[i - 1] * y) % N
        y_neg_n = [pow(y_inv, i + 1, N) for i in range(n)]
        zpow = [1] * (Q + 1)
        for i in range(1, Q + 1):
            zpow[i] = (zpow[i - 1] * z) % N
        z_1Q = zpow[1:Q + 1]

        def wsum(W):
            w = [0] * n
            for q in range(Q):
                zq = z_1Q[q]
                row = W[q]
                for i in range(n):
                    w[i] = (w[i] + zq * row[i]) % N
            return w

        wL, wR, wO = wsum(circuit.WL), wsum(circuit.WR), wsum(circuit.WO)
        k_yz = _dot(_hadamard(y_neg_n, wR), wL)
        delta_yz = (_dot(z_1Q, circuit.c) + k_yz) % N

        for i in (1, 3, 4, 5, 6):
            tr.append_point("T%d" % i, Ts[i])
        x = tr.challenge("x")
        xpow = [1] * 7
        for i in range(1, 7):
            xpow[i] = (xpow[i - 1] * x) % N

        tr.append_scalar("tau_x", tau_x)
        tr.append_scalar("mu", mu)
        tr.append_scalar("t_hat", t_hat)

        # check (a): unchanged (no AI/mu/V involved)
        lhs = curve.add(curve.mul(t_hat, g), curve.mul(tau_x, h))
        rhs = curve.mul((delta_yz * xpow[2]) % N, g)
        for i in (1, 3, 4, 5, 6):
            rhs = curve.add(rhs, curve.mul(xpow[i], Ts[i]))
        ok_poly = (lhs == rhs)
        checks.append(("poly", ok_poly))
        if not ok_poly:
            return False, checks

        # check (b): inner-product argument
        hp = [curve.mul(y_neg_n[i], hvec[i]) for i in range(n)]
        Pvef = curve.mul(x, AI)
        Pvef = curve.add(Pvef, curve.mul(xpow[2], AO))
        Pvef = curve.add(Pvef, curve.mul(xpow[3], S))
        Pvef = curve.add(Pvef, _vec_commit(_scale(x, _hadamard(y_neg_n, wR)), gvec))
        Pvef = curve.add(Pvef, _vec_commit([(N - 1)] * n, hvec))
        Pvef = curve.add(Pvef, _vec_commit(_add_vec(_scale(x, wL), wO), hp))
        # DELTA: supply the missing x*v*G_0 term from the public V.
        # (The x*gamma*H term is absorbed in the prover's mu.)
        Pvef = curve.add(Pvef, curve.mul(x, V))
        Pip = curve.add(Pvef, curve.mul((N - mu) % N, h))
        ipa_proof = ipa.from_jsonable(proof["ipa"])
        ok_ipa = bool(ipa.verify(gvec, hp, curve.U(), Pip, t_hat, ipa_proof, tr))
        checks.append(("ipa", ok_ipa))
        return ok_ipa, checks
    except (ValueError, KeyError, TypeError, AssertionError):
        return False, checks


def verify_committed(circuit, proof):
    """Verify a committed-input proof dict. Returns True/False."""
    valid, _ = verify_detail_committed(circuit, proof)
    return valid


def verify_detail(circuit, proof):
    """Verify a proof dict against the circuit.

    Returns (valid, checks) where checks is a list of
    (name, passed): "poly" (polynomial identity) and "ipa"
    (inner-product argument). A malformed proof yields
    (False, checks-so-far).
    """
    checks = []
    try:
        n, Q = circuit.n, circuit.q
        gvec = [curve.Gvec(i) for i in range(n)]
        hvec = [curve.Hvec(i) for i in range(n)]
        h = curve.H()
        g = curve.G

        AI = curve.decompress(bytes.fromhex(proof["AI"]))
        AO = curve.decompress(bytes.fromhex(proof["AO"]))
        S = curve.decompress(bytes.fromhex(proof["S"]))
        Ts = {i: curve.decompress(bytes.fromhex(proof["T"][str(i)]))
              for i in (1, 3, 4, 5, 6)}
        tau_x = int(proof["tau_x"], 16)
        mu = int(proof["mu"], 16)
        t_hat = int(proof["t_hat"], 16)
        for v in (tau_x, mu, t_hat):
            if not 0 <= v < N:
                return False

        tr = Transcript()
        tr.append_str("proto", "zkeq-circuit-v1")
        tr.append_int("n", n)
        tr.append_int("Q", Q)
        tr.append_str("eq", circuit.canonical)
        tr.append_point("AI", AI)
        tr.append_point("AO", AO)
        tr.append_point("S", S)
        y = tr.challenge("y")
        z = tr.challenge("z")

        y_inv = pow(y, N - 2, N)
        ypow = [1] * (n + 1)
        for i in range(1, n + 1):
            ypow[i] = (ypow[i - 1] * y) % N
        y_neg_n = [pow(y_inv, i + 1, N) for i in range(n)]
        zpow = [1] * (Q + 1)
        for i in range(1, Q + 1):
            zpow[i] = (zpow[i - 1] * z) % N
        z_1Q = zpow[1:Q + 1]

        def wsum(W):
            w = [0] * n
            for q in range(Q):
                zq = z_1Q[q]
                row = W[q]
                for i in range(n):
                    w[i] = (w[i] + zq * row[i]) % N
            return w

        wL, wR, wO = wsum(circuit.WL), wsum(circuit.WR), wsum(circuit.WO)
        k_yz = _dot(_hadamard(y_neg_n, wR), wL)
        delta_yz = (_dot(z_1Q, circuit.c) + k_yz) % N

        for i in (1, 3, 4, 5, 6):
            tr.append_point("T%d" % i, Ts[i])
        x = tr.challenge("x")
        xpow = [1] * 7
        for i in range(1, 7):
            xpow[i] = (xpow[i - 1] * x) % N

        tr.append_scalar("tau_x", tau_x)
        tr.append_scalar("mu", mu)
        tr.append_scalar("t_hat", t_hat)

        # check (a): g^{t_hat} h^{tau_x} =?= g^{delta*x^2} prod_{i!=2} T_i^{x^i}
        lhs = curve.add(curve.mul(t_hat, g), curve.mul(tau_x, h))
        rhs = curve.mul((delta_yz * xpow[2]) % N, g)
        for i in (1, 3, 4, 5, 6):
            rhs = curve.add(rhs, curve.mul(xpow[i], Ts[i]))
        ok_poly = (lhs == rhs)
        checks.append(("poly", ok_poly))
        if not ok_poly:
            return False, checks

        # check (b): inner-product argument
        hp = [curve.mul(y_neg_n[i], hvec[i]) for i in range(n)]
        Pvef = curve.mul(x, AI)
        Pvef = curve.add(Pvef, curve.mul(xpow[2], AO))
        Pvef = curve.add(Pvef, curve.mul(xpow[3], S))
        Pvef = curve.add(Pvef, _vec_commit(_scale(x, _hadamard(y_neg_n, wR)), gvec))
        Pvef = curve.add(Pvef, _vec_commit([(N - 1)] * n, hvec))
        Pvef = curve.add(Pvef, _vec_commit(_add_vec(_scale(x, wL), wO), hp))
        Pip = curve.add(Pvef, curve.mul((N - mu) % N, h))
        ipa_proof = ipa.from_jsonable(proof["ipa"])
        ok_ipa = bool(ipa.verify(gvec, hp, curve.U(), Pip, t_hat, ipa_proof, tr))
        checks.append(("ipa", ok_ipa))
        return ok_ipa, checks
    except (ValueError, KeyError, TypeError, AssertionError):
        return False, checks


def verify(circuit, proof):
    """Verify a proof dict against the circuit. Returns True/False."""
    valid, _ = verify_detail(circuit, proof)
    return valid
