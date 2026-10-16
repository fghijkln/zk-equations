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

from . import curve_backend as curve
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
    # MAX-SPEED: vartime for all (leaks secrets, experiment only)
    _mv = getattr(curve, 'msm_vartime', curve.msm)
    # Batched MSM: single FFI call instead of mul+add+msm
    AI = _mv([alpha] + aL + aR, [h] + gvec + hvec)
    AO = _mv([beta] + aO, [h] + gvec)
    S = _mv([rho] + sL + sR, [h] + gvec + hvec)
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
        Ts[i] = _mv([t[i], tau], [g, h])
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
    # hp[i] = y_neg_n[i] * hvec[i]; use batched mul (single FFI)
    # Pvef = x*AI + x^2*AO + x^3*S + <x*y^{-n}∘wR, gvec> + <-1^n, hvec> + <x*wL+wO, hp>
    hp = curve.batch_mul(y_neg_n, hvec)
    # MAX-SPEED EXPERIMENT: vartime Pippenger for prove too.
    # WARNING: leaks witness via timing. NOT for production.
    _msm_v = getattr(curve, 'msm_vartime', curve.msm)
    Pvef = _msm_v(
        [x, xpow[2], xpow[3]] +
        _scale(x, _hadamard(y_neg_n, wR)) +
        [(N - 1)] * n +
        _add_vec(_scale(x, wL), wO),
        [AI, AO, S] + gvec + hvec + hp)

    # ---- inner-product argument: <lx, rx> = t_hat ----
    # Pip = Pvef + (N-mu)*h
    Pip = _mv([1, (N - mu) % N], [Pvef, h])
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


def verify(circuit, proof):
    """Verify a proof dict against the circuit. Returns True/False."""
    try:
        # MAX-SPEED: vartime for the small MSMs too (all scalars public here)
        _msm_v = getattr(curve, 'msm_vartime', curve.msm)
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
        # Batched: lhs = MSM([t_hat, tau_x], [g, h])
        lhs = _msm_v([t_hat, tau_x], [g, h])
        # rhs = MSM([delta*x^2, x^1, x^3, x^4, x^5, x^6], [g, T1, T3, T4, T5, T6])
        rhs = _msm_v(
            [(delta_yz * xpow[2]) % N, xpow[1], xpow[3], xpow[4], xpow[5], xpow[6]],
            [g, Ts[1], Ts[3], Ts[4], Ts[5], Ts[6]])
        if lhs != rhs:
            return False

        # check (b): inner-product argument
        hp = curve.batch_mul(y_neg_n, hvec)
        # Pvef batched as in prove()
        # Use vartime Pippenger: all scalars public (challenges, proof elements)
        _msm_v = getattr(curve, 'msm_vartime', curve.msm)
        Pvef = _msm_v(
            [x, xpow[2], xpow[3]] +
            _scale(x, _hadamard(y_neg_n, wR)) +
            [(N - 1)] * n +
            _add_vec(_scale(x, wL), wO),
            [AI, AO, S] + gvec + hvec + hp)
        Pip = _msm_v([1, (N - mu) % N], [Pvef, h])
        ipa_proof = ipa.from_jsonable(proof["ipa"])
        return ipa.verify(gvec, hp, curve.U(), Pip, t_hat, ipa_proof, tr)
    except (ValueError, KeyError, TypeError, AssertionError):
        return False
