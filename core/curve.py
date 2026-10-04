"""secp256k1 group operations -- pure Python, standard library only.

Multi-scalar multiplication in Jacobian coordinates (one doubling chain
shared by all terms, single modular inverse at the end).
Scalars live in Z_n (n = curve order); field arithmetic is mod p.

Security note: this is a learning-grade implementation (no constant-time
guarantees). It is used here only for the ZK proof system, never for
signing or key exchange.
"""

import hashlib

# Field prime for secp256k1: p = 2**256 - 2**32 - 977
P = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2F
# Curve order (prime)
N = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141
# Standard generator, affine coordinates
G = (
    0x79BE667EF9DCBBAC55A06295CE870B07029BFCDB2DCE28D959F2815B16F81798,
    0x483ADA7726A3C4655DA4FBFC0E1108A8FD17B448A68554199C47D08FFB10D4B8,
)
# Point at infinity
INF = None


def _inv(a, m):
    """Modular inverse via Fermat (m prime)."""
    return pow(a, m - 2, m)


def is_inf(Pt):
    return Pt is INF


def neg(Pt):
    """Point negation: (x, y) -> (x, -y)."""
    if Pt is INF:
        return INF
    x, y = Pt
    return (x, (-y) % P)


def add(P1, P2):
    """Add two affine points (handles infinity and doubling)."""
    if P1 is INF:
        return P2
    if P2 is INF:
        return P1
    x1, y1 = P1
    x2, y2 = P2
    if x1 == x2:
        if (y1 + y2) % P == 0:
            return INF
        # doubling: s = 3*x1^2 / (2*y1)  (a = 0 for secp256k1)
        s = (3 * x1 * x1 * _inv((2 * y1) % P, P)) % P
    else:
        s = ((y2 - y1) * _inv((x2 - x1) % P, P)) % P
    x3 = (s * s - x1 - x2) % P
    y3 = (s * (x1 - x3) - y1) % P
    return (x3, y3)


def dbl(Pt):
    return add(Pt, Pt)


def mul(k, Pt):
    """Scalar multiplication via the Jacobian MSM (single inverse)."""
    return msm([k], [Pt])


def msm(scalars, points):
    """Multi-scalar multiplication: sum(scalars[i] * points[i]).

    Interleaved double-and-add in Jacobian coordinates: one doubling chain
    shared by all terms, then a single modular inverse at the end. Much
    faster than n separate affine scalar multiplications.
    """
    pairs = [(_to_jac(Pt), s % N) for s, Pt in zip(scalars, points)
             if s % N and Pt is not INF]
    if not pairs:
        return INF
    maxbits = max(s.bit_length() for _, s in pairs)
    acc = None
    for b in range(maxbits - 1, -1, -1):
        if acc is not None:
            acc = _jdbl(acc)
        for JPt, s in pairs:
            if (s >> b) & 1:
                acc = JPt if acc is None else _jadd(acc, JPt)
    return _from_jac(acc)


# --- Jacobian coordinates (X, Y, Z), infinity = None ---

def _to_jac(Pt):
    if Pt is INF:
        return None
    x, y = Pt
    return (x, y, 1)


def _from_jac(J):
    if J is None:
        return INF
    X, Y, Z = J
    zinv = _inv(Z, P)
    zinv2 = (zinv * zinv) % P
    return ((X * zinv2) % P, (Y * zinv2 % P * zinv) % P)


def _jdbl(J):
    # dbl-2009-l, a = 0
    X1, Y1, Z1 = J
    if Y1 == 0:
        return None
    XX = (X1 * X1) % P
    YY = (Y1 * Y1) % P
    YYYY = (YY * YY) % P
    S = (2 * (((X1 + YY) % P) ** 2 - XX - YYYY)) % P
    M = (3 * XX) % P
    T = (M * M - 2 * S) % P
    X3 = T
    Y3 = (M * (S - T) - 8 * YYYY) % P
    Z3 = (2 * Y1 * Z1) % P
    return (X3, Y3, Z3)


def _jadd(J1, J2):
    # add-2009-bl, a = 0
    if J1 is None:
        return J2
    if J2 is None:
        return J1
    X1, Y1, Z1 = J1
    X2, Y2, Z2 = J2
    Z1Z1 = (Z1 * Z1) % P
    Z2Z2 = (Z2 * Z2) % P
    U1 = (X1 * Z2Z2) % P
    U2 = (X2 * Z1Z1) % P
    S1 = (Y1 * Z2 % P * Z2Z2) % P
    S2 = (Y2 * Z1 % P * Z1Z1) % P
    if U1 == U2:
        if S1 != S2:
            return None
        return _jdbl(J1)
    H = (U2 - U1) % P
    I = ((2 * H) % P) ** 2 % P
    J = (H * I) % P
    r = (2 * (S2 - S1)) % P
    V = (U1 * I) % P
    X3 = (r * r - J - 2 * V) % P
    Y3 = (r * (V - X3) - 2 * S1 % P * J) % P
    Z3 = ((((Z1 + Z2) % P) ** 2 - Z1Z1 - Z2Z2) % P * H) % P
    return (X3, Y3, Z3)


def compress(Pt):
    """SEC1 compressed encoding: 33 bytes."""
    if Pt is INF:
        raise ValueError("cannot compress point at infinity")
    x, y = Pt
    return bytes([0x02 | (y & 1)]) + x.to_bytes(32, "big")


def decompress(data):
    """Decode SEC1 compressed point; raises ValueError if invalid."""
    if len(data) != 33 or data[0] not in (0x02, 0x03):
        raise ValueError("bad compressed point encoding")
    x = int.from_bytes(data[1:], "big")
    if x >= P:
        raise ValueError("x out of range")
    # y^2 = x^3 + 7 (mod p); p = 3 (mod 4) so sqrt = pow(y2, (p+1)//4)
    y2 = (pow(x, 3, P) + 7) % P
    y = pow(y2, (P + 1) // 4, P)
    if (y * y) % P != y2:
        raise ValueError("not on curve")
    if (y & 1) != (data[0] & 1):
        y = P - y
    return (x, y)


def hash_to_curve(tag):
    """Try-and-increment hash to a curve point (NUMS: discrete log vs G
    unknown to anyone). tag: bytes domain separator."""
    ctr = 0
    while True:
        x = int.from_bytes(hashlib.sha256(tag + b":" + bytes([ctr])).digest(), "big") % P
        y2 = (pow(x, 3, P) + 7) % P
        y = pow(y2, (P + 1) // 4, P)
        if (y * y) % P == y2:
            return (x, y)
        ctr += 1
        if ctr > 1000:
            raise RuntimeError("hash_to_curve failed to find point")


# --- NUMS generators (nothing-up-my-sleeve; discrete logs unknown) ---
# H  : single blinding generator for Pedersen commitments
# U  : blinding generator for the inner-product argument
# Gvec(i), Hvec(i): vector generators for the IPA / circuit protocol
_H = hash_to_curve(b"zkeq-H-v1")
_U = hash_to_curve(b"zkeq-U-v1")


def H():
    return _H


def U():
    return _U


_gvec_cache = {}
_hvec_cache = {}


def Gvec(i):
    if i not in _gvec_cache:
        _gvec_cache[i] = hash_to_curve(b"zkeq-Gvec-v1:" + str(i).encode())
    return _gvec_cache[i]


def Hvec(i):
    if i not in _hvec_cache:
        _hvec_cache[i] = hash_to_curve(b"zkeq-Hvec-v1:" + str(i).encode())
    return _hvec_cache[i]
