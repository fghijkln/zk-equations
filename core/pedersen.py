"""Pedersen commitments over secp256k1.

Commitment: C = v*G + r*H, with H a NUMS generator (discrete log of H
w.r.t. G unknown), hence computationally binding and perfectly hiding.
"""

from . import curve


def commit(value, blind):
    """Commit to scalar `value` with blinding `blind`."""
    return curve.add(curve.mul(value, curve.G), curve.mul(blind, curve.H()))


def commit_vec(values, gens):
    """Vector Pedersen commitment: sum(values[i]*gens[i])."""
    return curve.msm(values, gens)
