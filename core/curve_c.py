"""Curve API with native acceleration (debug/experimental).

Drop-in dispatcher: selects the best available backend
(rust > zig > c > python). See core/curve_backend.py.

All backends are byte-identical (verified by core/test_backends.py).
"""
from .curve_backend import (  # noqa: F401
    P, N, G, INF,
    is_inf, neg, add, dbl, mul, msm,
    compress, decompress,
    USING_C, BACKEND, available_backends, set_backend,
    hash_to_curve, H, U, Gvec, Hvec,
)
