"""Curve API with C acceleration (debug/experimental).

Tries the constant-time C extension; falls back to the pure-Python
implementation if the extension is not available. The two are
byte-identical (verified by core/cext/test_differential.py).

Import order:
  1. `import curve_ext` -- pip-installed (Chaquopy Android build).
  2. `from .cext import curve_ext` -- local `setup.py build_ext --inplace`.
  3. Pure Python fallback.

NUMS generators (H, U, Gvec, Hvec, hash_to_curve) stay in Python: they run
once at startup on public data and are not secret-dependent.
"""
from . import curve as _py

_c = None
try:
    import curve_ext as _c  # Chaquopy pip install (Android)
except ImportError:
    try:
        from .cext import curve_ext as _c  # local inplace build
    except ImportError:
        _c = None

USING_C = _c is not None

if USING_C:
    P = _py.P  # field prime (C module only exports N, G, INF)
    N = _c.N
    G = _c.G
    INF = _c.INF
    is_inf = _c.is_inf
    neg = _c.neg
    add = _c.add
    dbl = _c.dbl
    mul = _c.mul
    msm = _c.msm
    compress = _c.compress
    decompress = _c.decompress
else:
    P = _py.P
    N = _py.N
    G = _py.G
    INF = _py.INF
    is_inf = _py.is_inf
    neg = _py.neg
    add = _py.add
    dbl = _py.dbl
    mul = _py.mul
    msm = _py.msm
    compress = _py.compress
    decompress = _py.decompress

# NUMS generators: always Python (public, one-time)
hash_to_curve = _py.hash_to_curve
H = _py.H
U = _py.U
Gvec = _py.Gvec
Hvec = _py.Hvec
