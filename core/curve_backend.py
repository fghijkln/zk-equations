"""Unified curve backend dispatcher (debug/experimental).

Selects the fastest available native backend for secp256k1 group
operations, with pure-Python fallback. All backends are byte-identical
(verified by core/test_backends.py).

Backends (in priority order):
  rust   - core/rsext  (cdylib libcurve_rs.so,  C ABI, loaded via ctypes)
  zig    - core/zigext (cdylib libcurve_zig.so, C ABI, loaded via ctypes)
  c      - core/cext   (CPython extension module curve_ext)
  python - core/curve.py (pure Python reference)

The Rust/Zig backends expose a plain C ABI (see core/nbext/ABI.md) so no
Python headers are needed -- unlike the C extension, which needs pyconfig.h.

Backend selection:
  - env WITNESS_CURVE_BACKEND=rust|zig|c|python forces one (falls back
    to python with a warning if unavailable)
  - otherwise auto: first available in priority order
  - curve_backend.set_backend(name) switches at runtime

This module is a drop-in replacement for core/curve_c.py (same names).
"""
import ctypes
import os
import warnings

from . import curve as _py

# ---------------------------------------------------------------- constants

P = _py.P
N = _py.N
G = _py.G
INF = None

_MASK32 = 0xFFFFFFFF


def _int_to_limbs(v):
    v = int(v)
    if v < 0:
        raise ValueError("negative int to limbs")
    return [(v >> (32 * i)) & _MASK32 for i in range(8)]


def _limbs_to_int(limbs):
    v = 0
    for i in range(7, -1, -1):
        v = (v << 32) | int(limbs[i])
    return v


def _point_to_jac(p):
    """(x, y) tuple / None -> 24 limbs [X, Y, Z] Jacobian."""
    if p is None:
        return [0] * 24
    x, y = p
    return _int_to_limbs(x) + _int_to_limbs(y) + _int_to_limbs(1)


def _jac_to_point(j):
    """24 limbs -> (x, y) tuple / None. Uses backend to_affine."""
    raise NotImplementedError  # filled per-backend below


# ---------------------------------------------------------------- ctypes ABI

class _CABIBackend:
    """ctypes wrapper around a cdylib implementing core/nbext/ABI.md."""

    _FN_SIGS = {
        "wcurve_add": (None, "P32_P32P32"),
        "wcurve_dbl": (None, "P32_P32"),
        "wcurve_neg": (None, "P32_P32"),
        "wcurve_is_inf": (ctypes.c_int, "P32_"),
        "wcurve_msm": (None, "P32P32P32L"),
        "wcurve_from_affine": (None, "P32_P32P32"),
        "wcurve_to_affine": (None, "P32_P32P32PI"),
        "wcurve_compress": (None, "P32_P8"),
        "wcurve_decompress": (ctypes.c_int, "P8_P32"),
    }

    def __init__(self, lib):
        self._lib = lib
        u32p = ctypes.POINTER(ctypes.c_uint32)
        u8p = ctypes.POINTER(ctypes.c_uint8)
        intp = ctypes.POINTER(ctypes.c_int)
        _sigmap = {
            "P32_P32P32": [u32p, u32p, u32p],
            "P32_P32": [u32p, u32p],
            "P32_": [u32p],
            "P32P32P32L": [u32p, u32p, u32p, ctypes.c_size_t],
            "P32_P32P32PI": [u32p, u32p, u32p, intp],
            "P32_P8": [u32p, u8p],
            "P8_P32": [u8p, u32p],
        }
        for name, (restype, sig) in self._FN_SIGS.items():
            fn = getattr(lib, name)
            fn.restype = restype
            fn.argtypes = _sigmap[sig]

    # -- low-level helpers -------------------------------------------------
    @staticmethod
    def _u32buf(limbs):
        return (ctypes.c_uint32 * len(limbs))(*limbs)

    def _jac(self, p):
        return self._u32buf(_point_to_jac(p))

    def to_affine(self, jac_buf):
        """ctypes 24-limb Jacobian buffer -> (x, y) tuple / None."""
        x = (ctypes.c_uint32 * 8)()
        y = (ctypes.c_uint32 * 8)()
        is_inf = ctypes.c_int(0)
        self._lib.wcurve_to_affine(jac_buf, x, y, ctypes.byref(is_inf))
        if is_inf.value:
            return None
        return (_limbs_to_int(x), _limbs_to_int(y))

    # -- public API (mirrors curve_c) --------------------------------------
    def is_inf(self, p):
        if p is None:
            return True
        return bool(self._lib.wcurve_is_inf(self._jac(p)))

    def neg(self, p):
        if p is None:
            return None
        a = self._jac(p)
        r = (ctypes.c_uint32 * 24)()
        self._lib.wcurve_neg(a, r)
        return self.to_affine(r)

    def add(self, p1, p2):
        a, b = self._jac(p1), self._jac(p2)
        r = (ctypes.c_uint32 * 24)()
        self._lib.wcurve_add(a, b, r)
        return self.to_affine(r)

    def dbl(self, p):
        a = self._jac(p)
        r = (ctypes.c_uint32 * 24)()
        self._lib.wcurve_dbl(a, r)
        return self.to_affine(r)

    def mul(self, k, p):
        return self.msm([k], [p])

    def msm(self, scalars, points):
        n = len(scalars)
        if n != len(points):
            raise ValueError("length mismatch")
        if n == 0:
            return None
        if n > 100000:
            raise ValueError("too many terms")
        sbuf = (ctypes.c_uint32 * (8 * n))()
        pbuf = (ctypes.c_uint32 * (24 * n))()
        for i, (s, p) in enumerate(zip(scalars, points)):
            sm = int(s) % N
            for k in range(8):
                sbuf[i * 8 + k] = (sm >> (32 * k)) & _MASK32
            if sm == 0 or p is None:
                continue  # zero scalar / inf point -> stays all-zero (inf)
            x, y = p
            for k in range(8):
                pbuf[i * 24 + k] = (int(x) >> (32 * k)) & _MASK32
                pbuf[i * 24 + 8 + k] = (int(y) >> (32 * k)) & _MASK32
            pbuf[i * 24 + 16] = 1
        r = (ctypes.c_uint32 * 24)()
        self._lib.wcurve_msm(sbuf, pbuf, r, n)
        return self.to_affine(r)

    def compress(self, p):
        if p is None:
            raise ValueError("cannot compress infinity")
        a = self._jac(p)
        out = (ctypes.c_uint8 * 33)()
        self._lib.wcurve_compress(a, out)
        return bytes(out)

    def decompress(self, data):
        data = bytes(data)
        if len(data) != 33 or data[0] not in (0x02, 0x03):
            raise ValueError("bad compressed point encoding")
        buf = (ctypes.c_uint8 * 33).from_buffer_copy(data)
        r = (ctypes.c_uint32 * 24)()
        if self._lib.wcurve_decompress(buf, r) != 0:
            raise ValueError("decompress failed (x out of range or not on curve)")
        pt = self.to_affine(r)
        if pt is None:
            raise ValueError("decompress failed")
        return pt


def _find_cdylib(names):
    """Search likely locations for a cdylib; return path or None."""
    import sys
    here = os.path.dirname(os.path.abspath(__file__))
    candidates = []
    for nm in names:
        candidates += [
            os.path.join(here, nm),
            os.path.join(here, "nbext", nm),
            os.path.join(here, "rsext", "target", "release", nm),
            os.path.join(here, "zigext", nm),
        ]
    # Chaquopy Android: top-level python dir is on sys.path
    for sp in sys.path:
        if sp and os.path.isdir(sp):
            for nm in names:
                candidates.append(os.path.join(sp, nm))
    # cwd fallback
    for nm in names:
        candidates.append(os.path.join(os.getcwd(), nm))
    for p in candidates:
        if os.path.isfile(p):
            return p
    return None


def _load_cabibackend(names):
    path = _find_cdylib(names)
    if not path:
        return None
    try:
        return _CABIBackend(ctypes.CDLL(path))
    except OSError:
        return None


# ---------------------------------------------------------------- selection

def _load_c_ext():
    try:
        import curve_ext as _c  # Chaquopy pip install (Android)
        return _c
    except ImportError:
        pass
    try:
        from .cext import curve_ext as _c  # local inplace build
        return _c
    except ImportError:
        return None


_BACKENDS = {}   # name -> module-like object
_BACKEND_ORDER = ["rust", "zig", "c", "python"]


def _init_backends():
    _BACKENDS["python"] = _py
    _c = _load_c_ext()
    if _c is not None:
        _BACKENDS["c"] = _c
    rs = _load_cabibackend(["libcurve_rs.so", "curve_rs.so"])
    if rs is not None:
        _BACKENDS["rust"] = rs
    zg = _load_cabibackend(["libcurve_zig.so", "curve_zig.so"])
    if zg is not None:
        _BACKENDS["zig"] = zg


_init_backends()


def available_backends():
    return [b for b in _BACKEND_ORDER if b in _BACKENDS]


def set_backend(name):
    """Switch active backend; returns the backend name actually selected."""
    global BACKEND, USING_C, is_inf, neg, add, dbl, mul, msm, compress, decompress
    if name not in _BACKENDS:
        raise ValueError("backend '%s' not available (have: %s)"
                         % (name, available_backends()))
    mod = _BACKENDS[name]
    BACKEND = name
    USING_C = name in ("c", "rust", "zig")  # native (non-Python) backend active
    is_inf = mod.is_inf
    neg = mod.neg
    add = mod.add
    dbl = mod.dbl
    mul = mod.mul
    msm = mod.msm
    compress = mod.compress
    decompress = mod.decompress
    return BACKEND


# default selection
_want = os.environ.get("WITNESS_CURVE_BACKEND", "").strip().lower()
if _want:
    if _want in _BACKENDS:
        set_backend(_want)
    else:
        warnings.warn("WITNESS_CURVE_BACKEND=%s unavailable, falling back" % _want)
        set_backend(available_backends()[0])
else:
    set_backend(available_backends()[0])

# NUMS generators: always Python (public, one-time)
hash_to_curve = _py.hash_to_curve
H = _py.H
U = _py.U
Gvec = _py.Gvec
Hvec = _py.Hvec
