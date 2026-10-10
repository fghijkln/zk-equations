"""Numerical equation solver: find all real roots of an equation.

Pure Python, zero dependencies (must run on-device via Chaquopy).

The equation is always a polynomial here: user polynomials are
polynomials by construction, and sin/cos/exp/ln are replaced at parse
time by their fixed Chebyshev interpolation polynomials (see
circuit.py). So "solve the equation" = "find real roots of T(x) = 0"
for the normalized polynomial T (lhs - rhs), on the circuit domain
|x| <= R_DOM.

Pipeline (standard practice, see research notes):
  1. Aberth-Ehrlich simultaneous iteration (cubic convergence) finds
     ALL roots (complex) at once -- catches multiple and clustered
     roots that bracketing can miss.
  2. Grid scan for sign changes + Brent's method (bracketing,
     superlinear, guaranteed) as an independent safety net.
  3. Newton polish on every candidate.
  4. Keep real roots inside the domain, deduplicate, sort.

Accuracy: roots are polished in float (53-bit) to ~1e-12 residual;
more than enough for the 8-decimal display and the app-level
precision check |w - r| < 1/k (k = user-supplied denominator).
"""

import cmath
import math
from fractions import Fraction

from . import circuit as circuit_mod

R_DOM = float(circuit_mod.R_DOM)


# ---------------- polynomial utilities (float) ----------------

def _to_float_coeffs(poly):
    """{exp: Fraction or (re,im)} -> [c_n, ..., c_0] floats, highest degree first.
    For complex pairs, takes the real part (use _to_complex_coeffs for full)."""
    if not poly:
        return [0.0]
    n = max(poly.keys())
    out = []
    for e in range(n, -1, -1):
        c = poly.get(e, Fraction(0))
        if isinstance(c, tuple):
            out.append(float(c[0]))
        else:
            out.append(float(c))
    return out


def _to_complex_coeffs(poly):
    """{exp: (re, im)} -> [c_n, ..., c_0] complex, highest degree first."""
    if not poly:
        return [0j]
    n = max(poly.keys())
    out = []
    for e in range(n, -1, -1):
        c = poly.get(e, (Fraction(0), Fraction(0)))
        # c may be a pair (re, im) or a plain Fraction (real)
        if isinstance(c, tuple):
            out.append(complex(float(c[0]), float(c[1])))
        else:
            out.append(complex(float(c), 0.0))
    return out


def find_complex_roots(poly):
    """All complex roots of {exp: (re,im)} poly. List of complex."""
    coeffs = _to_complex_coeffs(poly)
    n = len(coeffs) - 1
    if n <= 0:
        return []
    roots = _aberth_ehrlich(coeffs)
    # Deduplicate (Aberth can return near-duplicates for multiple roots)
    uniq = []
    for z in roots:
        if not any(abs(z - u) < 1e-6 for u in uniq):
            uniq.append(z)
    return uniq


def _horner(coeffs, x):
    """Evaluate polynomial (highest-degree-first coeffs) at x."""
    v = 0.0
    for c in coeffs:
        v = v * x + c
    return v


def _horner2(coeffs, x):
    """Evaluate p(x) and p'(x) together."""
    p = 0.0
    dp = 0.0
    for c in coeffs:
        dp = dp * x + p
        p = p * x + c
    return p, dp


def _deriv_coeffs(coeffs):
    n = len(coeffs) - 1
    return [coeffs[i] * (n - i) for i in range(len(coeffs) - 1)] or [0.0]


# ---------------- Aberth-Ehrlich (all roots, simultaneous) ----------------

def _aberth_ehrlich(coeffs, max_iter=200, tol=1e-13):
    """Find all (complex) roots via the Aberth-Ehrlich iteration.

    z_k <- z_k - (p/p') / (1 - (p/p') * sum_{j != k} 1/(z_k - z_j))
    Cubic convergence for simple roots; linear (but reliable) for
    multiple roots. Initial guesses on a small circle.
    """
    n = len(coeffs) - 1
    if n <= 0:
        return []
    if n == 1:
        return [complex(-coeffs[1] / coeffs[0])]
    # initial guesses: n points on a circle of radius 0.4
    zs = [0.4 * cmath.exp(2j * cmath.pi * k / n) for k in range(n)]
    for _ in range(max_iter):
        max_step = 0.0
        for k in range(n):
            z = zs[k]
            p, dp = _horner2(coeffs, z)
            if abs(p) < tol:
                continue
            if abs(dp) < 1e-300:
                dp = 1e-300
            ratio = p / dp
            s = 0j
            for j in range(n):
                if j != k:
                    d = z - zs[j]
                    if abs(d) < 1e-300:
                        d = 1e-300
                    s += 1.0 / d
            denom = 1.0 - ratio * s
            if abs(denom) < 1e-300:
                denom = 1e-300
            step = ratio / denom
            zs[k] = z - step
            max_step = max(max_step, abs(step))
        if max_step < tol:
            break
    return zs


# ---------------- Brent's method (bracketing) ----------------

def _brent(f, a, b, tol=1e-13, max_iter=200):
    """Brent's root finder on a bracket [a, b] with f(a)*f(b) <= 0."""
    fa, fb = f(a), f(b)
    if abs(fa) < abs(fb):
        a, b, fa, fb = b, a, fb, fa
    c, fc = a, fa
    d = e = b - a
    for _ in range(max_iter):
        if abs(fc) < abs(fb):
            a, b, c = b, c, b
            fa, fb, fc = fb, fc, fb
        tol_b = 2.0 * 2.220446049250313e-16 * abs(b) + 0.5 * tol
        m = 0.5 * (c - b)
        if abs(m) <= tol_b or fb == 0.0:
            return b
        if abs(e) >= tol_b and abs(fa) > abs(fb):
            # inverse quadratic interpolation
            s = fb / fa
            if a == c:
                p = 2.0 * m * s
                q = 1.0 - s
            else:
                q = fa / fc
                r = fb / fc
                p = s * (2.0 * m * q * (q - r) - (b - a) * (r - 1.0))
                q = (q - 1.0) * (r - 1.0) * (s - 1.0)
            if p > 0:
                q = -q
            p = abs(p)
            if 2.0 * p < min(3.0 * m * q - abs(tol_b * q), abs(e * q)):
                e, d = d, p / q
            else:
                d = e = m
        else:
            d = e = m
        a, fa = b, fb
        if abs(d) > tol_b:
            b += d
        else:
            b += tol_b if m > 0 else -tol_b
        fb = f(b)
        if (fb > 0) == (fc > 0):
            c, fc, e, d = a, fa, b - a, b - a
    return b


def _newton(coeffs, x0, tol=1e-14, max_iter=200):
    """Newton polish from x0 (float). Returns refined root or None.

    Iterates until the step is tiny (not just the residual: for
    multiple roots the residual can be small while x is still off).
    """
    x = float(x0)
    dcoeffs = _deriv_coeffs(coeffs)
    for _ in range(max_iter):
        p, dp = _horner(coeffs, x), _horner(dcoeffs, x)
        if abs(dp) < 1e-300:
            break
        step = p / dp
        x -= step
        if abs(step) <= tol * max(1.0, abs(x)) and abs(p) < 1e-9:
            break
    else:
        pass
    return x if abs(_horner(coeffs, x)) < 1e-9 else None


# ---------------- main pipeline ----------------

def _meaningful_domain(eq_str, trans_used):
    """(lo, hi) to search. Transcendental functions are Chebyshev fits on
    fixed intervals; outside them the polynomial does not represent the
    true function, so we do not search there. Pure polynomials use the
    full circuit domain |x| <= R_DOM."""
    lo, hi = -R_DOM, R_DOM
    if not trans_used:
        return lo, hi
    s = eq_str.lower()
    if "sin" in s or "cos" in s:
        lo, hi = max(lo, -math.pi), min(hi, math.pi)
    if "exp" in s:
        lo, hi = max(lo, -2.0), min(hi, 2.0)
    if "ln" in s:
        lo, hi = max(lo, 1.0), min(hi, 2.0)
    return lo, hi


def find_roots(poly, lo=-5.0, hi=5.0, grid_n=2000):
    """All real roots of {exp: Fraction} poly in [lo, hi]. Sorted floats."""
    coeffs = _to_float_coeffs(poly)
    n = len(coeffs) - 1
    if n <= 0:
        return []
    lo = max(lo, -R_DOM)
    hi = min(hi, R_DOM)

    cands = []

    # 1) Aberth-Ehrlich: all roots at once (complex)
    for z in _aberth_ehrlich(coeffs):
        if abs(z.imag) < 1e-8 * max(1.0, abs(z.real)):
            r = _newton(coeffs, z.real)
            if r is not None and lo - 1e-9 <= r <= hi + 1e-9:
                cands.append(r)

    # 2) grid scan + Brent safety net (independent of Aberth)
    f = lambda x: _horner(coeffs, x)
    xs = [lo + (hi - lo) * i / grid_n for i in range(grid_n + 1)]
    fprev = f(xs[0])
    for i in range(1, grid_n + 1):
        fc = f(xs[i])
        if fprev == 0.0:
            cands.append(xs[i - 1])
        elif fc == 0.0:
            cands.append(xs[i])
        elif (fprev < 0) != (fc < 0):
            r = _brent(f, xs[i - 1], xs[i])
            r = _newton(coeffs, r) or r
            cands.append(r)
        elif abs(fc) < abs(fprev) and i + 1 <= grid_n:
            # possible even-multiplicity root: local minimum near zero
            fn = f(xs[i + 1]) if i + 1 <= grid_n else fc
            if abs(fc) <= abs(fn) and abs(fc) < 1e-6:
                r = _newton(coeffs, xs[i])
                if r is not None and abs(f(r)) < 1e-9:
                    cands.append(r)
        fprev = fc

    # 3) dedupe (1e-6: multiple roots found as a tight cluster merge;
    # distinct roots closer than that are indistinguishable anyway)
    cands.sort()
    roots = []
    for r in cands:
        r = max(lo, min(hi, r))
        if abs(f(r)) > 1e-6:
            continue
        if not roots or abs(r - roots[-1]) > 1e-6:
            roots.append(r)
    return roots


def _fmt_complex_root(z, decimals):
    """Format a complex root as "a+bi" string with given decimals."""
    q = 10.0 ** decimals
    def rnd(x):
        v = math.floor(x * q + 0.5) / q if x >= 0 else -math.floor(-x * q + 0.5) / q
        return 0.0 if v == 0 else v
    re = rnd(z.real)
    im = rnd(z.imag)
    # Format with decimals, strip trailing zeros
    def fmt(x):
        s = "%.*f" % (decimals, x)
        if "." in s:
            s = s.rstrip("0").rstrip(".")
        return s if s not in ("", "-0") else "0"
    re_s, im_s = fmt(re), fmt(abs(im))
    if im == 0:
        return re_s
    if re == 0:
        if im == 1:
            return "i"
        if im == -1:
            return "-i"
        return ("%s%s" % ("-" if im < 0 else "", im_s + "i"))
    # both nonzero
    if im == 1:
        return "%s+i" % re_s
    if im == -1:
        return "%s-i" % re_s
    sign = "+" if im > 0 else "-"
    return "%s%s%si" % (re_s, sign, im_s)


def _solve_complex(eq_str, poly, decimals, max_roots):
    """Solve a complex polynomial; return up to max_roots "a+bi" strings."""
    roots = find_complex_roots(poly)
    out, seen = [], set()
    for z in roots[:max_roots]:
        s = _fmt_complex_root(z, decimals)
        if s not in seen:
            seen.add(s)
            out.append(s)
    return out


def solve_equation(eq_str, decimals=8, max_roots=5):
    """Solve eq_str; return up to max_roots roots as strings, rounded to
    `decimals` places (for display). Real roots as "3" or "3.14";
    complex roots as "a+bi" (e.g. "i", "1+2i"). Raises ValueError on bad input.

    Note: transcendental functions are solved as their fixed Chebyshev
    polynomials (the same polynomials the proof is about); the roots
    agree with the true functions to ~1e-8 for well-conditioned cases.
    """
    if not (1 <= decimals <= 12):
        raise ValueError("decimals must be between 1 and 12")
    poly, trans_used, _, _, _ = circuit_mod.parse_polynomial(eq_str)
    if not poly:
        raise ValueError("empty equation")
    # Use complex root finder for all; real roots have im≈0.
    # For transcendental, restrict to the meaningful domain.
    if trans_used:
        lo, hi = _meaningful_domain(eq_str, trans_used)
        # For transcendental, use the real finder (complex not meaningful)
        roots = find_roots(poly, lo, hi)
        # Convert to complex for uniform formatting
        roots = [complex(r, 0) for r in roots]
    else:
        roots = find_complex_roots(poly)
        # Sort by real then imaginary for deterministic order
        roots.sort(key=lambda z: (z.real, z.imag))
    # Check if exact (integer/Gaussian integer) for integer display
    def _is_int_coeff(c):
        if isinstance(c, tuple):
            return c[0].denominator == 1 and c[1].denominator == 1
        return c.denominator == 1
    is_exact = (not trans_used
                and all(_is_int_coeff(c) for c in poly.values()))
    out, seen = [], set()
    eps_zero = 0.5 * 10**(-decimals)
    for z in roots[:max_roots]:
        # Snap tiny values to zero (numerical noise, esp. multiple roots)
        # Use 1e-6 threshold for the snap, then format to decimals.
        re = 0.0 if abs(z.real) < 1e-6 else z.real
        im = 0.0 if abs(z.imag) < 1e-6 else z.imag
        # Format: if imaginary part is negligible, show as real.
        if im == 0.0:
            # Real root
            v = re
            q = 10.0 ** decimals
            rv = math.floor(v * q + 0.5) / q if v >= 0 else -math.floor(-v * q + 0.5) / q
            if rv == 0:
                rv = 0.0
            if is_exact and abs(rv - round(rv)) < 1e-9:
                s = str(int(round(rv)))
            else:
                # Keep trailing zeros to match expected format (e.g. "1.41421356")
                s = "%.*f" % (decimals, rv)
        else:
            s = _fmt_complex_root(complex(re, im), decimals)
        if s not in seen:
            seen.add(s)
            out.append(s)
            if len(out) >= max_roots:
                break
    return out


def all_roots(eq_str):
    """All real roots (full float precision) -- for the app-level
    precision check |w - r| < 1/k (k = user-supplied denominator)."""
    poly, trans_used, _, _, _ = circuit_mod.parse_polynomial(eq_str)
    lo, hi = _meaningful_domain(eq_str, trans_used)
    return find_roots(poly, lo, hi)
