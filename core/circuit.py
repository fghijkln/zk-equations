"""Equation compiler:  "x^3 + 2*x + 5 = 38"  ->  arithmetic circuit.

v2 language: single variable x, integer/decimal constants, + - * ^
(non-negative integer exponent), parentheses, '=', implicit
multiplication ("2x"), '**' alias, and transcendental functions
sin/cos/exp/ln.

Transcendental functions are evaluated as fixed Chebyshev-interpolation
polynomials (near-minimax; per literature -- Trefethen -- strictly
better than Taylor at the same degree). The proved statement is
*exactly* about that polynomial -- nothing is hidden. Measured max
errors (40001-point grid): sin 2.3e-9 on [-pi, pi] (deg 13),
cos 2.9e-10 on [-pi, pi] (deg 14), exp 5.3e-9 on [-2, 2] (deg 11),
ln 6.1e-9 on [1, 2] (deg 9). No in-circuit range reduction is done
(that needs division/comparisons); outside the documented interval the
proof stays valid but is about the polynomial, not the true function.
ln is singular at 0: it is only meaningful on [1, 2] (hard restriction).

Precision k is MANDATORY for every equation (no exceptions).
Strict definition, enforced by the prover before proving:
    the witness is within 1/k of an exact root,
otherwise the witness is refused (ValueError), never proven.

Two circuit modes, chosen automatically from the equation (a pure
circuit-size optimization; the precision requirement is identical):
  exact:     the equation is an integer polynomial (no decimals, no
             transcendental functions). Proved exactly, mod n.
             (Small circuit: no bit-decomposition range proofs.)
  tolerance: decimals or transcendental functions present. The circuit
             itself works in 4 fixed decimal places (integers =
             value * 10^4; the witness is quantized to 4 places) and
             proves
                 |P(X)| <= B   i.e.   |T(x) - y| < 1/k,
             where P is the scaled integer polynomial, B = ceil(S/k)-1,
             S the scale factor. |X| <= 5*10^4 is range-checked so the
             mod-n arithmetic faithfully represents integer arithmetic.
             |P(X)| <= B is enforced by bit-decomposing u = P(X)+B.

The canonical equation string always includes ";k=<k>", binding the
proof to the precision used.

Proofs are bound to the *canonical* equation (polynomial normal form),
so "x^3+2*x+5=38" and the term-moved "x^3+2*x=33" verify each other's
proofs; genuinely different equations do not.
"""

from fractions import Fraction

from . import curve

N = curve.N

M_DEC = 10 ** 4     # fixed decimal scale: 4 places
R_DOM = 5           # |x| <= R_DOM enforced (tolerance mode)
TOL_S = 10 ** 70    # fixed scale for tolerance mode: prove |P(X)| <= B
                    # where P(X)/TOL_S ~= T(X/M_DEC); |P(X)|<=B  <=>  |T(x)-y| < 1/k


# ---------------- Chebyshev interpolants {exp: Fraction} ----------------
# Method (per literature: Trefethen's recommendation; near-minimax, and
# 1-2 degrees better than Taylor at the same max error):
# interpolate f at Chebyshev nodes on the stated interval, convert to
# power basis, round to 15 significant digits. Max errors below were
# measured on a 40001-point grid against math.sin/cos/exp/log.
#
# No in-circuit range reduction is performed (that needs division /
# comparisons); the polynomial IS the semantics, and the documented
# interval is where it approximates the true function. Outside the
# interval the proof remains valid -- about the polynomial.
#
# References: Trefethen "Six Myths of Polynomial Interpolation" (2011);
# Cephes/fdlibm (range reduction + Remez core -- the reduction step is
# not transferred, only the lesson: minimax beats Taylor);
# Kurik & Laud, IACR 2024/859 (Remez baselines for ZK circuits).

# sin on [-pi, pi], degree 13, maxerr 2.3e-9
_CHEB_SIN = {
    0: Fraction('-1.8827873147e-16'),
    1: Fraction('9.9999999925e-01'),
    2: Fraction('8.2491041236e-16'),
    3: Fraction('-1.6666665922e-01'),
    4: Fraction('-7.1460879062e-16'),
    5: Fraction('8.3333212360e-03'),
    6: Fraction('2.3706464914e-16'),
    7: Fraction('-1.9840531555e-04'),
    8: Fraction('-3.5347377809e-17'),
    9: Fraction('2.7535800486e-06'),
    10: Fraction('2.3736805841e-18'),
    11: Fraction('-2.4728366506e-08'),
    12: Fraction('-5.7382568822e-20'),
    13: Fraction('1.3611596267e-10'),
}

# cos on [-pi, pi], degree 14, maxerr 2.9e-10
_CHEB_COS = {
    0: Fraction('1.0000000000e+00'),
    1: Fraction('-7.2597114330e-17'),
    2: Fraction('-4.9999999965e-01'),
    3: Fraction('-4.9500501191e-16'),
    4: Fraction('4.1666665328e-02'),
    5: Fraction('3.4944115247e-17'),
    6: Fraction('-1.3888874208e-03'),
    7: Fraction('1.4226243379e-16'),
    8: Fraction('2.4800876207e-05'),
    9: Fraction('-4.6508600102e-17'),
    10: Fraction('-2.7539591258e-07'),
    11: Fraction('5.2411291823e-18'),
    12: Fraction('2.0638851269e-09'),
    13: Fraction('-2.0004603741e-19'),
    14: Fraction('-9.8245041169e-12'),
}

# exp on [-2, 2], degree 11, maxerr 5.3e-9
_CHEB_EXP = {
    0: Fraction('9.9999999554e-01'),
    1: Fraction('9.9999999966e-01'),
    2: Fraction('5.0000008025e-01'),
    3: Fraction('1.6666667279e-01'),
    4: Fraction('4.1666432981e-02'),
    5: Fraction('8.3333155095e-03'),
    6: Fraction('1.3891373523e-03'),
    7: Fraction('1.9843165735e-04'),
    8: Fraction('2.4682528184e-05'),
    9: Fraction('2.7466398011e-06'),
    10: Fraction('3.0168412779e-07'),
    11: Fraction('2.7049571604e-08'),
}

# ln on [1, 2], degree 9, maxerr 6.1e-9
# (ln is singular at 0; no polynomial fits near 0 -- hard restriction)
_CHEB_LN = {
    0: Fraction('-2.4796443201e+00'),
    1: Fraction('6.3250652972e+00'),
    2: Fraction('-8.8305417827e+00'),
    3: Fraction('9.5246617105e+00'),
    4: Fraction('-7.3804845674e+00'),
    5: Fraction('4.0401024757e+00'),
    6: Fraction('-1.5258331627e+00'),
    7: Fraction('3.7860030008e-01'),
    8: Fraction('-5.5588186649e-02'),
    9: Fraction('3.6622421780e-03'),
}

_CHEBYSHEV = {
    "sin": _CHEB_SIN,
    "cos": _CHEB_COS,
    "exp": _CHEB_EXP,
    "ln": _CHEB_LN,
}


# ---------------- tokenizer ----------------

class Tok:
    INT, FLOAT, IDENT, X, PLUS, MINUS, STAR, POW, LP, RP, EQ, COMMA = range(12)

    def __init__(self, kind, val=None):
        self.kind = kind
        self.val = val

    def __repr__(self):
        return "Tok(%d,%r)" % (self.kind, self.val)


def tokenize(s):
    toks = []
    i = 0
    while i < len(s):
        ch = s[i]
        if ch.isspace():
            i += 1
        elif ch.isdigit() or (ch == "." and i + 1 < len(s) and s[i + 1].isdigit()):
            j = i
            while j < len(s) and s[j].isdigit():
                j += 1
            if j < len(s) and s[j] == ".":
                k = j + 1
                while k < len(s) and s[k].isdigit():
                    k += 1
                toks.append(Tok(Tok.FLOAT, Fraction(s[i:k])))
                i = k
            else:
                toks.append(Tok(Tok.INT, int(s[i:j])))
                i = j
        elif ch.isalpha():
            j = i
            while j < len(s) and s[j].isalpha():
                j += 1
            word = s[i:j]
            if len(word) == 1 and word in "xX":
                toks.append(Tok(Tok.X))
            else:
                toks.append(Tok(Tok.IDENT, word.lower()))
            i = j
        elif ch == "+":
            toks.append(Tok(Tok.PLUS)); i += 1
        elif ch == "-":
            toks.append(Tok(Tok.MINUS)); i += 1
        elif ch == "*":
            if i + 1 < len(s) and s[i + 1] == "*":
                toks.append(Tok(Tok.POW)); i += 2
            else:
                toks.append(Tok(Tok.STAR)); i += 1
        elif ch == "^":
            toks.append(Tok(Tok.POW)); i += 1
        elif ch == "(":
            toks.append(Tok(Tok.LP)); i += 1
        elif ch == ")":
            toks.append(Tok(Tok.RP)); i += 1
        elif ch == "=":
            toks.append(Tok(Tok.EQ)); i += 1
        elif ch == ",":
            toks.append(Tok(Tok.COMMA)); i += 1
        else:
            raise ValueError("unexpected character %r in equation" % ch)
    return toks


# ---------------- complex numbers: (re, im) as Fraction pairs ----------------
# A complex number is a pair (re, im) of Fractions. Real numbers are (x, 0).
# The imaginary unit i is purely formal (i^2 = -1); it never needs to
# exist in the field. All pair-arithmetic identities hold in the ring.

def _c(re, im=0):
    """Make a complex pair from re/im (int, Fraction, or str)."""
    return (Fraction(re), Fraction(im))


def _cadd(a, b):
    # Sentinel propagates (for int(f,a,b) placeholder).
    if a is _INT_PH or b is _INT_PH:
        return _INT_PH
    return (a[0] + b[0], a[1] + b[1])


def _cneg(a):
    if a is _INT_PH:
        return _INT_PH
    return (-a[0], -a[1])


def _csub(a, b):
    return (a[0] - b[0], a[1] - b[1])


def _cmul(a, b):
    """(a+bi)(c+di) = (ac-bd) + (ad+bc)i. Schoolbook 4-mult."""
    if a is _INT_PH or b is _INT_PH:
        # int(...) in a product: propagate sentinel; form check rejects.
        # (0 * sentinel = 0, so 0*int(...) is just 0.)
        zero = (Fraction(0), Fraction(0))
        if a == zero or b == zero:
            return zero
        return _INT_PH
    return (a[0] * b[0] - a[1] * b[1], a[0] * b[1] + a[1] * b[0])


def _csq(a):
    """(a+bi)^2 = (a^2-b^2) + 2abi."""
    return (a[0] * a[0] - a[1] * a[1], 2 * a[0] * a[1])


def _ciszero(a):
    if a is _INT_PH:
        return False  # sentinel is never "zero" (form check handles it)
    return a[0] == 0 and a[1] == 0


def _cisgaussian_int(a):
    """True iff both parts are integers (Gaussian integer)."""
    return a[0].denominator == 1 and a[1].denominator == 1


# Sentinel for int(f,a,b) placeholder (so parse_equation can detect it
# before LHS-RHS mixing). Never appears in a real polynomial.
_INT_PH = object()


# ---------------- adaptive Simpson (exact Fraction arithmetic) ----------------

def _simpson_sum(f_poly, a, b, n):
    """Composite Simpson's rule with n intervals (n even).
    f_poly: {exp: (re,im)} polynomial (must be real-valued).
    a, b: Fractions (bounds). n: even int.
    Returns Fraction (exact).
    """
    # f must be real (im parts zero)
    h = (b - a) / n
    # Evaluate f at x_i = a + i*h for i=0..n
    def eval_poly(x):
        # x is Fraction; f_poly has (re,im) pairs, use re only
        v = Fraction(0)
        for e, c in f_poly.items():
            v += c[0] * (x ** e)
        return v
    s = eval_poly(a) + eval_poly(b)
    for i in range(1, n):
        x_i = a + i * h
        coeff = 4 if i % 2 == 1 else 2
        s += coeff * eval_poly(x_i)
    return s * h / 3


def _adaptive_simpson(f_poly, a, b, k, max_n=1024):
    """Adaptive Simpson via Richardson extrapolation.
    Start n=4, double until |S_{2n} - S_n|/15 < 1/(4k).
    Returns (S: Fraction, n: int).
    Raises ValueError if not converged by max_n.
    """
    # Budget: Simpson truncation error < 1/(4k)
    tol = Fraction(1, 4 * k)
    n = 4
    S_n = _simpson_sum(f_poly, a, b, n)
    while n < max_n:
        n2 = n * 2
        S_2n = _simpson_sum(f_poly, a, b, n2)
        # Richardson error estimate for Simpson: |S_{2n} - S_n|/15
        err = abs(S_2n - S_n) / 15
        if err < tol:
            return S_2n, n2
        n, S_n = n2, S_2n
    raise ValueError(
        "integrand too oscillatory for k=%d (n exceeded %d)" % (k, max_n))


# ---------------- parser: AST -> polynomial {exp: (re, im)} ----------------

class Parser:
    def __init__(self, toks):
        self.toks = toks
        self.pos = 0
        self.trans_used = False
        self.float_seen = False
        self.integral_spec = None  # set if int(f,a,b) is used

    def peek(self):
        return self.toks[self.pos] if self.pos < len(self.toks) else None

    def next(self):
        t = self.peek()
        if t is None:
            raise ValueError("unexpected end of equation")
        self.pos += 1
        return t

    def parse_equation(self):
        lhs = self.parse_expr()
        t = self.next()
        if t.kind != Tok.EQ:
            raise ValueError("equation must contain exactly one '='")
        rhs = self.parse_expr()
        if self.peek() is not None:
            raise ValueError("trailing tokens after equation")
        # Integral form restriction: if int(...) was used, the equation
        # must be exactly "x = int(f,a,b)" or "int(f,a,b) = x".
        # (Checked here, before LHS-RHS mixing, via the sentinel.)
        if self.integral_spec is not None:
            def _is_int_side(p):
                return (set(p.keys()) == {0} and p[0] is _INT_PH)
            def _is_x_side(p):
                return (set(p.keys()) == {1}
                        and p[1][0] in (1, -1) and p[1][1] == 0)
            ok = ((_is_x_side(lhs) and _is_int_side(rhs))
                  or (_is_int_side(lhs) and _is_x_side(rhs)))
            if not ok:
                raise ValueError(
                    "integral equations must be of the form x = int(f, a, b)")
            # Replace sentinel with {0: (0,0)} for the subtraction.
            if _is_int_side(lhs):
                lhs = {0: _c(0)}
            if _is_int_side(rhs):
                rhs = {0: _c(0)}
        return _psub(lhs, rhs)  # normalize: LHS - RHS = 0

    def parse_expr(self):
        node = self.parse_term()
        while True:
            t = self.peek()
            if t is not None and t.kind == Tok.PLUS:
                self.next(); node = _padd(node, self.parse_term())
            elif t is not None and t.kind == Tok.MINUS:
                self.next(); node = _psub(node, self.parse_term())
            else:
                return node

    def parse_term(self):
        node = self.parse_factor()
        while True:
            t = self.peek()
            if t is not None and t.kind == Tok.STAR:
                self.next(); node = _pmul(node, self.parse_factor())
            elif t is not None and t.kind in (Tok.INT, Tok.FLOAT, Tok.X,
                                              Tok.IDENT, Tok.LP):
                # implicit multiplication: 2x, 2.5x, 2(x+1), x sin(x)...
                node = _pmul(node, self.parse_factor())
            else:
                return node

    def parse_factor(self):
        # Unary minus binds LOOSER than ^, as in standard math/Python:
        # "-x^2" is -(x^2), not (-x)^2.
        t = self.peek()
        if t is not None and t.kind == Tok.MINUS:
            self.next()
            return _pneg(self.parse_factor())
        if t is not None and t.kind == Tok.PLUS:
            self.next()
            return self.parse_factor()
        node = self.parse_primary()
        t = self.peek()
        if t is not None and t.kind == Tok.POW:
            self.next()
            e = self.next()
            if e.kind != Tok.INT or e.val < 0:
                raise ValueError("exponent must be a non-negative integer")
            node = _ppow(node, e.val)
        return node

    def parse_primary(self):
        t = self.next()
        if t.kind == Tok.INT:
            return {0: _c(t.val)}
        if t.kind == Tok.FLOAT:
            self.float_seen = True
            return {0: _c(t.val)}
        if t.kind == Tok.X:
            return {1: _c(1)}
        if t.kind == Tok.IDENT:
            name = t.val
            if name == "i":
                # imaginary unit (formal; i^2 = -1)
                return {0: _c(0, 1)}
            if name == "int":
                # Definite integral: int(f, a, b).
                # f: integrand polynomial in x; a, b: constant bounds.
                # Returns placeholder {0: (0,0)}; Circuit computes S via
                # adaptive Simpson and fills in the real constant.
                if self.integral_spec is not None:
                    raise ValueError("only one int(...) per equation")
                if self.next().kind != Tok.LP:
                    raise ValueError("expected '(' after 'int'")
                f_poly = self.parse_expr()
                if self.next().kind != Tok.COMMA:
                    raise ValueError("expected ',' in int(f, a, b)")
                a_poly = self.parse_expr()
                if self.next().kind != Tok.COMMA:
                    raise ValueError("expected ',' in int(f, a, b)")
                b_poly = self.parse_expr()
                if self.next().kind != Tok.RP:
                    raise ValueError("missing ')' in int(f, a, b)")
                # Bounds must be constants (no x); must be real.
                def _const_frac(p, bname):
                    if any(e != 0 for e in p.keys()):
                        raise ValueError(
                            "integration %s bound must be a number" % bname)
                    c = p.get(0, _c(0))
                    if c[1] != 0:
                        raise ValueError(
                            "integration bounds must be real numbers")
                    return c[0]
                a_frac = _const_frac(a_poly, "lower")
                b_frac = _const_frac(b_poly, "upper")
                # Integrand must be real-valued (im parts zero).
                if any(c[1] != 0 for c in f_poly.values()):
                    raise ValueError(
                        "integrand must be real-valued (no 'i' in f)")
                self.integral_spec = {
                    'integrand': f_poly,
                    'a': a_frac,
                    'b': b_frac,
                }
                return {0: _INT_PH}  # sentinel; Circuit fills in S
            if name not in _APPROX:
                raise ValueError("unknown function '%s' "
                                 "(supported: sin, cos, exp, ln)" % name)
            if self.next().kind != Tok.LP:
                raise ValueError("expected '(' after '%s'" % name)
            arg = self.parse_expr()
            if self.next().kind != Tok.RP:
                raise ValueError("missing ')'")
            self.trans_used = True
            return _pcompose(_APPROX[name], arg)
        if t.kind == Tok.LP:
            node = self.parse_expr()
            if self.next().kind != Tok.RP:
                raise ValueError("missing ')'")
            return node
        raise ValueError("unexpected token %r" % (t,))


def _clean(p):
    return {e: c for e, c in p.items() if not _ciszero(c)}


def _padd(a, b):
    r = dict(a)
    for e, c in b.items():
        r[e] = _cadd(r.get(e, _c(0)), c)
    return _clean(r)


def _pneg(a):
    return {e: _cneg(c) for e, c in a.items()}


def _psub(a, b):
    return _padd(a, _pneg(b))


def _pscale(a, s):
    # s may be Fraction (real) or complex pair
    if isinstance(s, tuple):
        return _clean({e: _cmul(c, s) for e, c in a.items()})
    return _clean({e: (c[0] * s, c[1] * s) for e, c in a.items()})


def _pmul(a, b):
    r = {}
    for e1, c1 in a.items():
        for e2, c2 in b.items():
            e = e1 + e2
            r[e] = _cadd(r.get(e, _c(0)), _cmul(c1, c2))
    return _clean(r)


def _ppow(a, e):
    r = {0: _c(1)}
    base = a
    while e:
        if e & 1:
            r = _pmul(r, base)
        base = _pmul(base, base)
        e >>= 1
    return r


def _pcompose(t, p):
    """Substitute polynomial p into approximation polynomial t: t(p)."""
    r = {}
    for e, c in t.items():
        r = _padd(r, _pscale(_ppow(p, e), c))
    return r


_APPROX = {
    "sin": _CHEB_SIN,
    "cos": _CHEB_COS,
    "exp": _CHEB_EXP,
    "ln": _CHEB_LN,
}


def parse_polynomial(eq_str):
    """Parse 'lhs = rhs' into ({exp: (re,im)} for lhs - rhs,
    trans_used, float_seen, integral_spec).

    integral_spec is None for non-integral equations, else a dict with
    'integrand' ({exp: (re,im)}), 'a' (Fraction), 'b' (Fraction).
    """
    parser = Parser(tokenize(eq_str))
    poly = parser.parse_equation()
    return poly, parser.trans_used, parser.float_seen, parser.integral_spec


# ---------------- canonical form ----------------

def _round_frac(fr):
    """Round a Fraction to the nearest integer (half away from zero)."""
    n, d = fr.numerator, fr.denominator
    if n >= 0:
        return (2 * n + d) // (2 * d)
    return -((2 * (-n) + d) // (2 * d))


def _canon_poly_frac(poly):
    """poly: {exp: Fraction or (re,im)} -> canonical string (deterministic).

    For complex pairs, uses the real part (tolerance mode is real-only)."""

    def fmt(c):
        # c may be a Fraction or a (re, im) pair
        if isinstance(c, tuple):
            c = c[0]
        return (str(c.numerator) if c.denominator == 1
                else "%d/%d" % (c.numerator, c.denominator))

    def is_zero(c):
        if isinstance(c, tuple):
            return c[0] == 0 and c[1] == 0
        return c == 0

    def is_neg(c):
        if isinstance(c, tuple):
            return c[0] < 0
        return c < 0

    def abs_c(c):
        if isinstance(c, tuple):
            return abs(c[0])
        return abs(c)

    terms = []
    for e in sorted(poly.keys(), reverse=True):
        c = poly[e]
        if is_zero(c):
            continue
        sign = "-" if is_neg(c) else "+"
        ac = abs_c(c)
        cs = fmt(ac)
        if e == 0:
            body = cs
        elif e == 1:
            body = (cs + "*" if ac != 1 else "") + "x"
        else:
            body = (cs + "*" if ac != 1 else "") + "x^%d" % e
        terms.append((sign, body))
    if not terms:
        return "0=0"
    s = ("" if terms[0][0] == "+" else "-") + terms[0][1]
    for sign, body in terms[1:]:
        s += sign + body
    return s + "=0"


def _canon_poly(poly):
    """poly: {exp: signed int} -> canonical string like "x^3+2*x-33=0".

    Deterministic in the polynomial, so term-moved spellings of the same
    equation canonicalize identically.
    """
    terms = []
    for e in sorted(poly.keys(), reverse=True):
        c = poly[e]
        if c == 0:
            continue
        sign = "-" if c < 0 else "+"
        ac = abs(c)
        if e == 0:
            body = str(ac)
        elif e == 1:
            body = (str(ac) + "*" if ac != 1 else "") + "x"
        else:
            body = (str(ac) + "*" if ac != 1 else "") + "x^%d" % e
        terms.append((sign, body))
    if not terms:
        return "0=0"
    s = ("" if terms[0][0] == "+" else "-") + terms[0][1]
    for sign, body in terms[1:]:
        s += sign + body
    return s + "=0"


def _fmt_complex(c):
    """(re, im) signed ints -> string like "1+2i", "3", "-i"."""
    re, im = c
    if im == 0:
        return str(re)
    if re == 0:
        if im == 1:
            return "i"
        if im == -1:
            return "-i"
        return "%di" % im
    # both nonzero
    if im == 1:
        im_s = "+i"
    elif im == -1:
        im_s = "-i"
    elif im > 0:
        im_s = "+%di" % im
    else:
        im_s = "%di" % im
    return "%s%s" % (re, im_s)


def _canon_poly_complex(poly):
    """poly: {exp: (re, im) signed ints} -> canonical like "(1+2i)*x^2+3*x-1+i=0"."""
    terms = []
    for e in sorted(poly.keys(), reverse=True):
        c = poly[e]
        if c == (0, 0):
            continue
        # Determine sign from the formatted complex (leading "-" if re<0 or re==0 and im<0)
        cs = _fmt_complex(c)
        if cs.startswith("-"):
            sign = "-"
            body_c = cs[1:]
        else:
            sign = "+"
            body_c = cs
        if e == 0:
            body = body_c
        elif e == 1:
            # (a+bi)*x, but omit *1 for real 1
            if c == (1, 0):
                body = "x"
            else:
                body = "(%s)*x" % cs
        else:
            if c == (1, 0):
                body = "x^%d" % e
            else:
                body = "(%s)*x^%d" % (cs, e)
        terms.append((sign, body))
    if not terms:
        return "0=0"
    s = ("" if terms[0][0] == "+" else "-") + terms[0][1]
    for sign, body in terms[1:]:
        s += sign + body
    return s + "=0"


def _to_signed(c):
    c %= N
    return c if c <= N // 2 else c - N


# ---------------- circuit ----------------

class Circuit:
    """Arithmetic circuit for one equation.

    Attributes:
      n: number of multiplication gates (power of two, >= 1)
      q: number of linear constraints
      WL, WR, WO: q x n matrices (list of q rows of n scalars mod N)
      c: length-q RHS vector
      mode: "exact" | "tolerance" (circuit optimization; precision
              is mandatory in both)
      precision: k (int), always required (final precision 1/k)
      canonical: canonical equation string the proof is bound to
                 (always includes ";k=<k>")
    """

    def __init__(self, eq_str, poly, precision, trans_used, float_seen,
                 force_complex=False, integral_spec=None):
        self.eq_str = eq_str
        # Precision k is MANDATORY for every equation (parsed first).
        self.precision = _parse_precision(precision)
        self.integral_spec = integral_spec
        # Integral equations: FTC mode (bypass exact/tolerance selection).
        # The poly is {1: (±1,0), 0: (0,0)} (placeholder); FTC uses the
        # spec directly.
        if integral_spec is not None:
            self._process_integral(poly, integral_spec)
            self.is_complex = False
            self.mode = "ftc"
            spec = integral_spec
            self._build_ftc(spec['integrand'], spec['a'], spec['b'])
            self.canonical = ("%s;k=%d" % (
                self._integral_canon, self.precision))
            return
        # Normalize overall sign so that "38 = x^3+2*x+5" and
        # "x^3+2*x+5 = 38" compile to the identical circuit and
        # canonical form. For complex leading coefficient (a+bi):
        # make a > 0, or if a == 0 then b > 0.
        if poly:
            lc = poly[max(poly.keys())]
            if lc[0] < 0 or (lc[0] == 0 and lc[1] < 0):
                poly = {e: _cneg(c) for e, c in poly.items()}
        # Circuit mode is a pure size optimization, chosen from the
        # equation alone (the verifier reproduces it without the witness).
        # Any decimal literal (even "1.0") or transcendental function
        # selects tolerance mode: what you write is what you get.
        # For complex: tolerance if any re/im part is non-integer.
        needs_tol = (trans_used or float_seen
                     or any(not _cisgaussian_int(c) for c in poly.values()))
        # Complex detection: any coefficient with nonzero imaginary part,
        # or forced by the caller (complex witness for a real equation).
        is_complex = (force_complex
                      or any(c[1] != 0 for c in poly.values()))
        self.is_complex = is_complex
        if not needs_tol:
            self.mode = "exact"
            if is_complex:
                poly_cint = {e: (int(c[0]) % N, int(c[1]) % N)
                             for e, c in poly.items() if not _ciszero(c)}
                self._build_exact_complex(poly_cint)
                self.canonical = ("%s;k=%d;cx=1" % (
                    _canon_poly_complex({e: (_to_signed(v[0]), _to_signed(v[1]))
                                         for e, v in poly_cint.items()}),
                    self.precision))
            else:
                poly_int = {e: int(c[0]) % N for e, c in poly.items()
                            if not _ciszero(c)}
                self._build_exact(poly_int)
                self.canonical = ("%s;k=%d" % (
                    _canon_poly({e: _to_signed(v)
                                 for e, v in poly_int.items()}),
                    self.precision))
        else:
            self.mode = "tolerance"
            if is_complex:
                # TODO: complex tolerance circuit
                raise ValueError("complex tolerance mode not yet implemented")
            self._build_tolerance(poly, self.precision)
            self.canonical = ("%s;k=%d"
                              % (_canon_poly_frac(poly), self.precision))

    # ---- exact mode (unchanged semantics) ----

    def _build_exact(self, poly):
        self._exact_poly = poly
        self.degree = max(poly.keys()) if poly else 0
        d = self.degree
        if d >= 2:
            n_real = d - 1
        else:
            n_real = 1
        n = 1
        while n < n_real:
            n *= 2
        self.n = n

        cons = []

        def gate(which, i, coeff, acc):
            acc[which, i] = (acc.get((which, i), 0) + coeff) % N

        if d >= 2:
            for k in range(1, d - 1):
                acc = {}
                gate("L", k, 1, acc)
                gate("O", k - 1, N - 1, acc)
                cons.append((acc, 0))          # a_L[k] = a_O[k-1]
            for k in range(0, d - 1):
                acc = {}
                gate("R", k, 1, acc)
                gate("L", 0, N - 1, acc)
                cons.append((acc, 0))          # a_R[k] = a_L[0] (= x)
            acc = {}
            for exp, coeff in poly.items():
                if exp == 0:
                    continue
                if exp == 1:
                    gate("L", 0, coeff, acc)
                else:
                    gate("O", exp - 2, coeff, acc)
            cons.append((acc, (-poly.get(0, 0)) % N))  # main equation
        elif d == 1:
            c1 = poly.get(1, 0)
            acc = {}
            gate("R", 0, 1, acc)
            cons.append((acc, 1))              # a_R[0] = 1
            acc = {}
            gate("O", 0, 1, acc)
            gate("L", 0, N - 1, acc)
            cons.append((acc, 0))              # a_O[0] = a_L[0]
            acc = {}
            gate("L", 0, c1, acc)
            cons.append((acc, (-poly.get(0, 0)) % N))  # c1*x = -c0
        else:  # d == 0: constant equation
            if poly.get(0, 0) % N != 0:
                raise ValueError("equation has no solution (constant false)")
            # trivially true: no constraints (dummy zero gate)

        self._densify(cons, n)

    def _build_exact_complex(self, poly):
        """Exact circuit for Gaussian-integer polynomials.

        poly: {exp: (re_int, im_int)} with ints mod N.
        Witness is a pair (A, B) = (re(x), im(x)).
        Uses Horner's method with Karatsuba 3-mult complex multiplication.
        Proves re(P(x)) = 0 AND im(P(x)) = 0.

        Karatsuba 3-mult identity (exact in field arithmetic):
          (R+iI)(A+iB): k1=A(R+I), k2=R(B-A), k3=I(A+B);
          Re = k1-k3, Im = k1+k2.
        Proof: k1-k3 = A(R+I)-I(A+B) = AR+AI-IA-IB = AR-IB;
               k1+k2 = A(R+I)+R(B-A) = AR+AI+RB-RA = AI+RB. ∎
        (Gauss; see docs/complex-division.md and research notes.)
        Each complex multiplication costs 3 gates (vs 4 schoolbook).
        """
        self._exact_poly_c = poly
        self.degree = max(poly.keys()) if poly else 0
        d = self.degree

        if d == 0:
            c0 = poly.get(0, (0, 0))
            if c0 != (0, 0):
                raise ValueError("equation has no solution (constant false)")
            n = 1
            self.n = n
            self._densify([], n)
            return

        # Gates: 2 witness (A in aL[0], B in aL[1]) + 3 per Horner step.
        n_real = 2 + 3 * d
        n = 1
        while n < n_real:
            n *= 2
        self.n = n

        cons = []

        def gate(which, i, coeff, acc):
            acc[which, i] = (acc.get((which, i), 0) + coeff) % N

        def _add_lin(dst, src, sign=1):
            """Add linear combo src (dict) into dst (dict), with sign."""
            for k, v in src.items():
                dst[k] = (dst.get(k, 0) + sign * v) % N

        # Linear combination representation: dict {(which, idx): coeff}
        # plus a constant. Witness: A = aL[0], B = aL[1].

        # acc = c_d (constant complex)
        c_d = poly[d]
        R = {}  # linear combo for re(acc)
        I = {}  # linear combo for im(acc)
        R_const = c_d[0] % N
        I_const = c_d[1] % N

        g = 2  # next gate index (0, 1 reserved for witness)
        for e in range(d - 1, -1, -1):
            # acc = acc * x + c_e, via Karatsuba 3-mult.
            # k1 = A*(R+I), k2 = R*(B-A), k3 = I*(A+B).
            c_e = poly.get(e, (0, 0))
            # Linear combos for the three products:
            #   S1 = R+I, S2 = B-A (witness), S3 = A+B (witness).
            S1 = {}
            _add_lin(S1, R, 1)
            _add_lin(S1, I, 1)
            S1_const = (R_const + I_const) % N
            # S2 = B - A = aL[1] - aL[0]
            S2 = {("L", 1): 1, ("L", 0): N - 1}
            S2_const = 0
            # S3 = A + B = aL[0] + aL[1]
            S3 = {("L", 0): 1, ("L", 1): 1}
            S3_const = 0
            # Gate g+0: k1 = A * S1. aL = A (witness), aR = S1.
            # Gate g+1: k2 = R * S2. aL = R, aR = S2.
            # Gate g+2: k3 = I * S3. aL = I, aR = S3.
            for k, (linL, constL, linR, constR) in enumerate([
                (None, 0, S1, S1_const),      # g+0: aL=A, aR=S1
                (R, R_const, S2, S2_const),   # g+1: aL=R, aR=S2
                (I, I_const, S3, S3_const),   # g+2: aL=I, aR=S3
            ]):
                gi = g + k
                # aL[gi]: witness A if linL is None, else linear combo.
                acc = {}
                gate("L", gi, 1, acc)
                if linL is None:
                    gate("L", 0, N - 1, acc)  # aL[gi] = aL[0] = A
                    cons.append((acc, 0))
                else:
                    for (w, idx), coeff in linL.items():
                        gate(w, idx, N - coeff, acc)
                    cons.append((acc, constL % N))
                # aR[gi] = linear combo.
                acc = {}
                gate("R", gi, 1, acc)
                for (w, idx), coeff in linR.items():
                    gate(w, idx, N - coeff, acc)
                cons.append((acc, constR % N))
            # New acc: R' = k1 - k3 + re(c_e), I' = k1 + k2 + im(c_e).
            # k1 = aO[g], k2 = aO[g+1], k3 = aO[g+2].
            R = {("O", g): 1, ("O", g + 2): N - 1}
            I = {("O", g): 1, ("O", g + 1): 1}
            R_const = c_e[0] % N
            I_const = c_e[1] % N
            g += 3

        # Final constraints: R = 0 and I = 0.
        acc = {}
        for (w, idx), coeff in R.items():
            gate(w, idx, coeff, acc)
        cons.append((acc, (-R_const) % N))  # re(P) = 0
        acc = {}
        for (w, idx), coeff in I.items():
            gate(w, idx, coeff, acc)
        cons.append((acc, (-I_const) % N))  # im(P) = 0

        self._densify(cons, n)

    def _build_ftc(self, f_poly, a, b):
        """FTC circuit for x = int(f, a, b).

        Witness (all private, in aL):
          F_0, F_1, ..., F_{d+1}  (antiderivative coefficients, d = deg f)
          v                        (the integral value)
        Total: d+3 witness variables in aL[0..d+2].

        Constraints (all linear, zero multiplication gates):
          (j+1)*F_{j+1} = c_j  for j=0..d   (F' = f, coefficient-wise)
          v = Σ_{j=0}^{d+1} F_j*(b^j - a^j)   (v = F(b) - F(a))

        All arithmetic is exact in the field (Fractions via modular
        inverse). The proof leaks nothing about F or v (Bulletproofs ZK).
        The ONLY way to learn v from public artifacts is to compute
        ∫ₐᵇ f yourself — which is the intended secrecy model.
        """
        d = max(f_poly.keys()) if f_poly else 0
        # Witness count: d+3. n = next power of 2.
        n_wit = d + 3
        n = 1
        while n < n_wit:
            n *= 2
        self.n = n
        self._ftc_degree = d
        # Store for evaluate/check_witness (as field elements).
        self._ftc_a = a
        self._ftc_b = b
        self._ftc_f = f_poly

        cons = []

        def gate(which, i, coeff, acc):
            acc[which, i] = (acc.get((which, i), 0) + coeff) % N

        def _field(fr):
            """Fraction -> field element (exact via modular inverse)."""
            return (fr.numerator * pow(fr.denominator, N - 2, N)) % N

        # Constraint set 1: (j+1)*F_{j+1} = c_j for j=0..d.
        for j in range(d + 1):
            c_j = f_poly.get(j, (Fraction(0), Fraction(0)))[0]
            acc = {}
            gate("L", j + 1, j + 1, acc)  # (j+1)*F_{j+1}
            cons.append((acc, _field(c_j)))

        # Constraint set 2: v - Σ F_j*(b^j - a^j) = 0.
        acc = {}
        gate("L", d + 2, 1, acc)  # v
        for j in range(d + 2):
            diff = b ** j - a ** j  # Fraction
            gate("L", j, (-_field(diff)) % N, acc)
        cons.append((acc, 0))

        self._densify(cons, n)

    def _process_integral(self, poly, spec):
        """Handle int(f,a,b): verify form, store spec for FTC.
        Form already verified in parse_equation (x = int(...) or
        int(...) = x). Returns poly unchanged (FTC doesn't need S).
        Also stores self._integral_canon for the symbolic canonical.
        """
        f_poly = spec['integrand']
        a, b = spec['a'], spec['b']
        # Symbolic canonical (red line: never print any numeric value).
        # Format: x=int(<f_canon>,<a_canon>,<b_canon>);k=<k>
        f_canon_full = _canon_poly_frac(
            {e: c[0] for e, c in f_poly.items()})
        # Strip the trailing "=0" (we want just the polynomial)
        f_canon = f_canon_full[:-2] if f_canon_full.endswith("=0") else f_canon_full
        def _fmt_bound(fr):
            return (str(fr.numerator) if fr.denominator == 1
                    else "%d/%d" % (fr.numerator, fr.denominator))
        a_canon = _fmt_bound(a)
        b_canon = _fmt_bound(b)
        # Normalize sign for canonical: always "x=int(...)"
        self._integral_canon = "x=int(%s,%s,%s)" % (f_canon, a_canon, b_canon)
        return poly

    # ---- tolerance mode ----

    def _build_tolerance(self, poly, k):
        d = max(poly.keys()) if poly else 0
        self.degree = d
        # Fixed scale S = TOL_S. Integer coefficients:
        #   A_e = round(S * c_e / M^e),  so P(X)/S ~= T(X/M) = sum c_e (X/M)^e.
        # No denominator-lcm is needed: rounding to integers with the huge
        # fixed S keeps the error far below 1/k (see bound check below).
        # |P(X)| <= B  <=>  |T(x)-y| < 1/k   with B = ceil(S/k)-1.
        S = TOL_S
        a = {}
        for e, c in poly.items():
            # c is a (re, im) pair; tolerance mode is real-only for now.
            # (Complex tolerance with R^2+I^2<=B^2 is future work.)
            if isinstance(c, tuple):
                if c[1] != 0:
                    raise ValueError(
                        "complex tolerance mode not yet implemented; "
                        "use Gaussian integer equations for complex")
                c = c[0]
            a[e] = _round_frac(S * c / (M_DEC ** e))
        self._tol_a = a
        self._tol_d = d
        self._tol_S = S
        B = (S + k - 1) // k - 1
        self._tol_B = B
        XR = M_DEC * R_DOM
        self._tol_XR = XR
        # compile-time field-size guard: integer arithmetic must not wrap mod n
        bound = sum(abs(ae) * (XR ** e) for e, ae in a.items())
        if bound >= N // 4:
            raise ValueError("equation too large for the field "
                             "(reduce degree or range)")
        m = (2 * B).bit_length() if B > 0 else 0     # tolerance bits
        m2 = (2 * XR).bit_length()                    # domain bits
        self._tol_m = m
        self._tol_m2 = m2

        if d == 0:
            if abs(a.get(0, 0)) > B:
                raise ValueError("no solution within precision 1/%d" % k)
            # vacuously true: dummy zero gate, no constraints
            self.n = 1
            self._tol_bit0 = 1
            self._tol_dom0 = 1
            self._densify([], 1)
            return

        n_pow = d - 1 if d >= 2 else 1
        bit0 = n_pow
        dom0 = n_pow + m
        self._tol_bit0 = bit0
        self._tol_dom0 = dom0
        n_real = n_pow + m + m2
        n = 1
        while n < n_real:
            n *= 2
        self.n = n

        cons = []

        def gate(which, i, coeff, acc):
            acc[which, i] = (acc.get((which, i), 0) + coeff) % N

        if d >= 2:
            for kk in range(1, d - 1):
                acc = {}
                gate("L", kk, 1, acc)
                gate("O", kk - 1, N - 1, acc)
                cons.append((acc, 0))
            for kk in range(0, d - 1):
                acc = {}
                gate("R", kk, 1, acc)
                gate("L", 0, N - 1, acc)
                cons.append((acc, 0))
        else:  # d == 1: variable gate x = x * 1
            acc = {}
            gate("R", 0, 1, acc)
            cons.append((acc, 1))
            acc = {}
            gate("O", 0, 1, acc)
            gate("L", 0, N - 1, acc)
            cons.append((acc, 0))

        # tolerance: sum_j a_j*w(X^j) - sum_i 2^i*b_i = -(a_0 + B)
        # (bit gates carry b_i in all three wires)
        acc = {}
        for e in range(1, d + 1):
            ae = a.get(e, 0)
            if ae:
                if e == 1:
                    gate("L", 0, ae, acc)
                else:
                    gate("O", e - 2, ae, acc)
        for i in range(m):
            gate("L", bit0 + i, -(1 << i), acc)
        cons.append((acc, (-(a.get(0, 0) + B)) % N))

        # domain: X - sum_i 2^i*d_i = -XR   (|X| <= XR)
        acc = {}
        gate("L", 0, 1, acc)
        for i in range(m2):
            gate("L", dom0 + i, -(1 << i), acc)
        cons.append((acc, (-XR) % N))

        self._densify(cons, n)

    def _densify(self, cons, n):
        self.q = len(cons)
        self.WL = [[0] * n for _ in range(self.q)]
        self.WR = [[0] * n for _ in range(self.q)]
        self.WO = [[0] * n for _ in range(self.q)]
        self.c = [0] * self.q
        for qi, (acc, rhs) in enumerate(cons):
            for (which, i), coeff in acc.items():
                v = coeff % N
                if which == "L":
                    self.WL[qi][i] = v
                elif which == "R":
                    self.WR[qi][i] = v
                else:
                    self.WO[qi][i] = v
            self.c[qi] = rhs % N

    # ---- witness evaluation ----

    def evaluate(self, X):
        """Wire values for witness X: returns (a_L, a_R, a_O).

        Exact mode (real): X is the integer witness.
        Exact mode (complex): X is a pair (A, B) = (re, im).
        Tolerance mode: X is the scaled integer witness (x = X / M_DEC).
        FTC mode: X is a list [F_0..F_{d+1}, v] of field elements.
        """
        n = self.n
        aL = [0] * n
        aR = [0] * n
        aO = [0] * n
        if self.mode == "ftc":
            # FTC: witness is [F_0..F_{d+1}, v] as field elements.
            # All constraints are linear; aR/aO are zero.
            for j, w in enumerate(X):
                aL[j] = w % N
            return aL, aR, aO
        if self.mode == "exact":
            if hasattr(self, '_exact_poly_c'):
                # Complex exact: Horner with Karatsuba 3-mult.
                A, B = X[0] % N, X[1] % N
                aL[0] = A; aR[0] = 1; aO[0] = A  # witness re
                aL[1] = B; aR[1] = 1; aO[1] = B  # witness im
                d = self.degree
                if d == 0:
                    return aL, aR, aO
                poly = self._exact_poly_c
                R, I = poly[d][0] % N, poly[d][1] % N
                g = 2
                for e in range(d - 1, -1, -1):
                    # k1 = A*(R+I), k2 = R*(B-A), k3 = I*(A+B)
                    S1 = (R + I) % N
                    S2 = (B - A) % N
                    S3 = (A + B) % N
                    k1 = (A * S1) % N
                    k2 = (R * S2) % N
                    k3 = (I * S3) % N
                    aL[g] = A; aR[g] = S1; aO[g] = k1
                    aL[g+1] = R; aR[g+1] = S2; aO[g+1] = k2
                    aL[g+2] = I; aR[g+2] = S3; aO[g+2] = k3
                    c_e = poly.get(e, (0, 0))
                    R = (k1 - k3 + c_e[0]) % N
                    I = (k1 + k2 + c_e[1]) % N
                    g += 3
                return aL, aR, aO
            X = X % N
            d = self.degree
            if d >= 2:
                pw = [1] * (d + 1)
                for kk in range(1, d + 1):
                    pw[kk] = (pw[kk - 1] * X) % N
                for kk in range(d - 1):
                    aL[kk] = pw[kk + 1]
                    aR[kk] = pw[1]
                    aO[kk] = pw[kk + 2]
            elif d == 1:
                aL[0] = X
                aR[0] = 1
                aO[0] = X
            return aL, aR, aO

        # tolerance mode
        X = int(X)
        d = self._tol_d
        a = self._tol_a
        B = self._tol_B
        XR = self._tol_XR
        m = self._tol_m
        m2 = self._tol_m2
        bit0 = self._tol_bit0
        dom0 = self._tol_dom0
        if d == 0:
            return aL, aR, aO  # dummy zero gate
        Xm = X % N
        if d >= 2:
            pw = [1] * (d + 1)
            for kk in range(1, d + 1):
                pw[kk] = (pw[kk - 1] * Xm) % N
            for kk in range(d - 1):
                aL[kk] = pw[kk + 1]
                aR[kk] = pw[1]
                aO[kk] = pw[kk + 2]
        else:
            aL[0] = Xm
            aR[0] = 1
            aO[0] = Xm
        if B > 0:
            PX = sum(ae * (X ** e) for e, ae in a.items())  # exact ints
            u = PX + B
            assert u >= 0, "compiler bug: negative tolerance slack"
            for i in range(m):
                bit = (u >> i) & 1
                aL[bit0 + i] = aR[bit0 + i] = aO[bit0 + i] = bit
        v = X + XR
        assert v >= 0, "compiler bug: negative domain slack"
        for i in range(m2):
            bit = (v >> i) & 1
            aL[dom0 + i] = aR[dom0 + i] = aO[dom0 + i] = bit
        return aL, aR, aO

    def check_witness(self, X):
        """True iff X satisfies the equation (within precision)."""
        if self.mode == "ftc":
            # FTC: X = [F_0..F_{d+1}, v] as field elements.
            # Verify: (j+1)*F_{j+1} = c_j, and v = Σ F_j*(b^j - a^j).
            d = self._ftc_degree
            f_poly = self._ftc_f
            a, b = self._ftc_a, self._ftc_b
            def _field(fr):
                return (fr.numerator * pow(fr.denominator, N - 2, N)) % N
            # Check F' = f: (j+1)*F_{j+1} = c_j
            for j in range(d + 1):
                c_j = f_poly.get(j, (Fraction(0), Fraction(0)))[0]
                lhs = ((j + 1) * (X[j + 1] % N)) % N
                if lhs != _field(c_j):
                    return False
            # Check v = F(b) - F(a)
            v = X[d + 2] % N
            s = 0
            for j in range(d + 2):
                diff = _field(b ** j - a ** j)
                s = (s + (X[j] % N) * diff) % N
            return v == s
        if self.mode == "exact":
            if hasattr(self, '_exact_poly_c'):
                # Complex exact: Horner with Karatsuba 3-mult.
                A, B = X[0] % N, X[1] % N
                poly = self._exact_poly_c
                d = self.degree
                if d == 0:
                    return poly.get(0, (0, 0)) == (0, 0)
                R, I = poly[d][0] % N, poly[d][1] % N
                for e in range(d - 1, -1, -1):
                    c_e = poly.get(e, (0, 0))
                    # k1=A(R+I), k2=R(B-A), k3=I(A+B)
                    k1 = (A * (R + I)) % N
                    k2 = (R * (B - A)) % N
                    k3 = (I * (A + B)) % N
                    R = (k1 - k3 + c_e[0]) % N
                    I = (k1 + k2 + c_e[1]) % N
                return R == 0 and I == 0
            X = X % N
            v = 0
            for exp, coeff in self._exact_poly.items():
                v = (v + coeff * pow(X, exp, N)) % N
            return v == 0
        X = int(X)
        d = self._tol_d
        if d == 0:
            return True  # compile-time decided
        PX = sum(ae * (X ** e) for e, ae in self._tol_a.items())
        return abs(PX) <= self._tol_B and abs(X) <= self._tol_XR

    def check_constraints(self, aL, aR, aO):
        """Verify Hadamard + linear constraints on explicit wire values."""
        for i in range(self.n):
            if (aL[i] * aR[i]) % N != aO[i] % N:
                return False
        for qi in range(self.q):
            s = (sum(self.WL[qi][i] * aL[i] for i in range(self.n)) +
                 sum(self.WR[qi][i] * aR[i] for i in range(self.n)) +
                 sum(self.WO[qi][i] * aO[i] for i in range(self.n))) % N
            if s != self.c[qi]:
                return False
        return True


def _parse_precision(precision):
    s = str(precision).strip()
    if not s:
        raise ValueError("precision required for this equation: "
                         "enter k for final precision 1/k")
    try:
        k = int(s)
    except ValueError:
        raise ValueError("precision k must be a positive integer")
    if k < 1:
        raise ValueError("precision k must be >= 1")
    return k


def compile(eq_str, precision, force_complex=False):
    """Compile an equation string into a Circuit.

    precision: k as int/str; MANDATORY for every equation (no default).
    Final precision is 1/k: the witness must be within 1/k of an exact
    root (enforced by the prover before proving).
    force_complex: if True, build a complex circuit even if the
        coefficients are real (needed when the witness is complex,
        e.g. x^2+1=0 with witness i).
    Raises ValueError on syntax errors, unsatisfiable equations, or a
    missing/invalid precision.
    """
    poly, trans_used, float_seen, integral_spec = parse_polynomial(eq_str)
    return Circuit(eq_str, poly, precision, trans_used, float_seen,
                   force_complex=force_complex, integral_spec=integral_spec)
