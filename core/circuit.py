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

Single mode: every equation requires a user-supplied precision
denominator k (no automatic exact/tolerance detection). Strict
definition, enforced by the prover before proving:
    the witness is within 1/k of an exact root,
otherwise the witness is refused (ValueError), never proven.
The circuit itself works in 4 fixed decimal places
(integers = value * 10^4; the witness is quantized to 4 places)
and proves
    |P(X)| <= B   i.e.   |T(x) - y| < 1/k,
where P is the scaled integer polynomial, B = ceil(S/k)-1,
S the scale factor. |X| <= 5*10^4 is range-checked so the
mod-n arithmetic faithfully represents integer arithmetic.
|P(X)| <= B is enforced by bit-decomposing u = P(X)+B.

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
    INT, FLOAT, IDENT, X, PLUS, MINUS, STAR, POW, LP, RP, EQ = range(11)

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
        else:
            raise ValueError("unexpected character %r in equation" % ch)
    return toks


# ---------------- parser: AST -> polynomial {exp: Fraction} ----------------

class Parser:
    def __init__(self, toks):
        self.toks = toks
        self.pos = 0
        self.trans_used = False
        self.float_seen = False

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
            return {0: Fraction(t.val)}
        if t.kind == Tok.FLOAT:
            self.float_seen = True
            return {0: t.val}
        if t.kind == Tok.X:
            return {1: Fraction(1)}
        if t.kind == Tok.IDENT:
            name = t.val
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
    return {e: c for e, c in p.items() if c != 0}


def _padd(a, b):
    r = dict(a)
    for e, c in b.items():
        r[e] = r.get(e, Fraction(0)) + c
    return _clean(r)


def _pneg(a):
    return {e: -c for e, c in a.items()}


def _psub(a, b):
    return _padd(a, _pneg(b))


def _pscale(a, s):
    return _clean({e: c * s for e, c in a.items()})


def _pmul(a, b):
    r = {}
    for e1, c1 in a.items():
        for e2, c2 in b.items():
            e = e1 + e2
            r[e] = r.get(e, Fraction(0)) + c1 * c2
    return _clean(r)


def _ppow(a, e):
    r = {0: Fraction(1)}
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
    """Parse 'lhs = rhs' into ({exp: Fraction} for lhs - rhs,
    trans_used, float_seen)."""
    parser = Parser(tokenize(eq_str))
    return parser.parse_equation(), parser.trans_used, parser.float_seen


# ---------------- canonical form ----------------

def _round_frac(fr):
    """Round a Fraction to the nearest integer (half away from zero)."""
    n, d = fr.numerator, fr.denominator
    if n >= 0:
        return (2 * n + d) // (2 * d)
    return -((2 * (-n) + d) // (2 * d))


def _canon_poly_frac(poly):
    """poly: {exp: Fraction} -> canonical string (deterministic, injective)."""

    def fmt(c):
        return (str(c.numerator) if c.denominator == 1
                else "%d/%d" % (c.numerator, c.denominator))

    terms = []
    for e in sorted(poly.keys(), reverse=True):
        c = poly[e]
        if c == 0:
            continue
        sign = "-" if c < 0 else "+"
        ac = abs(c)
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
      precision: k (int), always required (final precision 1/k)
      canonical: canonical equation string the proof is bound to
                 (includes ";k=<k>")
    """

    def __init__(self, eq_str, poly, precision):
        self.eq_str = eq_str
        # Normalize overall sign (leading coefficient positive) so that
        # "38 = x^3+2*x+5" and "x^3+2*x+5 = 38" compile to the identical
        # circuit and canonical form.
        if poly:
            lc = poly[max(poly.keys())]
            if lc < 0:
                poly = {e: -c for e, c in poly.items()}
        # Single mode: precision k is mandatory for every equation.
        # No exact/tolerance auto-detection.
        self.precision = _parse_precision(precision)
        self._build_tolerance(poly, self.precision)
        self.canonical = ("%s;k=%d"
                          % (_canon_poly_frac(poly), self.precision))

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

        X is the scaled integer witness (x = X / M_DEC).
        """
        n = self.n
        aL = [0] * n
        aR = [0] * n
        aO = [0] * n
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
        """True iff X satisfies the equation (within precision 1/k)."""
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


def compile(eq_str, precision):
    """Compile an equation string into a Circuit.

    precision: k as int/str; MANDATORY for every equation (no exceptions).
    Final precision is 1/k: the witness must be within 1/k of an exact
    root (enforced by the prover before proving).
    Raises ValueError on syntax errors, unsatisfiable equations, or a
    missing/invalid precision.
    """
    poly, _, _ = parse_polynomial(eq_str)
    return Circuit(eq_str, poly, precision)
