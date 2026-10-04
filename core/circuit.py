"""Equation compiler:  "x^3 + 2*x + 5 = 38"  ->  arithmetic circuit.

v2 language: single variable x, integer/decimal constants, + - * ^
(non-negative integer exponent), parentheses, '=', implicit
multiplication ("2x"), '**' alias, and transcendental functions
sin/cos/exp/ln.

Transcendental functions are evaluated as fixed degree-12 Taylor
polynomials (sin/cos/exp around 0, ln around x=1 i.e. series in x-1).
The proved statement is *exactly* about that Taylor polynomial --
nothing is hidden: for |x| small it genuinely approximates the true
function. Accuracy guide (degree 12): sin/cos good for |x| <= 2*pi,
exp good for |x| <= 2, ln good for 0.5 <= x <= 1.5 (at x=2 the
degree-12 error is already ~0.04 -- the series converges slowly at the
interval edge, and the prover will honestly refuse a witness that does
not meet the requested precision).

Two modes, chosen automatically:
  exact:     the equation is an integer polynomial (no decimals, no
             transcendental functions). Proved exactly, mod n.
  tolerance: decimals or transcendental functions present. Values use
             4 fixed decimal places (integers = value * 10^4); the user
             supplies a precision denominator k and the circuit proves
                 |P(X)| <= B   i.e.   |T(x) - y| < 1/k,
             where P is the scaled integer polynomial, B = ceil(S/k)-1,
             S the scale factor. |X| <= 10*10^4 is range-checked so the
             mod-n arithmetic faithfully represents integer arithmetic.
             |P(X)| <= B is enforced by bit-decomposing u = P(X)+B.

Proofs are bound to the *canonical* equation (polynomial normal form),
so "x^3+2*x+5=38" and the term-moved "x^3+2*x=33" verify each other's
proofs; genuinely different equations do not.
"""

from fractions import Fraction
from math import gcd

from . import curve

N = curve.N

M_DEC = 10 ** 4     # fixed decimal scale: 4 places
R_DOM = 10          # |x| <= R_DOM enforced (tolerance mode)
TAYLOR_DEG = 12


# ---------------- Taylor polynomials {exp: Fraction} ----------------

def _fact(n):
    r = 1
    for i in range(2, n + 1):
        r *= i
    return r


def _taylor_sin():
    return {2 * j + 1: Fraction((-1) ** j, _fact(2 * j + 1))
            for j in range(TAYLOR_DEG // 2 + 1) if 2 * j + 1 <= TAYLOR_DEG}


def _taylor_cos():
    return {2 * j: Fraction((-1) ** j, _fact(2 * j))
            for j in range(TAYLOR_DEG // 2 + 1)}


def _taylor_exp():
    return {j: Fraction(1, _fact(j)) for j in range(TAYLOR_DEG + 1)}


def _taylor_ln():
    # ln(x) = sum_{j=1}^{12} (-1)^{j+1} (x-1)^j / j
    res = {}
    for j in range(1, TAYLOR_DEG + 1):
        p = _ppow({1: Fraction(1), 0: Fraction(-1)}, j)
        res = _padd(res, _pscale(p, Fraction((-1) ** (j + 1), j)))
    return res


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
            return {0: t.val}
        if t.kind == Tok.X:
            return {1: Fraction(1)}
        if t.kind == Tok.IDENT:
            name = t.val
            if name not in _TAYLORS:
                raise ValueError("unknown function '%s' "
                                 "(supported: sin, cos, exp, ln)" % name)
            if self.next().kind != Tok.LP:
                raise ValueError("expected '(' after '%s'" % name)
            arg = self.parse_expr()
            if self.next().kind != Tok.RP:
                raise ValueError("missing ')'")
            self.trans_used = True
            return _pcompose(_TAYLORS[name], arg)
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
    """Substitute polynomial p into Taylor polynomial t: t(p)."""
    r = {}
    for e, c in t.items():
        r = _padd(r, _pscale(_ppow(p, e), c))
    return r


_TAYLORS = {
    "sin": _taylor_sin(),
    "cos": _taylor_cos(),
    "exp": _taylor_exp(),
    "ln": _taylor_ln(),
}


def parse_polynomial(eq_str):
    """Parse 'lhs = rhs' into ({exp: Fraction} for lhs - rhs, trans_used)."""
    parser = Parser(tokenize(eq_str))
    return parser.parse_equation(), parser.trans_used


# ---------------- canonical form ----------------

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
      mode: "exact" | "tolerance"
      precision: k (int) in tolerance mode, else None
      canonical: canonical equation string the proof is bound to
    """

    def __init__(self, eq_str, poly, precision, trans_used):
        self.eq_str = eq_str
        # Normalize overall sign (leading coefficient positive) so that
        # "38 = x^3+2*x+5" and "x^3+2*x+5 = 38" compile to the identical
        # circuit and canonical form.
        if poly:
            lc = poly[max(poly.keys())]
            if lc < 0:
                poly = {e: -c for e, c in poly.items()}
        needs_tol = trans_used or any(c.denominator != 1
                                      for c in poly.values())
        if not needs_tol:
            self.mode = "exact"
            self.precision = None
            poly_int = {e: int(c) % N for e, c in poly.items() if c != 0}
            self._build_exact(poly_int)
            self.canonical = _canon_poly(
                {e: _to_signed(v) for e, v in poly_int.items()})
        else:
            self.mode = "tolerance"
            self.precision = _parse_precision(precision)
            self._build_tolerance(poly, self.precision)
            self.canonical = ("%s;k=%d"
                              % (_canon_poly(self._tol_a), self.precision))

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

    # ---- tolerance mode ----

    def _build_tolerance(self, poly, k):
        d = max(poly.keys()) if poly else 0
        self.degree = d
        # D = lcm of denominators; S = M_DEC^d * D; P(X) = S*T(X/M_DEC)
        D = 1
        for c in poly.values():
            D = D * c.denominator // gcd(D, c.denominator)
        S = (M_DEC ** d) * D
        a = {}
        for e, c in poly.items():
            a[e] = int(D * c) * (M_DEC ** (d - e))  # D*c is integral
        self._tol_a = a
        self._tol_d = d
        self._tol_S = S
        B = (S + k - 1) // k - 1   # ceil(S/k) - 1 ; |P(X)| <= B  <=>  |T(x)-y| < 1/k
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

        Exact mode: X is the integer witness.
        Tolerance mode: X is the scaled integer witness (x = X / M_DEC).
        """
        n = self.n
        aL = [0] * n
        aR = [0] * n
        aO = [0] * n
        if self.mode == "exact":
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
        if self.mode == "exact":
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


def compile(eq_str, precision=""):
    """Compile an equation string into a Circuit.

    precision: k as int/str; required iff the equation needs tolerance
    mode (decimals or transcendental functions). Final precision is 1/k.
    Raises ValueError on syntax errors, unsatisfiable equations, or a
    missing/invalid precision.
    """
    poly, trans_used = parse_polynomial(eq_str)
    return Circuit(eq_str, poly, precision, trans_used)
