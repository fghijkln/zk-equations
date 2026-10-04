"""Equation compiler:  "x^3 + 2*x + 5 = 42"  ->  arithmetic circuit.

v1 language: single variable x, integer constants, + - * ^ (non-negative
integer exponent), parentheses, '=', implicit multiplication ("2x").
'**' is accepted as an alias for '^'.

The equation is normalized to a univariate polynomial
    c_0 + c_1*x + ... + c_d*x^d = 0   (over Z_n, n = curve order)
and compiled to:
  * n multiplication gates  (constraints a_L o a_R = a_O),
  * Q linear constraints     (W_L*a_L + W_R*a_R + W_O*a_O = c).

Gate layout (d >= 2): gate k (0-based) computes x^{k+2} = x^{k+1} * x,
  i.e. a_L[k] = x^{k+1}, a_R[k] = x, a_O[k] = x^{k+2}.
For d <= 1, gate 0 is the "variable gate"  x = x * 1.
n is padded up to a power of two with zero gates (0*0 = 0); the inner
product argument (paper Section 3) requires n = 2^k.
"""

from . import curve

N = curve.N


# ---------------- tokenizer ----------------

class Tok:
    INT, X, PLUS, MINUS, STAR, POW, LP, RP, EQ = range(9)

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
        elif ch.isdigit():
            j = i
            while j < len(s) and s[j].isdigit():
                j += 1
            toks.append(Tok(Tok.INT, int(s[i:j])))
            i = j
        elif ch in "xX":
            toks.append(Tok(Tok.X))
            i += 1
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


# ---------------- parser: AST -> polynomial {exp: coeff} ----------------

class Parser:
    def __init__(self, toks):
        self.toks = toks
        self.pos = 0

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
            elif t is not None and t.kind in (Tok.INT, Tok.X, Tok.LP):
                # implicit multiplication: 2x, 2(x+1), x(x+1)
                node = _pmul(node, self.parse_factor())
            else:
                return node

    def parse_factor(self):
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
            return {0: t.val % N}
        if t.kind == Tok.X:
            return {1: 1}
        if t.kind == Tok.LP:
            node = self.parse_expr()
            if self.next().kind != Tok.RP:
                raise ValueError("missing ')'")
            return node
        if t.kind == Tok.MINUS:
            return _pneg(self.parse_primary())
        if t.kind == Tok.PLUS:
            return self.parse_primary()
        raise ValueError("unexpected token %r" % (t,))


def _clean(p):
    return {e: c % N for e, c in p.items() if c % N != 0}


def _padd(a, b):
    r = dict(a)
    for e, c in b.items():
        r[e] = (r.get(e, 0) + c) % N
    return _clean(r)


def _pneg(a):
    return {e: (-c) % N for e, c in a.items()}


def _psub(a, b):
    return _padd(a, _pneg(b))


def _pmul(a, b):
    r = {}
    for e1, c1 in a.items():
        for e2, c2 in b.items():
            e = e1 + e2
            r[e] = (r.get(e, 0) + c1 * c2) % N
    return _clean(r)


def _ppow(a, e):
    r = {0: 1}
    base = a
    while e:
        if e & 1:
            r = _pmul(r, base)
        base = _pmul(base, base)
        e >>= 1
    return r


def parse_polynomial(eq_str):
    """Parse 'lhs = rhs' into {exp: coeff} for lhs - rhs = 0 (mod N)."""
    return Parser(tokenize(eq_str)).parse_equation()


# ---------------- circuit ----------------

class Circuit:
    """Arithmetic circuit for one polynomial equation.

    Attributes:
      n: number of multiplication gates (power of two, >= 1)
      q: number of linear constraints
      WL, WR, WO: q x n matrices (list of q rows of n scalars)
      c: length-q RHS vector
      eq_str: original equation string
      degree: polynomial degree
    """

    def __init__(self, eq_str, poly):
        self.eq_str = eq_str
        self.poly = poly
        self.degree = max(poly.keys()) if poly else 0
        d = self.degree

        # --- gate layout ---
        if d >= 2:
            n_real = d - 1          # gates 0..d-2 compute x^2..x^d
        else:
            n_real = 1             # gate 0 is the variable gate / dummy
        n = 1
        while n < n_real:
            n *= 2
        self.n = n

        # --- linear constraints as (dict wire->coeff, rhs) then densified ---
        # wire ids: ("L", i), ("R", i), ("O", i)
        cons = []  # list of (dict, rhs)

        def gate_con(which, i, coeff, acc):
            acc[which, i] = (acc.get((which, i), 0) + coeff) % N

        if d >= 2:
            for k in range(1, d - 1):
                acc = {}
                gate_con("L", k, 1, acc)
                gate_con("O", k - 1, N - 1, acc)
                cons.append((acc, 0))          # a_L[k] = a_O[k-1]
            for k in range(0, d - 1):
                acc = {}
                gate_con("R", k, 1, acc)
                gate_con("L", 0, N - 1, acc)
                cons.append((acc, 0))          # a_R[k] = a_L[0] (= x)
            acc = {}
            for exp, coeff in poly.items():
                if exp == 0:
                    continue
                if exp == 1:
                    gate_con("L", 0, coeff, acc)
                else:
                    gate_con("O", exp - 2, coeff, acc)
            cons.append((acc, (-poly.get(0, 0)) % N))  # main equation
        elif d == 1:
            c1 = poly.get(1, 0)
            acc = {}
            gate_con("R", 0, 1, acc)
            cons.append((acc, 1))              # a_R[0] = 1
            acc = {}
            gate_con("O", 0, 1, acc)
            gate_con("L", 0, N - 1, acc)
            cons.append((acc, 0))              # a_O[0] = a_L[0]
            acc = {}
            gate_con("L", 0, c1, acc)
            cons.append((acc, (-poly.get(0, 0)) % N))  # c1*x = -c0
        else:  # d == 0: constant equation; gate 0 is a zero dummy gate
            if poly.get(0, 0) % N != 0:
                raise ValueError("equation has no solution (constant false)")
            # trivially true: no constraints

        self.q = len(cons)
        self.WL = [[0] * n for _ in range(self.q)]
        self.WR = [[0] * n for _ in range(self.q)]
        self.WO = [[0] * n for _ in range(self.q)]
        self.c = [0] * self.q
        for qi, (acc, rhs) in enumerate(cons):
            for (which, i), coeff in acc.items():
                if which == "L":
                    self.WL[qi][i] = coeff
                elif which == "R":
                    self.WR[qi][i] = coeff
                else:
                    self.WO[qi][i] = coeff
            self.c[qi] = rhs % N

    def evaluate(self, x):
        """Wire values for witness x: returns (a_L, a_R, a_O) length-n lists."""
        x = x % N
        d = self.degree
        n = self.n
        aL = [0] * n
        aR = [0] * n
        aO = [0] * n
        if d >= 2:
            pw = [1] * (d + 1)
            for k in range(1, d + 1):
                pw[k] = (pw[k - 1] * x) % N
            for k in range(d - 1):
                aL[k] = pw[k + 1]
                aR[k] = pw[1]
                aO[k] = pw[k + 2]
        elif d == 1:
            aL[0] = x
            aR[0] = 1
            aO[0] = x
        # d == 0: all-zero dummy gate; padding gates stay zero (0*0 = 0)
        return aL, aR, aO

    def check_witness(self, x):
        """True iff x satisfies the equation mod N."""
        v = 0
        for exp, coeff in self.poly.items():
            v = (v + coeff * pow(x % N, exp, N)) % N
        return v == 0

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


def compile(eq_str):
    """Compile an equation string into a Circuit. Raises ValueError on
    syntax errors or provably-unsatisfiable (constant-false) equations."""
    return Circuit(eq_str, parse_polynomial(eq_str))
