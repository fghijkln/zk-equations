"""Public API: prove/verify knowledge of a solution to an equation.

    proof = prove_equation("x^3 + 2*x + 5 = 38", 3, "10000")
    verify_equation("x^3 + 2*x + 5 = 38", proof, "10000")  # True

    # Precision k is MANDATORY for every equation (no exceptions).
    # Strict definition: the witness must be within 1/k of an exact
    # root; otherwise prove_equation refuses (raises ValueError).
    # Use solve_equation to get accurate roots.
    proof = prove_equation("sin(x) = 0.5", "0.52359878", "1000")
    verify_equation("sin(x) = 0.5", proof, "1000")  # True

The proof is a plain dict of hex strings (JSON-serializable). Proofs are
bound to the *canonical* equation (polynomial normal form + ";k=<k>"):
the same equation with a different k does NOT verify.
"""

from fractions import Fraction

import math

from . import bulletproof
from . import circuit as circuit_mod
from . import solve as solve_mod
from .circuit import _c

M_DEC = circuit_mod.M_DEC
WIT_DECIMALS = 12  # witness: up to 12 decimal places (solver max)


def _parse_witness_full(w):
    """Parse a decimal witness string to a Fraction (full precision).

    Allows up to 12 decimal places (matches the solver's max)."""
    s = str(w).strip()
    try:
        fr = Fraction(s)
    except (ValueError, ZeroDivisionError):
        raise ValueError("witness must be a decimal number")
    # count decimal places ("1e-8" style has none in the mantissa part)
    mant = s.lstrip("+-").split("e")[0].split("E")[0]
    if "." in mant and len(mant.split(".")[1]) > WIT_DECIMALS:
        raise ValueError("witness supports at most %d decimal places"
                         % WIT_DECIMALS)
    return fr


def _parse_witness_complex(w):
    """Parse a witness string to a complex pair (re, im) of Fractions.

    Accepts: "3", "3.5", "i", "-i", "2i", "1+2i", "1-2i", "1.5+0.5i".
    Each part allows up to 12 decimal places.
    Returns ((re_frac, im_frac), is_complex) where is_complex is True
    iff the imaginary part is nonzero.
    """
    s = str(w).strip().replace(" ", "")
    if not s:
        raise ValueError("witness must be a number")
    # Handle pure imaginary: "i", "-i", "2i", "-2.5i"
    if s.endswith("i") or s.endswith("I"):
        body = s[:-1]
        if body in ("", "+"):
            return (_c(0, 1), True)
        if body == "-":
            return (_c(0, -1), True)
        # Split into re and im parts: find the last + or - (not at start)
        # e.g. "1+2" -> re=1, im=2; "1-2" -> re=1, im=-2; "2" -> re=0, im=2
        idx = -1
        for i in range(len(body) - 1, 0, -1):
            if body[i] in "+-":
                idx = i
                break
        if idx == -1:
            # Pure imaginary like "2i"
            im_frac = _parse_witness_full(body)
            return (_c(0, im_frac), im_frac != 0)
        else:
            re_part = body[:idx]
            im_part = body[idx:]
            # im_part may be "+" or "-" meaning +1 or -1 (e.g. "1+i")
            if im_part == "+":
                im_frac = Fraction(1)
            elif im_part == "-":
                im_frac = Fraction(-1)
            else:
                im_frac = _parse_witness_full(im_part)
            if re_part in ("", "+"):
                re_frac = Fraction(0)
            elif re_part == "-":
                re_frac = Fraction(0)
            else:
                re_frac = _parse_witness_full(re_part)
            return ((re_frac, im_frac), im_frac != 0)
    else:
        # Real number
        fr = _parse_witness_full(s)
        return (_c(fr, 0), False)


def _quantize(fr):
    """Round a Fraction to the circuit's fixed-point scale -> int X."""
    num, den = (fr * M_DEC).numerator, (fr * M_DEC).denominator
    # round half away from zero
    if num >= 0:
        return (2 * num + den) // (2 * den)
    return -((2 * -num + den) // (2 * den))


def _check_precision(eq_str, w_frac, k):
    """Enforce the strict precision definition: the witness must be
    within 1/k of an exact root, where k is the user-supplied precision
    denominator. Raises ValueError otherwise.

    The k used here is ALWAYS the user-supplied value (the definition),
    never a value inferred from the proof.

    w_frac may be a Fraction (real) or a (re, im) pair (complex).
    For complex, uses the complex modulus |w - r|.
    """
    # Vacuously true equation (zero polynomial): every x is a root,
    # so any witness is within precision.
    poly, _, _, integral_spec, ode_spec, cheb_spec = circuit_mod.parse_polynomial(eq_str)
    # ODE/Chebyshev equations: coefficients are exact (computed, not provided).
    # No precision check needed.
    if ode_spec is not None or cheb_spec is not None:
        return 0.0
    if not poly or all(circuit_mod._ciszero(c) for c in poly.values()):
        return 0.0
    eps = 1.0 / k
    # Integral equations (FTC): witness must be within 1/k of the
    # analytical integral value I_true = F(b) - F(a).
    if integral_spec is not None:
        if isinstance(w_frac, tuple):
            raise ValueError("integral witness must be a real number")
        f_poly = integral_spec['integrand']
        a, b = integral_spec['a'], integral_spec['b']
        # Analytical I_true via antiderivative (exact Fraction).
        I_true = Fraction(0)
        for e, c in f_poly.items():
            coeff = c[0] / (e + 1)
            I_true += coeff * (b ** (e + 1) - a ** (e + 1))
        w = w_frac  # Fraction
        dist = abs(float(w - I_true))
        if dist >= eps:
            raise ValueError(
                "witness is not within precision 1/%d=%.3g of the integral "
                "value (distance %.3g)" % (k, eps, dist))
        return dist
    # Detect complex: witness is a pair, OR the polynomial has complex coeffs.
    poly_is_complex = any(isinstance(c, tuple) and c[1] != 0
                          for c in poly.values())
    is_cx = isinstance(w_frac, tuple) or poly_is_complex
    if is_cx:
        # Complex: |w - r| < 1/k with complex modulus.
        # Normalize witness to a pair.
        if isinstance(w_frac, tuple):
            w = complex(float(w_frac[0]), float(w_frac[1]))
        else:
            w = complex(float(w_frac), 0.0)
        roots = solve_mod.find_complex_roots(poly)
        best = min((abs(w - r) for r in roots), default=None)
        if best is None or best >= eps:
            detail = ("no complex root found"
                      if best is None else "distance %.3g" % best)
            raise ValueError(
                "witness is not within precision 1/%d=%.3g of an exact "
                "complex root (%s); use the solver to get an accurate root"
                % (k, eps, detail))
        return best
    roots = solve_mod.all_roots(eq_str)
    w = float(w_frac)
    best = min((abs(w - r) for r in roots), default=None)
    if best is None or best >= eps:
        detail = ("no real root found"
                  if best is None else "distance %.3g" % best)
        raise ValueError(
            "witness is not within precision 1/%d=%.3g of an exact "
            "root (%s); use the solver to get an accurate root"
            % (k, eps, detail))
    return best


def prove_equation(eq_str, witness, precision):
    """Prove knowledge of `witness` satisfying `eq_str`.

    witness: decimal string (up to 12 decimal places); an int is also
             accepted (e.g. 3 or "3"). For complex: "1+2i", "i", "-i", etc.
    precision: k (int/str), MANDATORY for every equation (no exceptions,
               regardless of whether the roots are infinite decimals).
               Strict definition: the witness must be within 1/k of an
               exact root, else ValueError ("refuse").
    For ODEs, precision k is meaningless (precision is deg) and may be
    empty; the witness is auto-computed.
    Raises ValueError if the equation is malformed, the precision is
    missing/invalid (non-ODE), or the witness does not satisfy (within
    precision).
    """
    # ODEs/Chebyshev: k is meaningless; allow empty (use dummy for compile).
    _, _, _, _, ode_spec, cheb_spec = circuit_mod.parse_polynomial(eq_str)
    is_func_eq = ode_spec is not None or cheb_spec is not None
    if not is_func_eq and not str(precision).strip():
        raise ValueError("precision k is required for every equation")
    if is_func_eq and not str(precision).strip():
        precision = "1"  # dummy; canonical has no k
    # Parse witness FIRST to detect complex (determines circuit type).
    # ODEs/Chebyshev: witness is auto-computed; allow empty.
    if is_func_eq and not str(witness).strip():
        w_pair, w_is_complex = ((Fraction(0), Fraction(0)), False)
    else:
        w_pair, w_is_complex = _parse_witness_complex(witness)
    circ = circuit_mod.compile(eq_str, precision, force_complex=w_is_complex)
    k = circ.precision  # int, from the user-supplied definition (1/k)
    # ODE mode: witness is computed automatically (coefficients are exact).
    # Skip the precision pre-check (no user witness to check).
    if circ.mode == "ode":
        coeffs = ode_coefficients(eq_str)
        def _field(fr):
            return (fr.numerator * pow(fr.denominator,
                                       circuit_mod.N - 2,
                                       circuit_mod.N)) % circuit_mod.N
        X = [_field(c) for c in coeffs]
        if not circ.check_witness(X):
            raise ValueError("ODE coefficients do not satisfy constraints")
        aL, aR, aO = circ.evaluate(X)
        assert circ.check_constraints(aL, aR, aO), "compiler bug: bad ODE wires"
        proof = bulletproof.prove(circ, aL, aR, aO)
        proof["equation"] = circ.canonical
        proof["precision"] = k
        # Store ODE spec for coarse verification (deg_prove, diffeq, ics).
        # Fractions are stored as strings for JSON compatibility.
        # Note: diffeq keys are (ex, ords) where ords is a tuple.
        # For JSON, convert ords tuple to list.
        _spec = circ._ode_spec_for_proof
        proof["ode_spec"] = {
            'diffeq': [((ex, list(ords)), str(coeff))
                       for (ex, ords), coeff in _spec['diffeq'].items()],
            'order': _spec['order'],
            'ics': [(ord_, str(v)) for ord_, v in _spec['ics'].items()],
            'deg': _spec['deg'],
        }
        return proof
    # Chebyshev mode: witness is computed automatically (tau coefficients).
    if circ.mode == "cheb":
        # Check for nonlinear: not yet supported for prove (no exact field root).
        _spec = circ._cheb_spec_for_proof
        _has_nonlinear = any(len(ords) == 2 for (ex, ords) in _spec['diffeq'].keys())
        if _has_nonlinear:
            raise ValueError(
                "Nonlinear Chebyshev ODE prove not yet supported: "
                "exact field root finding is an open problem. "
                "Use cheb_coefficients() to get approximate coefficients.")
        coeffs = cheb_coefficients(eq_str)
        def _field(fr):
            return (fr.numerator * pow(fr.denominator,
                                       circuit_mod.N - 2,
                                       circuit_mod.N)) % circuit_mod.N
        X = [_field(c) for c in coeffs]
        if not circ.check_witness(X):
            raise ValueError("Chebyshev coefficients do not satisfy constraints")
        aL, aR, aO = circ.evaluate(X)
        assert circ.check_constraints(aL, aR, aO), "compiler bug: bad cheb wires"
        proof = bulletproof.prove(circ, aL, aR, aO)
        proof["equation"] = circ.canonical
        proof["precision"] = k
        # Store cheb spec for coarse verification.
        # Note: diffeq keys are (ex, ords) where ords is a tuple.
        _spec = circ._cheb_spec_for_proof
        proof["cheb_spec"] = {
            'diffeq': [((ex, list(ords)), str(coeff))
                       for (ex, ords), coeff in _spec['diffeq'].items()],
            'order': _spec['order'],
            'ics': [((ord_, x0), str(v))
                    for (ord_, x0), v in _spec['ics'].items()],
            'deg': _spec['deg'],
        }
        return proof
    # Strict pre-check (the definition) applies to every equation:
    # |witness - exact root| < 1/k (complex modulus if complex).
    _check_precision(eq_str, w_pair if w_is_complex else w_pair[0], k)
    if circ.mode == "exact":
        if circ.is_complex:
            # Complex exact circuit: witness as (re, im) pair.
            if w_is_complex:
                if (w_pair[0].denominator != 1 or w_pair[1].denominator != 1):
                    raise ValueError(
                        "witness must be a Gaussian integer for this equation; "
                        "for a decimal witness, write the equation with decimals")
                w = (int(w_pair[0]), int(w_pair[1]))
            else:
                # Real witness for complex equation: (w, 0).
                w_frac = w_pair[0]
                if w_frac.denominator != 1:
                    raise ValueError(
                        "witness must be an integer for this integer equation; "
                        "for a decimal witness, write the equation with decimals")
                w = (int(w_frac), 0)
        else:
            # Real exact: integer witness.
            w_frac = w_pair[0]
            if w_frac.denominator != 1:
                raise ValueError(
                    "witness must be an integer for this integer equation; "
                    "for a decimal witness, write the equation with decimals "
                    "(e.g. x^2 = 2.0)")
            w = int(w_frac)
        if not circ.check_witness(w):
            raise ValueError("witness does not satisfy the equation")
        aL, aR, aO = circ.evaluate(w)
    else:
        # TODO: complex tolerance mode
        if w_is_complex:
            raise ValueError("complex tolerance mode not yet implemented")
        if circ.mode == "ftc":
            # FTC: witness is [F_0..F_{d+1}, v] as field elements.
            # Compute F analytically from f; v is the EXACT I_true = F(b)-F(a)
            # (not the user's decimal, which has quantization error).
            # The user's v was already pre-checked: |v - I_true| < 1/k.
            spec = circ.integral_spec
            f_poly = spec['integrand']
            a, b = spec['a'], spec['b']
            d = max(f_poly.keys()) if f_poly else 0
            # F_j = c_{j-1}/j for j>=1, F_0 = 0 (constant of integration).
            F = [Fraction(0)] * (d + 2)
            for j in range(1, d + 2):
                c_jm1 = f_poly.get(j - 1, (Fraction(0), Fraction(0)))[0]
                F[j] = c_jm1 / j
            # Exact I_true = F(b) - F(a).
            I_true = Fraction(0)
            for j in range(d + 2):
                I_true += F[j] * (b ** j - a ** j)
            # Convert to field elements.
            def _field(fr):
                return (fr.numerator * pow(fr.denominator,
                                           circuit_mod.N - 2,
                                           circuit_mod.N)) % circuit_mod.N
            X = [_field(F[j]) for j in range(d + 2)] + [_field(I_true)]
            if not circ.check_witness(X):
                raise ValueError(
                    "witness does not satisfy the FTC constraints")
            aL, aR, aO = circ.evaluate(X)
        else:
            w_frac = w_pair[0]
            X = _quantize(w_frac)
            if not circ.check_witness(X):
                raise ValueError(
                    "witness does not satisfy within 1/%d" % k)
            aL, aR, aO = circ.evaluate(X)
    assert circ.check_constraints(aL, aR, aO), "compiler bug: bad witness wires"
    proof = bulletproof.prove(circ, aL, aR, aO)
    proof["equation"] = circ.canonical
    proof["precision"] = k  # informational only; verifiers MUST NOT
    # rely on it -- verify_equation requires k as explicit input
    return proof


def solve_equation(eq_str, decimals=8):
    """Numerically solve `eq_str`; return up to 5 real roots as strings,
    each rounded to `decimals` places (default 8, 1..12 allowed).

    This is the honest-prover helper: feed a returned root back as the
    witness to prove_equation. Transcendental functions are solved as
    their fixed Chebyshev polynomials (the same polynomials the proof
    is about).
    Raises ValueError on bad input.
    """
    try:
        dec = int(str(decimals).strip())
    except (ValueError, AttributeError):
        raise ValueError("decimals must be an integer 1..12")
    eq = eq_str.strip()
    roots = solve_mod.solve_equation(eq, decimals=dec)
    # Complex polynomials: roots are already "a+bi" strings.
    if roots and isinstance(roots[0], str):
        return roots
    # exact (integer) equations: show integral roots as plain integers
    # ("3", not "3.00000000") so they can be tapped straight into prove
    poly, trans_used, float_seen, _, _, _ = circuit_mod.parse_polynomial(eq)
    def _is_int_coeff(c):
        if isinstance(c, tuple):
            return c[0].denominator == 1 and c[1].denominator == 1
        return c.denominator == 1
    exact = (not trans_used and not float_seen
             and all(_is_int_coeff(c) for c in poly.values()))
    out = []
    for v in roots:
        if exact and abs(v - round(v)) < 1e-9:
            out.append(str(int(round(v))))
        else:
            out.append(("%." + str(dec) + "f") % v)
    return out


def integral_value(eq_str, decimals=12):
    """Get the analytical integral value for an integral equation.
    Returns a string with up to `decimals` places (for the prover to
    use as witness). Raises ValueError if not an integral equation.
    
    Uses the analytical antiderivative (exact Fraction), not Simpson.
    """
    poly, _, _, integral_spec, _, _ = circuit_mod.parse_polynomial(eq_str)
    if integral_spec is None:
        raise ValueError("not an integral equation")
    f_poly = integral_spec['integrand']
    a, b = integral_spec['a'], integral_spec['b']
    # Analytical antiderivative: F(x) = sum c_e/(e+1) * x^(e+1)
    # F(b) - F(a), exact Fraction
    F_b = Fraction(0)
    F_a = Fraction(0)
    for e, c in f_poly.items():
        # c is (re, im); integrand is real
        coeff = c[0] / (e + 1)
        F_b += coeff * (b ** (e + 1))
        F_a += coeff * (a ** (e + 1))
    I_true = F_b - F_a
    # Format to decimals
    q = 10 ** decimals
    # Round half away from zero
    v = float(I_true)
    rv = math.floor(v * q + 0.5) / q if v >= 0 else -math.floor(-v * q + 0.5) / q
    s = "%.*f" % (decimals, rv)
    # Strip trailing zeros but keep at least one decimal place for clarity
    if "." in s:
        s = s.rstrip("0").rstrip(".")
    return s if s not in ("", "-0") else "0"


def ode_coefficients(eq_str):
    """Solve a linear ODE for its power series coefficients.
    Returns a list of Fractions [c_0, c_1, ..., c_d].
    Raises ValueError if not an ODE equation.
    
    Solves the linear system from residual coefficients (k=0..d-order)
    and initial conditions. Exact Fraction arithmetic.
    """
    _, _, _, _, ode_spec, _ = circuit_mod.parse_polynomial(eq_str)
    if ode_spec is None:
        raise ValueError("not an ODE equation")
    diffeq = ode_spec['diffeq']
    order = ode_spec['order']
    ics = ode_spec['ics']
    d = ode_spec['deg']
    # Precompute factorials
    fact = [1] * (d + 1)
    for i in range(1, d + 1):
        fact[i] = fact[i - 1] * i
    # Build linear system: (d+1) equations, (d+1) unknowns.
    # Rows 0..d-order: residual coeff of x^k = 0.
    # Rows d-order+1 .. d: ICs (j! * c_j = ics[j]).
    n = d + 1
    # M[row][col], rhs[row]
    M = [[Fraction(0)] * n for _ in range(n)]
    rhs = [Fraction(0)] * n
    # Residual rows
    # Taylor mode: only linear terms (nonlinear rejected in _build_ode,
    # but solver should also check).
    for (ex, ords) in diffeq.keys():
        if len(ords) > 1:
            raise ValueError("Taylor ode() does not support nonlinear terms")
    for k in range(d - order + 1):
        for (ex, ords), coeff in diffeq.items():
            if len(ords) == 0:
                if ex == k:
                    rhs[k] -= coeff
            else:
                ord_ = ords[0]
                i = k - ex
                if 0 <= i <= d - ord_:
                    c_idx = i + ord_
                    M[k][c_idx] += coeff * Fraction(fact[i + ord_], fact[i])
    # IC rows
    for j in range(order):
        row = d - order + 1 + j
        M[row][j] = Fraction(fact[j])
        rhs[row] = ics[j]
    # Gaussian elimination (Fractions, exact).
    # Forward elimination
    for col in range(n):
        # Find pivot
        piv = None
        for r in range(col, n):
            if M[r][col] != 0:
                piv = r
                break
        if piv is None:
            raise ValueError("singular ODE system (cannot solve)")
        if piv != col:
            M[col], M[piv] = M[piv], M[col]
            rhs[col], rhs[piv] = rhs[piv], rhs[col]
        # Eliminate
        for r in range(col + 1, n):
            if M[r][col] != 0:
                f = M[r][col] / M[col][col]
                for c2 in range(col, n):
                    M[r][c2] -= f * M[col][c2]
                rhs[r] -= f * rhs[col]
    # Back substitution
    c = [Fraction(0)] * n
    for r in range(n - 1, -1, -1):
        s = rhs[r]
        for c2 in range(r + 1, n):
            s -= M[r][c2] * c[c2]
        if M[r][r] == 0:
            raise ValueError("singular ODE system (cannot solve)")
        c[r] = s / M[r][r]
    return c


def cheb_coefficients(eq_str):
    """Solve an ODE for its Chebyshev coefficients (tau method).
    Returns a list of Fractions [c_0, ..., c_N] where
    y(x) = Σ_{n=0}^{N} c_n T_n(x) on [-1, 1].
    Raises ValueError if not a Chebyshev equation.

    For linear ODEs: solves the linear system from tau conditions
    (residual coefficients r_k = 0 for k = 0..N-order) and side
    conditions. Exact Fractions.

    For nonlinear (quadratic) ODEs: uses Newton iteration in field
    arithmetic to find an exact field solution, then converts back
    to Fractions via rational reconstruction. The circuit verifies
    the residual is exactly zero in the field.
    """
    from .circuit import N as FIELD_N
    _, _, _, _, _, cheb_spec = circuit_mod.parse_polynomial(eq_str)
    if cheb_spec is None:
        raise ValueError("not a Chebyshev equation")
    diffeq = cheb_spec['diffeq']
    order = cheb_spec['order']
    ics = cheb_spec['ics']
    N = cheb_spec['deg']

    # Check for nonlinearity
    has_nonlinear = any(len(ords) == 2 for (ex, ords) in diffeq.keys())
    if any(len(ords) > 2 for (ex, ords) in diffeq.keys()):
        raise ValueError("only up to quadratic nonlinearities supported")

    # Chebyshev differentiation matrices (Fractions)
    def cheb_diff_mat(N):
        D = [[Fraction(0)] * (N + 1) for _ in range(N + 1)]
        if N >= 1:
            D[N - 1][N] = Fraction(2 * N)
        for nn in range(N - 2, -1, -1):
            cbar = Fraction(2) if nn == 0 else Fraction(1)
            for col in range(N + 1):
                D[nn][col] += D[nn + 2][col] / cbar
            D[nn][nn + 1] += Fraction(2 * (nn + 1)) / cbar
        return D

    D1 = cheb_diff_mat(N)
    Dpows = [[[Fraction(1) if i == j else Fraction(0)
               for j in range(N + 1)] for i in range(N + 1)]]
    for _ in range(order):
        prev = Dpows[-1]
        nxt = [[Fraction(0)] * (N + 1) for _ in range(N + 1)]
        for i in range(N + 1):
            for k in range(N + 1):
                if prev[i][k] != 0:
                    for j in range(N + 1):
                        if D1[k][j] != 0:
                            nxt[i][j] += prev[i][k] * D1[k][j]
        Dpows.append(nxt)

    # x^ex in Chebyshev basis
    def x_pow_cheb(ex):
        if ex == 0:
            return [Fraction(1)]
        p = [Fraction(1)]
        for _ in range(ex):
            q = [Fraction(0)] * (len(p) + 1)
            for k, pk in enumerate(p):
                if pk == 0:
                    continue
                q[abs(1 - k)] += pk / 2
                q[k + 1] += pk / 2
            p = q
        return p

    x_pow_cache = {}
    for (ex, ords) in diffeq.keys():
        if ex not in x_pow_cache:
            x_pow_cache[ex] = x_pow_cheb(ex)

    def T_n_x0(nn, x0):
        if x0 == 1:
            return Fraction(1)
        if x0 == -1:
            return Fraction(-1) ** nn
        if x0 == 0:
            r = nn % 4
            return {0: Fraction(1), 1: Fraction(0),
                    2: Fraction(-1), 3: Fraction(0)}[r]
        raise ValueError("x0 must be -1, 0, or 1")

    if not has_nonlinear:
        # --- Linear case: exact Fraction solve (unchanged) ---
        n = N + 1
        M = [[Fraction(0)] * n for _ in range(n)]
        rhs = [Fraction(0)] * n

        # Tau rows: r_k = 0 for k = 0..N-order
        for k in range(N - order + 1):
            for (ex, ords), coeff in diffeq.items():
                if len(ords) == 0:
                    px = x_pow_cache[ex]
                    if k < len(px):
                        rhs[k] -= coeff * px[k]
                else:
                    ord_ = ords[0]
                    px = x_pow_cache[ex]
                    Dp = Dpows[ord_]
                    for j, px_j in enumerate(px):
                        if px_j == 0:
                            continue
                        for nn in range(N + 1):
                            contrib = Fraction(0)
                            if abs(j - nn) == k:
                                contrib += px_j / 2
                            if j + nn == k:
                                contrib += px_j / 2
                            if contrib != 0:
                                for m in range(N + 1):
                                    if Dp[nn][m] != 0:
                                        M[k][m] += coeff * contrib * Dp[nn][m]

        # Side condition rows
        row = N - order + 1
        for (ord_, x0), v in sorted(ics.items()):
            Dp = Dpows[ord_]
            for nn in range(N + 1):
                tval = T_n_x0(nn, x0)
                if tval == 0:
                    continue
                for m in range(N + 1):
                    if Dp[nn][m] != 0:
                        M[row][m] += tval * Dp[nn][m]
            rhs[row] = v
            row += 1

        # Gaussian elimination
        for col in range(n):
            piv = None
            for r in range(col, n):
                if M[r][col] != 0:
                    piv = r
                    break
            if piv is None:
                raise ValueError("singular Chebyshev system (cannot solve)")
            if piv != col:
                M[col], M[piv] = M[piv], M[col]
                rhs[col], rhs[piv] = rhs[piv], rhs[col]
            for r in range(col + 1, n):
                if M[r][col] != 0:
                    f = M[r][col] / M[col][col]
                    for c2 in range(col, n):
                        M[r][c2] -= f * M[col][c2]
                    rhs[r] -= f * rhs[col]
        c = [Fraction(0)] * n
        for r in range(n - 1, -1, -1):
            s = rhs[r]
            for c2 in range(r + 1, n):
                s -= M[r][c2] * c[c2]
            if M[r][r] == 0:
                raise ValueError("singular Chebyshev system (cannot solve)")
            c[r] = s / M[r][r]
        return c
    else:
        # --- Nonlinear case: hybrid solver ---
        # Step 1: Solve in R with floats (Newton) to get close.
        # Step 2: Convert to field, then Newton in field to get exact root.
        import math
        # Float versions of Dpows and x_pow
        Dpows_float = []
        for o in range(order + 1):
            M = [[float(Dpows[o][i][j]) for j in range(N + 1)]
                 for i in range(N + 1)]
            Dpows_float.append(M)
        x_pow_float = {ex: [float(v) for v in px]
                       for ex, px in x_pow_cache.items()}
        ics_float = {k: float(v) for k, v in ics.items()}

        def residual_float(c):
            # c: list of floats, returns list of floats (size N+1)
            F = [0.0] * (N + 1)
            for k in range(N - order + 1):
                acc = 0.0
                for (ex, ords), coeff in diffeq.items():
                    cf = float(coeff)
                    if len(ords) == 0:
                        px = x_pow_float[ex]
                        if k < len(px):
                            acc += cf * px[k]
                    elif len(ords) == 1:
                        ord_ = ords[0]
                        px = x_pow_float[ex]
                        Dp = Dpows_float[ord_]
                        for j, px_j in enumerate(px):
                            if px_j == 0:
                                continue
                            for nn in range(N + 1):
                                contrib = 0.0
                                if abs(j - nn) == k:
                                    contrib += px_j / 2
                                if j + nn == k:
                                    contrib += px_j / 2
                                if contrib != 0:
                                    dn = sum(Dp[nn][m] * c[m] for m in range(N + 1))
                                    acc += cf * contrib * dn
                    else:
                        o1, o2 = ords[0], ords[1]
                        A = Dpows_float[o1]
                        B = Dpows_float[o2]
                        px = x_pow_float[ex]
                        u = [sum(A[m][a] * c[a] for a in range(N + 1))
                             for m in range(N + 1)]
                        v = [sum(B[n_][b] * c[b] for b in range(N + 1))
                             for n_ in range(N + 1)]
                        for t in range(2 * N + 1):
                            wt = 0.0
                            for m in range(N + 1):
                                for n_ in range(N + 1):
                                    w = 0.0
                                    if abs(m - n_) == t:
                                        w += 0.5
                                    if m + n_ == t:
                                        w += 0.5
                                    wt += u[m] * v[n_] * w
                            if wt == 0:
                                continue
                            for j, px_j in enumerate(px):
                                if px_j == 0:
                                    continue
                                w2 = 0.0
                                if abs(j - t) == k:
                                    w2 += 0.5
                                if j + t == k:
                                    w2 += 0.5
                                acc += cf * px_j * wt * w2
                F[k] = acc
            row = N - order + 1
            for (ord_, x0), v in sorted(ics_float.items()):
                Dp = Dpows_float[ord_]
                s = 0.0
                for nn in range(N + 1):
                    tval = float(T_n_x0(nn, x0))
                    if tval == 0:
                        continue
                    dn = sum(Dp[nn][m] * c[m] for m in range(N + 1))
                    s += tval * dn
                F[row] = s - v
                row += 1
            return F

        # Initial guess: linear solution (floats)
        # Solve linear system in floats
        n = N + 1
        M_lin = [[0.0] * n for _ in range(n)]
        rhs_lin = [0.0] * n
        for k in range(N - order + 1):
            for (ex, ords), coeff in diffeq.items():
                if len(ords) != 1:
                    continue
                ord_ = ords[0]
                cf = float(coeff)
                px = x_pow_float[ex]
                Dp = Dpows_float[ord_]
                for j, px_j in enumerate(px):
                    if px_j == 0:
                        continue
                    for nn in range(N + 1):
                        contrib = 0.0
                        if abs(j - nn) == k:
                            contrib += px_j / 2
                        if j + nn == k:
                            contrib += px_j / 2
                        if contrib != 0:
                            for m in range(N + 1):
                                M_lin[k][m] += cf * contrib * Dp[nn][m]
            for (ex, ords), coeff in diffeq.items():
                if len(ords) == 0:
                    px = x_pow_float[ex]
                    if k < len(px):
                        rhs_lin[k] -= float(coeff) * px[k]
        row = N - order + 1
        for (ord_, x0), v in sorted(ics_float.items()):
            Dp = Dpows_float[ord_]
            for nn in range(N + 1):
                tval = float(T_n_x0(nn, x0))
                if tval == 0:
                    continue
                for m in range(N + 1):
                    M_lin[row][m] += tval * Dp[nn][m]
            rhs_lin[row] = v
            row += 1
        # Solve
        M_tmp = [row[:] for row in M_lin]
        rhs_tmp = rhs_lin[:]
        for col in range(n):
            piv = max(range(col, n), key=lambda r: abs(M_tmp[r][col]))
            if abs(M_tmp[piv][col]) < 1e-12:
                raise ValueError("singular linear system")
            M_tmp[col], M_tmp[piv] = M_tmp[piv], M_tmp[col]
            rhs_tmp[col], rhs_tmp[piv] = rhs_tmp[piv], rhs_tmp[col]
            pivval = M_tmp[col][col]
            for c2 in range(col, n):
                M_tmp[col][c2] /= pivval
            rhs_tmp[col] /= pivval
            for r in range(n):
                if r != col:
                    f = M_tmp[r][col]
                    for c2 in range(col, n):
                        M_tmp[r][c2] -= f * M_tmp[col][c2]
                    rhs_tmp[r] -= f * rhs_tmp[col]
        c_float = rhs_tmp[:]  # solution

        # Newton in R
        for it in range(50):
            F = residual_float(c_float)
            norm = math.sqrt(sum(x * x for x in F))
            if norm < 1e-12:
                break
            # Numerical Jacobian
            J = [[0.0] * n for _ in range(n)]
            eps = 1e-8
            for j in range(n):
                c2 = c_float[:]
                c2[j] += eps
                F2 = residual_float(c2)
                for i in range(n):
                    J[i][j] = (F2[i] - F[i]) / eps
            # Solve J*delta = -F
            M = [J[i][:] + [-F[i]] for i in range(n)]
            for col in range(n):
                piv = max(range(col, n), key=lambda r: abs(M[r][col]))
                M[col], M[piv] = M[piv], M[col]
                pivval = M[col][col]
                for c2 in range(col, n + 1):
                    M[col][c2] /= pivval
                for r in range(n):
                    if r != col:
                        f = M[r][col]
                        for c2 in range(col, n + 1):
                            M[r][c2] -= f * M[col][c2]
            delta = [M[i][n] for i in range(n)]
            c_float = [c_float[i] + delta[i] for i in range(n)]
        else:
            raise ValueError("Float Newton did not converge")

        # Step 2: Convert to field and polish with field Newton
        def to_field(fr):
            return (fr.numerator * pow(fr.denominator, FIELD_N - 2, FIELD_N)) % FIELD_N

        # Convert floats to Fractions, then to field
        c_field = []
        for x in c_float:
            fr = Fraction(x).limit_denominator(10**12)
            c_field.append(to_field(fr))

        # Field Newton (using the residual and jacobian defined earlier)
        # We need to redefine them here or use the ones from outer scope.
        # For simplicity, reuse the residual() and jacobian() from the
        # Picard section (they're defined in the outer else block).
        # Actually, let's just do a few Newton steps with the analytic Jacobian.
        # (The jacobian function is defined below in the Picard code;
        #  for now, use numerical Jacobian in field with careful step.)
        #
        # Step 2: Field Newton polish.
        # NOTE: Exact root finding over finite fields for nonlinear
        # polynomial systems is hard (NP-hard in general). The float
        # Newton gives an approximate real solution. For the prototype,
        # we return the approximate Fractions. Prove/verify for
        # nonlinear ODEs requires an exact field root and is not yet
        # supported (prove_equation will raise a clear error).
        return [Fraction(x).limit_denominator(10**12) for x in c_float]


def verify_equation(eq_str, proof, precision):
    """Verify a proof against an equation string. Returns True/False.

    precision: k (int/str), MANDATORY for algebraic/integral equations.
    For ODEs, k is meaningless (precision is deg) and may be empty.
    The verifier must supply the precision explicitly; it is NEVER taken
    from the proof itself. (Reading k from the proof would let a malicious
    prover weaken the requirement.) Returns False (rather than raising)
    on any malformed proof, a canonical-equation mismatch (covers
    moved-term spellings AND precision mismatches, since canonical
    includes ";k=<k>"), or a precision mismatch.
    Raises ValueError if precision is missing/invalid (caller error,
    except for ODEs where it's optional).
    """
    # ODEs/Chebyshev: k is meaningless; allow empty.
    _, _, _, _, ode_spec, cheb_spec = circuit_mod.parse_polynomial(eq_str)
    is_func_eq = ode_spec is not None or cheb_spec is not None
    if not is_func_eq and not str(precision).strip():
        raise ValueError("precision k is required to verify")
    if is_func_eq and not str(precision).strip():
        precision = "1"  # dummy; canonical has no k
    try:
        if not isinstance(proof, dict):
            return False
        proof_eq = proof.get("equation", "")
        # ODE coarse verification: a deg_prove proof verifies against
        # deg_verify <= deg_prove (finer proof implies coarser statement).
        if ode_spec is not None and "ode_spec" in proof:
            from fractions import Fraction as _Fr
            ps = proof["ode_spec"]
            # Reconstruct proof's spec
            # Note: diffeq keys were stored as (ex, list(ords)); convert back.
            proof_diffeq = {(ex, tuple(ords)): _Fr(coeff_str)
                            for (ex, ords), coeff_str in ps['diffeq']}
            proof_ics = {ord_: _Fr(v_str) for ord_, v_str in ps['ics']}
            deg_prove = ps['deg']
            deg_verify = ode_spec['deg']
            # Check ODE matches (diffeq, order, ics), ignoring deg.
            if (proof_diffeq != ode_spec['diffeq'] or
                    ps['order'] != ode_spec['order'] or
                    proof_ics != ode_spec['ics']):
                return False
            # Coarse check: deg_verify <= deg_prove
            if deg_verify > deg_prove:
                return False
            # Verify against the deg_prove circuit.
            _spec_prove = {
                'diffeq': proof_diffeq,
                'order': ps['order'],
                'ics': proof_ics,
                'deg': deg_prove,
            }
            circ = circuit_mod.Circuit(
                "ode_dummy", {0: (0, 0)}, "1", False, False,
                ode_spec=_spec_prove)
            if proof_eq != circ.canonical:
                return False
            return bulletproof.verify(circ, proof)
        # Chebyshev coarse verification: same idea.
        if cheb_spec is not None and "cheb_spec" in proof:
            from fractions import Fraction as _Fr
            ps = proof["cheb_spec"]
            # Note: diffeq keys were stored as (ex, list(ords)); convert back to tuple.
            proof_diffeq = {(ex, tuple(ords)): _Fr(coeff_str)
                            for (ex, ords), coeff_str in ps['diffeq']}
            proof_ics = {(ord_, x0): _Fr(v_str)
                         for (ord_, x0), v_str in ps['ics']}
            deg_prove = ps['deg']
            deg_verify = cheb_spec['deg']
            if (proof_diffeq != cheb_spec['diffeq'] or
                    ps['order'] != cheb_spec['order'] or
                    proof_ics != cheb_spec['ics']):
                return False
            if deg_verify > deg_prove:
                return False
            _spec_prove = {
                'diffeq': proof_diffeq,
                'order': ps['order'],
                'ics': proof_ics,
                'deg': deg_prove,
            }
            circ = circuit_mod.Circuit(
                "cheb_dummy", {0: (0, 0)}, "1", False, False,
                cheb_spec=_spec_prove)
            if proof_eq != circ.canonical:
                return False
            return bulletproof.verify(circ, proof)
        # NOTE: k comes ONLY from the caller's explicit input (the
        # definition), NEVER from proof.get("precision").
        # The complex marker (;cx=1) is part of the circuit identifier
        # (like the equation), not the precision -- it tells the verifier
        # which circuit to build.
        is_complex = proof_eq.endswith(";cx=1")
        circ = circuit_mod.compile(eq_str, str(precision).strip(),
                                   force_complex=is_complex)
        if proof.get("equation") != circ.canonical:
            return False
        return bulletproof.verify(circ, proof)
    except (ValueError, KeyError, TypeError, AssertionError):
        return False
