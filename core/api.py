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
    poly, _, _ = circuit_mod.parse_polynomial(eq_str)
    if not poly or all(circuit_mod._ciszero(c) for c in poly.values()):
        return 0.0
    eps = 1.0 / k
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
    Raises ValueError if the equation is malformed, the precision is
    missing/invalid, or the witness does not satisfy (within precision).
    """
    if not str(precision).strip():
        raise ValueError("precision k is required for every equation")
    # Parse witness FIRST to detect complex (determines circuit type).
    w_pair, w_is_complex = _parse_witness_complex(witness)
    circ = circuit_mod.compile(eq_str, precision, force_complex=w_is_complex)
    k = circ.precision  # int, from the user-supplied definition (1/k)
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
    poly, trans_used, float_seen = circuit_mod.parse_polynomial(eq)
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


def verify_equation(eq_str, proof, precision):
    """Verify a proof against an equation string. Returns True/False.

    precision: k (int/str), MANDATORY. The verifier must supply the
    precision explicitly; it is NEVER taken from the proof itself.
    (Reading k from the proof would let a malicious prover weaken the
    requirement.) Returns False (rather than raising) on any malformed
    proof, a canonical-equation mismatch (covers moved-term spellings
    AND precision mismatches, since canonical includes ";k=<k>"),
    or a precision mismatch.
    Raises ValueError if precision is missing/invalid (caller error).
    """
    if not str(precision).strip():
        raise ValueError("precision k is required to verify")
    try:
        if not isinstance(proof, dict):
            return False
        # NOTE: k comes ONLY from the caller's explicit input (the
        # definition), NEVER from proof.get("precision").
        # The complex marker (;cx=1) is part of the circuit identifier
        # (like the equation), not the precision -- it tells the verifier
        # which circuit to build.
        proof_eq = proof.get("equation", "")
        is_complex = proof_eq.endswith(";cx=1")
        circ = circuit_mod.compile(eq_str, str(precision).strip(),
                                   force_complex=is_complex)
        if proof.get("equation") != circ.canonical:
            return False
        return bulletproof.verify(circ, proof)
    except (ValueError, KeyError, TypeError, AssertionError):
        return False
