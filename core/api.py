"""Public API: prove/verify knowledge of a solution to an equation.

    proof = prove_equation("x^3 + 2*x + 5 = 38", 3)
    verify_equation("x^3 + 2*x + 5 = 38", proof)   # True, reveals nothing about x

    # transcendental / decimal equations need a precision denominator k.
    # Strict definition: the witness must be within 1/(k*10^4) of an
    # exact root; otherwise prove_equation refuses (raises ValueError).
    # Use solve_equation to get accurate roots (8 decimals by default).
    proof = prove_equation("sin(x) = 0.5", "0.52359878", "1000")
    verify_equation("sin(x) = 0.5", proof, "1000")  # True

The proof is a plain dict of hex strings (JSON-serializable). Proofs are
bound to the *canonical* equation (polynomial normal form): term-moved
spellings of the same equation verify each other's proofs, genuinely
different equations do not.
"""

from fractions import Fraction

from . import bulletproof
from . import circuit as circuit_mod
from . import solve as solve_mod

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


def _quantize(fr):
    """Round a Fraction to the circuit's fixed-point scale -> int X."""
    num, den = (fr * M_DEC).numerator, (fr * M_DEC).denominator
    # round half away from zero
    if num >= 0:
        return (2 * num + den) // (2 * den)
    return -((2 * -num + den) // (2 * den))


def _check_precision(eq_str, w_frac, k):
    """Enforce the strict precision definition: the witness must be
    within 1/(k*10^4) of an exact root. Raises ValueError otherwise."""
    eps = 1.0 / (k * 10 ** 4)
    roots = solve_mod.all_roots(eq_str)
    w = float(w_frac)
    best = min((abs(w - r) for r in roots), default=None)
    if best is None or best >= eps:
        detail = ("no real root found"
                  if best is None else "distance %.3g" % best)
        raise ValueError(
            "witness is not within precision 1/(%d*10^4)=%.3g of an exact "
            "root (%s); use the solver to get an accurate root"
            % (k, eps, detail))
    return best


def prove_equation(eq_str, witness, precision=""):
    """Prove knowledge of `witness` satisfying `eq_str`.

    witness: int (exact mode) or decimal string (tolerance mode, up to
             12 decimal places; an int is also accepted).
    precision: k (int/str), required in tolerance mode. Strict
               definition: the witness must be within 1/(k*10^4) of an
               exact root, else ValueError ("refuse"). Ignored in exact
               mode.
    Raises ValueError if the equation is malformed, the precision is
    missing/invalid, or the witness does not satisfy (within precision).
    """
    circ = circuit_mod.compile(eq_str, precision)
    if circ.mode == "exact":
        if isinstance(witness, int):
            w = witness
        else:
            s = str(witness).strip()
            try:
                fr = Fraction(s)
            except (ValueError, ZeroDivisionError):
                fr = None
            if fr is not None and fr.denominator == 1:
                # "3", "3.0", "3.0000000000" (e.g. tapped from the
                # solver) are all the integer 3
                w = int(fr)
            else:
                raise ValueError(
                    "witness must be an integer for this exact equation; "
                    "for a decimal witness, write the equation with decimals "
                    "(e.g. x^5 = 1.0) and set precision k")
        if not circ.check_witness(w):
            raise ValueError("witness does not satisfy the equation")
        aL, aR, aO = circ.evaluate(w)
    else:
        w_frac = _parse_witness_full(witness)
        _check_precision(eq_str, w_frac, circ.precision)
        X = _quantize(w_frac)
        if not circ.check_witness(X):
            raise ValueError(
                "witness does not satisfy within 1/%d" % circ.precision)
        aL, aR, aO = circ.evaluate(X)
    assert circ.check_constraints(aL, aR, aO), "compiler bug: bad witness wires"
    proof = bulletproof.prove(circ, aL, aR, aO)
    proof["equation"] = circ.canonical
    proof["precision"] = circ.precision  # int k, or None in exact mode
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
    # exact (integer) equations: show integral roots as plain integers
    # ("3", not "3.00000000") so they can be tapped straight into prove
    poly, trans_used, float_seen = circuit_mod.parse_polynomial(eq)
    exact = (not trans_used and not float_seen
             and all(c.denominator == 1 for c in poly.values()))
    out = []
    for v in roots:
        if exact and abs(v - round(v)) < 1e-9:
            out.append(str(int(round(v))))
        else:
            out.append(("%." + str(dec) + "f") % v)
    return out


def verify_equation(eq_str, proof, precision=""):
    """Verify a proof against an equation string. Returns True/False.

    precision: k for tolerance equations; if empty, it is taken from the
    proof itself. Returns False (rather than raising) on any malformed
    proof, a canonical-equation mismatch (covers moved-term spellings),
    or a precision mismatch.
    """
    try:
        if not isinstance(proof, dict):
            return False
        k = str(precision).strip() or proof.get("precision")
        circ = circuit_mod.compile(eq_str, k if k is not None else "")
        if proof.get("equation") != circ.canonical:
            return False
        return bulletproof.verify(circ, proof)
    except (ValueError, KeyError, TypeError, AssertionError):
        return False
