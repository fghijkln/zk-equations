"""Public API: prove/verify knowledge of a solution to an equation.

    proof = prove_equation("x^3 + 2*x + 5 = 38", 3)
    verify_equation("x^3 + 2*x + 5 = 38", proof)   # True, reveals nothing about x

    # transcendental / decimal equations need a precision denominator k
    # (final precision 1/k):
    proof = prove_equation("sin(x) = 0.5", "0.5236", "1000")
    verify_equation("sin(x) = 0.5", proof, "1000")  # True

The proof is a plain dict of hex strings (JSON-serializable). Proofs are
bound to the *canonical* equation (polynomial normal form): term-moved
spellings of the same equation verify each other's proofs, genuinely
different equations do not.
"""

from fractions import Fraction

from . import bulletproof
from . import circuit as circuit_mod

M_DEC = circuit_mod.M_DEC


def _parse_decimal_witness(w):
    """Parse a decimal witness string to the scaled integer X = round(x*M)."""
    s = str(w).strip()
    try:
        fr = Fraction(s)
    except (ValueError, ZeroDivisionError):
        raise ValueError("witness must be a decimal number")
    num, den = (fr * M_DEC).numerator, (fr * M_DEC).denominator
    # round half away from zero
    if num >= 0:
        return (2 * num + den) // (2 * den)
    return -((2 * -num + den) // (2 * den))


def prove_equation(eq_str, witness, precision=""):
    """Prove knowledge of `witness` satisfying `eq_str`.

    witness: int (exact mode) or decimal string (tolerance mode; an int
             is also accepted and treated as an exact decimal).
    precision: k (int/str), required in tolerance mode; final precision
               is 1/k. Ignored in exact mode.
    Raises ValueError if the equation is malformed, the precision is
    missing/invalid, or the witness does not satisfy (within precision).
    """
    circ = circuit_mod.compile(eq_str, precision)
    if circ.mode == "exact":
        if isinstance(witness, int):
            w = witness
        else:
            try:
                w = int(str(witness).strip())
            except (ValueError, AttributeError):
                raise ValueError(
                    "witness must be an integer for this exact equation; "
                    "for a decimal witness, write the equation with decimals "
                    "(e.g. x^5 = 1.0) and set precision k")
        if not circ.check_witness(w):
            raise ValueError("witness does not satisfy the equation")
        aL, aR, aO = circ.evaluate(w)
    else:
        X = _parse_decimal_witness(witness)
        if not circ.check_witness(X):
            raise ValueError(
                "witness does not satisfy within precision 1/%d" % circ.precision)
        aL, aR, aO = circ.evaluate(X)
    assert circ.check_constraints(aL, aR, aO), "compiler bug: bad witness wires"
    proof = bulletproof.prove(circ, aL, aR, aO)
    proof["equation"] = circ.canonical
    proof["precision"] = circ.precision  # int k, or None in exact mode
    return proof


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
