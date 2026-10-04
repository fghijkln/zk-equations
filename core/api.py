"""Public API: prove/verify knowledge of a solution to an arithmetic equation.

    proof = prove_equation("x^3 + 2*x + 5 = 38", 3)
    verify_equation("x^3 + 2*x + 5 = 38", proof)   # True, reveals nothing about x

The proof is a plain dict of hex strings (JSON-serializable).
"""

from . import bulletproof
from . import circuit as circuit_mod


def prove_equation(eq_str, witness):
    """Prove knowledge of an integer `witness` satisfying `eq_str`.

    Raises ValueError if the equation is malformed or the witness does
    not satisfy it.
    """
    circ = circuit_mod.compile(eq_str)
    if not circ.check_witness(witness):
        raise ValueError("witness does not satisfy the equation")
    aL, aR, aO = circ.evaluate(witness)
    assert circ.check_constraints(aL, aR, aO), "compiler bug: bad witness wires"
    proof = bulletproof.prove(circ, aL, aR, aO)
    proof["equation"] = eq_str
    return proof


def verify_equation(eq_str, proof):
    """Verify a proof against an equation string. Returns True/False.

    Returns False (rather than raising) on any malformed proof, and if
    the proof's embedded equation differs from `eq_str`.
    """
    try:
        if proof.get("equation") != eq_str:
            return False
        circ = circuit_mod.compile(eq_str)
        return bulletproof.verify(circ, proof)
    except (ValueError, KeyError, TypeError):
        return False
