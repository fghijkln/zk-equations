"""End-to-end demo: prove knowledge of x with x^3 + 2*x + 5 = 38 (x = 3).

NOTE on the equation: an early draft of this demo used "= 42", but
3^3 + 2*3 + 5 = 38, not 42 -- x = 3 is not a witness for "= 42"
(and the prover correctly refuses it; see the negative check below).
Run with:  python3 -m core.demo
"""

import json
import time

from . import api


def main():
    eq = "x^3 + 2*x + 5 = 38"
    x = 3
    print("equation :", eq)
    print("witness  : x = %d  (kept secret from the verifier)" % x)

    t0 = time.time()
    proof = api.prove_equation(eq, x)
    t_prove = time.time() - t0
    print("prove    : %.2fs" % t_prove)

    size = len(json.dumps(proof))
    print("proof size: %d bytes JSON" % size)

    t0 = time.time()
    ok = api.verify_equation(eq, proof)
    t_verify = time.time() - t0
    print("verify   : %.2fs -> %s" % (t_verify, ok))
    assert ok, "honest proof must verify"

    # A wrong witness must be refused at prove time...
    try:
        api.prove_equation("x^3 + 2*x + 5 = 42", 3)
        print("ERROR: bad witness was accepted!")
    except ValueError as e:
        print("bad witness correctly refused: %s" % e)

    # ...and a proof for one equation must not verify against another.
    ok2 = api.verify_equation("x^3 + 2*x + 5 = 42", proof)
    print("proof replayed against '= 42' verifies: %s (must be False)" % ok2)
    assert not ok2


if __name__ == "__main__":
    main()
