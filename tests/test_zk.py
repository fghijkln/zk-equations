"""Correctness + tamper-resistance tests for the ZK core.

Run:  python3 -m tests.test_zk
(plain asserts, no third-party test framework needed)
"""

import copy
import json
import time

from core import api, bulletproof
from core import circuit as circuit_mod

PASS = []
FAIL = []


def check(name, cond):
    (PASS if cond else FAIL).append(name)
    print(("PASS " if cond else "FAIL ") + name)


def expect_raises(name, fn):
    try:
        fn()
    except (ValueError, AssertionError):
        check(name, True)
    except Exception:
        check(name, False)
    else:
        check(name, False)


# ---------------- 1. honest proofs verify ----------------

CASES = [
    ("x^3 + 2*x + 5 = 38", 3),
    ("x + 5 = 8", 3),
    ("2x + 1 = 7", 3),            # implicit multiplication
    ("x^2 = 16", 4),
    ("x^2 = 16", -4),             # negative witness (mod n)
    ("x^2 = 0", 0),               # zero witness
    ("(x + 1)^2 = 16", 3),        # parentheses
    ("x**2 + 3*x = 18", 3),       # ** alias
    ("x^4 - 5*x^2 + 4 = 0", 2),   # degree 4
    ("x^5 = 32", 2),              # degree 5
    ("3 = 3", 0),                 # constant-true equation
    ("x - x = 0", 12345),         # degenerate: any witness works
    ("-x^2 = -4", 2),             # unary minus binds looser than ^
    ("-x^2 = -4", -2),
    ("-(x + 1)^2 = -9", 2),
]

for eq, w in CASES:
    t0 = time.time()
    proof = api.prove_equation(eq, w)
    tp = time.time() - t0
    t0 = time.time()
    ok = api.verify_equation(eq, proof)
    tv = time.time() - t0
    # JSON round-trip must preserve verifiability
    ok2 = api.verify_equation(eq, json.loads(json.dumps(proof)))
    check("honest %-28s x=%-6d verifies (prove %.2fs, verify %.2fs, json-ok)"
          % (eq, w, tp, tv), ok and ok2)

# ---------------- 2. wrong witnesses fail ----------------

for eq, w in [("x^3 + 2*x + 5 = 38", 4), ("x^2 = 16", 5), ("x + 5 = 8", 0)]:
    expect_raises("bad witness refused: %s x=%d" % (eq, w),
                  lambda: api.prove_equation(eq, w))

# unsatifiable constant equation is a compile error
expect_raises("constant-false equation rejected",
              lambda: api.prove_equation("3 = 4", 0))

# malformed equations are compile errors
for bad in ["x^2 = ", "= 5", "x + = 3", "2y + 1 = 5", "x^(-1) = 2", "x^2 == 4"]:
    expect_raises("malformed rejected: %r" % bad,
                  lambda: api.prove_equation(bad, 1))

# ---------------- 3. soundness at the protocol level ----------------
# Bypass the API's witness check: feed bad wires straight into the prover.
# A cheating prover must still fail verification.

eq, w = "x^3 + 2*x + 5 = 38", 3
circ = circuit_mod.compile(eq)
aL, aR, aO = circ.evaluate(w)

# 3a. break the Hadamard constraint (aO[0] is x^2 = 9 -> 10)
bL, bR, bO = list(aL), list(aR), list(aO)
bO[0] = (bO[0] + 1) % circuit_mod.N
pf = bulletproof.prove(circ, bL, bR, bO)
check("broken Hadamard gate -> verify False",
      bulletproof.verify(circ, pf) is False)

# 3b. break a linear constraint (use x=4's wires but keep equation)
aL4, aR4, aO4 = circ.evaluate(4)
pf = bulletproof.prove(circ, aL4, aR4, aO4)
check("wires for x=4 against '=38' circuit -> verify False",
      bulletproof.verify(circ, pf) is False)

# 3c. tamper every top-level proof field
good = api.prove_equation(eq, w)


def tampered(path, mutate):
    p = copy.deepcopy(good)
    node = p
    for key in path[:-1]:
        node = node[key]
    node[path[-1]] = mutate(node[path[-1]])
    return p


def flip_hex(h):
    return ("0" if h[0] != "0" else "1") + h[1:]


tamper_cases = [
    (["AI"], flip_hex),
    (["AO"], flip_hex),
    (["S"], flip_hex),
    (["T", "1"], flip_hex),
    (["T", "3"], flip_hex),
    (["tau_x"], flip_hex),
    (["mu"], flip_hex),
    (["t_hat"], flip_hex),
    (["ipa", "a"], flip_hex),
    (["ipa", "b"], flip_hex),
    (["ipa", "Ls", 0], flip_hex),
    (["ipa", "Rs", 0], flip_hex),
]
for path, mut in tamper_cases:
    p = tampered(path, mut)
    check("tamper %s -> verify False" % ".".join(map(str, path)),
          api.verify_equation(eq, p) is False)

# garbage / truncated proofs must not raise, just return False
check("empty dict -> False", api.verify_equation(eq, {}) is False)
check("garbage hex -> False",
      api.verify_equation(eq, tampered(["AI"], lambda h: "zz")) is False)
p = copy.deepcopy(good)
del p["ipa"]
check("missing ipa section -> False", api.verify_equation(eq, p) is False)

# proof bound to its equation
check("equation mismatch -> False",
      api.verify_equation("x^3 + 2*x + 5 = 39", good) is False)

# ---------------- 4. zero-knowledge sanity ----------------
# proofs are randomized: two proofs of the same statement differ
p1 = api.prove_equation(eq, w)
p2 = api.prove_equation(eq, w)
check("proofs are randomized (AI differs)",
      p1["AI"] != p2["AI"] and p1["S"] != p2["S"])
check("both randomized proofs verify",
      api.verify_equation(eq, p1) and api.verify_equation(eq, p2))

# ---------------- 5. scaling spot check ----------------
eq_big = "x^9 + x^7 - 3*x^4 + 2*x = 596"   # x=2: 512+128-48+4 = 596
t0 = time.time()
pb = api.prove_equation(eq_big, 2)
tp = time.time() - t0
t0 = time.time()
ok = api.verify_equation(eq_big, pb)
tv = time.time() - t0
check("degree-9 equation verifies (prove %.2fs, verify %.2fs)" % (tp, tv), ok)

# ---------------- 6. canonical binding: moved terms verify ----------------
# (regression test for the "moved term shows as failed" bug: proofs are
# bound to the polynomial normal form, not the raw string)

pm = api.prove_equation("x^3 + 2*x + 5 = 38", 3)
check("moved term verifies", api.verify_equation("x^3 + 2*x = 33", pm))
check("reordered verifies",
      api.verify_equation("38 = x^3 + 2*x + 5", pm))
check("different equation still fails",
      api.verify_equation("x^3 + 2*x + 5 = 39", pm) is False)
check("scaled (different poly) still fails",
      api.verify_equation("2*x^3 + 4*x + 10 = 76", pm) is False)

# ---------------- 7. tolerance mode: decimals & transcendental ----------------

t0 = time.time()
# witnesses come from the solver (8 decimals), as the strict precision
# 1/(k*10^4) requires
w_sin = api.solve_equation("sin(x) = 0.5", decimals=8)[0]
pt = api.prove_equation("sin(x) = 0.5", w_sin, "1000")
tpt = time.time() - t0
t0 = time.time()
okt = api.verify_equation("sin(x) = 0.5", pt, "1000")
tvt = time.time() - t0
check("sin(x)=0.5 @1/1000 verifies (prove %.1fs, verify %.1fs)" % (tpt, tvt),
      okt and api.verify_equation("sin(x) = 0.5", pt))  # k from proof
check("sin moved term verifies",
      api.verify_equation("sin(x) - 0.5 = 0", pt, "1000"))
check("sin wrong k fails",
      api.verify_equation("sin(x) = 0.5", pt, "100") is False)
expect_raises("sin bad witness refused",
              lambda: api.prove_equation("sin(x) = 0.5", "1.0", "1000"))
expect_raises("sin missing precision refused",
              lambda: api.prove_equation("sin(x) = 0.5", "0.5236", ""))

pdec = api.prove_equation("x = 0.5", "0.5", "1000")
check("decimal x=0.5 verifies", api.verify_equation("x = 0.5", pdec, "1000"))
expect_raises("decimal bad witness refused",
              lambda: api.prove_equation("x = 0.5", "0.6", "1000"))

pexp = api.prove_equation("exp(x) = 2.7183",
                         api.solve_equation("exp(x) = 2.7183", decimals=8)[0],
                         "1000")
check("exp(x)=2.7183 verifies",
      api.verify_equation("exp(x) = 2.7183", pexp, "1000"))
pln = api.prove_equation("ln(x) = 0.6931",
                        api.solve_equation("ln(x) = 0.6931", decimals=8)[0],
                        "1000")
check("ln(x)=0.6931 @1/1000 verifies",
      api.verify_equation("ln(x) = 0.6931", pln, "1000"))
# Chebyshev on [1,2] is accurate at x=2 (Taylor-12 was ~0.04 off there)
pcheb = api.prove_equation("cos(x) = 0.5",
                           api.solve_equation("cos(x) = 0.5", decimals=8)[1],
                           "1000")
check("cos(x)=0.5 @1/1000 verifies",
      api.verify_equation("cos(x) = 0.5", pcheb, "1000"))

# strict precision: witness must be within 1/(k*10^4) of an exact root
expect_raises("4-decimal witness now refused (not within 1/(k*10^4))",
              lambda: api.prove_equation("sin(x) = 0.5", "0.5236", "1000"))
expect_raises("13-decimal witness refused",
              lambda: api.prove_equation("sin(x) = 0.5", "0.5235987755983",
                                        "1000"))
# 10-decimal witness (from solver) is accepted
w10 = api.solve_equation("sin(x) = 0.5", decimals=10)[0]
check("10-decimal solver root accepted",
      api.verify_equation("sin(x) = 0.5",
                          api.prove_equation("sin(x) = 0.5", w10, "1000"),
                          "1000"))
# solver: accuracy + display
sr = api.solve_equation("x^2 = 2", decimals=8)
check("solver x^2=2 -> +-1.41421356",
      sr == ["-1.41421356", "1.41421356"])
sr5 = api.solve_equation("(x-1)*(x-2)*(x-3)*(x-4)*(x-5)*(x-6) = 0",
                         decimals=8)
check("solver shows at most 5 roots", len(sr5) == 5 and sr5[0] == "1")
check("solver no real roots -> []",
      api.solve_equation("x^2 + 1 = 0", decimals=8) == [])
check("solver double root -> single '0'",
      api.solve_equation("x^2 = 0", decimals=8) == ["0"])
# exact equation: solver returns plain integer; tapping it proves fine
check("solver exact eq returns '3'",
      api.solve_equation("x^3 + 2*x + 5 = 38", decimals=10) == ["3"])
pex = api.prove_equation("x^3 + 2*x + 5 = 38", "3")
check("tapped integer witness proves",
      api.verify_equation("x^3 + 2*x + 5 = 38", pex))
pex2 = api.prove_equation("x^3 + 2*x + 5 = 38", "3.0000000000")
check("decimal-form integer witness proves",
      api.verify_equation("x^3 + 2*x + 5 = 38", pex2))
expect_raises("non-integer witness refused in exact mode",
              lambda: api.prove_equation("x^3 + 2*x + 5 = 38", "3.5"))

# tampering with the embedded precision must fail (k falls back to the
# proof's own value, so the canonical form mismatches)
pt2 = copy.deepcopy(pt)
pt2["precision"] = 100
check("tampered precision -> False",
      api.verify_equation("sin(x) = 0.5", pt2) is False)

# ---------------- 8. exact/tolerance mode agreement ----------------
# exact equation + decimal witness -> friendly error, not a raw int() crash
try:
    api.prove_equation("x^5=1", "1.001", "1000")
    check("decimal witness on exact equation refused", False)
except ValueError as e:
    check("decimal witness on exact equation refused",
          "integer" in str(e) and "invalid literal" not in str(e))
# exact equation with a (now ignored) precision still works
pe = api.prove_equation("x^5=1", 1, "1000")
check("exact equation ignores precision",
      api.verify_equation("x^5=1", pe))
# decimal equation forces tolerance mode and requires k
expect_raises("decimal equation without precision refused",
              lambda: api.prove_equation("x^5 = 1.0", "1.0", ""))
pt5 = api.prove_equation("x^5 = 1.0", "1.0", "1000")
check("x^5=1.0 @1/1000 verifies",
      api.verify_equation("x^5 = 1.0", pt5, "1000"))

# ---------------- 9. proof certificate (describe_proof) ----------------
# valid proof -> valid cert with QED
pc = api.prove_equation("x^3 + 2*x + 5 = 38", 3)
cert = api.describe_proof(pc)
check("certificate of valid proof: valid=True",
      cert["valid"] is True)
check("certificate of valid proof: qed=True",
      cert["qed"] is True)
check("certificate checks all pass",
      [c["ok"] for c in cert["checks"]] == [True, True])
check("certificate proposition shows user input",
      "x^3 + 2*x + 5 = 38" in cert["proposition"]["zh"])
# JSON-serializable
json.dumps(cert)
check("certificate is JSON-serializable", True)

# tampered t_hat -> invalid, no QED, short-circuit [False, None]
ptc = copy.deepcopy(pc)
ptc["t_hat"] = "%064x" % ((int(ptc["t_hat"], 16) + 1) % 2**256)
cert_t = api.describe_proof(ptc)
check("tampered t_hat -> valid=False",
      cert_t["valid"] is False)
check("tampered t_hat -> qed=False (no QED)",
      cert_t["qed"] is False)
check("tampered t_hat -> checks [False, None] (short-circuit)",
      [c["ok"] for c in cert_t["checks"]] == [False, None])

# input swapped to another equation -> binding fails -> valid=False
pbc = copy.deepcopy(pc)
pbc["input"] = "x + 5 = 8"
cert_b = api.describe_proof(pbc)
check("swapped input -> valid=False (binding)",
      cert_b["valid"] is False)

# legacy proof (no input field) -> valid=None, not False
plc = copy.deepcopy(pc)
del plc["input"]
cert_l = api.describe_proof(plc)
check("legacy proof (no input) -> valid=None",
      cert_l["valid"] is None)

# tolerance mode: proposition shows user input, canonical as footnote
pct = api.prove_equation("sin(x) = 0.5", "0.52359878", "1000")
cert_tol = api.describe_proof(pct)
check("tolerance cert valid=True",
      cert_tol["valid"] is True)
check("tolerance proposition shows user input",
      "sin(x) = 0.5" in cert_tol["proposition"]["zh"])
check("tolerance canonical kept as binding footnote",
      "13611596267" in cert_tol["binding_note"]["zh"])

# malformed proof still raises
expect_raises("describe_proof rejects malformed proof",
              lambda: api.describe_proof({"equation": "x=1"}))

print()
print("passed %d/%d" % (len(PASS), len(PASS) + len(FAIL)))
if FAIL:
    print("FAILURES:")
    for f in FAIL:
        print("  - " + f)
    raise SystemExit(1)
