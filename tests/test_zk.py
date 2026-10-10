"""Correctness + tamper-resistance tests for the ZK core.

Run:  python3 -m tests.test_zk
(plain asserts, no third-party test framework needed)
"""

import copy
import json
import time
from fractions import Fraction

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
# Precision k is mandatory for every equation (final precision 1/k).
K = "1000"

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
    proof = api.prove_equation(eq, w, K)
    tp = time.time() - t0
    t0 = time.time()
    ok = api.verify_equation(eq, proof, K)
    tv = time.time() - t0
    # JSON round-trip must preserve verifiability
    ok2 = api.verify_equation(eq, json.loads(json.dumps(proof)), K)
    check("honest %-28s x=%-6d verifies (prove %.2fs, verify %.2fs, json-ok)"
          % (eq, w, tp, tv), ok and ok2)

# ---------------- 2. wrong witnesses fail ----------------

for eq, w in [("x^3 + 2*x + 5 = 38", 4), ("x^2 = 16", 5), ("x + 5 = 8", 0)]:
    expect_raises("bad witness refused: %s x=%d" % (eq, w),
                  lambda: api.prove_equation(eq, w, K))

# unsatifiable constant equation is a compile error
expect_raises("constant-false equation rejected",
              lambda: api.prove_equation("3 = 4", 0, K))

# malformed equations are compile errors
for bad in ["x^2 = ", "= 5", "x + = 3", "2y + 1 = 5", "x^(-1) = 2", "x^2 == 4"]:
    expect_raises("malformed rejected: %r" % bad,
                  lambda: api.prove_equation(bad, 1, K))

# precision is mandatory for every equation (no exceptions)
expect_raises("prove without precision refused",
              lambda: api.prove_equation("x^3 + 2*x + 5 = 38", 3, ""))
expect_raises("prove with blank precision refused",
              lambda: api.prove_equation("x^2 = 16", 4, "   "))

# ---------------- 3. soundness at the protocol level ----------------
# Bypass the API's witness check: feed bad wires straight into the prover.
# A cheating prover must still fail verification.

eq, w = "x^3 + 2*x + 5 = 38", 3
circ = circuit_mod.compile(eq, K)
assert circ.mode == "exact"  # integer polynomial -> small exact circuit
aL, aR, aO = circ.evaluate(w)  # exact mode: integer witness directly

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
good = api.prove_equation(eq, w, K)


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
          api.verify_equation(eq, p, K) is False)

# garbage / truncated proofs must not raise, just return False
check("empty dict -> False", api.verify_equation(eq, {}, K) is False)
check("garbage hex -> False",
      api.verify_equation(eq, tampered(["AI"], lambda h: "zz"), K) is False)
p = copy.deepcopy(good)
del p["ipa"]
check("missing ipa section -> False", api.verify_equation(eq, p, K) is False)

# proof bound to its equation
check("equation mismatch -> False",
      api.verify_equation("x^3 + 2*x + 5 = 39", good, K) is False)

# ---------------- 4. zero-knowledge sanity ----------------
# proofs are randomized: two proofs of the same statement differ
p1 = api.prove_equation(eq, w, K)
p2 = api.prove_equation(eq, w, K)
check("proofs are randomized (AI differs)",
      p1["AI"] != p2["AI"] and p1["S"] != p2["S"])
check("both randomized proofs verify",
      api.verify_equation(eq, p1, K) and api.verify_equation(eq, p2, K))

# ---------------- 5. scaling spot check ----------------
eq_big = "x^9 + x^7 - 3*x^4 + 2*x = 596"   # x=2: 512+128-48+4 = 596
t0 = time.time()
pb = api.prove_equation(eq_big, 2, K)
tp = time.time() - t0
t0 = time.time()
ok = api.verify_equation(eq_big, pb, K)
tv = time.time() - t0
check("degree-9 equation verifies (prove %.2fs, verify %.2fs)" % (tp, tv), ok)

# ---------------- 6. canonical binding: moved terms verify ----------------
# (regression test for the "moved term shows as failed" bug: proofs are
# bound to the polynomial normal form, not the raw string)

pm = api.prove_equation("x^3 + 2*x + 5 = 38", 3, K)
check("moved term verifies", api.verify_equation("x^3 + 2*x = 33", pm, K))
check("reordered verifies",
      api.verify_equation("38 = x^3 + 2*x + 5", pm, K))
check("different equation still fails",
      api.verify_equation("x^3 + 2*x + 5 = 39", pm, K) is False)
check("scaled (different poly) still fails",
      api.verify_equation("2*x^3 + 4*x + 10 = 76", pm, K) is False)

# ---------------- 7. precision: mandatory, by definition ----------------
# Final precision is 1/k (k = user-supplied denominator). The witness
# must be within 1/k of an exact root. The verifier MUST supply k
# explicitly; it is NEVER taken from the proof.

t0 = time.time()
w_sin = api.solve_equation("sin(x) = 0.5", decimals=8)[0]
pt = api.prove_equation("sin(x) = 0.5", w_sin, "1000")
tpt = time.time() - t0
t0 = time.time()
okt = api.verify_equation("sin(x) = 0.5", pt, "1000")
tvt = time.time() - t0
check("sin(x)=0.5 @1/1000 verifies (prove %.1fs, verify %.1fs)" % (tpt, tvt),
      okt)
# verify WITHOUT k must raise (never fall back to the proof's k)
expect_raises("verify without precision raises",
              lambda: api.verify_equation("sin(x) = 0.5", pt, ""))
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

# strict precision: witness must be within 1/k of an exact root.
# sin(x)=0.5 root is ~0.5235987756; "0.53" is ~0.0014 away > 1/1000.
expect_raises("witness beyond 1/k refused",
              lambda: api.prove_equation("sin(x) = 0.5", "0.53", "1000"))
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
check("solver no real roots -> complex roots",
      api.solve_equation("x^2 + 1 = 0", decimals=8) == ["-i", "i"])
check("solver double root -> single '0'",
      api.solve_equation("x^2 = 0", decimals=8) == ["0"])
# exact equation: solver returns plain integer for the real root;
# complex roots are also shown (v0.4.3)
sr_ex = api.solve_equation("x^3 + 2*x + 5 = 38", decimals=10)
check("solver exact eq returns '3' for real root",
      "3" in sr_ex)
pex = api.prove_equation("x^3 + 2*x + 5 = 38", "3", K)
check("tapped integer witness proves",
      api.verify_equation("x^3 + 2*x + 5 = 38", pex, K))
pex2 = api.prove_equation("x^3 + 2*x + 5 = 38", "3.0000000000", K)
check("decimal-form integer witness proves",
      api.verify_equation("x^3 + 2*x + 5 = 38", pex2, K))
# non-integer witness for an integer equation is refused (exact circuit
# needs an integer witness; write the equation with decimals for that)
expect_raises("non-integer witness refused for integer equation",
              lambda: api.prove_equation("x^3 + 2*x + 5 = 38", "3.5", K))

# the proof's "precision" field is informational only: tampering with it
# must NOT affect verification, which uses the caller-supplied k
# (the definition), never the proof's value
pt2 = copy.deepcopy(pt)
pt2["precision"] = 100
check("tampered proof precision field ignored (explicit k verifies)",
      api.verify_equation("sin(x) = 0.5", pt2, "1000"))
check("tampered proof precision + wrong k still fails",
      api.verify_equation("sin(x) = 0.5", pt2, "100") is False)

# ---------------- 8. precision mandatory in both modes ----------------
# k is required for every equation; the circuit mode (exact/tolerance)
# is a pure size optimization chosen from the equation alone.
expect_raises("integer equation without precision refused",
              lambda: api.prove_equation("x^5=1", 1, ""))
expect_raises("verify without precision raises",
              lambda: api.verify_equation(
                  "x^5=1", api.prove_equation("x^5=1", 1, "1000"), ""))
pt5 = api.prove_equation("x^5 = 1.0", "1.0", "1000")
check("x^5=1.0 @1/1000 verifies",
      api.verify_equation("x^5 = 1.0", pt5, "1000"))
# integer equation with explicit precision works (exact circuit + 1/k check)
pe = api.prove_equation("x^5=1", 1, "1000")
check("integer equation with precision verifies",
      api.verify_equation("x^5=1", pe, "1000"))
# decimal witness for integer equation: helpful error, not a crash
expect_raises("decimal witness for integer equation refused",
              lambda: api.prove_equation("x^2 = 2", "1.41421356", "1000"))

# ---------------- 9. complex numbers (v0.4.3) ----------------
# Gaussian integer equations with complex witnesses.
pc1 = api.prove_equation("x^2 + 1 = 0", "i", "1000")
check("x^2+1=0 witness i verifies",
      api.verify_equation("x^2 + 1 = 0", pc1, "1000"))
check("x^2+1=0 witness i wrong k fails",
      api.verify_equation("x^2 + 1 = 0", pc1, "999") is False)
pc2 = api.prove_equation("x^2 + 1 = 0", "-i", "1000")
check("x^2+1=0 witness -i verifies",
      api.verify_equation("x^2 + 1 = 0", pc2, "1000"))
# Complex coefficients
pc3 = api.prove_equation("x = 1+i", "1+i", "1000")
check("x=1+i verifies",
      api.verify_equation("x = 1+i", pc3, "1000"))
pc4 = api.prove_equation("(1+2i)*x = 1+2i", "1", "1000")
check("(1+2i)*x=1+2i witness 1 verifies",
      api.verify_equation("(1+2i)*x = 1+2i", pc4, "1000"))
# Complex solver
cr = api.solve_equation("x^2 + 1 = 0")
check("solver finds i and -i for x^2+1=0",
      "i" in cr and "-i" in cr)
# Wrong witness refused
expect_raises("wrong complex witness refused",
              lambda: api.prove_equation("x^2 + 1 = 0", "1+i", "1000"))
# Non-Gaussian-integer complex witness for exact equation refused
expect_raises("non-Gaussian witness refused",
              lambda: api.prove_equation("x^2 + 1 = 0", "0.5i", "1000"))

# ---- integral equations: x = int(f, a, b) ----
# integral_value helper
iv = api.integral_value("x = int(x^2, 0, 1)")
check("integral_value(x^2,0,1) ~ 1/3", abs(float(iv) - 1/3) < 1e-9)
# Prove and verify
pi = api.prove_equation("x = int(x^2, 0, 1)", iv, "1000")
check("integral x=int(x^2,0,1) proves",
      api.verify_equation("x = int(x^2, 0, 1)", pi, "1000"))
# Red line: no numeric value in the equation field
check("integral equation hides numeric value",
      "0.333" not in pi["equation"] and pi["equation"].startswith("x=int("))
# Wrong k does not verify
check("integral wrong k rejected",
      not api.verify_equation("x = int(x^2, 0, 1)", pi, "999"))
# Reversed form
pi2 = api.prove_equation("int(x^2, 0, 1) = x", iv, "1000")
check("integral int=x form proves",
      api.verify_equation("int(x^2, 0, 1) = x", pi2, "1000"))
# Wrong witness refused
expect_raises("integral wrong witness refused",
              lambda: api.prove_equation("x = int(x^2, 0, 1)", "0.5", "1000"))
# Form restriction: only x = int(...) allowed
expect_raises("integral x+int rejected",
              lambda: api.prove_equation("x + int(x^2, 0, 1) = 1", "0.6", "100"))
expect_raises("integral 2x=int rejected",
              lambda: api.prove_equation("2*x = int(x^2, 0, 1)", "0.6", "100"))

# FTC secrecy: circuit contains no trace of the integral value.
# (The value 1/3 = inv(3) mod N must not appear in constraints.)
from core import circuit as _cm
_fc = _cm.compile("x = int(x^2, 0, 1)", "1000")
_inv3 = pow(3, _cm.N - 2, _cm.N)
_flat = []
for qi in range(_fc.q):
    _flat.extend(_fc.WL[qi])
    _flat.extend(_fc.WR[qi])
    _flat.extend(_fc.WO[qi])
    _flat.append(_fc.c[qi])
check("FTC circuit leaks no value", _inv3 not in _flat and _inv3 not in _fc.c)
check("FTC mode is ftc", _fc.mode == "ftc")

# ---- ODE: general linear ODEs, any order ----
# First-order, general form (y' - y = 0, not forced y' = ...)
oc = api.ode_coefficients("ode(y' - y = 0, y(0) = 1, deg = 5)")
check("ode general form coeffs", oc == [Fraction(1), Fraction(1), Fraction(1,2),
      Fraction(1,6), Fraction(1,24), Fraction(1,120)])
po = api.prove_equation("ode(y' - y = 0, y(0) = 1, deg = 5)", "", "")
check("ode general form proves",
      api.verify_equation("ode(y' - y = 0, y(0) = 1, deg = 5)", po, ""))
# Canonical is symbolic, no k
check("ode canonical symbolic, no k",
      po["equation"].startswith("ode(") and ";k=" not in po["equation"])
# k is truly meaningless: any k verifies
check("ode k meaningless",
      api.verify_equation("ode(y' - y = 0, y(0) = 1, deg = 5)", po, "999"))
# Coarse deg verifies: deg=5 proof against deg=3 statement
check("ode coarse deg verifies",
      api.verify_equation("ode(y' - y = 0, y(0) = 1, deg = 3)", po, ""))
# Finer deg does NOT verify
check("ode finer deg rejected",
      not api.verify_equation("ode(y' - y = 0, y(0) = 1, deg = 6)", po, ""))
# Different ODE does NOT verify
check("ode different ODE rejected",
      not api.verify_equation("ode(y' - 2*y = 0, y(0) = 1, deg = 3)", po, ""))
# Second-order: y'' + y = 0, y(0)=0, y'(0)=1 -> sin(x)
oc2 = api.ode_coefficients("ode(y'' + y = 0, y(0) = 0, y'(0) = 1, deg = 6)")
check("ode 2nd order coeffs (sin)", oc2 == [Fraction(0), Fraction(1),
      Fraction(0), Fraction(-1,6), Fraction(0), Fraction(1,120), Fraction(0)])
po2 = api.prove_equation("ode(y'' + y = 0, y(0) = 0, y'(0) = 1, deg = 6)", "", "")
check("ode 2nd order proves",
      api.verify_equation("ode(y'' + y = 0, y(0) = 0, y'(0) = 1, deg = 6)", po2, ""))
check("ode 2nd order coarse",
      api.verify_equation("ode(y'' + y = 0, y(0) = 0, y'(0) = 1, deg = 4)", po2, ""))

print()
print("passed %d/%d" % (len(PASS), len(PASS) + len(FAIL)))
if FAIL:
    print("FAILURES:")
    for f in FAIL:
        print("  - " + f)
    raise SystemExit(1)
