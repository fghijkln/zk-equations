"""Simpson tolerance measurement for the Witness integral-proof design.

Question: for `x = int(f, a, b)` compiled as the tolerance-mode equation
    T(x) = x - S_n  (|T(x)| < 1/k in-circuit),
where S_n is the compile-time Simpson sum (exact Fractions) of the
integrand polynomial f on [a,b] with n subintervals, how large must n be
so an honest prover (witness within 1/(k*10^4) of the TRUE integral
I_true) always passes?

Acceptance rule used here: e_n = |S_n - I_true| < 1/(2k)
(half the 1/k budget; the other half covers witness distance 1/(k*10^4)
plus quantization, both negligible in comparison).

I_true is the analytic integral of the polynomial the proof is actually
about (for sin/cos/exp/ln: the fixed Chebyshev interpolant -- the same
honesty note as the main language). All arithmetic in Fractions: the
measurement replicates compile-time semantics exactly.

Battery: polynomials of several degrees + each transcendental +
one composition stress case (exp(x^2), degree 22 after composition).
"""
import sys
import time
from fractions import Fraction

sys.path.insert(0, "/home/hatch/workspace/zk-equations")
from core import circuit as C

BATTERY = [
    ("x^2", "0", "1"),                    # user's example: I = 1/3
    ("x^3+2*x+5", "0", "0.5"),            # I = 2.765625
    ("x^10", "0", "1"),                   # high-degree stress: I = 1/11
    ("sin(x)", "0", "3.14159265"),        # I ~= 2 (of the Chebyshev poly)
    ("exp(x)", "0", "1"),                # I ~= e-1
    ("cos(x)", "0", "1.57079632"),       # I ~= 1
    ("ln(x)", "1", "2"),                 # I = 2ln2-1 ~= 0.386
    ("exp(x^2)", "0", "1"),              # composition stress, deg 22
]

NS = [2, 4, 8, 16, 32, 64, 128, 256]
KS = [10, 100, 1000, 10000]


def parse_integrand(expr):
    p = C.Parser(C.tokenize(expr))
    poly = p.parse_expr()
    if p.peek() is not None:
        raise ValueError("trailing tokens in integrand")
    return poly


def horner(poly, x):
    # exact Fraction Horner, highest-degree-first
    n = max(poly.keys()) if poly else 0
    v = Fraction(0)
    for e in range(n, -1, -1):
        v = v * x + poly.get(e, Fraction(0))
    return v


def true_integral(poly, a, b):
    # analytic: sum c_e (b^{e+1} - a^{e+1})/(e+1), exact Fractions
    t = Fraction(0)
    for e, c in poly.items():
        t += c * (b ** (e + 1) - a ** (e + 1)) / (e + 1)
    return t


def simpson(poly, a, b, n):
    # exact Fraction Simpson with n subintervals (n even)
    h = (b - a) / n
    s = horner(poly, a) + horner(poly, b)
    for i in range(1, n):
        w = 4 if i % 2 else 2
        s += w * horner(poly, a + i * h)
    return s * h / 3


def main():
    print("Simpson truncation error |S_n - I_true| (exact Fractions)")
    print("=" * 78)
    results = {}
    for expr, sa, sb in BATTERY:
        poly = parse_integrand(expr)
        deg = max(poly.keys()) if poly else 0
        a, b = Fraction(sa), Fraction(sb)
        I = true_integral(poly, a, b)
        assert abs(I) <= 5, "integral out of witness domain |x|<=5"
        t0 = time.time()
        errs = {}
        for n in NS:
            S = simpson(poly, a, b, n)
            errs[n] = abs(float(S - I))
        dt = time.time() - t0
        results[expr] = errs
        line = "f=%-14s deg=%2d I~%10.6f | " % (expr, deg, float(I))
        line += " ".join("n=%-3d %.1e" % (n, errs[n]) for n in NS)
        print(line)
        print("    (timing for all n: %.2fs)" % dt)
    print()
    print("Minimal n with |S_n - I_true| < 1/(2k)  (half-budget rule)")
    print("=" * 78)
    hdr = "%-14s | " % "integrand" + " ".join("k=%-5d" % k for k in KS)
    print(hdr)
    for expr, _, _ in BATTERY:
        row = "%-14s | " % expr
        for k in KS:
            bound = 1.0 / (2 * k)
            pick = next((n for n in NS if results[expr][n] < bound), None)
            row += "%-7s " % ("n=%d" % pick if pick else ">256")
        print(row)
    print()
    print("Worst-case compile-time cost check: largest recommended n, "
          "exact Fractions.")
    # time the single most expensive recommended evaluation
    expr, sa, sb = "exp(x^2)", "0", "1"
    poly = parse_integrand(expr)
    a, b = Fraction(sa), Fraction(sb)
    t0 = time.time()
    S = simpson(poly, a, b, 256)
    print("exp(x^2) n=256: %.2fs, |err|=%.2e" %
          (time.time() - t0, abs(float(S - true_integral(poly, a, b)))))


if __name__ == "__main__":
    main()
