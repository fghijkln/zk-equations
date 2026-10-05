"""Differential test: C extension vs pure-Python curve implementation."""
import random
import sys

sys.path.insert(0, ".")

from core import curve as py
from core.cext import curve_ext as c

random.seed(1234)
N = py.N

fails = 0

def check(name, a, b):
    global fails
    if a != b:
        fails += 1
        print("FAIL", name, str(a)[:60], str(b)[:60])

# field-level via group ops; random points
pts = [py.G]
for _ in range(20):
    k = random.randrange(1, N)
    pts.append(py.mul(k, py.G))

# add / dbl / neg / mul / is_inf
for _ in range(300):
    p1, p2 = random.choice(pts), random.choice(pts)
    check("add", c.add(p1, p2), py.add(p1, p2))
    check("dbl", c.dbl(p1), py.dbl(p1))
    check("neg", c.neg(p1), py.neg(p1))
    check("is_inf", c.is_inf(p1), py.is_inf(p1))
check("add inf", c.add(None, pts[0]), py.add(None, pts[0]))
check("add inf2", c.add(pts[0], None), py.add(pts[0], None))
check("is_inf none", c.is_inf(None), True)
check("neg none", c.neg(None), None)
for _ in range(100):
    k = random.randrange(0, N)
    p = random.choice(pts)
    check("mul", c.mul(k, p), py.mul(k, p))
check("mul 0", c.mul(0, pts[0]), py.mul(0, pts[0]))
check("mul N", c.mul(N, pts[0]), py.mul(N, pts[0]))
check("mul neg scalar", c.mul(-5, pts[1]), py.mul(-5, pts[1]))

# msm
for trial in range(30):
    n = random.randint(1, 12)
    ss = [random.randrange(-N, 2 * N) for _ in range(n)]
    ps = [random.choice(pts) for _ in range(n)]
    check("msm n=%d" % n, c.msm(ss, ps), py.msm(ss, ps))
# msm with zero scalar / inf point / empty
check("msm zero scalar", c.msm([0, 5], [pts[0], pts[1]]), py.msm([0, 5], [pts[0], pts[1]]))
check("msm inf point", c.msm([3, 5], [None, pts[1]]), py.msm([3, 5], [None, pts[1]]))
check("msm empty", c.msm([], []), py.msm([], []))

# compress / decompress roundtrip
for _ in range(100):
    p = py.mul(random.randrange(1, N), py.G)
    check("compress", c.compress(p), py.compress(p))
    check("decompress", c.decompress(c.compress(p)), p)
# decompress known vector
check("decompress G", c.decompress(py.compress(py.G)), py.G)

print("fails:", fails)
sys.exit(1 if fails else 0)
