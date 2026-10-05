"""Differential test: rust / zig / c / python backends must agree bit-for-bit.

Run:  python3 -m core.test_backends
"""
import random
import sys

from . import curve as py_ref
from . import curve_backend as cb


def _rand_point(rng):
    # random multiple of G (always on curve)
    k = rng.randrange(1, py_ref.N)
    return py_ref.mul(k, py_ref.G)


def check(name, got, want):
    if got != want:
        print("FAIL %s:\n  got  %r\n  want %r" % (name, got, want))
        return False
    return True


def run_suite(mod, tag, rng):
    ok = True
    pts = [_rand_point(rng) for _ in range(8)] + [None, py_ref.G]

    for p in pts:
        ok &= check("%s is_inf" % tag, mod.is_inf(p), py_ref.is_inf(p))
        ok &= check("%s dbl" % tag, mod.dbl(p), py_ref.dbl(p))
        ok &= check("%s neg" % tag, mod.neg(p), py_ref.neg(p))
    for a in pts:
        for b in pts:
            ok &= check("%s add" % tag, mod.add(a, b), py_ref.add(a, b))

    for k in [0, 1, 2, 123456789, py_ref.N - 1, py_ref.N, py_ref.N + 5,
              (1 << 256) - 1]:
        for p in pts:
            ok &= check("%s mul k=%d" % (tag, k), mod.mul(k, p),
                        py_ref.mul(k, p))

    for n in [1, 2, 5, 16, 64]:
        ss = [rng.randrange(0, py_ref.N * 2) for _ in range(n)]
        ps = [_rand_point(rng) if rng.random() > 0.15 else None
              for _ in range(n)]
        if rng.random() < 0.3:
            ss[rng.randrange(n)] = 0
        ok &= check("%s msm n=%d" % (tag, n), mod.msm(ss, ps),
                    py_ref.msm(ss, ps))

    # compress / decompress
    for p in pts:
        if p is None:
            continue
        c = mod.compress(p)
        ok &= check("%s compress" % tag, c, py_ref.compress(p))
        ok &= check("%s decompress" % tag, mod.decompress(c), p)
    # bad decompress inputs must raise on all backends
    for bad in [b"\x00" * 33, b"\x04" + b"\x00" * 32,
                b"\x02" + b"\xff" * 32, b"\x02" + b"\x00" * 31]:
        try:
            mod.decompress(bad)
            print("FAIL %s decompress should reject %r" % (tag, bad[:4]))
            ok = False
        except ValueError:
            pass
        except Exception as e:  # noqa: BLE001
            print("FAIL %s decompress wrong exc %r" % (tag, e))
            ok = False
    return ok


def main():
    rng = random.Random(0xC0FFEE)
    all_ok = True
    print("available backends:", cb.available_backends())
    for name in cb.available_backends():
        mod = cb._BACKENDS[name]
        # python backend: self-check skipped (it IS the reference)
        if name == "python":
            print("backend python: reference (self-check skipped)")
            continue
        ok = run_suite(mod, name, rng)
        print("backend %-6s: %s" % (name, "OK" if ok else "MISMATCH"))
        all_ok &= ok
    print("ALL BACKENDS AGREE" if all_ok else "DIFFERENCES FOUND")
    return 0 if all_ok else 1


if __name__ == "__main__":
    sys.exit(main())
