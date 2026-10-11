"""Benchmark: C extension vs pure Python for MSM and batch operations.

Tests:
1. msm() - multi-scalar multiplication (the hot path)
2. batch_mul() - n scalar-point multiplications
3. Full prove()/verify() with current backend

Run: python3 bench.py
"""
import time
import sys
import random

sys.path.insert(0, '.')


def bench_msm(n):
    """Benchmark MSM of size n: Python vs C."""
    from core import curve as curve_py
    
    # Generate random scalars and points
    random.seed(42)
    scalars = [random.getrandbits(256) for _ in range(n)]
    points = [curve_py.Gvec(i % 100) for i in range(n)]  # reuse cached gens
    
    # Python MSM
    t0 = time.perf_counter()
    r_py = curve_py.msm(scalars, points)
    t_py = time.perf_counter() - t0
    
    # C MSM (if available)
    t_c = None
    try:
        from core.cext import curve_ext
        # Convert points to C format (tuples work)
        t0 = time.perf_counter()
        r_c = curve_ext.msm(scalars, points)
        t_c = time.perf_counter() - t0
        # Verify correctness
        assert r_py == r_c or (r_py is None and r_c is None), "MSM mismatch!"
    except ImportError:
        pass
    
    return t_py, t_c


def bench_batch_mul(n):
    """Benchmark batch_mul of size n."""
    from core import curve as curve_py
    
    random.seed(42)
    scalars = [random.getrandbits(256) for _ in range(n)]
    points = [curve_py.Gvec(i % 100) for i in range(n)]
    
    # Python (loop)
    t0 = time.perf_counter()
    r_py = curve_py.batch_mul(scalars, points)
    t_py = time.perf_counter() - t0
    
    # C (if available)
    t_c = None
    try:
        from core.cext import curve_ext
        if hasattr(curve_ext, 'batch_mul'):
            t0 = time.perf_counter()
            r_c = curve_ext.batch_mul(scalars, points)
            t_c = time.perf_counter() - t0
            assert r_py == r_c, "batch_mul mismatch!"
    except ImportError:
        pass
    
    return t_py, t_c


def bench_prove_verify():
    """Benchmark full prove/verify (uses current backend)."""
    from core import api
    
    t0 = time.perf_counter()
    proof = api.prove_equation('x^2 = 4', '2', '1000')
    t_prove = time.perf_counter() - t0
    
    t0 = time.perf_counter()
    ok = api.verify_equation('x^2 = 4', proof, '1000')
    t_verify = time.perf_counter() - t0
    
    assert ok
    return t_prove, t_verify


def main():
    print("=" * 65)
    print("Benchmark: Python vs C extension")
    print("=" * 65)
    
    try:
        from core.cext import curve_ext
        has_c = True
        has_batch = hasattr(curve_ext, 'batch_mul')
        print(f"C extension: YES (batch_mul={has_batch}, batch_fold={hasattr(curve_ext, 'batch_fold')})")
    except ImportError:
        has_c = False
        print("C extension: NO")
    print()
    
    # MSM benchmarks
    print("MSM (multi-scalar multiplication):")
    print("-" * 65)
    print(f"{'n':>6} | {'Python':>12} | {'C':>12} | {'Speedup':>8}")
    print("-" * 65)
    for n in [16, 64, 256]:
        t_py, t_c = bench_msm(n)
        if t_c:
            speedup = t_py / t_c
            print(f"{n:>6} | {t_py*1000:>10.1f}ms | {t_c*1000:>10.1f}ms | {speedup:>7.1f}x")
        else:
            print(f"{n:>6} | {t_py*1000:>10.1f}ms | {'N/A':>12} | {'N/A':>8}")
    print()
    
    # batch_mul benchmarks
    if has_c and has_batch:
        print("batch_mul (n scalar-point multiplications):")
        print("-" * 65)
        print(f"{'n':>6} | {'Python':>12} | {'C':>12} | {'Speedup':>8}")
        print("-" * 65)
        for n in [16, 64, 256]:
            t_py, t_c = bench_batch_mul(n)
            if t_c:
                speedup = t_py / t_c
                print(f"{n:>6} | {t_py*1000:>10.1f}ms | {t_c*1000:>10.1f}ms | {speedup:>7.1f}x")
        print()
    
    # Full prove/verify
    print("Full prove/verify (x^2 = 4):")
    print("-" * 65)
    t_prove, t_verify = bench_prove_verify()
    print(f"  prove:  {t_prove*1000:.1f} ms")
    print(f"  verify: {t_verify*1000:.1f} ms")


if __name__ == "__main__":
    main()
