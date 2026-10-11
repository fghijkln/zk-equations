"""Quick benchmark: Rust vs Python vs C for MSM."""
import time
import sys
import ctypes
import random

sys.path.insert(0, '.')

def bench_rust_msm(n):
    """Benchmark Rust MSM via ctypes."""
    lib = ctypes.CDLL('/home/hatch/workspace/zk-equations/core/rsext/target/release/libcurve_rs.so')
    
    # Setup argtypes (from curve_backend.py)
    # wcurve_msm(scalars: *const u32, points: *const u32, r: *mut u32, n: i64)
    lib.wcurve_msm.argtypes = [
        ctypes.POINTER(ctypes.c_uint32),
        ctypes.POINTER(ctypes.c_uint32),
        ctypes.POINTER(ctypes.c_uint32),
        ctypes.c_int64
    ]
    lib.wcurve_msm.restype = None
    
    # Helper to convert int to 8 limbs
    def to_limbs(v):
        return [(v >> (32*i)) & 0xFFFFFFFF for i in range(8)]
    
    # Generate test data
    random.seed(42)
    from core import curve as curve_py
    scalars = [random.getrandbits(256) for _ in range(n)]
    points = [curve_py.Gvec(i % 100) for i in range(n)]
    
    # Convert to limb arrays
    # scalars: n*8 u32
    # points: n*24 u32 (Jacobian X,Y,Z)
    import array
    scalars_arr = array.array('I')
    points_arr = array.array('I')
    for s in scalars:
        scalars_arr.extend(to_limbs(s % curve_py.N))
    for (x, y) in points:
        points_arr.extend(to_limbs(x))
        points_arr.extend(to_limbs(y))
        points_arr.extend(to_limbs(1))  # Z=1
    
    # Result buffer: 24 u32
    result_arr = array.array('I', [0]*24)
    
    # Get pointers
    s_ptr = (ctypes.c_uint32 * len(scalars_arr)).from_buffer(scalars_arr)
    p_ptr = (ctypes.c_uint32 * len(points_arr)).from_buffer(points_arr)
    r_ptr = (ctypes.c_uint32 * 24).from_buffer(result_arr)
    
    t0 = time.perf_counter()
    lib.wcurve_msm(s_ptr, p_ptr, r_ptr, n)
    t_rust = time.perf_counter() - t0
    
    return t_rust

def main():
    from core import curve as curve_py
    
    print("MSM Benchmark: Python vs C vs Rust")
    print("=" * 50)
    
    for n in [16, 64]:
        print(f"\nn={n}:")
        
        # Python
        random.seed(42)
        scalars = [random.getrandbits(256) for _ in range(n)]
        points = [curve_py.Gvec(i % 100) for i in range(n)]
        t0 = time.perf_counter()
        curve_py.msm(scalars, points)
        t_py = time.perf_counter() - t0
        print(f"  Python: {t_py*1000:.1f} ms")
        
        # C
        try:
            from core.cext import curve_ext
            t0 = time.perf_counter()
            curve_ext.msm(scalars, points)
            t_c = time.perf_counter() - t0
            print(f"  C:      {t_c*1000:.1f} ms ({t_py/t_c:.1f}x)")
        except ImportError:
            print(f"  C:      N/A")
        
        # Rust
        try:
            t_rust = bench_rust_msm(n)
            print(f"  Rust:   {t_rust*1000:.1f} ms ({t_py/t_rust:.1f}x)")
        except Exception as e:
            print(f"  Rust:   ERROR {e}")

if __name__ == "__main__":
    main()
