# Constant-time curve operations (debug^2, experimental)

Branch: `exp/ct-curve`. This is a **debug/experimental** implementation, not
audited, not for production.

## What it is

A C extension (`core/cext/curve_ext`) that reimplements `core/curve.py`'s
group operations with a **fixed operation sequence** — no secret-dependent
branches in the hot path (MSM / scalar multiplication).

This answers the reviewer's "pure Python can't be constant-time" point:
Python's big ints have data-dependent branches (Karatsuba thresholds,
early exits); the C code does the same limb operations regardless of
secret scalar values.

## What's constant-time

- **MSM**: fixed 256 iterations (not `maxbits`), masked conditional adds.
  Scalar bits are read but never branched on.
- **Field ops** (`fe_add`, `fe_sub`, `fe_mul`): straight-line limb
  arithmetic + masked selects. No branches on values.
- **Group ops** (`gej_add`, `gej_dbl`): branch-free complete formulas.
  Infinity and doubling cases selected via masks, not branches.
- **Inverse / sqrt**: square-and-multiply-always (Fermat).

## What's NOT claimed

- **Not audited.** This is a from-scratch implementation. It is
  byte-identical to the Python implementation on the differential test
  suite, but that is correctness, not a security audit.
- **Compiler-introduced branches.** The C source has no secret-dependent
  branches, but the compiler (`-O2`) may introduce them (e.g., for
  `__uint128_t` division, or cmov lowering). Verify the assembly if this
  matters.
- **Python-level overhead.** Point serialization, memory allocation, and
  the Python↔C boundary are not constant-time. The guarantee covers the
  cryptographic operation sequence, not the entire process.
- **Side channels beyond timing.** Cache, power, EM are out of scope.

## Performance

~2x faster than pure Python on 128-term MSM (0.41s vs 0.56s). The win is
constant-time, not speed — Python's big ints are already C-optimized.

## API

`core/curve_c.py` is a drop-in replacement for `core/curve.py`:
- Tries `core.cext.curve_ext`, falls back to pure Python if not built.
- `USING_C` flag indicates which backend is active.
- NUMS generators (`H`, `U`, `Gvec`, `Hvec`) stay in Python (public,
  one-time setup).

Build: `cd core/cext && python3 setup.py build_ext --inplace`

Test: `python3 core/cext/test_differential.py` (must show `fails: 0`)

## Files

- `core/cext/curve_ext.c` — the C implementation (~700 lines)
- `core/cext/setup.py` — build script
- `core/cext/test_differential.py` — C vs Python differential test
- `core/curve_c.py` — shim with fallback

## Bugs found during development

1. `gej_dbl` reused `YYYY` as a temp, corrupting `8*YYYY`.
2. `fe_sqrt` exponent `(p+1)/4` limbs were wrong.
3. `decompress` had reversed byte order in x parsing.
4. `fe_mul` accumulator overflow: 8 terms × 2^64 needs 2^67, but
   `uint64_t` holds 2^64. Fixed with `__uint128_t`.
5. `fe_mul` dropped residual `t[8]` after folds.
6. `gej_dbl`/`gej_add` not in-place safe (`gej_dbl(&acc, &acc)` corrupted
   `a->Y` after writing `r->Y`). Fixed by reading inputs before writes.

All caught by differential testing against the Python implementation.
