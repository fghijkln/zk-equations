# Native backends C ABI (v1)

Shared C ABI implemented by the Rust (`curve_rs`) and Zig (`curve_zig`)
cdylib backends. Semantics must match `core/cext/curve_ext.c` exactly.

## Types

- `fe`: 8 x `uint32_t` limbs, little-endian, fully reduced mod
  p = 2^256 - 2^32 - 977 after every operation.
- `scalar`: 8 x `uint32_t` limbs, little-endian. **Caller** reduces mod
  N (curve order) before calling; native code assumes 0 <= s < 2^256.
- `point`: 24 x `uint32_t` = X[8], Y[8], Z[8] (Jacobian).
  All-zero (Z == 0) = point at infinity.

## Functions

```c
void wcurve_add(const uint32_t *a, const uint32_t *b, uint32_t *r);
// r = a + b. Branch-free complete formula: infinity / doubling cases
// selected via masks (same as gej_add in curve_ext.c: add-2009-bl main
// formula + dbl-2009-l + explicit infinity handling).

void wcurve_dbl(const uint32_t *a, uint32_t *r);
// r = 2*a. dbl-2009-l. Z=0 -> infinity, no branch. In-place safe.

void wcurve_neg(const uint32_t *a, uint32_t *r);
// r = -a: negate Y in Jacobian (X, -Y mod p, Z). In-place safe.

int wcurve_is_inf(const uint32_t *a);
// 1 if Z == 0 (branch-free compare), else 0.

void wcurve_msm(const uint32_t *scalars, const uint32_t *points,
                uint32_t *r, size_t n);
// r = sum_{i=0}^{n-1} scalars[i] * points[i].
// Fixed 256 iterations (bit 255..0), masked conditional adds:
//   bit = (scalars[i][b>>5] >> (b&31)) & 1   // public loop index, never branched on
//   sel = points[i] & (0 - bit)              // masked point (all-zero = inf if bit=0)
//   acc = acc + sel                          // via wcurve_add (complete formula)
// acc starts at infinity. n == 0 -> r = infinity. In-place safe on r.

void wcurve_from_affine(const uint32_t *x, const uint32_t *y, uint32_t *r);
// r = (x, y, 1). Assumes x, y < p (caller guarantees).

void wcurve_to_affine(const uint32_t *a, uint32_t *x, uint32_t *y, int *is_inf);
// If Z == 0: *is_inf = 1 (x, y untouched). Else *is_inf = 0 and
// x = X/Z^2, y = Y/Z^3 mod p (via Fermat inverse, square-and-multiply-always).

void wcurve_compress(const uint32_t *a, uint8_t *out /*33*/);
// out[0] = 0x02 | (y & 1); out[1..33] = x big-endian.
// Caller guarantees a is not infinity. (Affine conversion inside.)

int wcurve_decompress(const uint8_t *in /*33*/, uint32_t *r /*24*/);
// Returns 0 on success, nonzero on error. Checks:
//   - in[0] == 0x02 or 0x03 (else error)
//   - x = big-endian in[1..33] satisfies x < p (else error)
//   - y2 = x^3 + 7; y = y2^((p+1)/4); y^2 == y2 (else "not on curve" error)
//   - if (y[0] & 1) != (in[0] & 1): y = -y
//   - r = (x, y, 1)
// (p+1)/4 limbs (little-endian): BFFFFF0C FFFFFFFF FFFFFFFF FFFFFFFF
//                                FFFFFFFF FFFFFFFF FFFFFFFF 3FFFFFFF
```

## Constant-time rules

- No branches on secret data (scalar bits, point coordinates). Loop
  bounds are public constants. Bit extraction uses public loop index.
- All selections via masks / cmov. No table lookups.
- Same caveats as the C version: compiler may introduce branches;
  verify assembly if it matters. Not audited.
