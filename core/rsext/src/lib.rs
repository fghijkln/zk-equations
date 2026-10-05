/* Constant-time secp256k1 group operations (debug/experimental) — Rust backend.
 *
 * Line-by-line translation of core/cext/curve_ext.c into Rust, exposing the
 * C ABI specified in core/nbext/ABI.md as `wcurve_*`.
 *
 * Purpose: answer the "pure Python can't be constant-time" limitation with a
 * native implementation where the hot path (MSM / scalar mult) runs a fixed
 * operation sequence with no secret-dependent branches.
 *
 * Design (mirrors the C version):
 *  - Field: 8 x 32-bit limbs, fully reduced mod p after every op.
 *    All field ops are straight-line limb arithmetic + masked selects.
 *    Rust's u128 is used for the 128-bit multiply accumulator
 *    (8 products per limb need ~67 bits; a u64 accumulator would overflow --
 *    this was a real bug in the C version's first draft).
 *  - Group: Jacobian coordinates. Addition uses a branch-free complete
 *    formula (all cases -- infinity, doubling -- selected via masks).
 *    Doubling is dbl-2009-l (handles Z=0 correctly, no branch).
 *  - MSM: fixed 256 iterations, masked conditional adds. The scalar bits
 *    are read (public loop index) but never branched on.
 *  - Scalars are 256-bit (8 x 32-bit limbs).
 *
 * Rust notes: the C version writes through `fe *r` out-params and relies on
 * read-before-write for in-place calls like `fe_mul(y2, y2, x)`. Rust's
 * borrow checker rejects aliasing a `&mut` with a `&`, so the internal
 * functions return `Fe`/`Gej` by value instead (`y2 = fe_mul(&y2, &x)`).
 * Each C statement still maps to exactly one Rust statement; the operation
 * sequence and constant-time discipline are unchanged. In-place safety is
 * then trivially guaranteed.
 *
 * Constant-time discipline: no branches on secret data, all selections via
 * masks / cmov, fixed loop bounds, no table lookups.
 *
 * What is NOT claimed: this has not been independently audited. The
 * constant-time property holds for the operation sequence; it does not
 * cover FFI-level overhead, or compiler-introduced branches
 * (verify the assembly if this matters to you).
 *
 * Zero external dependencies. `unsafe` is used only at the FFI boundary
 * (raw pointer <-> owned value copies); the arithmetic itself is safe Rust
 * using wrapping operations so it can never panic.
 */

/* ================= field arithmetic mod p = 2^256 - 2^32 - 977 ================= */

type Fe = [u32; 8];

const FE_P: Fe = [
    0xFFFFFC2F, 0xFFFFFFFE, 0xFFFFFFFF, 0xFFFFFFFF,
    0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF,
];
/* p - 2 (for Fermat inverse) */
const FE_PM2: Fe = [
    0xFFFFFC2D, 0xFFFFFFFE, 0xFFFFFFFF, 0xFFFFFFFF,
    0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF,
];
/* (p + 1) / 4 (for sqrt; p = 3 mod 4).
 * Verified limb by limb: (p+1)/4 = 2^254 - 2^30 - 244, so limb0 =
 * 0 - 0x40000000 - 0xF4 = 0xBFFFFF0C (borrow 1), limbs 1..6 = 0xFFFFFFFF,
 * limb7 = 0x40000000 - 1 = 0x3FFFFFFF. (Wrong limbs here were a real bug
 * in the C version's first draft.) */
const FE_SQRT_E: Fe = [
    0xBFFFFF0C, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF,
    0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0x3FFFFFFF,
];

fn fe_set_int(v: u32) -> Fe {
    let mut r = [0u32; 8];
    r[0] = v;
    r
}

/* 1 if a == 0 else 0, branch-free */
fn fe_is_zero(a: &Fe) -> u32 {
    let mut acc: u32 = 0;
    for i in 0..8 {
        acc |= a[i];
    }
    ((acc | acc.wrapping_neg()) >> 31) ^ 1
}

/* 1 if a == b else 0, branch-free */
fn fe_equal(a: &Fe, b: &Fe) -> u32 {
    let mut acc: u32 = 0;
    for i in 0..8 {
        acc |= a[i] ^ b[i];
    }
    ((acc | acc.wrapping_neg()) >> 31) ^ 1
}

/* flag ? a : r, flag in {0,1}, branch-free */
fn fe_cmov(r: &Fe, a: &Fe, flag: u32) -> Fe {
    let mask = flag.wrapping_neg();
    let mut out = [0u32; 8];
    for i in 0..8 {
        out[i] = (a[i] & mask) | (r[i] & !mask);
    }
    out
}

/* (a + b) mod p */
fn fe_add(a: &Fe, b: &Fe) -> Fe {
    let mut t = [0u32; 8];
    let mut c: u64 = 0;
    for i in 0..8 {
        let v = a[i] as u64 + b[i] as u64 + c;
        t[i] = v as u32;
        c = v >> 32;
    }
    let mut d = [0u32; 8];
    let mut br: u64 = 0;
    for i in 0..8 {
        let v = (t[i] as u64)
            .wrapping_sub(FE_P[i] as u64)
            .wrapping_sub(br);
        d[i] = v as u32;
        br = (v >> 32) & 1;
    }
    /* subtract p iff c==1 or t>=p (br==0) */
    let need = (c | (br ^ 1)) as u32;
    fe_cmov(&t, &d, need)
}

/* (a - b) mod p */
fn fe_sub(a: &Fe, b: &Fe) -> Fe {
    let mut t = [0u32; 8];
    let mut br: u64 = 0;
    for i in 0..8 {
        let v = (a[i] as u64)
            .wrapping_sub(b[i] as u64)
            .wrapping_sub(br);
        t[i] = v as u32;
        br = (v >> 32) & 1;
    }
    let mut d = [0u32; 8];
    let mut c: u64 = 0;
    for i in 0..8 {
        let v = t[i] as u64 + FE_P[i] as u64 + c;
        d[i] = v as u32;
        c = v >> 32;
    }
    let need = br as u32; /* borrow => add p back */
    fe_cmov(&t, &d, need)
}

/* (-a) mod p */
fn fe_neg(a: &Fe) -> Fe {
    fe_sub(&[0u32; 8], a)
}

/* (a * b) mod p. 5 fixed fold rounds. */
fn fe_mul(a: &Fe, b: &Fe) -> Fe {
    let mut acc = [0u128; 16];
    for i in 0..8 {
        for j in 0..8 {
            acc[i + j] += a[i] as u128 * b[j] as u128;
        }
    }

    let mut t = [0u32; 17];
    {
        let mut c: u64 = 0;
        for i in 0..16 {
            let v: u128 = acc[i] + c as u128;
            t[i] = v as u32;
            c = (v >> 32) as u64;
        }
        t[16] = c as u32;
    }
    /* Each fold: value V = lo + hi*2^256 -> lo + hi*(2^32+977).
     * V0 < 2^512; V1 < 2^290; V2 < 2^257; V3 < 2^256+2^34; V4 < 2^256+2^33. */
    for _iter in 0..5 {
        let mut u = [0u64; 11];
        for i in 0..8 {
            u[i] = t[i] as u64;
        }
        for i in 8..=16 {
            u[i - 8] = u[i - 8].wrapping_add(t[i] as u64 * 977);
            u[i - 7] = u[i - 7].wrapping_add(t[i] as u64);
        }
        let mut c: u64 = 0;
        for i in 0..11 {
            let v = u[i].wrapping_add(c);
            t[i] = v as u32;
            c = v >> 32;
        }
        t[11] = c as u32;
        for i in 12..=16 {
            t[i] = 0;
        }
    }
    let mut r = [0u32; 8];
    for i in 0..8 {
        r[i] = t[i];
    }
    /* Fold residual t[8] (< 8 after 5 folds; t[9..16] are 0).
     * r_true = r + t[8]*2^256 = r + t[8]*(2^32+977) (mod p).
     * (Dropping t[8] here was a real bug in the C version's first draft.) */
    {
        let t8 = t[8] as u64;
        let v0 = r[0] as u64 + t8 * 977;
        let v1 = r[1] as u64 + t8 + (v0 >> 32);
        r[0] = v0 as u32;
        r[1] = v1 as u32;
        let mut c = v1 >> 32;
        for i in 2..8 {
            let v = r[i] as u64 + c;
            r[i] = v as u32;
            c = v >> 32;
        }
        /* c < 8: fold c*2^256 once more (fixed, constant-time) */
        let v0 = r[0] as u64 + c * 977;
        let v1 = r[1] as u64 + c + (v0 >> 32);
        r[0] = v0 as u32;
        r[1] = v1 as u32;
        c = v1 >> 32;
        for i in 2..8 {
            let v = r[i] as u64 + c;
            r[i] = v as u32;
            c = v >> 32;
        }
        /* c is now 0; value < 2^256 + 2^34 */
    }
    /* single conditional subtract */
    let mut d = [0u32; 8];
    let mut br: u64 = 0;
    for i in 0..8 {
        let v = (r[i] as u64)
            .wrapping_sub(FE_P[i] as u64)
            .wrapping_sub(br);
        d[i] = v as u32;
        br = (v >> 32) & 1;
    }
    fe_cmov(&r, &d, (br ^ 1) as u32)
}

fn fe_sqr(a: &Fe) -> Fe {
    fe_mul(a, a)
}

/* a^(p-2) mod p (Fermat inverse), square-and-multiply-always */
fn fe_inv(a: &Fe) -> Fe {
    let mut result = fe_set_int(1);
    for b in (0..256).rev() {
        let tmp = fe_sqr(&result);
        let r2 = fe_mul(&tmp, a);
        let bit = (FE_PM2[b >> 5] >> (b & 31)) & 1;
        let tmp = fe_cmov(&tmp, &r2, bit); /* tmp = bit ? r2 : tmp */
        result = tmp;
    }
    result
}

/* sqrt(a) mod p for p = 3 mod 4: a^((p+1)/4), square-and-multiply-always */
fn fe_sqrt(a: &Fe) -> Fe {
    let mut result = fe_set_int(1);
    for b in (0..256).rev() {
        let tmp = fe_sqr(&result);
        let r2 = fe_mul(&tmp, a);
        let bit = (FE_SQRT_E[b >> 5] >> (b & 31)) & 1;
        let tmp = fe_cmov(&tmp, &r2, bit);
        result = tmp;
    }
    result
}

/* ================= group operations (Jacobian) ================= */

#[derive(Clone, Copy)]
struct Gej {
    x: Fe,
    y: Fe,
    z: Fe,
}

impl Gej {
    const fn zero() -> Gej {
        Gej {
            x: [0; 8],
            y: [0; 8],
            z: [0; 8],
        }
    }
}

/* 2*a. dbl-2009-l. Handles Z=0 (infinity) -> infinity, no branch. */
fn gej_dbl(a: &Gej) -> Gej {
    let xx = fe_sqr(&a.x);
    let yy = fe_sqr(&a.y);
    let mut yyyy = fe_sqr(&yy);
    let mut tmp = fe_add(&a.x, &yy);
    tmp = fe_sqr(&tmp);
    tmp = fe_sub(&tmp, &xx);
    tmp = fe_sub(&tmp, &yyyy);
    let s = fe_add(&tmp, &tmp);
    let mut m = fe_add(&xx, &xx);
    m = fe_add(&m, &xx);
    let mut t = fe_sqr(&m);
    tmp = fe_add(&s, &s);
    t = fe_sub(&t, &tmp);
    let y1z1 = fe_mul(&a.y, &a.z);
    let x_out = t;
    tmp = fe_sub(&s, &t);
    tmp = fe_mul(&m, &tmp);
    /* 8*YYYY. NOTE: yyyy must stay untouched until here -- it holds (Y^2)^2
     * needed for Y3 = M*(S-T) - 8*YYYY. Reusing it as scratch earlier was
     * a real bug in the C version's first draft. */
    yyyy = fe_add(&yyyy, &yyyy);
    yyyy = fe_add(&yyyy, &yyyy);
    yyyy = fe_add(&yyyy, &yyyy);
    let y_out = fe_sub(&tmp, &yyyy);
    let z_out = fe_add(&y1z1, &y1z1);
    Gej {
        x: x_out,
        y: y_out,
        z: z_out,
    }
}

/* a + b. Branch-free complete: infinity/doubling cases selected by masks. */
fn gej_add(a: &Gej, b: &Gej) -> Gej {
    let a_inf = fe_is_zero(&a.z);
    let b_inf = fe_is_zero(&b.z);

    let z1z1 = fe_sqr(&a.z);
    let z2z2 = fe_sqr(&b.z);
    let u1 = fe_mul(&a.x, &z2z2);
    let u2 = fe_mul(&b.x, &z1z1);
    let mut tmp = fe_mul(&a.y, &b.z);
    let s1 = fe_mul(&tmp, &z2z2);
    tmp = fe_mul(&b.y, &a.z);
    let s2 = fe_mul(&tmp, &z1z1);

    let u_eq = fe_equal(&u1, &u2);
    let s_eq = fe_equal(&s1, &s2);

    /* main formula (add-2009-bl) */
    let hh = fe_sub(&u2, &u1);
    tmp = fe_add(&hh, &hh);
    let ii = fe_sqr(&tmp);
    let jj = fe_mul(&hh, &ii);
    tmp = fe_sub(&s2, &s1);
    let rr = fe_add(&tmp, &tmp);
    let vv = fe_mul(&u1, &ii);
    tmp = fe_sqr(&rr);
    tmp = fe_sub(&tmp, &jj);
    let mut tmp2 = fe_add(&vv, &vv);
    let x3 = fe_sub(&tmp, &tmp2);
    tmp = fe_sub(&vv, &x3);
    tmp = fe_mul(&rr, &tmp);
    tmp2 = fe_mul(&s1, &jj);
    tmp2 = fe_add(&tmp2, &tmp2);
    let y3 = fe_sub(&tmp, &tmp2);
    tmp = fe_add(&a.z, &b.z);
    tmp = fe_sqr(&tmp);
    tmp = fe_sub(&tmp, &z1z1);
    tmp = fe_sub(&tmp, &z2z2);
    let z3 = fe_mul(&tmp, &hh);

    let dbl_a = gej_dbl(a);

    /* mutually exclusive selectors */
    let na = a_inf ^ 1;
    let nb = b_inf ^ 1;
    let m_a = a_inf; /* r = b */
    let m_b = b_inf & na; /* r = a */
    let m_d = u_eq & s_eq & na & nb; /* r = 2a */
    let m_i = u_eq & (s_eq ^ 1) & na & nb; /* r = INF */
    /* default: main formula (m_m implicit) */

    let mut rout = Gej { x: x3, y: y3, z: z3 };
    rout.x = fe_cmov(&rout.x, &dbl_a.x, m_d);
    rout.y = fe_cmov(&rout.y, &dbl_a.y, m_d);
    rout.z = fe_cmov(&rout.z, &dbl_a.z, m_d);
    {
        /* INF = (0,0,0) */
        let mask = m_i.wrapping_neg();
        for i in 0..8 {
            rout.x[i] &= !mask;
            rout.y[i] &= !mask;
            rout.z[i] &= !mask;
        }
    }
    rout.x = fe_cmov(&rout.x, &a.x, m_b);
    rout.y = fe_cmov(&rout.y, &a.y, m_b);
    rout.z = fe_cmov(&rout.z, &a.z, m_b);
    rout.x = fe_cmov(&rout.x, &b.x, m_a);
    rout.y = fe_cmov(&rout.y, &b.y, m_a);
    rout.z = fe_cmov(&rout.z, &b.z, m_a);
    rout
}

/* sum_i scalars[i] * points[i]. Fixed 256 iterations, masked adds.
 * The caller guarantees scalars.len() == points.len(). */
fn gej_msm(scalars: &[[u32; 8]], points: &[Gej]) -> Gej {
    let n = scalars.len();
    let mut acc = Gej::zero(); /* INF = (0,0,0) */
    for b in (0..256).rev() {
        acc = gej_dbl(&acc);
        for i in 0..n {
            /* public loop index b: read the bit, never branch on it */
            let bit = (scalars[i][b >> 5] >> (b & 31)) & 1;
            let mask = bit.wrapping_neg();
            let mut sel = Gej::zero();
            for k in 0..8 {
                sel.x[k] = points[i].x[k] & mask;
                sel.y[k] = points[i].y[k] & mask;
                sel.z[k] = points[i].z[k] & mask;
            }
            acc = gej_add(&acc, &sel);
        }
    }
    acc
}

/* ================= FFI boundary =================
 * The only `unsafe` in this crate: copying between raw pointers and owned
 * values. Every exported function reads its inputs into owned locals first,
 * so overlapping input/output pointers (in-place calls) are always safe. */

unsafe fn read_fe(p: *const u32) -> Fe {
    let mut f = [0u32; 8];
    for i in 0..8 {
        f[i] = *p.add(i);
    }
    f
}

unsafe fn read_gej(p: *const u32) -> Gej {
    Gej {
        x: read_fe(p),
        y: read_fe(p.add(8)),
        z: read_fe(p.add(16)),
    }
}

unsafe fn write_fe(p: *mut u32, f: &Fe) {
    for i in 0..8 {
        *p.add(i) = f[i];
    }
}

unsafe fn write_gej(p: *mut u32, g: &Gej) {
    write_fe(p, &g.x);
    write_fe(p.add(8), &g.y);
    write_fe(p.add(16), &g.z);
}

/* x = X/Z^2, y = Y/Z^3 mod p (via Fermat inverse). Caller: Z != 0. */
fn affine_of(a: &Gej) -> (Fe, Fe) {
    let zinv = fe_inv(&a.z);
    let zinv2 = fe_sqr(&zinv);
    let x = fe_mul(&a.x, &zinv2);
    let tmp = fe_mul(&zinv2, &zinv);
    let y = fe_mul(&a.y, &tmp);
    (x, y)
}

/* ================= exported C ABI (see core/nbext/ABI.md) ================= */

/// r = a + b. Branch-free complete formula.
#[no_mangle]
pub unsafe extern "C" fn wcurve_add(a: *const u32, b: *const u32, r: *mut u32) {
    let ga = read_gej(a);
    let gb = read_gej(b);
    let gr = gej_add(&ga, &gb);
    write_gej(r, &gr);
}

/// r = 2*a. dbl-2009-l. In-place safe.
#[no_mangle]
pub unsafe extern "C" fn wcurve_dbl(a: *const u32, r: *mut u32) {
    let ga = read_gej(a);
    let gr = gej_dbl(&ga);
    write_gej(r, &gr);
}

/// r = -a: negate Y in Jacobian (X, -Y mod p, Z). In-place safe.
#[no_mangle]
pub unsafe extern "C" fn wcurve_neg(a: *const u32, r: *mut u32) {
    let ga = read_gej(a);
    let gr = Gej {
        x: ga.x,
        y: fe_neg(&ga.y),
        z: ga.z,
    };
    write_gej(r, &gr);
}

/// 1 if Z == 0 (branch-free compare), else 0.
#[no_mangle]
pub unsafe extern "C" fn wcurve_is_inf(a: *const u32) -> i32 {
    let ga = read_gej(a);
    fe_is_zero(&ga.z) as i32
}

/// r = sum_{i=0}^{n-1} scalars[i] * points[i].
/// Fixed 256 iterations, masked conditional adds. n == 0 -> r = infinity.
/// In-place safe on r (inputs are copied before computing).
#[no_mangle]
pub unsafe extern "C" fn wcurve_msm(
    scalars: *const u32,
    points: *const u32,
    r: *mut u32,
    n: usize,
) {
    if n == 0 {
        /* do not touch the (possibly null) input pointers */
        write_gej(r, &Gej::zero());
        return;
    }
    let mut sc: Vec<[u32; 8]> = Vec::with_capacity(n);
    let mut pt: Vec<Gej> = Vec::with_capacity(n);
    for i in 0..n {
        sc.push(read_fe(scalars.add(i * 8)));
        pt.push(read_gej(points.add(i * 24)));
    }
    let gr = gej_msm(&sc, &pt);
    write_gej(r, &gr);
}

/// r = (x, y, 1). Assumes x, y < p (caller guarantees).
#[no_mangle]
pub unsafe extern "C" fn wcurve_from_affine(x: *const u32, y: *const u32, r: *mut u32) {
    let gr = Gej {
        x: read_fe(x),
        y: read_fe(y),
        z: [1, 0, 0, 0, 0, 0, 0, 0],
    };
    write_gej(r, &gr);
}

/// If Z == 0: *is_inf = 1 (x, y untouched). Else *is_inf = 0 and
/// x = X/Z^2, y = Y/Z^3 mod p (via Fermat inverse).
#[no_mangle]
pub unsafe extern "C" fn wcurve_to_affine(
    a: *const u32,
    x: *mut u32,
    y: *mut u32,
    is_inf: *mut i32,
) {
    let ga = read_gej(a);
    if fe_is_zero(&ga.z) == 1 {
        *is_inf = 1;
    } else {
        *is_inf = 0;
        let (xo, yo) = affine_of(&ga);
        write_fe(x, &xo);
        write_fe(y, &yo);
    }
}

/// out[0] = 0x02 | (y & 1); out[1..33] = x big-endian.
/// Caller guarantees a is not infinity. (Affine conversion inside.)
#[no_mangle]
pub unsafe extern "C" fn wcurve_compress(a: *const u32, out: *mut u8) {
    let ga = read_gej(a);
    let (x, y) = affine_of(&ga);
    *out = 0x02 | (y[0] & 1) as u8;
    for i in 0..8 {
        /* x is < p < 2^256; big-endian: most significant limb first.
         * (Byte order here was a real bug in the C version's first draft.) */
        let w = x[7 - i];
        *out.add(1 + i * 4 + 0) = (w >> 24) as u8;
        *out.add(1 + i * 4 + 1) = (w >> 16) as u8;
        *out.add(1 + i * 4 + 2) = (w >> 8) as u8;
        *out.add(1 + i * 4 + 3) = w as u8;
    }
}

/// Returns 0 on success, nonzero on error.
/// Checks: prefix byte, x < p, y^2 == x^3 + 7, parity match. r = (x, y, 1).
#[no_mangle]
pub unsafe extern "C" fn wcurve_decompress(inp: *const u8, r: *mut u32) -> i32 {
    let prefix = *inp;
    if prefix != 0x02 && prefix != 0x03 {
        return 1;
    }
    /* in[1..33] is x big-endian; x[0] is the least-significant limb */
    let mut x = [0u32; 8];
    for i in 0..8 {
        let mut v: u32 = 0;
        for j in 0..4 {
            v |= (*inp.add(1 + (7 - i) * 4 + j) as u32) << (8 * (3 - j));
        }
        x[i] = v;
    }
    /* check x < p: borrow out of (x - p) means x < p */
    {
        let mut br: u64 = 0;
        for i in 0..8 {
            let v = (x[i] as u64)
                .wrapping_sub(FE_P[i] as u64)
                .wrapping_sub(br);
            br = (v >> 32) & 1;
        }
        if br == 0 {
            return 2; /* x >= p */
        }
    }
    /* y2 = x^3 + 7 */
    let mut y2 = fe_sqr(&x);
    y2 = fe_mul(&y2, &x);
    let t7 = fe_set_int(7);
    y2 = fe_add(&y2, &t7);
    /* y = y2^((p+1)/4); verify y^2 == y2 */
    let mut y = fe_sqrt(&y2);
    let t = fe_sqr(&y);
    if fe_equal(&t, &y2) != 1 {
        return 3; /* not on curve */
    }
    if (y[0] & 1) != (prefix as u32 & 1) {
        y = fe_neg(&y);
    }
    let gr = Gej {
        x,
        y,
        z: [1, 0, 0, 0, 0, 0, 0, 0],
    };
    write_gej(r, &gr);
    0
}
