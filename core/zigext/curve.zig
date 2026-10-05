// Zig translation of core/cext/curve_ext.c (constant-time secp256k1).
//
// This is a line-by-line port of the C implementation's semantics:
//  - Field: 8 x 32-bit limbs, fully reduced mod p = 2^256 - 2^32 - 977.
//  - Group: Jacobian; dbl-2009-l; branch-free complete add (add-2009-bl
//    main formula + masked infinity/doubling selection).
//  - MSM: fixed 256 iterations, masked conditional adds.
//  - Inverse: Fermat, square-and-multiply-always.
//  - sqrt: exponent (p+1)/4, square-and-multiply-always (p = 3 mod 4).
//
// C ABI: see core/nbext/ABI.md. Points are 24 x u32 = X[8], Y[8], Z[8];
// all-zero (Z == 0) is infinity.
//
// Constant-time discipline: no secret-dependent branches, masked selects,
// fixed loop bounds. All limb arithmetic uses wrapping operators so the
// code cannot panic in any build mode.
//
// NOT audited. Same caveats as the C version: the compiler may introduce
// branches; verify the assembly if it matters.

const Fe = [8]u32;

const Gej = struct {
    x: Fe,
    y: Fe,
    z: Fe,
};

const FE_P: [8]u32 = .{
    0xFFFFFC2F, 0xFFFFFFFE, 0xFFFFFFFF, 0xFFFFFFFF,
    0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF,
};
// p - 2 (Fermat inverse exponent)
const FE_PM2: [8]u32 = .{
    0xFFFFFC2D, 0xFFFFFFFE, 0xFFFFFFFF, 0xFFFFFFFF,
    0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF,
};
// (p + 1) / 4 (sqrt exponent, p = 3 mod 4)
const FE_SQRT_E: [8]u32 = .{
    0xBFFFFF0C, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF,
    0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0x3FFFFFFF,
};

// ================= field arithmetic =================

fn fe_copy(r: *Fe, a: *const Fe) void {
    for (0..8) |i| r[i] = a[i];
}

fn fe_set_int(r: *Fe, v: u32) void {
    r[0] = v;
    for (1..8) |i| r[i] = 0;
}

// 1 if a == 0 else 0, branch-free
fn fe_is_zero(a: *const Fe) u32 {
    var acc: u32 = 0;
    for (0..8) |i| acc |= a[i];
    return (((acc | (0 -% acc)) >> 31) ^ 1);
}

// 1 if a == b else 0, branch-free
fn fe_equal(a: *const Fe, b: *const Fe) u32 {
    var acc: u32 = 0;
    for (0..8) |i| acc |= a[i] ^ b[i];
    return (((acc | (0 -% acc)) >> 31) ^ 1);
}

// r = flag ? a : r, flag in {0,1}, branch-free
fn fe_cmov(r: *Fe, a: *const Fe, flag: u32) void {
    const mask: u32 = 0 -% flag;
    for (0..8) |i| r[i] = (a[i] & mask) | (r[i] & ~mask);
}

// r = (a + b) mod p
fn fe_add(r: *Fe, a: *const Fe, b: *const Fe) void {
    var t: Fe = undefined;
    var d: Fe = undefined;
    var c: u64 = 0;
    for (0..8) |i| {
        const v: u64 = @as(u64, a[i]) +% @as(u64, b[i]) +% c;
        t[i] = @as(u32, @truncate(v));
        c = v >> 32;
    }
    var br: u64 = 0;
    for (0..8) |i| {
        const v: u64 = @as(u64, t[i]) -% @as(u64, FE_P[i]) -% br;
        d[i] = @as(u32, @truncate(v));
        br = (v >> 32) & 1;
    }
    // subtract p iff c==1 or t>=p (br==0)
    const need: u32 = @as(u32, @truncate(c)) | @as(u32, @truncate(br ^ 1));
    const m: u32 = 0 -% need;
    for (0..8) |i| r[i] = (d[i] & m) | (t[i] & ~m);
}

// r = (a - b) mod p
fn fe_sub(r: *Fe, a: *const Fe, b: *const Fe) void {
    var t: Fe = undefined;
    var d: Fe = undefined;
    var br: u64 = 0;
    for (0..8) |i| {
        const v: u64 = @as(u64, a[i]) -% @as(u64, b[i]) -% br;
        t[i] = @as(u32, @truncate(v));
        br = (v >> 32) & 1;
    }
    var c: u64 = 0;
    for (0..8) |i| {
        const v: u64 = @as(u64, t[i]) +% @as(u64, FE_P[i]) +% c;
        d[i] = @as(u32, @truncate(v));
        c = v >> 32;
    }
    const need: u32 = @as(u32, @truncate(br)); // borrow => add p back
    const m: u32 = 0 -% need;
    for (0..8) |i| r[i] = (d[i] & m) | (t[i] & ~m);
}

// r = (-a) mod p
fn fe_neg(r: *Fe, a: *const Fe) void {
    const zero: Fe = [_]u32{0} ** 8;
    fe_sub(r, &zero, a);
}

// r = (a * b) mod p. 5 fixed fold rounds.
// acc needs 128 bits: 8 products per limb -> up to 2^67.
fn fe_mul(r: *Fe, a: *const Fe, b: *const Fe) void {
    var acc: [16]u128 = [_]u128{0} ** 16;
    for (0..8) |i| {
        for (0..8) |j| {
            acc[i + j] +%= @as(u128, a[i]) * @as(u128, b[j]);
        }
    }

    var t: [17]u32 = [_]u32{0} ** 17;
    {
        var c: u64 = 0;
        for (0..16) |i| {
            const v: u128 = acc[i] +% @as(u128, c);
            t[i] = @as(u32, @truncate(v));
            c = @as(u64, @truncate(v >> 32));
        }
        t[16] = @as(u32, @truncate(c));
    }
    // Each fold: value V = lo + hi*2^256 -> lo + hi*(2^32+977).
    // V0 < 2^512; V1 < 2^290; V2 < 2^257; V3 < 2^256+2^34; V4 < 2^256+2^33.
    var iter: usize = 0;
    while (iter < 5) : (iter += 1) {
        var u: [11]u64 = [_]u64{0} ** 11;
        for (0..8) |i| u[i] = @as(u64, t[i]);
        for (8..17) |i| {
            u[i - 8] +%= @as(u64, t[i]) * 977;
            u[i - 7] +%= @as(u64, t[i]);
        }
        var c: u64 = 0;
        for (0..11) |i| {
            const v: u64 = u[i] +% c;
            t[i] = @as(u32, @truncate(v));
            c = v >> 32;
        }
        t[11] = @as(u32, @truncate(c));
        for (12..17) |i| t[i] = 0;
    }
    for (0..8) |i| r[i] = t[i];
    // Fold residual t[8] (< 8 after 5 folds; t[9..16] are 0).
    // r_true = r + t[8]*2^256 = r + t[8]*(2^32+977) (mod p).
    {
        const t8: u32 = t[8];
        var v0: u64 = @as(u64, r[0]) +% @as(u64, t8) *% 977;
        var v1: u64 = @as(u64, r[1]) +% @as(u64, t8) +% (v0 >> 32);
        r[0] = @as(u32, @truncate(v0));
        r[1] = @as(u32, @truncate(v1));
        var c: u64 = v1 >> 32;
        for (2..8) |i| {
            const v: u64 = @as(u64, r[i]) +% c;
            r[i] = @as(u32, @truncate(v));
            c = v >> 32;
        }
        // fold c*2^256 once more (fixed, constant-time)
        v0 = @as(u64, r[0]) +% c *% 977;
        v1 = @as(u64, r[1]) +% c +% (v0 >> 32);
        r[0] = @as(u32, @truncate(v0));
        r[1] = @as(u32, @truncate(v1));
        c = v1 >> 32;
        for (2..8) |i| {
            const v: u64 = @as(u64, r[i]) +% c;
            r[i] = @as(u32, @truncate(v));
            c = v >> 32;
        }
        // value < 2^256 + 2^34
    }
    // single conditional subtract
    var d: Fe = undefined;
    var br: u64 = 0;
    for (0..8) |i| {
        const v: u64 = @as(u64, r[i]) -% @as(u64, FE_P[i]) -% br;
        d[i] = @as(u32, @truncate(v));
        br = (v >> 32) & 1;
    }
    fe_cmov(r, &d, @as(u32, @truncate(br ^ 1)));
}

fn fe_sqr(r: *Fe, a: *const Fe) void {
    fe_mul(r, a, a);
}

// r = a^(p-2) mod p (Fermat inverse), square-and-multiply-always
fn fe_inv(r: *Fe, a: *const Fe) void {
    var result: Fe = undefined;
    var tmp: Fe = undefined;
    var r2: Fe = undefined;
    fe_set_int(&result, 1);
    var b: usize = 256;
    while (b > 0) {
        b -= 1;
        fe_sqr(&tmp, &result);
        fe_mul(&r2, &tmp, a);
        const bit: u32 = (FE_PM2[b >> 5] >> @as(u5, @intCast(b & 31))) & 1;
        fe_cmov(&tmp, &r2, bit); // tmp = bit ? r2 : tmp
        fe_copy(&result, &tmp);
    }
    fe_copy(r, &result);
}

// r = sqrt(a) mod p for p = 3 mod 4: a^((p+1)/4), square-and-multiply-always
fn fe_sqrt(r: *Fe, a: *const Fe) void {
    var result: Fe = undefined;
    var tmp: Fe = undefined;
    var r2: Fe = undefined;
    fe_set_int(&result, 1);
    var b: usize = 256;
    while (b > 0) {
        b -= 1;
        fe_sqr(&tmp, &result);
        fe_mul(&r2, &tmp, a);
        const bit: u32 = (FE_SQRT_E[b >> 5] >> @as(u5, @intCast(b & 31))) & 1;
        fe_cmov(&tmp, &r2, bit);
        fe_copy(&result, &tmp);
    }
    fe_copy(r, &result);
}

// ================= group operations (Jacobian) =================

fn load_gej(p: [*]const u32) Gej {
    var g: Gej = undefined;
    for (0..8) |i| {
        g.x[i] = p[i];
        g.y[i] = p[8 + i];
        g.z[i] = p[16 + i];
    }
    return g;
}

fn store_gej(p: [*]u32, g: *const Gej) void {
    for (0..8) |i| {
        p[i] = g.x[i];
        p[8 + i] = g.y[i];
        p[16 + i] = g.z[i];
    }
}

// r = 2*a. dbl-2009-l. Handles Z=0 (infinity) -> infinity, no branch.
// In-place safe (r may equal a): a is fully read before any write to r.
fn gej_dbl(r: *Gej, a: *const Gej) void {
    var xx: Fe = undefined;
    var yy: Fe = undefined;
    var yyyy: Fe = undefined;
    var s: Fe = undefined;
    var m: Fe = undefined;
    var t: Fe = undefined;
    var tmp: Fe = undefined;
    var y1z1: Fe = undefined;
    fe_sqr(&xx, &a.x);
    fe_sqr(&yy, &a.y);
    fe_sqr(&yyyy, &yy);
    fe_add(&tmp, &a.x, &yy);
    fe_sqr(&tmp, &tmp);
    fe_sub(&tmp, &tmp, &xx);
    fe_sub(&tmp, &tmp, &yyyy);
    fe_add(&s, &tmp, &tmp);
    fe_add(&m, &xx, &xx);
    fe_add(&m, &m, &xx);
    fe_sqr(&t, &m);
    fe_add(&tmp, &s, &s);
    fe_sub(&t, &t, &tmp);
    fe_mul(&y1z1, &a.y, &a.z); // read a before any write to r
    fe_copy(&r.x, &t);
    fe_sub(&tmp, &s, &t);
    fe_mul(&tmp, &m, &tmp);
    // 8*yyyy (yyyy untouched until here)
    fe_add(&yyyy, &yyyy, &yyyy);
    fe_add(&yyyy, &yyyy, &yyyy);
    fe_add(&yyyy, &yyyy, &yyyy);
    fe_sub(&r.y, &tmp, &yyyy);
    fe_add(&r.z, &y1z1, &y1z1);
}

// r = a + b. Branch-free complete: infinity/doubling selected by masks.
// In-place safe (r may equal a or b): writes go to rout first.
fn gej_add(r: *Gej, a: *const Gej, b: *const Gej) void {
    const a_inf: u32 = fe_is_zero(&a.z);
    const b_inf: u32 = fe_is_zero(&b.z);

    var z1z1: Fe = undefined;
    var z2z2: Fe = undefined;
    var uu1: Fe = undefined;
    var uu2: Fe = undefined;
    var s1: Fe = undefined;
    var s2: Fe = undefined;
    var hh: Fe = undefined;
    var ii: Fe = undefined;
    var jj: Fe = undefined;
    var rr: Fe = undefined;
    var vv: Fe = undefined;
    var x3: Fe = undefined;
    var y3: Fe = undefined;
    var z3: Fe = undefined;
    var tmp: Fe = undefined;
    var tmp2: Fe = undefined;
    fe_sqr(&z1z1, &a.z);
    fe_sqr(&z2z2, &b.z);
    fe_mul(&uu1, &a.x, &z2z2);
    fe_mul(&uu2, &b.x, &z1z1);
    fe_mul(&tmp, &a.y, &b.z);
    fe_mul(&s1, &tmp, &z2z2);
    fe_mul(&tmp, &b.y, &a.z);
    fe_mul(&s2, &tmp, &z1z1);

    const u_eq: u32 = fe_equal(&uu1, &uu2);
    const s_eq: u32 = fe_equal(&s1, &s2);

    // main formula (add-2009-bl)
    fe_sub(&hh, &uu2, &uu1);
    fe_add(&tmp, &hh, &hh);
    fe_sqr(&ii, &tmp);
    fe_mul(&jj, &hh, &ii);
    fe_sub(&tmp, &s2, &s1);
    fe_add(&rr, &tmp, &tmp);
    fe_mul(&vv, &uu1, &ii);
    fe_sqr(&tmp, &rr);
    fe_sub(&tmp, &tmp, &jj);
    fe_add(&tmp2, &vv, &vv);
    fe_sub(&x3, &tmp, &tmp2);
    fe_sub(&tmp, &vv, &x3);
    fe_mul(&tmp, &rr, &tmp);
    fe_mul(&tmp2, &s1, &jj);
    fe_add(&tmp2, &tmp2, &tmp2);
    fe_sub(&y3, &tmp, &tmp2);
    fe_add(&tmp, &a.z, &b.z);
    fe_sqr(&tmp, &tmp);
    fe_sub(&tmp, &tmp, &z1z1);
    fe_sub(&tmp, &tmp, &z2z2);
    fe_mul(&z3, &tmp, &hh);

    var dbl_a: Gej = undefined;
    gej_dbl(&dbl_a, a);

    // mutually exclusive selectors
    const na: u32 = a_inf ^ 1;
    const nb: u32 = b_inf ^ 1;
    const m_a: u32 = a_inf; // r = b
    const m_b: u32 = b_inf & na; // r = a
    const m_d: u32 = u_eq & s_eq & na & nb; // r = 2a
    const m_i: u32 = u_eq & (s_eq ^ 1) & na & nb; // r = INF
    // default: main formula (m_m implicit)

    // write to temp first for in-place safety (r may equal a or b)
    var rout: Gej = undefined;
    fe_copy(&rout.x, &x3);
    fe_copy(&rout.y, &y3);
    fe_copy(&rout.z, &z3);
    fe_cmov(&rout.x, &dbl_a.x, m_d);
    fe_cmov(&rout.y, &dbl_a.y, m_d);
    fe_cmov(&rout.z, &dbl_a.z, m_d);
    { // INF = (0,0,0)
        const mask: u32 = 0 -% m_i;
        for (0..8) |i| {
            rout.x[i] &= ~mask;
            rout.y[i] &= ~mask;
            rout.z[i] &= ~mask;
        }
    }
    fe_cmov(&rout.x, &a.x, m_b);
    fe_cmov(&rout.y, &a.y, m_b);
    fe_cmov(&rout.z, &a.z, m_b);
    fe_cmov(&rout.x, &b.x, m_a);
    fe_cmov(&rout.y, &b.y, m_a);
    fe_cmov(&rout.z, &b.z, m_a);
    r.* = rout;
}

// affine conversion shared by compress/to_affine (caller: not infinity)
fn gej_to_affine(a: *const Gej, x: *Fe, y: *Fe) void {
    var zinv: Fe = undefined;
    var zinv2: Fe = undefined;
    var t: Fe = undefined;
    fe_inv(&zinv, &a.z);
    fe_sqr(&zinv2, &zinv);
    fe_mul(x, &a.x, &zinv2);
    fe_mul(&t, &zinv2, &zinv);
    fe_mul(y, &a.y, &t);
}

// ================= exported C ABI (see core/nbext/ABI.md) =================

export fn wcurve_add(a: [*]const u32, b: [*]const u32, r: [*]u32) void {
    const pa = load_gej(a);
    const pb = load_gej(b);
    var pr: Gej = undefined;
    gej_add(&pr, &pa, &pb);
    store_gej(r, &pr);
}

export fn wcurve_dbl(a: [*]const u32, r: [*]u32) void {
    const pa = load_gej(a);
    var pr: Gej = undefined;
    gej_dbl(&pr, &pa);
    store_gej(r, &pr);
}

export fn wcurve_neg(a: [*]const u32, r: [*]u32) void {
    var pa = load_gej(a);
    var ny: Fe = undefined;
    fe_neg(&ny, &pa.y);
    pa.y = ny;
    store_gej(r, &pa);
}

export fn wcurve_is_inf(a: [*]const u32) c_int {
    var z: Fe = undefined;
    for (0..8) |i| z[i] = a[16 + i];
    return @as(c_int, @intCast(fe_is_zero(&z)));
}

// r = sum_i scalars[i] * points[i]. Fixed 256 iterations, masked adds.
// n == 0 -> r = infinity. In-place safe on r (written once at the end).
export fn wcurve_msm(scalars: [*]const u32, points: [*]const u32, r: [*]u32, n: usize) void {
    var acc: Gej = .{
        .x = [_]u32{0} ** 8,
        .y = [_]u32{0} ** 8,
        .z = [_]u32{0} ** 8,
    };
    var b: usize = 256;
    while (b > 0) {
        b -= 1;
        gej_dbl(&acc, &acc);
        var i: usize = 0;
        while (i < n) : (i += 1) {
            // public loop index, never branched on
            const bit: u32 = (scalars[i *% 8 +% (b >> 5)] >> @as(u5, @intCast(b & 31))) & 1;
            const mask: u32 = 0 -% bit;
            var sel: Gej = undefined;
            const base = i *% 24;
            for (0..8) |k| {
                sel.x[k] = points[base +% k] & mask;
                sel.y[k] = points[base +% 8 +% k] & mask;
                sel.z[k] = points[base +% 16 +% k] & mask;
            }
            var new_acc: Gej = undefined;
            gej_add(&new_acc, &acc, &sel);
            acc = new_acc;
        }
    }
    store_gej(r, &acc);
}

export fn wcurve_from_affine(x: [*]const u32, y: [*]const u32, r: [*]u32) void {
    // copy via locals so r may alias x or y
    var lx: Fe = undefined;
    var ly: Fe = undefined;
    for (0..8) |i| {
        lx[i] = x[i];
        ly[i] = y[i];
    }
    for (0..8) |i| {
        r[i] = lx[i];
        r[8 + i] = ly[i];
        r[16 + i] = 0;
    }
    r[16] = 1; // Z = 1
}

export fn wcurve_to_affine(a: [*]const u32, x: [*]u32, y: [*]u32, is_inf: *c_int) void {
    const pa = load_gej(a);
    var xa: Fe = undefined;
    var ya: Fe = undefined;
    gej_to_affine(&pa, &xa, &ya); // safe on Z=0: yields (0,0), masked out below
    const inf: u32 = fe_is_zero(&pa.z);
    is_inf.* = @as(c_int, @intCast(inf));
    // masked write: x, y untouched when infinity (branch-free)
    const mask: u32 = 0 -% (inf ^ 1);
    for (0..8) |i| {
        x[i] = (xa[i] & mask) | (x[i] & ~mask);
        y[i] = (ya[i] & mask) | (y[i] & ~mask);
    }
}

export fn wcurve_compress(a: [*]const u32, out: [*]u8) void {
    const pa = load_gej(a);
    var xa: Fe = undefined;
    var ya: Fe = undefined;
    gej_to_affine(&pa, &xa, &ya);
    out[0] = 0x02 | @as(u8, @truncate(ya[0] & 1));
    for (0..8) |i| {
        // x big-endian: most significant limb first
        const w: u32 = xa[7 - i];
        out[1 + i * 4 + 0] = @as(u8, @truncate(w >> 24));
        out[1 + i * 4 + 1] = @as(u8, @truncate(w >> 16));
        out[1 + i * 4 + 2] = @as(u8, @truncate(w >> 8));
        out[1 + i * 4 + 3] = @as(u8, @truncate(w));
    }
}

// Returns 0 on success, nonzero on error.
export fn wcurve_decompress(in: [*]const u8, r: [*]u32) c_int {
    const prefix: u8 = in[0];
    if (prefix != 0x02 and prefix != 0x03) return 1;
    var x: Fe = undefined;
    for (0..8) |i| {
        // in[1..33] is x big-endian; x[0] is the least-significant limb
        var v: u32 = 0;
        for (0..4) |j| {
            v |= @as(u32, in[1 +% (7 - i) *% 4 +% j]) << @as(u5, @intCast(8 * (3 - j)));
        }
        x[i] = v;
    }
    // check x < p
    {
        var br: u64 = 0;
        for (0..8) |i| {
            const v: u64 = @as(u64, x[i]) -% @as(u64, FE_P[i]) -% br;
            br = (v >> 32) & 1;
        }
        if (br == 0) return 2; // x >= p
    }
    var y2: Fe = undefined;
    var y: Fe = undefined;
    var t: Fe = undefined;
    fe_sqr(&y2, &x);
    fe_mul(&y2, &y2, &x);
    fe_set_int(&t, 7);
    fe_add(&y2, &y2, &t); // y2 = x^3 + 7
    fe_sqrt(&y, &y2);
    fe_sqr(&t, &y);
    if (fe_equal(&t, &y2) == 0) return 3; // not on curve
    // parity adjust, branch-free
    const flip: u32 = (y[0] & 1) ^ @as(u32, prefix & 1);
    var ny: Fe = undefined;
    fe_neg(&ny, &y);
    fe_cmov(&y, &ny, flip);
    // r = (x, y, 1)
    for (0..8) |i| {
        r[i] = x[i];
        r[8 + i] = y[i];
        r[16 + i] = 0;
    }
    r[16] = 1;
    return 0;
}
