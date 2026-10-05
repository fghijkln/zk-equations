/* Constant-time secp256k1 group operations (debug/experimental).
 *
 * Purpose: answer the "pure Python can't be constant-time" limitation with a
 * C implementation where the hot path (MSM / scalar mult) runs a fixed
 * operation sequence with no secret-dependent branches.
 *
 * Design:
 *  - Field: 8 x 32-bit limbs, fully reduced mod p after every op.
 *    All field ops are straight-line limb arithmetic + masked selects.
 *  - Group: Jacobian coordinates. Addition uses a branch-free complete
 *    formula (all cases -- infinity, doubling -- selected via masks).
 *    Doubling is dbl-2009-l (handles Z=0 correctly, no branch).
 *  - MSM: fixed 256 iterations, masked conditional adds. The scalar bits
 *    are read (public loop index) but never branched on.
 *  - Scalars are 256-bit (8 x 32-bit limbs).
 *
 * What is NOT claimed: this has not been independently audited. The
 * constant-time property holds for the operation sequence; it does not
 * cover Python-level overhead, memory-access patterns in table lookups
 * (there are none -- no tables), or compiler-introduced branches
 * (build with -O2 and verify the assembly if this matters to you).
 */

#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <stdint.h>
#include <string.h>


/* ================= field arithmetic mod p = 2^256 - 2^32 - 977 ================= */

typedef uint32_t fe[8];

static const uint32_t FE_P[8] = {
    0xFFFFFC2Fu, 0xFFFFFFFEu, 0xFFFFFFFFu, 0xFFFFFFFFu,
    0xFFFFFFFFu, 0xFFFFFFFFu, 0xFFFFFFFFu, 0xFFFFFFFFu
};
/* p - 2 (for Fermat inverse) */
static const uint32_t FE_PM2[8] = {
    0xFFFFFC2Du, 0xFFFFFFFEu, 0xFFFFFFFFu, 0xFFFFFFFFu,
    0xFFFFFFFFu, 0xFFFFFFFFu, 0xFFFFFFFFu, 0xFFFFFFFFu
};

static void fe_copy(fe r, const fe a) {
    for (int i = 0; i < 8; i++) r[i] = a[i];
}

static void fe_set_int(fe r, uint32_t v) {
    r[0] = v;
    for (int i = 1; i < 8; i++) r[i] = 0;
}

/* 1 if a == 0 else 0, branch-free */
static uint32_t fe_is_zero(const fe a) {
    uint32_t acc = 0;
    for (int i = 0; i < 8; i++) acc |= a[i];
    return (((acc | (0u - acc)) >> 31) ^ 1u);
}

/* 1 if a == b else 0, branch-free */
static uint32_t fe_equal(const fe a, const fe b) {
    uint32_t acc = 0;
    for (int i = 0; i < 8; i++) acc |= a[i] ^ b[i];
    return (((acc | (0u - acc)) >> 31) ^ 1u);
}

/* r = flag ? a : r, flag in {0,1}, branch-free */
static void fe_cmov(fe r, const fe a, uint32_t flag) {
    uint32_t mask = 0u - flag;
    for (int i = 0; i < 8; i++) r[i] = (a[i] & mask) | (r[i] & ~mask);
}

/* r = (a + b) mod p */
static void fe_add(fe r, const fe a, const fe b) {
    uint32_t t[8], d[8];
    uint64_t c = 0;
    for (int i = 0; i < 8; i++) {
        uint64_t v = (uint64_t)a[i] + b[i] + c;
        t[i] = (uint32_t)v; c = v >> 32;
    }
    uint64_t br = 0;
    for (int i = 0; i < 8; i++) {
        uint64_t v = (uint64_t)t[i] - FE_P[i] - br;
        d[i] = (uint32_t)v; br = (v >> 32) & 1u;
    }
    /* subtract p iff c==1 or t>=p (br==0) */
    uint32_t need = (uint32_t)(c | (br ^ 1u));
    for (int i = 0; i < 8; i++) r[i] = (d[i] & (0u - need)) | (t[i] & ~(0u - need));
}

/* r = (a - b) mod p */
static void fe_sub(fe r, const fe a, const fe b) {
    uint32_t t[8], d[8];
    uint64_t br = 0;
    for (int i = 0; i < 8; i++) {
        uint64_t v = (uint64_t)a[i] - b[i] - br;
        t[i] = (uint32_t)v; br = (v >> 32) & 1u;
    }
    uint64_t c = 0;
    for (int i = 0; i < 8; i++) {
        uint64_t v = (uint64_t)t[i] + FE_P[i] + c;
        d[i] = (uint32_t)v; c = v >> 32;
    }
    uint32_t need = (uint32_t)br; /* borrow => add p back */
    for (int i = 0; i < 8; i++) r[i] = (d[i] & (0u - need)) | (t[i] & ~(0u - need));
}

/* r = (-a) mod p */
static void fe_neg(fe r, const fe a) {
    static const fe zero = {0,0,0,0,0,0,0,0};
    fe_sub(r, zero, a);
}

/* r = (a * b) mod p. 5 fixed fold rounds; see docs for the bound. */
static void fe_mul(fe r, const fe a, const fe b) {
    /* 8 products per limb -> up to 2^67, needs 128-bit accumulator */
    unsigned __int128 acc[16] = {0};
    for (int i = 0; i < 8; i++)
        for (int j = 0; j < 8; j++)
            acc[i+j] += (unsigned __int128)a[i] * b[j];

    uint32_t t[17];
    {
        uint64_t c = 0;
        for (int i = 0; i < 16; i++) {
            unsigned __int128 v = acc[i] + c;
            t[i] = (uint32_t)v; c = (uint64_t)(v >> 32);
        }
        t[16] = (uint32_t)c;
    }
    /* Each fold: value V = lo + hi*2^256 -> lo + hi*(2^32+977).
     * V0 < 2^512; V1 < 2^290; V2 < 2^257; V3 < 2^256+2^34; V4 < 2^256+2^33. */
    for (int iter = 0; iter < 5; iter++) {
        uint64_t u[11] = {0};
        for (int i = 0; i < 8; i++) u[i] = t[i];
        for (int i = 8; i <= 16; i++) {
            u[i-8] += (uint64_t)t[i] * 977u;
            u[i-7] += (uint64_t)t[i];
        }
        uint64_t c = 0;
        for (int i = 0; i < 11; i++) {
            uint64_t v = u[i] + c;
            t[i] = (uint32_t)v; c = v >> 32;
        }
        t[11] = (uint32_t)c;
        for (int i = 12; i <= 16; i++) t[i] = 0;
    }
    for (int i = 0; i < 8; i++) r[i] = t[i];
    /* Fold residual t[8] (< 8 after 5 folds; t[9..16] are 0).
     * r_true = r + t[8]*2^256 = r + t[8]*(2^32+977) (mod p). */
    {
        uint32_t t8 = t[8];
        uint64_t v0 = (uint64_t)r[0] + (uint64_t)t8 * 977u;
        uint64_t v1 = (uint64_t)r[1] + (uint64_t)t8 + (v0 >> 32);
        r[0] = (uint32_t)v0;
        r[1] = (uint32_t)v1;
        uint64_t c = v1 >> 32;
        for (int i = 2; i < 8; i++) {
            uint64_t v = (uint64_t)r[i] + c;
            r[i] = (uint32_t)v;
            c = v >> 32;
        }
        /* c < 8: fold c*2^256 once more (fixed, constant-time) */
        v0 = (uint64_t)r[0] + c * 977u;
        v1 = (uint64_t)r[1] + c + (v0 >> 32);
        r[0] = (uint32_t)v0;
        r[1] = (uint32_t)v1;
        c = v1 >> 32;
        for (int i = 2; i < 8; i++) {
            uint64_t v = (uint64_t)r[i] + c;
            r[i] = (uint32_t)v;
            c = v >> 32;
        }
        /* c is now 0; value < 2^256 + 2^34 */
    }
    /* single conditional subtract */
    uint32_t d[8];
    uint64_t br = 0;
    for (int i = 0; i < 8; i++) {
        uint64_t v = (uint64_t)r[i] - FE_P[i] - br;
        d[i] = (uint32_t)v; br = (v >> 32) & 1u;
    }
    fe_cmov(r, d, (uint32_t)(br ^ 1u));
}

static void fe_sqr(fe r, const fe a) { fe_mul(r, a, a); }

/* r = a^(p-2) mod p (Fermat inverse), square-and-multiply-always */
static void fe_inv(fe r, const fe a) {
    fe result, tmp, r2;
    fe_set_int(result, 1);
    for (int b = 255; b >= 0; b--) {
        fe_sqr(tmp, result);
        fe_mul(r2, tmp, a);
        uint32_t bit = (FE_PM2[b >> 5] >> (b & 31)) & 1u;
        fe_cmov(tmp, r2, bit); /* tmp = bit ? r2 : tmp */
        fe_copy(result, tmp);
    }
    fe_copy(r, result);
}

/* ================= group operations (Jacobian) ================= */

typedef struct { fe X, Y, Z; } gej;

static void gej_copy(gej *r, const gej *a) {
    fe_copy(r->X, a->X); fe_copy(r->Y, a->Y); fe_copy(r->Z, a->Z);
}

/* r = 2*a. dbl-2009-l. Handles Z=0 (infinity) -> infinity, no branch. */
/* r = 2*a. In-place safe (r may equal a). */
static void gej_dbl(gej *r, const gej *a) {
    fe XX, YY, YYYY, S, M, T, tmp, Y1Z1;
    fe_sqr(XX, a->X);
    fe_sqr(YY, a->Y);
    fe_sqr(YYYY, YY);
    fe_add(tmp, a->X, YY); fe_sqr(tmp, tmp);
    fe_sub(tmp, tmp, XX); fe_sub(tmp, tmp, YYYY);
    fe_add(S, tmp, tmp);
    fe_add(M, XX, XX); fe_add(M, M, XX);
    fe_sqr(T, M);
    fe_add(tmp, S, S); fe_sub(T, T, tmp);
    fe_mul(Y1Z1, a->Y, a->Z);   /* read a before any write to r */
    fe_copy(r->X, T);
    fe_sub(tmp, S, T); fe_mul(tmp, M, tmp);
    /* 8*YYYY (YYYY untouched until here) */
    fe_add(YYYY, YYYY, YYYY);
    fe_add(YYYY, YYYY, YYYY);
    fe_add(YYYY, YYYY, YYYY);
    fe_sub(r->Y, tmp, YYYY);
    fe_add(r->Z, Y1Z1, Y1Z1);
}

/* r = a + b. Branch-free complete: infinity/doubling selected by masks. */
static void gej_add(gej *r, const gej *a, const gej *b) {
    uint32_t a_inf = fe_is_zero(a->Z);
    uint32_t b_inf = fe_is_zero(b->Z);

    fe Z1Z1, Z2Z2, U1, U2, S1, S2, Hh, Ii, Jj, Rr, Vv, X3, Y3, Z3, tmp, tmp2;
    fe_sqr(Z1Z1, a->Z); fe_sqr(Z2Z2, b->Z);
    fe_mul(U1, a->X, Z2Z2); fe_mul(U2, b->X, Z1Z1);
    fe_mul(tmp, a->Y, b->Z); fe_mul(S1, tmp, Z2Z2);
    fe_mul(tmp, b->Y, a->Z); fe_mul(S2, tmp, Z1Z1);

    uint32_t u_eq = fe_equal(U1, U2);
    uint32_t s_eq = fe_equal(S1, S2);

    /* main formula (add-2009-bl) */
    fe_sub(Hh, U2, U1);
    fe_add(tmp, Hh, Hh); fe_sqr(Ii, tmp);
    fe_mul(Jj, Hh, Ii);
    fe_sub(tmp, S2, S1); fe_add(Rr, tmp, tmp);
    fe_mul(Vv, U1, Ii);
    fe_sqr(tmp, Rr);
    fe_sub(tmp, tmp, Jj);
    fe_add(tmp2, Vv, Vv);
    fe_sub(X3, tmp, tmp2);
    fe_sub(tmp, Vv, X3);
    fe_mul(tmp, Rr, tmp);
    fe_mul(tmp2, S1, Jj);
    fe_add(tmp2, tmp2, tmp2);
    fe_sub(Y3, tmp, tmp2);
    fe_add(tmp, a->Z, b->Z);
    fe_sqr(tmp, tmp);
    fe_sub(tmp, tmp, Z1Z1);
    fe_sub(tmp, tmp, Z2Z2);
    fe_mul(Z3, tmp, Hh);

    gej dbl_a;
    gej_dbl(&dbl_a, a);

    /* mutually exclusive selectors */
    uint32_t na = a_inf ^ 1u, nb = b_inf ^ 1u;
    uint32_t m_a = a_inf;                 /* r = b */
    uint32_t m_b = b_inf & na;            /* r = a */
    uint32_t m_d = u_eq & s_eq & na & nb; /* r = 2a */
    uint32_t m_i = u_eq & (s_eq ^ 1u) & na & nb; /* r = INF */
    /* default: main formula (m_m implicit) */

    /* write to temp first for in-place safety (r may equal a or b) */
    gej rout;
    fe_copy(rout.X, X3); fe_copy(rout.Y, Y3); fe_copy(rout.Z, Z3);
    fe_cmov(rout.X, dbl_a.X, m_d); fe_cmov(rout.Y, dbl_a.Y, m_d); fe_cmov(rout.Z, dbl_a.Z, m_d);
    { /* INF = (0,0,0) */
        uint32_t mask = 0u - m_i;
        for (int i = 0; i < 8; i++) {
            rout.X[i] &= ~mask; rout.Y[i] &= ~mask; rout.Z[i] &= ~mask;
        }
    }
    fe_cmov(rout.X, a->X, m_b); fe_cmov(rout.Y, a->Y, m_b); fe_cmov(rout.Z, a->Z, m_b);
    fe_cmov(rout.X, b->X, m_a); fe_cmov(rout.Y, b->Y, m_a); fe_cmov(rout.Z, b->Z, m_a);
    gej_copy(r, &rout);
}

/* r = sum_i scalars[i] * points[i]. Fixed 256 iterations, masked adds. */
static void gej_msm(gej *r, const uint32_t scalars[][8], const gej *points, int n) {
    gej acc;
    memset(&acc, 0, sizeof(acc)); /* INF = (0,0,0) */
    for (int b = 255; b >= 0; b--) {
        gej_dbl(&acc, &acc);
        for (int i = 0; i < n; i++) {
            uint32_t bit = (scalars[i][b >> 5] >> (b & 31)) & 1u;
            gej sel, new_acc;
            uint32_t mask = 0u - bit;
            for (int k = 0; k < 8; k++) {
                sel.X[k] = points[i].X[k] & mask;
                sel.Y[k] = points[i].Y[k] & mask;
                sel.Z[k] = points[i].Z[k] & mask;
            }
            gej_add(&new_acc, &acc, &sel);
            gej_copy(&acc, &new_acc);
        }
    }
    gej_copy(r, &acc);
}

/* ================= Python bindings ================= */

/* Python (x, y) tuple -> fe; returns 0 on success. */
static int py_to_fe(fe r, PyObject *o) {
    if (!PyLong_Check(o)) return -1;
    /* convert via 32-bit limbs */
    PyObject *tmp = PyNumber_Long(o);
    if (!tmp) return -1;
    for (int i = 0; i < 8; i++) {
        PyObject *shift = PyLong_FromLong(i * 32);
        PyObject *sh = PyNumber_Rshift(tmp, shift);
        Py_DECREF(shift);
        PyObject *mask = PyLong_FromLong(0xFFFFFFFF);
        PyObject *lo = PyNumber_And(sh, mask);
        Py_DECREF(mask);
        unsigned long v = PyLong_AsUnsignedLong(lo);
        Py_DECREF(sh); Py_DECREF(lo);
        if (PyErr_Occurred()) { Py_DECREF(tmp); return -1; }
        r[i] = (uint32_t)v;
    }
    Py_DECREF(tmp);
    return 0;
}

/* fe -> Python int */
static PyObject *fe_to_py(const fe a) {
    PyObject *r = PyLong_FromLong(0);
    PyObject *c32 = PyLong_FromLong(32);
    for (int i = 7; i >= 0; i--) {
        PyObject *sh = PyNumber_Lshift(r, c32);
        PyObject *nw = PyNumber_Or(sh, PyLong_FromUnsignedLong(a[i]));
        Py_DECREF(sh); Py_SETREF(r, nw);
    }
    Py_DECREF(c32);
    return r;
}

/* Python point ((x,y) tuple or None) -> gej. Returns 0 ok, -1 error. */
static int py_to_gej(gej *r, PyObject *o) {
    if (o == Py_None) {
        memset(r, 0, sizeof(*r));
        return 0;
    }
    if (!PyTuple_Check(o) || PyTuple_Size(o) != 2) return -1;
    if (py_to_fe(r->X, PyTuple_GetItem(o, 0))) return -1;
    if (py_to_fe(r->Y, PyTuple_GetItem(o, 1))) return -1;
    fe_set_int(r->Z, 1);
    return 0;
}

/* gej -> Python ((x,y) tuple or None). Does the affine conversion. */
static PyObject *gej_to_py(const gej *a) {
    if (fe_is_zero(a->Z)) Py_RETURN_NONE;
    fe zinv, zinv2, x, y, tmp;
    fe_inv(zinv, a->Z);
    fe_sqr(zinv2, zinv);
    fe_mul(x, a->X, zinv2);
    fe_mul(tmp, zinv2, zinv);
    fe_mul(y, a->Y, tmp);
    PyObject *px = fe_to_py(x), *py = fe_to_py(y);
    PyObject *t = PyTuple_Pack(2, px, py);
    Py_DECREF(px); Py_DECREF(py);
    return t;
}

static PyObject *meth_add(PyObject *self, PyObject *args) {
    PyObject *p1, *p2;
    if (!PyArg_ParseTuple(args, "OO", &p1, &p2)) return NULL;
    gej a, b, r;
    if (py_to_gej(&a, p1) || py_to_gej(&b, p2)) {
        PyErr_SetString(PyExc_ValueError, "bad point");
        return NULL;
    }
    gej_add(&r, &a, &b);
    return gej_to_py(&r);
}

static PyObject *meth_dbl(PyObject *self, PyObject *args) {
    PyObject *p;
    if (!PyArg_ParseTuple(args, "O", &p)) return NULL;
    gej a, r;
    if (py_to_gej(&a, p)) {
        PyErr_SetString(PyExc_ValueError, "bad point");
        return NULL;
    }
    gej_dbl(&r, &a);
    return gej_to_py(&r);
}

static PyObject *meth_neg(PyObject *self, PyObject *args) {
    PyObject *p;
    if (!PyArg_ParseTuple(args, "O", &p)) return NULL;
    if (p == Py_None) Py_RETURN_NONE;
    gej a;
    if (py_to_gej(&a, p)) {
        PyErr_SetString(PyExc_ValueError, "bad point");
        return NULL;
    }
    fe_neg(a.Y, a.Y);
    return gej_to_py(&a);
}

static PyObject *meth_is_inf(PyObject *self, PyObject *args) {
    PyObject *p;
    if (!PyArg_ParseTuple(args, "O", &p)) return NULL;
    if (p == Py_None) Py_RETURN_TRUE;
    gej a;
    if (py_to_gej(&a, p)) {
        PyErr_SetString(PyExc_ValueError, "bad point");
        return NULL;
    }
    /* affine point is never infinity */
    Py_RETURN_FALSE;
}

/* msm(scalars: [int], points: [(x,y)] ) -> point or None */
static PyObject *meth_msm(PyObject *self, PyObject *args) {
    PyObject *ss, *ps;
    if (!PyArg_ParseTuple(args, "OO", &ss, &ps)) return NULL;
    Py_ssize_t n = PySequence_Size(ss);
    if (n != PySequence_Size(ps) || n < 0) {
        PyErr_SetString(PyExc_ValueError, "length mismatch");
        return NULL;
    }
    if (n == 0) Py_RETURN_NONE;
    if (n > 100000) {
        PyErr_SetString(PyExc_ValueError, "too many terms");
        return NULL;
    }
    uint32_t (*scalars)[8] = PyMem_Malloc(n * sizeof(*scalars));
    gej *points = PyMem_Malloc(n * sizeof(*points));
    if (!scalars || !points) {
        PyMem_Free(scalars); PyMem_Free(points);
        return PyErr_NoMemory();
    }
    /* N = curve order; reduce scalars mod N (matches Python msm) */
    static const uint32_t N_LIMBS[8] = {
        0xD0364141u, 0xBFD25E8Cu, 0xAF48A03Bu, 0xBAAEDCE6u,
        0xFFFFFFFEu, 0xFFFFFFFFu, 0xFFFFFFFFu, 0xFFFFFFFFu
    };
    int fail = 0;
    PyObject *Npy = PyLong_FromString(
        "115792089237316195423570985008687907852837564279074904382605163141518161494337", NULL, 10);
    for (Py_ssize_t i = 0; i < n && !fail; i++) {
        PyObject *so = PySequence_GetItem(ss, i);
        PyObject *po = PySequence_GetItem(ps, i);
        if (!so || !po) { fail = 1; }
        else {
            /* scalar -> 8 limbs mod N. Use Python for the mod (one-time). */
            PyObject *sm = PyNumber_Remainder(so, Npy);
            if (!sm || py_to_fe(scalars[i], sm)) fail = 1;
            Py_XDECREF(sm);
            if (po == Py_None) {
                memset(&points[i], 0, sizeof(points[i]));
            } else if (py_to_gej(&points[i], po)) {
                fail = 1;
            } else {
                /* skip zero scalars / inf points: set scalar 0 (add becomes no-op via mask) */
                uint32_t allz = 1;
                for (int k = 0; k < 8; k++) allz &= (scalars[i][k] == 0);
                if (allz) memset(&points[i], 0, sizeof(points[i]));
            }
        }
        Py_XDECREF(so); Py_XDECREF(po);
    }
    Py_DECREF(Npy);
    PyObject *ret = NULL;
    if (!fail) {
        gej r;
        gej_msm(&r, (const uint32_t(*)[8])scalars, points, (int)n);
        ret = gej_to_py(&r);
    } else {
        PyErr_SetString(PyExc_ValueError, "bad msm input");
    }
    PyMem_Free(scalars); PyMem_Free(points);
    return ret;
}

static PyObject *meth_mul(PyObject *self, PyObject *args) {
    PyObject *k, *p;
    if (!PyArg_ParseTuple(args, "OO", &k, &p)) return NULL;
    PyObject *ss = PyTuple_Pack(1, k);
    PyObject *ps = PyTuple_Pack(1, p);
    PyObject *a = PyTuple_Pack(2, ss, ps);
    Py_DECREF(ss); Py_DECREF(ps);
    PyObject *r = meth_msm(self, a);
    Py_DECREF(a);
    return r;
}

static PyObject *meth_compress(PyObject *self, PyObject *args) {
    PyObject *p;
    if (!PyArg_ParseTuple(args, "O", &p)) return NULL;
    gej a;
    if (py_to_gej(&a, p) || fe_is_zero(a.Z)) {
        PyErr_SetString(PyExc_ValueError, "cannot compress infinity");
        return NULL;
    }
    /* to affine */
    PyObject *ap = gej_to_py(&a);
    PyObject *x = PyTuple_GetItem(ap, 0), *y = PyTuple_GetItem(ap, 1);
    unsigned long ybit = PyLong_AsUnsignedLong(PyNumber_And(y, PyLong_FromLong(1)));
    fe xf;
    if (py_to_fe(xf, x)) { Py_DECREF(ap); PyErr_SetString(PyExc_ValueError, "x convert"); return NULL; }
    char out[33];
    out[0] = (char)(0x02 | ybit);
    for (int i = 0; i < 8; i++) {
        /* x is < p < 2^256; big-endian */
        uint32_t w = xf[7 - i];
        out[1 + i*4 + 0] = (char)(w >> 24);
        out[1 + i*4 + 1] = (char)(w >> 16);
        out[1 + i*4 + 2] = (char)(w >> 8);
        out[1 + i*4 + 3] = (char)(w);
    }
    Py_DECREF(ap);
    return PyBytes_FromStringAndSize(out, 33);
}

/* modular sqrt for p = 3 mod 4: y = y2^((p+1)/4), constant-time */
static void fe_sqrt(fe r, const fe a) {
    static const uint32_t E[8] = { /* (p+1)/4 */
        0xBFFFFF0Cu, 0xFFFFFFFFu, 0xFFFFFFFFu, 0xFFFFFFFFu,
        0xFFFFFFFFu, 0xFFFFFFFFu, 0xFFFFFFFFu, 0x3FFFFFFFu
    };
    fe result, tmp, r2;
    fe_set_int(result, 1);
    for (int b = 255; b >= 0; b--) {
        fe_sqr(tmp, result);
        fe_mul(r2, tmp, a);
        uint32_t bit = (E[b >> 5] >> (b & 31)) & 1u;
        fe_cmov(tmp, r2, bit);
        fe_copy(result, tmp);
    }
    fe_copy(r, result);
}

static PyObject *meth_decompress(PyObject *self, PyObject *args) {
    const char *data;
    Py_ssize_t len;
    if (!PyArg_ParseTuple(args, "y#", &data, &len)) return NULL;
    if (len != 33 || (data[0] != 0x02 && data[0] != 0x03)) {
        PyErr_SetString(PyExc_ValueError, "bad compressed point encoding");
        return NULL;
    }
    fe x, y2, y, t;
    for (int i = 0; i < 8; i++) {
        /* data[1..32] is x big-endian; x[0] is the least-significant limb */
        uint32_t v = 0;
        for (int j = 0; j < 4; j++)
            v |= ((uint32_t)(unsigned char)data[1 + (7-i)*4 + j]) << (8*(3-j));
        x[i] = v;
    }
    /* check x < p */
    {
        uint32_t d[8]; uint64_t br = 0;
        for (int i = 0; i < 8; i++) {
            uint64_t v = (uint64_t)x[i] - FE_P[i] - br;
            d[i] = (uint32_t)v; br = (v >> 32) & 1u;
        }
        if (br == 0) { PyErr_SetString(PyExc_ValueError, "x out of range"); return NULL; }
    }
    fe_sqr(y2, x);
    fe_mul(y2, y2, x);
    fe_set_int(t, 7);
    fe_add(y2, y2, t);
    fe_sqrt(y, y2);
    /* verify y^2 == y2 */
    fe_sqr(t, y);
    if (!fe_equal(t, y2)) {
        PyErr_SetString(PyExc_ValueError, "not on curve");
        return NULL;
    }
    if ((y[0] & 1u) != ((unsigned)data[0] & 1u)) fe_neg(y, y);
    PyObject *px = fe_to_py(x), *py = fe_to_py(y);
    PyObject *tup = PyTuple_Pack(2, px, py);
    Py_DECREF(px); Py_DECREF(py);
    return tup;
}

static PyMethodDef methods[] = {
    {"add", meth_add, METH_VARARGS, "add(p1, p2)"},
    {"dbl", meth_dbl, METH_VARARGS, "dbl(p)"},
    {"neg", meth_neg, METH_VARARGS, "neg(p)"},
    {"is_inf", meth_is_inf, METH_VARARGS, "is_inf(p)"},
    {"mul", meth_mul, METH_VARARGS, "mul(k, p)"},
    {"msm", meth_msm, METH_VARARGS, "msm(scalars, points)"},
    {"compress", meth_compress, METH_VARARGS, "compress(p)"},
    {"decompress", meth_decompress, METH_VARARGS, "decompress(bytes)"},
    {NULL, NULL, 0, NULL}
};

static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT, "curve_ext", NULL, -1, methods
};

PyMODINIT_FUNC PyInit_curve_ext(void) {
    PyObject *m = PyModule_Create(&moduledef);
    if (!m) return NULL;
    PyObject *p = PyLong_FromString(
        "115792089237316195423570985008687907852837564279074904382605163141518161494337", NULL, 10);
    PyModule_AddObject(m, "N", p);
    /* G as tuple */
    PyObject *gx = PyLong_FromString(
        "55066263022277343669578718895168534326250603453777594175500187360389116729240", NULL, 10);
    PyObject *gy = PyLong_FromString(
        "32670510020758816978083085130507043184471273380659243275938904335757337482424", NULL, 10);
    PyObject *g = PyTuple_Pack(2, gx, gy);
    Py_DECREF(gx); Py_DECREF(gy);
    PyModule_AddObject(m, "G", g);
    PyModule_AddObject(m, "INF", Py_None);
    return m;
}
