/*
 * fi_u128.c — two-limb Montgomery arithmetic and Pollard–Brent rho for
 * 2^64 <= n < 2^126.
 *
 * GMP's mpz_mul + mpz_mod costs two calls, a normalisation and a quotient
 * per step on a two-limb modulus; a hand-written word-by-word REDC with
 * R = 2^128 is roughly an order of magnitude cheaper, and this range is
 * exactly where Automatic's cheap rho tier and the PollardRho method spend
 * their time.  n < 2^126 keeps every intermediate (t < n^2 + n R < 2^255,
 * result < 2n < 2^127) inside the four product limbs plus one carry.
 */
#include "fi_internal.h"

typedef struct { uint64_t lo, hi; } U2;

/* hi:lo = a*b + c + d (cannot overflow 128 bits). */
static inline void mul_add2(uint64_t a, uint64_t b, uint64_t c, uint64_t d,
                            uint64_t* hi, uint64_t* lo) {
    uint64_t h, l;
    fi_mul64(a, b, &h, &l);
    l += c; h += (l < c);
    l += d; h += (l < d);
    *hi = h; *lo = l;
}

static inline int geq(U2 a, U2 b) { return a.hi > b.hi || (a.hi == b.hi && a.lo >= b.lo); }
static inline U2 sub2(U2 a, U2 b) {
    U2 r;
    r.lo = a.lo - b.lo;
    r.hi = a.hi - b.hi - (a.lo < b.lo);
    return r;
}
static inline U2 addmod2(U2 a, U2 b, U2 n) {      /* a, b < n < 2^126 */
    U2 r;
    r.lo = a.lo + b.lo;
    r.hi = a.hi + b.hi + (r.lo < a.lo);
    if (geq(r, n)) r = sub2(r, n);
    return r;
}
static inline U2 absdiff2(U2 a, U2 b) { return geq(a, b) ? sub2(a, b) : sub2(b, a); }

typedef struct { U2 n; uint64_t ninv; U2 one, r2; } Mont2;   /* ninv = -n^-1 mod 2^64 */

static inline U2 mulredc2(U2 a, U2 b, const Mont2* M) {
    uint64_t t0, t1, t2, t3, t4 = 0, c, h;
    /* 256-bit product. */
    mul_add2(a.lo, b.lo, 0, 0, &c, &t0);
    mul_add2(a.lo, b.hi, c, 0, &t2, &t1);
    mul_add2(a.hi, b.lo, t1, 0, &c, &t1);
    mul_add2(a.hi, b.hi, t2, c, &t3, &t2);
    /* REDC, limb 0. */
    uint64_t m = t0 * M->ninv;
    mul_add2(m, M->n.lo, t0, 0, &c, &h);            /* low word becomes 0 */
    mul_add2(m, M->n.hi, t1, c, &c, &t1);
    t2 += c; c = (t2 < c); t3 += c; t4 += (t3 < c);
    /* REDC, limb 1. */
    m = t1 * M->ninv;
    mul_add2(m, M->n.lo, t1, 0, &c, &h);
    mul_add2(m, M->n.hi, t2, c, &c, &t2);
    t3 += c; t4 += (t3 < c);
    U2 r = { t2, t3 };
    if (t4 || geq(r, M->n)) r = sub2(r, M->n);
    (void)h;
    return r;
}

static void mont2_setup(Mont2* M, const mpz_t n) {
    M->n.lo = mpz_getlimbn(n, 0);
    M->n.hi = mpz_getlimbn(n, 1);
    uint64_t x = M->n.lo;
    for (int i = 0; i < 5; i++) x *= 2 - M->n.lo * x;
    M->ninv = 0 - x;
    /* one = 2^128 mod n, r2 = 2^256 mod n, via GMP once. */
    mpz_t t;
    mpz_init(t);
    mpz_setbit(t, 128);
    mpz_mod(t, t, n);
    M->one.lo = mpz_getlimbn(t, 0);
    M->one.hi = mpz_size(t) > 1 ? mpz_getlimbn(t, 1) : 0;
    mpz_set_ui(t, 0);
    mpz_setbit(t, 256);
    mpz_mod(t, t, n);
    M->r2.lo = mpz_getlimbn(t, 0);
    M->r2.hi = mpz_size(t) > 1 ? mpz_getlimbn(t, 1) : 0;
    mpz_clear(t);
}

static inline U2 to_mont2(uint64_t a, const Mont2* M) {
    U2 x = { a, 0 };
    return mulredc2(x, M->r2, M);
}

static void u2_to_mpz(mpz_t z, U2 a) {
    mth_mpz_set_u64(z, (uint64_t)a.hi);
    mpz_mul_2exp(z, z, 64);
    mth_mpz_add_u64(z, z, (uint64_t)a.lo);
}

int fi_u128_rho(mpz_t d, const mpz_t n, uint64_t c, uint64_t max_iters) {
    size_t bits = mpz_sizeinbase(n, 2);
    if (bits <= 64 || bits > 126 || mpz_even_p(n)) return 0;
    Mont2 M;
    mont2_setup(&M, n);
    U2 cm = to_mont2(c ? c : 1, &M);
    U2 y = to_mont2(2, &M), x = y, ys = y, q = M.one;
    mpz_t g, qz;
    mpz_init_set_ui(g, 1);
    mpz_init(qz);
    const uint64_t m = 128;
    uint64_t done = 0;
    for (uint64_t r = 1; mpz_cmp_ui(g, 1) == 0 && done < max_iters; r <<= 1) {
        x = y;
        for (uint64_t i = 0; i < r; i++) y = addmod2(mulredc2(y, y, &M), cm, M.n);
        done += r;
        for (uint64_t k = 0; k < r && mpz_cmp_ui(g, 1) == 0; k += m) {
            ys = y;
            uint64_t lim = (r - k < m) ? r - k : m;
            for (uint64_t i = 0; i < lim; i++) {
                y = addmod2(mulredc2(y, y, &M), cm, M.n);
                q = mulredc2(q, absdiff2(x, y), &M);
            }
            done += lim;
            u2_to_mpz(qz, q);
            mpz_gcd(g, qz, n);
        }
    }
    if (mpz_cmp(g, n) == 0 || mpz_sgn(g) == 0) {    /* backtrack the last batch */
        mpz_set_ui(g, 1);
        for (uint64_t i = 0; i < 4 * m && mpz_cmp_ui(g, 1) == 0; i++) {
            ys = addmod2(mulredc2(ys, ys, &M), cm, M.n);
            u2_to_mpz(qz, absdiff2(x, ys));
            mpz_gcd(g, qz, n);
        }
    }
    int ok = mpz_cmp_ui(g, 1) > 0 && mpz_cmp(g, n) < 0;
    if (ok) mpz_set(d, g);
    mpz_clears(g, qz, NULL);
    return ok;
}

/* ------------------------------------------------------------------ */
/* Two-limb micro-ECM                                                  */
/* ------------------------------------------------------------------ */

/* The word engine's micro-ECM (fi_u64.c) on two limbs: Montgomery curves
 * with Suyama's parametrisation, a prime-power ladder for stage 1 and a
 * baby-step/giant-step stage 2 over q = m D +- j, D = 210.  For a 64..126
 * bit composite this finds a 10..20 digit factor in tens of microseconds per
 * curve, the gap between cheap rho and the quadratic sieve. */

static inline U2 submod2(U2 a, U2 b, U2 n) {
    if (geq(a, b)) return sub2(a, b);
    U2 r = sub2(n, sub2(b, a));
    return r;
}

typedef struct { U2 X, Z; } P2;

static inline P2 e2_dbl(P2 P, U2 a24, const Mont2* M) {
    U2 s = addmod2(P.X, P.Z, M->n), d = submod2(P.X, P.Z, M->n);
    U2 s2 = mulredc2(s, s, M), d2 = mulredc2(d, d, M);
    U2 t = submod2(s2, d2, M->n);
    P2 R;
    R.X = mulredc2(s2, d2, M);
    R.Z = mulredc2(t, addmod2(d2, mulredc2(a24, t, M), M->n), M);
    return R;
}

static inline P2 e2_add(P2 P, P2 Q, P2 D, const Mont2* M) {
    U2 u = mulredc2(submod2(P.X, P.Z, M->n), addmod2(Q.X, Q.Z, M->n), M);
    U2 v = mulredc2(addmod2(P.X, P.Z, M->n), submod2(Q.X, Q.Z, M->n), M);
    U2 sp = addmod2(u, v, M->n), sm = submod2(u, v, M->n);
    P2 R;
    R.X = mulredc2(D.Z, mulredc2(sp, sp, M), M);
    R.Z = mulredc2(D.X, mulredc2(sm, sm, M), M);
    return R;
}

static P2 e2_mul(P2 P, uint64_t k, U2 a24, const Mont2* M) {
    if (k <= 1) return P;
    P2 R0 = P, R1 = e2_dbl(P, a24, M);
    int top = 63;
    while (!((k >> top) & 1)) top--;
    for (int b = top - 1; b >= 0; b--) {
        if ((k >> b) & 1) { R0 = e2_add(R1, R0, P, M); R1 = e2_dbl(R1, a24, M); }
        else              { R1 = e2_add(R1, R0, P, M); R0 = e2_dbl(R0, a24, M); }
    }
    return R0;
}

static void mpz_from_u2(mpz_t z, U2 a) { u2_to_mpz(z, a); }
static U2 u2_from_mpz(const mpz_t z) {
    U2 r = { mpz_getlimbn(z, 0), mpz_size(z) > 1 ? mpz_getlimbn(z, 1) : 0 };
    return r;
}

#define T2_D 210
static const uint8_t T2_J[] = {
    1, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73,
    79, 83, 89, 97, 101, 103
};
#define T2_NJ (sizeof T2_J / sizeof T2_J[0])

/* One curve; returns 1 with g set to gcd(.., n) when it is not 1. */
static int t2_curve(mpz_t g, const mpz_t n, uint64_t sigma, uint32_t B1, uint32_t B2,
                    const Mont2* M, const unsigned char* isp, mpz_t t) {
    /* Suyama: u = s^2 - 5, v = 4 s, x0 = u^3, z0 = v^3,
     * (A+2)/4 = (v-u)^3 (3u+v) / (16 u^3 v); the inverse is taken once in
     * GMP (a failed inverse is itself a factor). */
    mpz_t u, v, x, z, num, den;
    mpz_inits(u, v, x, z, num, den, NULL);
    mth_mpz_set_u64(u, (uint64_t)sigma); mpz_mul(u, u, u); mth_mpz_sub_u64(u, u, 5); mpz_mod(u, u, n);
    mth_mpz_set_u64(v, (uint64_t)sigma); mth_mpz_mul_u64(v, v, 4); mpz_mod(v, v, n);
    mpz_powm_ui(x, u, 3, n);
    mpz_powm_ui(z, v, 3, n);
    mpz_sub(num, v, u); mpz_powm_ui(num, num, 3, n);
    mpz_mul_ui(t, u, 3); mpz_add(t, t, v); mpz_mul(num, num, t); mpz_mod(num, num, n);
    mpz_mul_ui(den, x, 16); mpz_mul(den, den, v); mpz_mod(den, den, n);
    int hit = 0;
    if (!mpz_invert(t, den, n)) {
        mpz_gcd(g, den, n);
        hit = mpz_cmp_ui(g, 1) != 0;
        mpz_clears(u, v, x, z, num, den, NULL);
        return hit;
    }
    mpz_mul(num, num, t); mpz_mod(num, num, n);
    /* Into Montgomery form: a * R mod n via a * r2 * R^-1. */
    U2 a24 = mulredc2(u2_from_mpz(num), M->r2, M);
    P2 P;
    P.X = mulredc2(u2_from_mpz(x), M->r2, M);
    P.Z = mulredc2(u2_from_mpz(z), M->r2, M);
    mpz_clears(u, v, x, z, num, den, NULL);

    uint64_t acc = 1;
    for (uint32_t p = 2; p <= B1; p++) {
        if (!isp[p]) continue;
        uint64_t q = p;
        while (q * p <= B1) q *= p;
        uint64_t hi, lo;
        fi_mul64(acc, q, &hi, &lo);
        if (hi) { P = e2_mul(P, acc, a24, M); acc = q; }
        else acc = lo;
    }
    P = e2_mul(P, acc, a24, M);
    mpz_from_u2(t, P.Z);
    mpz_gcd(g, t, n);
    if (mpz_cmp_ui(g, 1) != 0) return 1;

    P2 baby[T2_NJ];
    P2 Pd = e2_dbl(P, a24, M);
    {
        P2 prev = P, cur = e2_add(Pd, P, P, M);
        size_t bi = 0;
        baby[bi++] = P;
        for (uint32_t j = 3; j < T2_D / 2 && bi < T2_NJ; j += 2) {
            if (j == T2_J[bi]) baby[bi++] = cur;
            P2 nx = e2_add(cur, Pd, prev, M);
            prev = cur; cur = nx;
        }
    }
    P2 G = e2_mul(P, T2_D, a24, M);
    uint32_t m0 = (B1 + T2_D / 2) / T2_D;
    if (m0 < 1) m0 = 1;
    P2 Mc = e2_mul(G, m0, a24, M);
    P2 Mp = m0 > 1 ? e2_mul(G, m0 - 1, a24, M) : G;
    U2 prod = M->one;
    for (uint32_t m = m0; (uint64_t)m * T2_D - T2_D / 2 <= B2; m++) {
        for (size_t k = 0; k < T2_NJ; k++) {
            uint32_t j = T2_J[k], q1 = m * T2_D - j, q2 = m * T2_D + j;
            int h = (q1 > B1 && q1 <= B2 && isp[q1]) || (q2 > B1 && q2 <= B2 && isp[q2]);
            if (h) prod = mulredc2(prod, submod2(mulredc2(Mc.X, baby[k].Z, M), mulredc2(baby[k].X, Mc.Z, M), M->n), M);
        }
        P2 Mn = (m == 1) ? e2_dbl(Mc, a24, M) : e2_add(Mc, G, Mp, M);
        Mp = Mc; Mc = Mn;
    }
    mpz_from_u2(t, prod);
    mpz_gcd(g, t, n);
    return mpz_cmp_ui(g, 1) != 0;
}

int fi_u128_ecm(mpz_t d, const mpz_t n, uint64_t* rng, uint32_t B1, int max_curves) {
    size_t bits = mpz_sizeinbase(n, 2);
    if (bits <= 64 || bits > 126 || mpz_even_p(n)) return 0;
    uint32_t lim;
    const unsigned char* isp = fi_prime_bitmap(&lim);
    uint32_t B2 = 25 * B1;
    if (B2 > lim - T2_D) B2 = lim - T2_D;
    if (B1 > B2) B1 = B2;
    Mont2 M;
    mont2_setup(&M, n);
    mpz_t g, t;
    mpz_inits(g, t, NULL);
    int ok = 0;
    for (int c = 0; c < max_curves && !ok; c++) {
        uint64_t x = *rng;
        x ^= x << 13; x ^= x >> 7; x ^= x << 17;
        *rng = x;
        uint64_t sigma = 6 + (x % (UINT64_C(1) << 30));
        if (t2_curve(g, n, sigma, B1, B2, &M, isp, t) && mpz_cmp(g, n) < 0) { mpz_set(d, g); ok = 1; }
    }
    mpz_clears(g, t, NULL);
    return ok;
}
