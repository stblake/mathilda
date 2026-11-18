/* dirichlet_l.c -- DirichletL[k, j, s], the Dirichlet L-function L(chi, s) of
 * the character chi = DirichletCharacter[k, j, .].
 *
 * Attributes: Listable | NHoldAll | Protected. Not NumericFunction: NumericQ
 * is False and 1. + DirichletL[4, 2, 3] stays unevaluated, while an explicit N
 * keeps k and j exact and numericalizes s (numericalize_rec in src/numeric.c).
 *
 * Evaluation order:
 *   1. arity (argrx); an NDArray / packed s -> the buffer loop below.
 *   2. (k, j) validation via the character engine (intp / invl messages).
 *   3. k = 1: L = Zeta[s] for numeric s (and Interval s).
 *   4. Interval s, real character: a rigorous arb enclosure.
 *   5. machine s: the raw-double kernel sf_dirichlet_l_tab_* -- the same one
 *      the NDArray loop and Compile[] use (src/ndkernels.c ndk_DirichletL_n).
 *   6. exact s: the closed forms Mathematica gives -- s a non-positive integer
 *      (generalized Bernoulli numbers, all characters) and s = 1 for a
 *      non-principal character -- else unevaluated.
 *   7. inexact s beyond machine (or a machine decline): FLINT acb_dirichlet_l on
 *      the primitive character plus Euler factors, at the precision of s; the
 *      fallback build sums HurwitzZeta templates.
 *
 * Exact values live in the cyclotomic field Q(zeta_O), O = order of chi; they
 * are computed as coefficient vectors in Q[x]/(x^O - 1), reduced modulo the
 * cyclotomic polynomial (so the representation is canonical), and rendered as
 * sums of rational multiples of roots of unity in the DirichletCharacter normal
 * form (unit * E^(I Pi r)).
 */
#include "dirichlet_l.h"
#include "dirichlet.h"
#include "sf_machine_dirichlet.h"
#include "sf_machine_zeta.h"
#include "sf_scalar.h"
#include "symtab.h"
#include "attr.h"
#include "eval.h"
#include "common.h"
#include "message.h"
#include "sym_names.h"
#include "ndarray.h"
#include "arithmetic.h"
#include "numeric.h"
#include "flint_num_bridge.h"
#include <gmp.h>
#ifdef USE_MPFR
#include <mpfr.h>
#endif
#include <math.h>
#include <stdlib.h>
#include <string.h>

/* Resource bounds for the exact closed forms (cost and expression size grow
 * with the conductor); beyond them the call stays unevaluated, never wrong. */
#define DL_EXACT_FMAX   ((uint64_t)1 << 20)   /* s <= 0: O(f m) big-integer work */
#define DL_S1_NMAX      4096                  /* s = 1: lcm(order, f) cyclotomic */
#define DL_REDUCE_NMAX  2048                  /* reduce modulo Phi_n up to here  */

static bool is_inexact_number(const Expr* e);

static Expr* mk1(const char* h, Expr* a) {
    return expr_new_function(expr_new_symbol(h), (Expr*[]){ a }, 1);
}
static Expr* mk2(const char* h, Expr* a, Expr* b) {
    return expr_new_function(expr_new_symbol(h), (Expr*[]){ a, b }, 2);
}

/* ------------------------------------------------------------------ */
/*  Cyclotomic numbers: Q[x]/(x^n - 1)                                 */
/* ------------------------------------------------------------------ */

typedef struct { uint64_t n; mpq_t* c; } Cyc;

static void cyc_init(Cyc* z, uint64_t n) {
    z->n = n;
    z->c = (mpq_t*)malloc(n * sizeof(mpq_t));
    for (uint64_t i = 0; i < n; i++) mpq_init(z->c[i]);
}
static void cyc_clear(Cyc* z) {
    for (uint64_t i = 0; i < z->n; i++) mpq_clear(z->c[i]);
    free(z->c);
}
static void cyc_mul(Cyc* out, const Cyc* a, const Cyc* b) {   /* out != a, b */
    uint64_t n = a->n;
    mpq_t t; mpq_init(t);
    for (uint64_t i = 0; i < n; i++) mpq_set_ui(out->c[i], 0, 1);
    for (uint64_t i = 0; i < n; i++) {
        if (mpq_sgn(a->c[i]) == 0) continue;
        for (uint64_t j = 0; j < n; j++) {
            if (mpq_sgn(b->c[j]) == 0) continue;
            mpq_mul(t, a->c[i], b->c[j]);
            mpq_add(out->c[(i + j) % n], out->c[(i + j) % n], t);
        }
    }
    mpq_clear(t);
}
static void cyc_mul_inplace(Cyc* a, const Cyc* b) {
    Cyc t; cyc_init(&t, a->n);
    cyc_mul(&t, a, b);
    for (uint64_t i = 0; i < a->n; i++) mpq_swap(a->c[i], t.c[i]);
    cyc_clear(&t);
}
static void cyc_scale(Cyc* a, const mpq_t q) {
    for (uint64_t i = 0; i < a->n; i++) mpq_mul(a->c[i], a->c[i], q);
}

/* Integer polynomial helpers for Phi_n (coefficients low -> high). */
static void poly_mul_xd_minus1(mpz_t* p, uint64_t* deg, uint64_t d) {
    /* p * (x^d - 1) */
    for (uint64_t i = *deg + d + 1; i-- > 0;) {
        if (i >= d && i - d <= *deg) {
            if (i <= *deg) mpz_sub(p[i], p[i - d], p[i]);
            else mpz_set(p[i], p[i - d]);
        } else if (i <= *deg) {
            mpz_neg(p[i], p[i]);
        }
    }
    *deg += d;
}
static void poly_div_xd_minus1(mpz_t* p, uint64_t* deg, uint64_t d) {
    /* p / (x^d - 1), exact: q_i = q_{i+d} - p_{i+d} going downward */
    uint64_t qd = *deg - d;
    mpz_t* q = (mpz_t*)malloc((qd + 1) * sizeof(mpz_t));
    for (uint64_t i = 0; i <= qd; i++) mpz_init(q[i]);
    for (uint64_t i = qd + 1; i-- > 0;) {
        /* coefficient of x^(i+d) in p = q_i - q_{i+d} */
        mpz_set(q[i], p[i + d]);
        if (i + d <= qd) mpz_add(q[i], q[i], q[i + d]);
    }
    for (uint64_t i = 0; i <= *deg; i++) mpz_set_ui(p[i], 0);
    for (uint64_t i = 0; i <= qd; i++) { mpz_set(p[i], q[i]); mpz_clear(q[i]); }
    free(q);
    *deg = qd;
}
static int moebius(uint64_t n) {
    int mu = 1;
    for (uint64_t p = 2; p * p <= n; p++) {
        if (n % p) continue;
        n /= p;
        if (n % p == 0) return 0;
        mu = -mu;
    }
    if (n > 1) mu = -mu;
    return mu;
}

/* Reduce z modulo Phi_n (n = z->n), leaving a canonical vector supported on
 * 0..phi(n)-1. Skipped (left unreduced, still correct) for n > DL_REDUCE_NMAX. */
static void cyc_reduce(Cyc* z) {
    uint64_t n = z->n;
    if (n <= 1 || n > DL_REDUCE_NMAX) return;
    /* Phi_n = prod_{d|n} (x^d - 1)^mu(n/d): multiply the mu = +1 factors, then
     * divide out the mu = -1 ones. Degree stays <= sum of divisors. */
    uint64_t cap = 1;
    for (uint64_t d = 1; d <= n; d++) if (n % d == 0) cap += d;
    mpz_t* p = (mpz_t*)malloc((cap + 1) * sizeof(mpz_t));
    for (uint64_t i = 0; i <= cap; i++) mpz_init(p[i]);
    mpz_set_ui(p[0], 1);
    uint64_t deg = 0;
    for (uint64_t d = 1; d <= n; d++)
        if (n % d == 0 && moebius(n / d) == 1) poly_mul_xd_minus1(p, &deg, d);
    for (uint64_t d = 1; d <= n; d++)
        if (n % d == 0 && moebius(n / d) == -1) poly_div_xd_minus1(p, &deg, d);
    /* Phi_n is monic of degree deg = phi(n): reduce from the top. */
    mpq_t t; mpq_init(t);
    for (uint64_t i = n; i-- > deg;) {
        if (mpq_sgn(z->c[i]) == 0) continue;
        for (uint64_t k = 0; k < deg; k++) {
            if (mpz_sgn(p[k]) == 0) continue;
            mpq_set_z(t, p[k]);
            mpq_mul(t, t, z->c[i]);
            mpq_sub(z->c[i - deg + k], z->c[i - deg + k], t);
        }
        mpq_set_ui(z->c[i], 0, 1);
    }
    mpq_clear(t);
    for (uint64_t i = 0; i <= cap; i++) mpz_clear(p[i]);
    free(p);
}

static Expr* mpq_expr(const mpq_t q) {
    return make_rational_mpz(mpq_numref(q), mpq_denref(q));
}

/* Render sum_m c_m zeta_n^m (after reduction) as an evaluated Expr. */
static Expr* cyc_render(Cyc* z) {
    cyc_reduce(z);
    uint64_t nz = 0;
    for (uint64_t m = 0; m < z->n; m++) if (mpq_sgn(z->c[m])) nz++;
    if (nz == 0) return expr_new_integer(0);
    Expr** terms = (Expr**)malloc(nz * sizeof(Expr*));
    uint64_t t = 0;
    for (uint64_t m = 0; m < z->n; m++) {
        if (!mpq_sgn(z->c[m])) continue;
        Expr* root = dch_root_of_unity_expr(0, m, z->n);
        terms[t++] = mk2(SYM_Times, mpq_expr(z->c[m]), root);
    }
    Expr* sum = expr_new_function(expr_new_symbol(SYM_Plus), terms, (size_t)nz);
    free(terms);
    return eval_and_free(sum);
}

/* An exact Integer / BigInt / Rational Expr -> mpq. */
static bool expr_to_mpq(const Expr* e, mpq_t out) {
    if (e->type == EXPR_INTEGER) { mth_mpq_set_i64(out, e->data.integer, 1); return true; }
    if (e->type == EXPR_BIGINT) { mpq_set_z(out, e->data.bigint); return true; }
    int64_t n, d;
    if (is_rational((Expr*)e, &n, &d)) { mth_mpq_set_i64(out, n, (uint64_t)d); mpq_canonicalize(out); return true; }
    if (e->type == EXPR_FUNCTION && e->data.function.arg_count == 2
        && e->data.function.head->type == EXPR_SYMBOL
        && e->data.function.head->data.symbol.name == SYM_Rational) {
        const Expr* a = e->data.function.args[0];
        const Expr* b = e->data.function.args[1];
        if (a->type == EXPR_BIGINT || b->type == EXPR_BIGINT) {
            mpz_t x, y; mpz_init(x); mpz_init(y);
            expr_to_mpz(a, x); expr_to_mpz(b, y);
            mpq_set_num(out, x); mpq_set_den(out, y); mpq_canonicalize(out);
            mpz_clear(x); mpz_clear(y);
            return true;
        }
    }
    return false;
}

/* ------------------------------------------------------------------ */
/*  Exact values                                                       */
/* ------------------------------------------------------------------ */

/* chi*(a) exponent over the order (a coprime to f); uses the dense table. */
static int64_t star_e(const DchTable* T, uint64_t a) {
    return T->estar[a % T->f];
}

/* L(chi, 1 - m), m >= 1: -B_{m,chi*}/m * prod (1 - chi*(p) p^(m-1)), with
 * B_{m,chi} = sum_{i=0}^m C(m,i) B_i f^(i-1) sum_{a=1}^f chi(a) a^(m-i). */
static Expr* dl_exact_nonpos(const DchTable* T, int64_t m) {
    uint64_t f = T->f, O = T->order;
    if (!T->estar || f > DL_EXACT_FMAX) return NULL;
    /* Bernoulli numbers B_0..B_m (BernoulliB[1] = -1/2). */
    mpq_t* B = (mpq_t*)malloc((size_t)(m + 1) * sizeof(mpq_t));
    bool ok = true;
    for (int64_t i = 0; i <= m; i++) {
        mpq_init(B[i]);
        if (!ok) continue;
        Expr* b = eval_and_free(mk1(SYM_BernoulliB, expr_new_integer(i)));
        if (!b || !expr_to_mpq(b, B[i])) ok = false;
        if (b) expr_free(b);
    }
    Cyc acc; cyc_init(&acc, O);
    if (ok) {
        /* Power sums P_t = sum_a chi(a) a^t, t = 0..m, as cyclotomic vectors. */
        Cyc* P = (Cyc*)malloc((size_t)(m + 1) * sizeof(Cyc));
        for (int64_t t = 0; t <= m; t++) cyc_init(&P[t], O);
        mpz_t pw; mpz_init(pw);
        for (uint64_t a = 1; a <= f; a++) {
            int64_t e = star_e(T, a);
            if (e < 0) continue;
            mpz_set_ui(pw, 1);
            for (int64_t t = 0; t <= m; t++) {
                /* integer entries: the denominator stays 1 */
                mpz_add(mpq_numref(P[t].c[e]), mpq_numref(P[t].c[e]), pw);
                if (t < m) mpz_mul_ui(pw, pw, (unsigned long)a);
            }
        }
        mpz_clear(pw);
        /* acc = sum_i C(m,i) B_i f^(i-1) P_{m-i} */
        mpq_t coef, fp; mpq_init(coef); mpq_init(fp);
        mpz_t bin, fz; mpz_init(bin); mpz_init(fz);
        for (int64_t i = 0; i <= m; i++) {
            if (mpq_sgn(B[i]) == 0) continue;
            mpz_bin_uiui(bin, (unsigned long)m, (unsigned long)i);
            mpq_set_z(coef, bin);
            mpq_mul(coef, coef, B[i]);
            /* f^(i-1) */
            mpz_ui_pow_ui(fz, (unsigned long)f, (unsigned long)(i > 0 ? i - 1 : 1));
            mpq_set_z(fp, fz);
            if (i == 0) mpq_inv(fp, fp);
            mpq_mul(coef, coef, fp);
            const Cyc* src = &P[m - i];
            for (uint64_t r = 0; r < O; r++) {
                if (mpq_sgn(src->c[r]) == 0) continue;
                mpq_t t2; mpq_init(t2);
                mpq_mul(t2, coef, src->c[r]);
                mpq_add(acc.c[r], acc.c[r], t2);
                mpq_clear(t2);
            }
        }
        mpq_clear(coef); mpq_clear(fp); mpz_clear(bin); mpz_clear(fz);
        for (int64_t t = 0; t <= m; t++) cyc_clear(&P[t]);
        free(P);
        /* * (-1/m) */
        mpq_t s; mpq_init(s); mpq_set_si(s, -1, (unsigned long)m); mpq_canonicalize(s);
        cyc_scale(&acc, s);
        mpq_clear(s);
        /* Euler factors (1 - chi*(p) p^(m-1)) */
        for (int i = 0; i < T->neul; i++) {
            Cyc F; cyc_init(&F, O);
            mpq_set_ui(F.c[0], 1, 1);
            mpz_t pp; mpz_init(pp);
            mpz_ui_pow_ui(pp, (unsigned long)T->eul_p[i], (unsigned long)(m - 1));
            mpq_t q; mpq_init(q); mpq_set_z(q, pp);
            mpq_sub(F.c[T->eul_e[i] % O], F.c[T->eul_e[i] % O], q);
            mpq_clear(q); mpz_clear(pp);
            cyc_mul_inplace(&acc, &F);
            cyc_clear(&F);
        }
    }
    for (int64_t i = 0; i <= m; i++) mpq_clear(B[i]);
    free(B);
    Expr* out = ok ? cyc_render(&acc) : NULL;
    cyc_clear(&acc);
    return out;
}

/* prod (1 - chi*(p)/p) as a cyclotomic vector in Q(zeta_N), N a multiple of
 * the order. */
static void dl_euler_at1(const DchTable* T, Cyc* acc) {
    uint64_t N = acc->n, O = T->order;
    for (int i = 0; i < T->neul; i++) {
        Cyc F; cyc_init(&F, N);
        mpq_set_ui(F.c[0], 1, 1);
        mpq_t q; mpq_init(q); mpq_set_ui(q, 1, (unsigned long)T->eul_p[i]);
        uint64_t idx = (T->eul_e[i] % O) * (N / O) % N;
        mpq_sub(F.c[idx], F.c[idx], q);
        mpq_clear(q);
        cyc_mul_inplace(acc, &F);
        cyc_clear(&F);
    }
}

/* L(chi, 1) for a non-principal chi: the classical closed forms for the
 * primitive chi* mod f times the Euler factors prod (1 - chi*(p)/p).
 *   odd chi*:   L = Pi I tau(chi*)/f^2  sum_a conj(chi*(a)) a
 *   even chi*:  L = -2 tau(chi*)/f sum_{a<f/2} conj(chi*(a)) Log[Sin[Pi a/f]]
 * For a real chi*, tau = Sqrt[f] (even) or I Sqrt[f] (odd). */
static Expr* dl_exact_one(const DchTable* T) {
    uint64_t f = T->f, O = T->order;
    if (T->principal || !T->estar) return NULL;
    uint64_t N = O / dch_gcd(O, f) * f;
    if (N > DL_S1_NMAX) return NULL;
    Expr* res = NULL;
    if (T->real) {
        /* Euler factors, rational. */
        mpq_t E; mpq_init(E); mpq_set_ui(E, 1, 1);
        for (int i = 0; i < T->neul; i++) {
            mpq_t q; mpq_init(q);
            mpq_set_si(q, (T->eul_e[i] % O == 0) ? -1 : 1, (unsigned long)T->eul_p[i]);
            mpq_canonicalize(q);
            mpq_t one; mpq_init(one); mpq_set_ui(one, 1, 1);
            mpq_add(q, one, q);                 /* 1 - chi(p)/p */
            mpq_mul(E, E, q);
            mpq_clear(q); mpq_clear(one);
        }
        Expr* sqrt_f = mk2(SYM_Power, expr_new_integer((int64_t)f), make_rational(-1, 2));
        if (T->odd) {
            int64_t S = 0;
            for (uint64_t a = 1; a < f; a++) {
                int64_t e = T->estar[a];
                if (e >= 0) S += (e == 0) ? (int64_t)a : -(int64_t)a;
            }
            res = expr_new_function(expr_new_symbol(SYM_Times),
                      (Expr*[]){ make_rational(-S, (int64_t)f), expr_new_symbol(SYM_Pi),
                                 sqrt_f, mpq_expr(E) }, 4);
        } else {
            uint64_t half = (f - 1) / 2, cnt = 0;
            Expr** terms = (Expr**)malloc((half ? half : 1) * sizeof(Expr*));
            for (uint64_t a = 1; a <= half; a++) {
                int64_t e = T->estar[a];
                if (e < 0) continue;
                Expr* lg = mk1(SYM_Log, mk1(SYM_Sin,
                               mk2(SYM_Times, make_rational((int64_t)a, (int64_t)f),
                                   expr_new_symbol(SYM_Pi))));
                terms[cnt++] = (e == 0) ? lg : mk2(SYM_Times, expr_new_integer(-1), lg);
            }
            Expr* sum = expr_new_function(expr_new_symbol(SYM_Plus), terms, (size_t)cnt);
            free(terms);
            res = expr_new_function(expr_new_symbol(SYM_Times),
                      (Expr*[]){ expr_new_integer(-2), sqrt_f, sum, mpq_expr(E) }, 4);
        }
        mpq_clear(E);
        return eval_and_free(res);
    }
    /* Complex character: tau(chi*) and the sums in Q(zeta_N). */
    Cyc tau; cyc_init(&tau, N);
    for (uint64_t a = 1; a < f; a++) {
        int64_t e = T->estar[a];
        if (e < 0) continue;
        uint64_t idx = ((uint64_t)e * (N / O) + a * (N / f)) % N;
        mpq_t one; mpq_init(one); mpq_set_ui(one, 1, 1);
        mpq_add(tau.c[idx], tau.c[idx], one);
        mpq_clear(one);
    }
    dl_euler_at1(T, &tau);
    if (T->odd) {
        Cyc S; cyc_init(&S, N);
        for (uint64_t a = 1; a < f; a++) {
            int64_t e = T->estar[a];
            if (e < 0) continue;
            uint64_t idx = ((O - (uint64_t)e) % O) * (N / O) % N;   /* conj */
            mpq_t q; mpq_init(q); mpq_set_ui(q, (unsigned long)a, 1);
            mpq_add(S.c[idx], S.c[idx], q);
            mpq_clear(q);
        }
        cyc_mul_inplace(&tau, &S);
        cyc_clear(&S);
        Expr* alg = cyc_render(&tau);
        res = expr_new_function(expr_new_symbol(SYM_Times),
                  (Expr*[]){ make_complex(expr_new_integer(0), make_rational(1, (int64_t)(f * f))),
                             expr_new_symbol(SYM_Pi), alg }, 3);
    } else {
        Expr* alg = cyc_render(&tau);
        uint64_t half = (f - 1) / 2, cnt = 0;
        Expr** terms = (Expr**)malloc((half ? half : 1) * sizeof(Expr*));
        for (uint64_t a = 1; a <= half; a++) {
            int64_t e = T->estar[a];
            if (e < 0) continue;
            Expr* c = dch_root_of_unity_expr(0, (O - (uint64_t)e) % O, O);
            Expr* lg = mk1(SYM_Log, mk1(SYM_Sin,
                           mk2(SYM_Times, make_rational((int64_t)a, (int64_t)f),
                               expr_new_symbol(SYM_Pi))));
            terms[cnt++] = mk2(SYM_Times, c, lg);
        }
        Expr* sum = expr_new_function(expr_new_symbol(SYM_Plus), terms, (size_t)cnt);
        free(terms);
        res = expr_new_function(expr_new_symbol(SYM_Times),
                  (Expr*[]){ make_rational(-2, (int64_t)f), alg, sum }, 3);
    }
    cyc_clear(&tau);
    return eval_and_free(res);
}

/* Exact s: the values Mathematica generates; NULL leaves the call alone. */
static Expr* dl_exact(const DchTable* T, const Expr* s) {
    if (s->type != EXPR_INTEGER) return NULL;
    int64_t sv = s->data.integer;
    if (sv <= 0) return dl_exact_nonpos(T, 1 - sv);
    if (sv == 1) return dl_exact_one(T);
    return NULL;
}

/* ------------------------------------------------------------------ */
/*  Arbitrary precision                                                */
/* ------------------------------------------------------------------ */

static bool dl_chi_cb(const void* ctx, uint64_t n, uint64_t* num) {
    return dch_star_exp((const DchTable*)ctx, n, num);
}

/* Fallback (no FLINT): f^-s sum_a chi*(a) HurwitzZeta[s, a/f] times the Euler
 * factors, evaluated at the precision of s. Needs f <= the table limit. */
static Expr* dl_numeric_template(const DchTable* T, const Expr* s) {
    if (!T->estar || T->f > 4096) return NULL;
    uint64_t f = T->f, O = T->order;
    Expr** terms = (Expr**)malloc(f * sizeof(Expr*));
    size_t cnt = 0;
    for (uint64_t a = 1; a <= f; a++) {
        int64_t e = T->estar[a % f];
        if (e < 0) continue;
        Expr* hz = mk2(SYM_HurwitzZeta, expr_copy((Expr*)s),
                       make_rational((int64_t)a, (int64_t)f));
        terms[cnt++] = mk2(SYM_Times, dch_root_of_unity_expr(0, (uint64_t)e, O), hz);
    }
    Expr* sum = expr_new_function(expr_new_symbol(SYM_Plus), terms, cnt);
    free(terms);
    Expr* prod = mk2(SYM_Times,
                     mk2(SYM_Power, expr_new_integer((int64_t)f),
                         mk2(SYM_Times, expr_new_integer(-1), expr_copy((Expr*)s))),
                     sum);
    for (int i = 0; i < T->neul; i++) {
        Expr* ps = mk2(SYM_Power, expr_new_integer((int64_t)T->eul_p[i]),
                       mk2(SYM_Times, expr_new_integer(-1), expr_copy((Expr*)s)));
        Expr* c = dch_root_of_unity_expr(0, T->eul_e[i], O);
        Expr* fac = mk2(SYM_Plus, expr_new_integer(1),
                        expr_new_function(expr_new_symbol(SYM_Times),
                            (Expr*[]){ expr_new_integer(-1), c, ps }, 3));
        prod = mk2(SYM_Times, prod, fac);
    }
    Expr* v = eval_and_free(prod);
    if (v && !is_inexact_number(v)) { expr_free(v); return NULL; }
    return v;
}

/* Demote an MPFR-leafed numeric result to machine precision. */
static Expr* dl_to_machine(Expr* v) {
    if (!v) return NULL;
    return eval_and_free(mk1(SYM_N, v));
}

static Expr* dl_numeric(const DchTable* T, Expr* s, bool machine) {
    Expr* v = flint_num_dirichlet_l(s, T->f, T->order, dl_chi_cb, T,
                                    T->neul, T->eul_p, T->eul_e);
    if (!v) v = dl_numeric_template(T, s);
    if (v && machine) v = dl_to_machine(v);
    return v;
}

/* ------------------------------------------------------------------ */
/*  Interval                                                           */
/* ------------------------------------------------------------------ */

/* An Interval endpoint as a double, rounded OUTWARD (down for a lower endpoint,
 * up for an upper one) so the double interval contains the exact one: a machine
 * Real as is, an exact Integer / BigInt / Rational (or a wider MPFR) to the
 * nearest double on the outer side. */
static bool dl_leaf_d(const Expr* e, bool up, double* d) {
    if (e->type == EXPR_REAL) { *d = e->data.real; return true; }
#ifdef USE_MPFR
    if (dch_is_mpfr(e)) {
        *d = mpfr_get_d(e->data.mpfr, up ? MPFR_RNDU : MPFR_RNDD);
        return isfinite(*d);
    }
#endif
    if (e->type != EXPR_INTEGER && e->type != EXPR_BIGINT && !is_rational_like(e))
        return false;
    mpq_t q, t;
    mpq_init(q); mpq_init(t);
    mpz_t zn;                                  /* expr_to_mpz initialises its out */
    if (e->type == EXPR_FUNCTION) {            /* Rational[n, d] */
        mpz_t zd;
        expr_to_mpz(e->data.function.args[0], zn);
        expr_to_mpz(e->data.function.args[1], zd);
        mpq_set_num(q, zn);
        mpq_set_den(q, zd);
        mpq_canonicalize(q);
        mpz_clear(zd);
    } else {
        expr_to_mpz(e, zn);
        mpq_set_z(q, zn);
    }
    mpz_clear(zn);
    double v = mpq_get_d(q);                   /* truncated toward zero */
    bool ok = isfinite(v);
    if (ok) {
        mpq_set_d(t, v);
        int c = mpq_cmp(t, q);
        if (up && c < 0) v = nextafter(v, INFINITY);
        else if (!up && c > 0) v = nextafter(v, -INFINITY);
        *d = v;
        ok = isfinite(v);
    }
    mpq_clear(q); mpq_clear(t);
    return ok;
}

/* Interval[{a1, b1}, ...] for a real character: rigorous enclosure per piece. */
static Expr* dl_interval(const DchTable* T, const Expr* iv) {
    if (!T->real || !T->estar) return NULL;
    size_t np = iv->data.function.arg_count;
    if (np == 0) return NULL;
    Expr** pieces = (Expr**)malloc(np * sizeof(Expr*));
    bool ok = true;
    size_t built = 0;
    for (size_t i = 0; i < np && ok; i++) {
        const Expr* pr = iv->data.function.args[i];
        double lo, hi, ylo, yhi;
        ok = pr->type == EXPR_FUNCTION && pr->data.function.arg_count == 2
             && pr->data.function.head->type == EXPR_SYMBOL
             && pr->data.function.head->data.symbol.name == SYM_List
             && dl_leaf_d(pr->data.function.args[0], false, &lo)
             && dl_leaf_d(pr->data.function.args[1], true, &hi)
             && flint_dirichlet_l_interval(lo, hi, T->f, T->order, dl_chi_cb, T,
                                           T->neul, T->eul_p, T->eul_e, &ylo, &yhi);
        if (ok) pieces[built++] = mk2(SYM_List, expr_new_real(ylo), expr_new_real(yhi));
    }
    Expr* out = NULL;
    if (ok) out = eval_and_free(expr_new_function(expr_new_symbol(SYM_Interval), pieces, built));
    else for (size_t i = 0; i < built; i++) expr_free(pieces[i]);
    free(pieces);
    return out;
}

/* ------------------------------------------------------------------ */
/*  NDArray / packed s                                                 */
/* ------------------------------------------------------------------ */

/* Thread the machine kernel over an NDArray s with scalar exact (k, j). A real
 * character over a real buffer stays real; anything else is complex. NULL on
 * any decline -> the caller delists and threads the accurate scalar path. */
static Expr* dirichletl_ndarray(Expr* res) {
    Expr** a = res->data.function.args;
    if (is_ndarray(a[0]) || is_ndarray(a[1]) || !is_ndarray(a[2])) return NULL;
    if (a[0]->type != EXPR_INTEGER || a[1]->type != EXPR_INTEGER) return NULL;
    const Expr* nd = a[2];
    NDType dt = nd->data.ndarray.dtype;
    if (dt == NDT_INT64 || dt == NDT_BOOL) return NULL;   /* exact s: per element */
    DchTable* owned = NULL;
    const DchTable* T = dch_get(a[0]->data.integer, a[1]->data.integer, true, &owned, NULL);
    if (!T || !T->estar) { dch_free(owned); return NULL; }
    bool real_out = T->real && !ndt_is_complex(dt);
    NDType odt = real_out ? dt : ndt_as_complex(dt);
    size_t n = ndarray_size(nd);
    void* out = malloc(ndt_elem_size(odt) * (n ? n : 1));
    bool ok = out != NULL;
    for (size_t i = 0; i < n && ok; i++) {
        double sr, si, vr, vi = 0.0;
        ndt_get(nd->data.ndarray.data, i, dt, &sr, &si);
        if (T->k == 1) {
            ok = real_out ? sf_zeta_r(sr, &vr) : sf_zeta_c(sr, si, &vr, &vi);
        } else if (real_out) {
            ok = sf_dirichlet_l_tab_r(T, sr, &vr);
        } else {
            ok = sf_dirichlet_l_tab_c(T, sr, si, &vr, &vi);
        }
        if (ok) ndt_set(out, i, odt, vr, vi);
    }
    dch_free(owned);
    if (!ok) { free(out); return NULL; }
    return expr_new_ndarray_like(nd, nd->data.ndarray.rank, nd->data.ndarray.dims, out, odt);
}

/* ------------------------------------------------------------------ */
/*  Builtin                                                            */
/* ------------------------------------------------------------------ */

static bool is_interval(const Expr* e) {
    return e->type == EXPR_FUNCTION && e->data.function.head->type == EXPR_SYMBOL
        && e->data.function.head->data.symbol.name == SYM_Interval;
}

static bool has_mpfr(const Expr* e) {
    Expr *re, *im;
    if (dch_is_mpfr(e)) return true;
    if (is_complex((Expr*)e, &re, &im)) return dch_is_mpfr(re) || dch_is_mpfr(im);
    return false;
}

static bool is_inexact_number(const Expr* e) {
    Expr *re, *im;
    if (e->type == EXPR_REAL || dch_is_mpfr(e)) return true;
    if (is_complex((Expr*)e, &re, &im))
        return re->type == EXPR_REAL || dch_is_mpfr(re)
            || im->type == EXPR_REAL || dch_is_mpfr(im);
    return false;
}

static bool is_exact_number(const Expr* e) {
    Expr *re, *im;
    if (e->type == EXPR_INTEGER || e->type == EXPR_BIGINT || is_rational_like(e)) return true;
    if (is_complex((Expr*)e, &re, &im)) return !is_inexact_number(e);
    return false;
}

Expr* builtin_dirichletl(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    if (argc != 3) return builtin_arg_error("DirichletL", argc, 3, 3);
    Expr** a = res->data.function.args;

    if (is_ndarray(a[0]) || is_ndarray(a[1]) || is_ndarray(a[2])) {
        Expr* out = dirichletl_ndarray(res);
        if (out) return out;
        return ndarray_delist_and_reeval(res);
    }

    DchTable* owned = NULL;
    const DchTable* T = dch_validate_kj("DirichletL", res, true, &owned);
    if (!T) return NULL;
    Expr* s = a[2];
    Expr* out = NULL;

    if (T->k == 1) {
        /* L(chi_1, s) = Zeta[s]: delegate numbers and intervals wholesale, so
         * DirichletL[1, 1, s] is bit-identical to Zeta[s] on every surface. */
        if (s->type == EXPR_INTEGER || is_inexact_number(s) || is_interval(s))
            out = eval_and_free(mk1(SYM_Zeta, expr_copy(s)));
        goto done;
    }

    if (is_interval(s)) { out = dl_interval(T, s); goto done; }

    {
        SFMachineArgs ma;
        if (sf_machine_scalar_args(res, &ma)) {
            double sr = ma.re[2], si = ma.im[2], vr, vi;
            if (T->principal && sr == 1.0 && si == 0.0) {
                out = expr_new_symbol(SYM_ComplexInfinity);
                goto done;
            }
            if (!ma.is_complex[2] && T->real) {
                if (sf_dirichlet_l_tab_r(T, sr, &vr)) { out = expr_new_real(vr); goto done; }
            } else if (sf_dirichlet_l_tab_c(T, sr, si, &vr, &vi)) {
                out = sf_machine_complex_result(vr, vi);
                goto done;
            }
            out = dl_numeric(T, s, true);       /* machine decline: arb, rounded */
            goto done;
        }
    }

    if (is_exact_number(s)) { out = dl_exact(T, s); goto done; }

    if (has_mpfr(s)) {
#ifdef USE_MPFR
        if (T->principal && dch_is_mpfr(s) && mpfr_cmp_ui(s->data.mpfr, 1) == 0) {
            out = expr_new_symbol(SYM_ComplexInfinity);
            goto done;
        }
#endif
        out = dl_numeric(T, s, false);
    }

done:
    dch_free(owned);
    return out;
}

void dirichletl_init(void) {
    symtab_add_builtin("DirichletL", builtin_dirichletl);
    symtab_get_def("DirichletL")->attributes |=
        (ATTR_LISTABLE | ATTR_NHOLDALL | ATTR_PROTECTED);
    /* Docstring lives in info.c (info_init). */
}
