/* simp_kummer.c — denesting radicals over a Kummer field, and Root objects of
 * composed binomials, in native C.
 *
 * 1. Denesting.  Let r be a constant built from rationals and real radicals
 *    of positive rationals, and K = Q(p_1^(1/m_1), ..., p_k^(1/m_k)) the
 *    Kummer field its radicals generate (primes p_i, so the monomials
 *    prod p_i^(e_i/m_i), 0 <= e_i < m_i, are a Q-basis of K).  A radical
 *    r^(1/n) with r > 0 denests in a real radical extension of K exactly when
 *        r^(1/n) = c^(1/n) u,     c in Q_{>0},  u in K
 *    (Ramanujan's (2^(1/3) - 1)^(1/3) = (1 - 2^(1/3) + 4^(1/3))/9^(1/3);
 *    Sqrt[28^(1/3) - 3] = (98^(1/3) - 28^(1/3) - 1)/3).  The candidates for
 *    c are the products of a few small primes -- those of the generators, of
 *    the m_i and n, and of r's coordinates -- with exponents below n.  For
 *    each, u = r^(1/n)/c^(1/n) is found as an integer relation between its
 *    256-bit value and the basis values (FLINT fmpz_lll), then VERIFIED
 *    exactly: c u^n == r in K (rational arithmetic on the basis).  Nothing is
 *    returned without that proof.
 *
 * 2. Root objects.  Root[q(x^m) &, k] with q of degree <= 2, or a cubic that
 *    a Tschirnhaus shift makes binomial, is solved in radicals: y from q, then
 *    the real m-th root, matched to the Root's numeric value.  Ramanujan's
 *    Cos[2 Pi/7]^(1/3) - (-Cos[4 Pi/7])^(1/3) - (-Cos[8 Pi/7])^(1/3) is the
 *    Root of 4 x^9 - 30 x^6 + 75 x^3 + 32, i.e. -(3 7^(1/3)/2 - 5/2)^(1/3).
 */
#include "simp.h"
#include "simp_internal.h"
#include "arithmetic.h"
#include "eval.h"
#include "core.h"        /* leaf_count_internal */
#include "expr.h"
#include "numeric.h"
#include "sym_names.h"
#include <gmp.h>
#include <math.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#ifdef USE_MPFR
#include <mpfr.h>
#endif
#ifdef USE_FLINT
#include <flint/fmpz.h>
#include <flint/fmpz_mat.h>
#include <flint/fmpz_lll.h>
#endif

#define KF_MAX_PRIMES 4
#define KF_MAX_DIM    48
#define KF_MAX_N      12
#define KF_BITS       320

/* ---------------------------------------------------------------- field -- */

typedef struct {
    int     k;                       /* generator primes                */
    int64_t p[KF_MAX_PRIMES];
    int     m[KF_MAX_PRIMES];        /* radical degree per prime        */
    int     D;                       /* basis size = prod m             */
} KField;

typedef struct { mpq_t* c; int D; } KElt;    /* coordinates on the basis  */

static void ke_init(KElt* a, int D) {
    a->D = D; a->c = malloc(sizeof(mpq_t) * (size_t)D);
    for (int i = 0; i < D; i++) mpq_init(a->c[i]);
}
static void ke_clear(KElt* a) {
    for (int i = 0; i < a->D; i++) mpq_clear(a->c[i]);
    free(a->c); a->c = NULL;
}
static void ke_set(KElt* d, const KElt* s) { for (int i = 0; i < d->D; i++) mpq_set(d->c[i], s->c[i]); }
static bool ke_is_zero(const KElt* a) {
    for (int i = 0; i < a->D; i++) if (mpq_sgn(a->c[i])) return false;
    return true;
}

/* basis index <-> exponent vector (mixed radix) */
static void kf_exps(const KField* F, int idx, int* e) {
    for (int i = 0; i < F->k; i++) { e[i] = idx % F->m[i]; idx /= F->m[i]; }
}
static int kf_index(const KField* F, const int* e) {
    int idx = 0, mul = 1;
    for (int i = 0; i < F->k; i++) { idx += e[i] * mul; mul *= F->m[i]; }
    return idx;
}

/* d = a * b */
static void ke_mul(const KField* F, KElt* d, const KElt* a, const KElt* b) {
    KElt t; ke_init(&t, F->D);
    mpq_t q; mpq_init(q);
    int ea[KF_MAX_PRIMES], eb[KF_MAX_PRIMES], es[KF_MAX_PRIMES];
    for (int i = 0; i < F->D; i++) {
        if (!mpq_sgn(a->c[i])) continue;
        kf_exps(F, i, ea);
        for (int j = 0; j < F->D; j++) {
            if (!mpq_sgn(b->c[j])) continue;
            kf_exps(F, j, eb);
            mpq_mul(q, a->c[i], b->c[j]);
            for (int s = 0; s < F->k; s++) {
                es[s] = ea[s] + eb[s];
                if (es[s] >= F->m[s]) {                /* p^(m/m) = p */
                    es[s] -= F->m[s];
                    mpz_mul_si(mpq_numref(q), mpq_numref(q), (long)F->p[s]);
                }
            }
            mpq_canonicalize(q);
            int idx = kf_index(F, es);
            mpq_add(t.c[idx], t.c[idx], q);
        }
    }
    ke_set(d, &t);
    ke_clear(&t); mpq_clear(q);
}

/* d = a^-1 by solving (multiplication-by-a) x = 1; false when a == 0 */
static bool ke_inv(const KField* F, KElt* d, const KElt* a) {
    int D = F->D;
    if (ke_is_zero(a)) return false;
    /* column j = a * basis_j */
    mpq_t* M = malloc(sizeof(mpq_t) * (size_t)(D * (D + 1)));
    for (int i = 0; i < D * (D + 1); i++) mpq_init(M[i]);
    KElt bj, col; ke_init(&bj, D); ke_init(&col, D);
    for (int j = 0; j < D; j++) {
        for (int i = 0; i < D; i++) mpq_set_ui(bj.c[i], i == j, 1);
        ke_mul(F, &col, a, &bj);
        for (int i = 0; i < D; i++) mpq_set(M[i * (D + 1) + j], col.c[i]);
    }
    for (int i = 0; i < D; i++) mpq_set_ui(M[i * (D + 1) + D], i == 0, 1);
    ke_clear(&bj); ke_clear(&col);
    bool ok = true;
    mpq_t f; mpq_init(f);
    for (int c = 0; c < D && ok; c++) {
        int piv = -1;
        for (int r = c; r < D; r++) if (mpq_sgn(M[r * (D + 1) + c])) { piv = r; break; }
        if (piv < 0) { ok = false; break; }
        if (piv != c)
            for (int t = 0; t <= D; t++) mpq_swap(M[piv * (D + 1) + t], M[c * (D + 1) + t]);
        for (int r = 0; r < D; r++) {
            if (r == c || !mpq_sgn(M[r * (D + 1) + c])) continue;
            mpq_div(f, M[r * (D + 1) + c], M[c * (D + 1) + c]);
            for (int t = c; t <= D; t++) {
                mpq_t g; mpq_init(g);
                mpq_mul(g, f, M[c * (D + 1) + t]);
                mpq_sub(M[r * (D + 1) + t], M[r * (D + 1) + t], g);
                mpq_clear(g);
            }
        }
    }
    if (ok)
        for (int i = 0; i < D; i++) mpq_div(d->c[i], M[i * (D + 1) + D], M[i * (D + 1) + i]);
    mpq_clear(f);
    for (int i = 0; i < D * (D + 1); i++) mpq_clear(M[i]);
    free(M);
    return ok;
}

static void ke_pow(const KField* F, KElt* d, const KElt* a, int n) {
    KElt r, b; ke_init(&r, F->D); ke_init(&b, F->D);
    mpq_set_ui(r.c[0], 1, 1);
    ke_set(&b, a);
    while (n > 0) {
        if (n & 1) ke_mul(F, &r, &r, &b);
        n >>= 1;
        if (n) ke_mul(F, &b, &b, &b);
    }
    ke_set(d, &r);
    ke_clear(&r); ke_clear(&b);
}

/* N_{K/Q}(a) = det(multiplication by a) */
static void ke_norm(const KField* F, mpq_t out, const KElt* a) {
    int D = F->D;
    mpq_t* M = malloc(sizeof(mpq_t) * (size_t)(D * D));
    for (int i = 0; i < D * D; i++) mpq_init(M[i]);
    KElt bj, col; ke_init(&bj, D); ke_init(&col, D);
    for (int j = 0; j < D; j++) {
        for (int i = 0; i < D; i++) mpq_set_ui(bj.c[i], i == j, 1);
        ke_mul(F, &col, a, &bj);
        for (int i = 0; i < D; i++) mpq_set(M[i * D + j], col.c[i]);
    }
    ke_clear(&bj); ke_clear(&col);
    mpq_set_ui(out, 1, 1);
    mpq_t f, g; mpq_init(f); mpq_init(g);
    for (int c = 0; c < D; c++) {
        int piv = -1;
        for (int r = c; r < D; r++) if (mpq_sgn(M[r * D + c])) { piv = r; break; }
        if (piv < 0) { mpq_set_ui(out, 0, 1); break; }
        if (piv != c) { for (int t = 0; t < D; t++) mpq_swap(M[piv * D + t], M[c * D + t]); mpq_neg(out, out); }
        mpq_mul(out, out, M[c * D + c]);
        for (int r = c + 1; r < D; r++) {
            if (!mpq_sgn(M[r * D + c])) continue;
            mpq_div(f, M[r * D + c], M[c * D + c]);
            for (int t = c; t < D; t++) { mpq_mul(g, f, M[c * D + t]); mpq_sub(M[r * D + t], M[r * D + t], g); }
        }
    }
    mpq_clear(f); mpq_clear(g);
    for (int i = 0; i < D * D; i++) mpq_clear(M[i]);
    free(M);
}

/* v_p of a nonzero rational */
static int kf_val(const mpq_t q, int64_t p) {
    int v = 0;
    mpz_t t; mpz_init_set(t, mpq_numref(q));
    while (mpz_sgn(t) && mpz_divisible_ui_p(t, (unsigned long)p)) { mpz_divexact_ui(t, t, (unsigned long)p); v++; }
    mpz_set(t, mpq_denref(q));
    while (mpz_divisible_ui_p(t, (unsigned long)p)) { mpz_divexact_ui(t, t, (unsigned long)p); v--; }
    mpz_clear(t);
    return v;
}

/* ------------------------------------------------------- reading an Expr -- */

static bool kf_head(const Expr* e, const char* h) {
    return e && e->type == EXPR_FUNCTION && e->data.function.head->type == EXPR_SYMBOL
        && e->data.function.head->data.symbol.name == h;
}

static bool kf_rational(const Expr* e, mpq_t out) {
    if (e->type == EXPR_INTEGER) { mpq_set_si(out, (long)e->data.integer, 1); return true; }
    if (e->type == EXPR_BIGINT)  { mpq_set_z(out, e->data.bigint); return true; }
    if (kf_head(e, SYM_Rational) && e->data.function.arg_count == 2) {
        mpq_t a, b; mpq_init(a); mpq_init(b);
        bool ok = kf_rational(e->data.function.args[0], a) && kf_rational(e->data.function.args[1], b)
               && mpq_sgn(b);
        if (ok) mpq_div(out, a, b);
        mpq_clear(a); mpq_clear(b);
        return ok;
    }
    return false;
}

/* Factor |z| over small primes; false if a cofactor > 1 remains. */
static bool kf_factor_small(const mpz_t z, int64_t* ps, int* vs, int* np, int cap) {
    mpz_t t; mpz_init_set(t, z); mpz_abs(t, t);
    *np = 0;
    for (int64_t p = 2; p < 1000 && mpz_cmp_ui(t, 1) > 0; p++) {
        if (p > 2 && p % 2 == 0) continue;
        int v = 0;
        while (mpz_divisible_ui_p(t, (unsigned long)p)) { mpz_divexact_ui(t, t, (unsigned long)p); v++; }
        if (v) {
            if (*np >= cap) { mpz_clear(t); return false; }
            ps[*np] = p; vs[*np] = v; (*np)++;
        }
    }
    bool ok = mpz_cmp_ui(t, 1) == 0;
    mpz_clear(t);
    return ok;
}

/* A radical atom q^(a/b), q a positive rational, a/b a non-integer rational. */
static bool kf_radical_atom(const Expr* e, mpq_t q, int64_t* a, int64_t* b) {
    if (!kf_head(e, SYM_Power) || e->data.function.arg_count != 2) return false;
    const Expr* ex = e->data.function.args[1];
    if (!kf_head(ex, SYM_Rational)) return false;
    if (ex->data.function.args[0]->type != EXPR_INTEGER || ex->data.function.args[1]->type != EXPR_INTEGER)
        return false;
    if (!kf_rational(e->data.function.args[0], q) || mpq_sgn(q) <= 0) return false;
    *a = ex->data.function.args[0]->data.integer;
    *b = ex->data.function.args[1]->data.integer;
    return *b > 1;
}

static int64_t kf_gcd(int64_t a, int64_t b) { if (a < 0) a = -a; while (b) { int64_t t = a % b; a = b; b = t; } return a; }

/* Collect the generator primes and degrees of e's radicals. */
static bool kf_collect(const Expr* e, KField* F) {
    mpq_t q; mpq_init(q);
    int64_t a, b;
    bool ok = true;
    if (kf_radical_atom(e, q, &a, &b)) {
        for (int part = 0; part < 2 && ok; part++) {
            int64_t ps[8]; int vs[8], np;
            if (!kf_factor_small(part ? mpq_denref(q) : mpq_numref(q), ps, vs, &np, 8)) { ok = false; break; }
            for (int i = 0; i < np && ok; i++) {
                int64_t num = (int64_t)vs[i] * a, g = kf_gcd(num, b), den = b / g;
                if (den == 1) continue;
                int s = 0;
                while (s < F->k && F->p[s] != ps[i]) s++;
                if (s == F->k) {
                    if (F->k == KF_MAX_PRIMES) { ok = false; break; }
                    F->p[F->k] = ps[i]; F->m[F->k] = 1; F->k++;
                }
                int64_t l = F->m[s] / kf_gcd(F->m[s], den) * den;
                if (l > 64) { ok = false; break; }
                F->m[s] = (int)l;
            }
        }
    } else if (e->type == EXPR_FUNCTION) {
        if (kf_head(e, SYM_Plus) || kf_head(e, SYM_Times) || kf_head(e, SYM_Power)
            || kf_head(e, SYM_Rational)) {
            for (size_t i = 0; i < e->data.function.arg_count && ok; i++)
                ok = kf_collect(e->data.function.args[i], F);
        } else ok = false;
    } else if (e->type != EXPR_INTEGER && e->type != EXPR_BIGINT) ok = false;
    mpq_clear(q);
    return ok;
}

/* The monomial q^(a/b) as an element: coefficient times one basis vector. */
static bool kf_atom_elt(const KField* F, const mpq_t q, int64_t a, int64_t b, KElt* out) {
    for (int i = 0; i < F->D; i++) mpq_set_ui(out->c[i], 0, 1);
    int e[KF_MAX_PRIMES] = {0};
    mpq_t coef; mpq_init(coef); mpq_set_ui(coef, 1, 1);
    bool ok = true;
    for (int part = 0; part < 2 && ok; part++) {
        int64_t ps[8]; int vs[8], np;
        if (!kf_factor_small(part ? mpq_denref(q) : mpq_numref(q), ps, vs, &np, 8)) { ok = false; break; }
        for (int i = 0; i < np && ok; i++) {
            int s = 0;
            while (s < F->k && F->p[s] != ps[i]) s++;
            /* exponent of ps[i]: sign * vs * a / b = whole + j/m */
            int64_t num = (part ? -1 : 1) * (int64_t)vs[i] * a;
            int64_t mden = s < F->k ? F->m[s] : 1;
            if ((num * mden) % b) { ok = false; break; }
            int64_t t = num * mden / b;                  /* exponent * m */
            int64_t whole = t >= 0 ? t / mden : -((-t + mden - 1) / mden);
            int64_t frac = t - whole * mden;
            if (s == F->k && frac) { ok = false; break; }
            if (s < F->k) e[s] = (e[s] + (int)frac) % (int)mden;
            mpz_t pw; mpz_init(pw);
            mpz_ui_pow_ui(pw, (unsigned long)ps[i], (unsigned long)(whole >= 0 ? whole : -whole));
            if (whole >= 0) mpz_mul(mpq_numref(coef), mpq_numref(coef), pw);
            else            mpz_mul(mpq_denref(coef), mpq_denref(coef), pw);
            mpz_clear(pw);
            mpq_canonicalize(coef);
        }
    }
    if (ok) mpq_set(out->c[kf_index(F, e)], coef);
    mpq_clear(coef);
    return ok;
}

/* Read e as an element of K. */
static bool kf_read(const KField* F, const Expr* e, KElt* out) {
    mpq_t q; mpq_init(q);
    int64_t a, b;
    bool ok = true;
    if (kf_rational(e, q)) {
        for (int i = 0; i < F->D; i++) mpq_set_ui(out->c[i], 0, 1);
        mpq_set(out->c[0], q);
    } else if (kf_radical_atom(e, q, &a, &b)) {
        ok = kf_atom_elt(F, q, a, b, out);
    } else if (kf_head(e, SYM_Plus) || kf_head(e, SYM_Times)) {
        bool plus = kf_head(e, SYM_Plus);
        KElt t; ke_init(&t, F->D);
        for (int i = 0; i < F->D; i++) mpq_set_ui(out->c[i], 0, 1);
        if (!plus) mpq_set_ui(out->c[0], 1, 1);
        for (size_t i = 0; i < e->data.function.arg_count && ok; i++) {
            ok = kf_read(F, e->data.function.args[i], &t);
            if (!ok) break;
            if (plus) for (int j = 0; j < F->D; j++) mpq_add(out->c[j], out->c[j], t.c[j]);
            else ke_mul(F, out, out, &t);
        }
        ke_clear(&t);
    } else if (kf_head(e, SYM_Power) && e->data.function.arg_count == 2
               && e->data.function.args[1]->type == EXPR_INTEGER) {
        int64_t n = e->data.function.args[1]->data.integer;
        if (n > 64 || n < -64) ok = false;
        else {
            KElt t; ke_init(&t, F->D);
            ok = kf_read(F, e->data.function.args[0], &t);
            if (ok && n < 0) ok = ke_inv(F, &t, &t);
            if (ok) ke_pow(F, out, &t, (int)(n < 0 ? -n : n));
            ke_clear(&t);
        }
    } else ok = false;
    mpq_clear(q);
    return ok;
}

#if defined(USE_MPFR) && defined(USE_FLINT)

/* basis values: prod p^(e/m) */
static void kf_basis_values(const KField* F, mpfr_t* bv) {
    int e[KF_MAX_PRIMES];
    mpfr_t t; mpfr_init2(t, KF_BITS);
    for (int i = 0; i < F->D; i++) {
        kf_exps(F, i, e);
        mpfr_set_ui(bv[i], 1, MPFR_RNDN);
        for (int s = 0; s < F->k; s++) {
            if (!e[s]) continue;
            mpfr_set_si(t, (long)F->p[s], MPFR_RNDN);
            mpfr_rootn_ui(t, t, (unsigned long)F->m[s], MPFR_RNDN);
            mpfr_pow_ui(t, t, (unsigned long)e[s], MPFR_RNDN);
            mpfr_mul(bv[i], bv[i], t, MPFR_RNDN);
        }
    }
    mpfr_clear(t);
}

static void ke_value(const KField* F, const KElt* a, mpfr_t* bv, mpfr_t out) {
    mpfr_t t; mpfr_init2(t, KF_BITS);
    mpfr_set_ui(out, 0, MPFR_RNDN);
    for (int i = 0; i < F->D; i++) {
        if (!mpq_sgn(a->c[i])) continue;
        mpfr_mul_q(t, bv[i], a->c[i], MPFR_RNDN);
        mpfr_add(out, out, t, MPFR_RNDN);
    }
    mpfr_clear(t);
}

/* An integer relation a0 u + sum a_j b_j = 0 (a0 != 0) by LLL; u's
 * coordinates to *out.  false when the reduced basis has no such row. */
static bool kf_recognise(const KField* F, mpfr_t u, mpfr_t* bv, KElt* out) {
    int D = F->D, n = D + 1;
    fmpz_mat_t M; fmpz_mat_init(M, n, n + 1);
    mpfr_t t; mpfr_init2(t, KF_BITS);
    mpz_t z; mpz_init(z);
    const long scale = KF_BITS - 64;
    for (int i = 0; i < n; i++) {
        fmpz_one(fmpz_mat_entry(M, i, i));
        mpfr_mul_2si(t, i == 0 ? u : bv[i - 1], scale, MPFR_RNDN);
        mpfr_get_z(z, t, MPFR_RNDN);
        fmpz_set_mpz(fmpz_mat_entry(M, i, n), z);
    }
    fmpz_lll_t fl;
    fmpz_lll_context_init(fl, 0.99, 0.51, Z_BASIS, APPROX);
    fmpz_lll(M, NULL, fl);
    bool found = false;
    /* the first row with a nonzero u-coefficient and a small residual */
    for (int r = 0; r < n && !found; r++) {
        if (fmpz_is_zero(fmpz_mat_entry(M, r, 0))) continue;
        if (fmpz_bits(fmpz_mat_entry(M, r, n)) > 40) continue;
        bool small = true;
        for (int j = 0; j < n && small; j++) small = fmpz_bits(fmpz_mat_entry(M, r, j)) < 60;
        if (!small) continue;
        mpz_t a0; mpz_init(a0);
        fmpz_get_mpz(a0, fmpz_mat_entry(M, r, 0));
        for (int j = 0; j < D; j++) {
            fmpz_get_mpz(z, fmpz_mat_entry(M, r, j + 1));
            mpz_neg(z, z);
            mpq_set_num(out->c[j], z);
            mpq_set_den(out->c[j], a0);
            if (mpz_sgn(a0) < 0) { mpq_neg(out->c[j], out->c[j]); mpz_neg(mpq_denref(out->c[j]), mpq_denref(out->c[j])); }
            mpq_canonicalize(out->c[j]);
        }
        mpz_clear(a0);
        found = true;
    }
    mpz_clear(z); mpfr_clear(t);
    fmpz_mat_clear(M);
    return found;
}

/* ------------------------------------------------------ building the answer */

static Expr* kf_fn2(const char* h, Expr* a, Expr* b) {
    return expr_new_function(expr_new_symbol(h), (Expr*[]){ a, b }, 2);
}
/* an integer node, machine-sized when it fits (a BigInt 2 is not the
 * Integer 2: Sqrt[BigInt 2] would never cancel against Sqrt[2]) */
static Expr* kf_z_expr(const mpz_t z) { return expr_bigint_normalize(expr_new_bigint_from_mpz(z)); }
static Expr* kf_q_expr(const mpq_t q) {
    Expr* n = kf_z_expr(mpq_numref(q));
    if (mpz_cmp_ui(mpq_denref(q), 1) == 0) return n;
    return kf_fn2(SYM_Times, n, kf_fn2(SYM_Power, kf_z_expr(mpq_denref(q)), expr_new_integer(-1)));
}

/* c^(1/n) * u, term by term: each term a rational times ONE radical of an
 * integer (exponents per prime summed, whole parts pulled out). */
static Expr* kf_build(const KField* F, const int64_t* cp, const int* cv, int nc, int n, const KElt* u) {
    size_t cap = (size_t)F->D, w = 0;
    Expr** terms = malloc(sizeof(Expr*) * cap);
    int e[KF_MAX_PRIMES];
    for (int i = 0; i < F->D; i++) {
        if (!mpq_sgn(u->c[i])) continue;
        kf_exps(F, i, e);
        /* exponent per prime as num/den, den = lcm */
        int64_t P[2 * KF_MAX_PRIMES]; int64_t num[2 * KF_MAX_PRIMES], den[2 * KF_MAX_PRIMES];
        int np = 0;
        for (int s = 0; s < F->k; s++) { P[np] = F->p[s]; num[np] = e[s]; den[np] = F->m[s]; np++; }
        for (int j = 0; j < nc; j++) {
            int s = 0;
            while (s < np && P[s] != cp[j]) s++;
            if (s == np) { P[np] = cp[j]; num[np] = 0; den[np] = 1; np++; }
            /* num/den + cv/n */
            num[s] = num[s] * n + (int64_t)cv[j] * den[s];
            den[s] = den[s] * n;
            int64_t g = kf_gcd(num[s], den[s]); if (g) { num[s] /= g; den[s] /= g; }
        }
        mpq_t coef; mpq_init(coef); mpq_set(coef, u->c[i]);
        int64_t L = 1;
        for (int s = 0; s < np; s++) {
            int64_t whole = num[s] / den[s];
            num[s] -= whole * den[s];
            if (whole) {
                mpz_t pw; mpz_init(pw); mpz_ui_pow_ui(pw, (unsigned long)P[s], (unsigned long)whole);
                mpz_mul(mpq_numref(coef), mpq_numref(coef), pw); mpq_canonicalize(coef); mpz_clear(pw);
            }
            if (num[s]) L = L / kf_gcd(L, den[s]) * den[s];
        }
        mpz_t R; mpz_init_set_ui(R, 1);
        for (int s = 0; s < np; s++) {
            if (!num[s]) continue;
            mpz_t pw; mpz_init(pw);
            mpz_ui_pow_ui(pw, (unsigned long)P[s], (unsigned long)(num[s] * (L / den[s])));
            mpz_mul(R, R, pw); mpz_clear(pw);
        }
        Expr* t = kf_q_expr(coef);
        if (L > 1)
            t = kf_fn2(SYM_Times, t, kf_fn2(SYM_Power, kf_z_expr(R),
                                            make_rational(1, L)));
        terms[w++] = t;
        mpz_clear(R); mpq_clear(coef);
    }
    Expr* s = expr_new_function(expr_new_symbol(SYM_Plus), terms, w);
    free(terms);
    Expr* r = evaluate(s);
    expr_free(s);
    return r;
}

/* A short quotient u = A/B: LLL on the lattice {(B, u B) : B in Z-span of the
 * basis}, scaled to integers, returns a B with u B also short
 * (u = (1 + 5^(1/4))/(5^(1/4) - 1) is 4 terms as a sum, 2 + 2 as a quotient).
 * Fills B (integer coordinates) and A = u B; false if nothing shorter. */
static bool kf_short_quotient(const KField* F, const KElt* u, KElt* A, KElt* B) {
    int D = F->D;
    if (D < 2) return false;
    KElt bj, col; ke_init(&bj, D); ke_init(&col, D);
    mpq_t* cols = malloc(sizeof(mpq_t) * (size_t)(D * D));
    mpz_t L; mpz_init_set_ui(L, 1);
    for (int j = 0; j < D; j++) {
        for (int i = 0; i < D; i++) mpq_set_ui(bj.c[i], i == j, 1);
        ke_mul(F, &col, u, &bj);
        for (int i = 0; i < D; i++) {
            mpq_init(cols[j * D + i]); mpq_set(cols[j * D + i], col.c[i]);
            mpz_lcm(L, L, mpq_denref(col.c[i]));
        }
    }
    fmpz_mat_t M; fmpz_mat_init(M, D, 2 * D);
    mpz_t z; mpz_init(z);
    mpq_t t; mpq_init(t);
    for (int j = 0; j < D; j++) {
        fmpz_set_ui(fmpz_mat_entry(M, j, j), 1);
        for (int i = 0; i < D; i++) {
            mpq_set_z(t, L); mpq_mul(t, t, cols[j * D + i]);      /* integer */
            mpz_set(z, mpq_numref(t));
            fmpz_set_mpz(fmpz_mat_entry(M, j, D + i), z);
        }
    }
    fmpz_lll_t fl;
    fmpz_lll_context_init(fl, 0.99, 0.51, Z_BASIS, APPROX);
    fmpz_lll(M, NULL, fl);
    /* the reduced row with the fewest nonzero entries in (B, u B) */
    int best = -1, bsup = 1 << 30;
    for (int r = 0; r < D; r++) {
        int sup = 0, bnz = 0;
        for (int i = 0; i < 2 * D; i++)
            if (!fmpz_is_zero(fmpz_mat_entry(M, r, i))) { sup++; if (i < D) bnz++; }
        if (bnz && sup < bsup) { bsup = sup; best = r; }
    }
    bool ok = best >= 0;
    if (ok) {
        for (int i = 0; i < D; i++) {
            fmpz_get_mpz(z, fmpz_mat_entry(M, best, i));
            mpq_set_z(B->c[i], z);
        }
        ke_mul(F, A, u, B);
        /* cancel the common monomial: the smallest exponent per prime over
         * both supports (Sqrt[5] (1 + 5^(1/4)) / (Sqrt[5] (5^(1/4) - 1))) */
        int mn[KF_MAX_PRIMES], e[KF_MAX_PRIMES];
        for (int s = 0; s < F->k; s++) mn[s] = F->m[s];
        for (int pass = 0; pass < 2; pass++) {
            const KElt* X = pass ? B : A;
            for (int i = 0; i < D; i++) {
                if (!mpq_sgn(X->c[i])) continue;
                kf_exps(F, i, e);
                for (int s = 0; s < F->k; s++) if (e[s] < mn[s]) mn[s] = e[s];
            }
        }
        bool shift = false;
        for (int s = 0; s < F->k; s++) shift = shift || mn[s] > 0;
        if (shift) {
            KElt T; ke_init(&T, D);
            for (int pass = 0; pass < 2; pass++) {
                KElt* X = pass ? B : A;
                for (int i = 0; i < D; i++) mpq_set_ui(T.c[i], 0, 1);
                for (int i = 0; i < D; i++) {
                    if (!mpq_sgn(X->c[i])) continue;
                    kf_exps(F, i, e);
                    for (int s = 0; s < F->k; s++) e[s] -= mn[s];
                    mpq_set(T.c[kf_index(F, e)], X->c[i]);
                }
                ke_set(X, &T);
            }
            ke_clear(&T);
        }
    }
    mpq_clear(t); mpz_clear(z); mpz_clear(L);
    fmpz_mat_clear(M);
    for (int i = 0; i < D * D; i++) mpq_clear(cols[i]);
    free(cols);
    ke_clear(&bj); ke_clear(&col);
    return ok;
}

static int kf_support(const KElt* a) {
    int n = 0;
    for (int i = 0; i < a->D; i++) n += mpq_sgn(a->c[i]) != 0;
    return n;
}

/* r^(1/n), r > 0 in K: the denested form, or NULL. */
static Expr* kf_denest_root(const Expr* r, int n) {
    KField F; memset(&F, 0, sizeof F);
    if (!kf_collect(r, &F) || F.k == 0) return NULL;
    F.D = 1;
    for (int i = 0; i < F.k; i++) { F.D *= F.m[i]; if (F.D > KF_MAX_DIM) return NULL; }
    KElt R; ke_init(&R, F.D);
    if (!kf_read(&F, r, &R) || ke_is_zero(&R)) { ke_clear(&R); return NULL; }
    mpfr_t* bv = malloc(sizeof(mpfr_t) * (size_t)F.D);
    for (int i = 0; i < F.D; i++) mpfr_init2(bv[i], KF_BITS);
    kf_basis_values(&F, bv);
    mpfr_t rv, ev, cv, uv; mpfr_inits2(KF_BITS, rv, ev, cv, uv, (mpfr_ptr)0);
    ke_value(&F, &R, bv, rv);
    Expr* out = NULL;
    if (mpfr_sgn(rv) > 0) {
        mpfr_rootn_ui(ev, rv, (unsigned long)n, MPFR_RNDN);
        /* candidate primes for c */
        int64_t cp[KF_MAX_PRIMES + 4]; int ncp = 0;
        for (int i = 0; i < F.k; i++) cp[ncp++] = F.p[i];
        int64_t extra[] = { 2, 3, 5 };
        int64_t mm = n;
        for (int i = 0; i < F.k; i++) mm *= F.m[i];
        for (size_t t = 0; t < 3 && ncp < KF_MAX_PRIMES + 1; t++) {
            if (mm % extra[t]) continue;
            bool have = false;
            for (int i = 0; i < ncp; i++) have = have || cp[i] == extra[t];
            if (!have) cp[ncp++] = extra[t];
        }
        int total = 1;
        for (int i = 0; i < ncp; i++) total *= n;
        KElt U, T; ke_init(&U, F.D); ke_init(&T, F.D);
        int cvv[KF_MAX_PRIMES + 4];
        /* Norm prune: c u^n = r gives N(r) = c^D N(u)^n, so for every
         * candidate prime v_p(N(r)) - D v_p(c) == 0 (mod n). */
        mpq_t nr; mpq_init(nr);
        ke_norm(&F, nr, &R);
        int vn[KF_MAX_PRIMES + 4];
        for (int i = 0; i < ncp; i++) vn[i] = kf_val(nr, cp[i]);
        mpq_clear(nr);
        for (int code = 0; code < total && !out; code++) {
            int x = code;
            for (int i = 0; i < ncp; i++) { cvv[i] = x % n; x /= n; }
            bool admissible = true;
            for (int i = 0; i < ncp && admissible; i++)
                admissible = ((vn[i] - F.D * cvv[i]) % n + n) % n == 0;
            if (!admissible) continue;
            /* c^(1/n) */
            mpfr_set_ui(cv, 1, MPFR_RNDN);
            for (int i = 0; i < ncp; i++) {
                if (!cvv[i]) continue;
                mpfr_t pw; mpfr_init2(pw, KF_BITS);
                mpfr_set_si(pw, (long)cp[i], MPFR_RNDN);
                mpfr_pow_ui(pw, pw, (unsigned long)cvv[i], MPFR_RNDN);
                mpfr_rootn_ui(pw, pw, (unsigned long)n, MPFR_RNDN);
                mpfr_mul(cv, cv, pw, MPFR_RNDN);
                mpfr_clear(pw);
            }
            mpfr_div(uv, ev, cv, MPFR_RNDN);
            if (!kf_recognise(&F, uv, bv, &U)) continue;
            /* exact proof: c u^n == r */
            ke_pow(&F, &T, &U, n);
            mpq_t c; mpq_init(c); mpq_set_ui(c, 1, 1);
            for (int i = 0; i < ncp; i++) {
                mpz_t pw; mpz_init(pw); mpz_ui_pow_ui(pw, (unsigned long)cp[i], (unsigned long)cvv[i]);
                mpz_mul(mpq_numref(c), mpq_numref(c), pw); mpz_clear(pw);
            }
            bool eq = true;
            for (int i = 0; i < F.D && eq; i++) {
                mpq_t a; mpq_init(a); mpq_mul(a, T.c[i], c);
                eq = mpq_equal(a, R.c[i]);
                mpq_clear(a);
            }
            mpq_clear(c);
            if (eq) {
                out = kf_build(&F, cp, cvv, ncp, n, &U);
                /* the quotient form when it has fewer terms */
                KElt QA, QB; ke_init(&QA, F.D); ke_init(&QB, F.D);
                if (kf_short_quotient(&F, &U, &QA, &QB)
                    && kf_support(&QA) + kf_support(&QB) <= kf_support(&U)) {
                    Expr* num = kf_build(&F, cp, cvv, ncp, n, &QA);
                    Expr* den = kf_build(&F, NULL, NULL, 0, n, &QB);
                    Expr* qt = kf_fn2(SYM_Times, num, kf_fn2(SYM_Power, den, expr_new_integer(-1)));
                    Expr* qv = evaluate(qt);
                    expr_free(qt);
                    if (leaf_count_internal(qv, true) < leaf_count_internal(out, true)) {
                        expr_free(out);
                        out = qv;
                    } else expr_free(qv);
                }
                ke_clear(&QA); ke_clear(&QB);
            }
        }
        ke_clear(&U); ke_clear(&T);
    }
    mpfr_clears(rv, ev, cv, uv, (mpfr_ptr)0);
    for (int i = 0; i < F.D; i++) mpfr_clear(bv[i]);
    free(bv);
    ke_clear(&R);
    return out;
}

#endif /* USE_MPFR && USE_FLINT */

/* ------------------------------------------------------------ the walker -- */

/* Does e contain a radical (a Power with a non-integer rational exponent)? */
static bool kf_has_radical(const Expr* e) {
    if (!e || e->type != EXPR_FUNCTION) return false;
    if (kf_head(e, SYM_Power) && e->data.function.arg_count == 2
        && kf_head(e->data.function.args[1], SYM_Rational)) return true;
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (kf_has_radical(e->data.function.args[i])) return true;
    return false;
}

static Expr* kf_walk(const Expr* e, bool* changed) {
    if (!e || e->type != EXPR_FUNCTION) return expr_copy((Expr*)e);
    size_t n = e->data.function.arg_count;
    Expr** a = malloc(sizeof(Expr*) * (n ? n : 1));
    bool sub = false;
    for (size_t i = 0; i < n; i++) {
        a[i] = kf_walk(e->data.function.args[i], &sub);
    }
    Expr* node = expr_new_function(expr_copy(e->data.function.head), a, n);
    free(a);
    if (sub) { Expr* ev = evaluate(node); expr_free(node); node = ev; *changed = true; }
#if defined(USE_MPFR) && defined(USE_FLINT)
    /* Power[r, p/q] with r a radical expression */
    if (kf_head(node, SYM_Power) && node->data.function.arg_count == 2
        && kf_head(node->data.function.args[1], SYM_Rational)
        && kf_has_radical(node->data.function.args[0])) {
        const Expr* ex = node->data.function.args[1];
        if (ex->data.function.args[0]->type == EXPR_INTEGER && ex->data.function.args[1]->type == EXPR_INTEGER) {
            int64_t pn = ex->data.function.args[0]->data.integer, q = ex->data.function.args[1]->data.integer;
            if (q >= 2 && q <= KF_MAX_N && (pn == 1 || pn == -1)) {
                Expr* d = kf_denest_root(node->data.function.args[0], (int)q);
                if (d) {
                    if (pn == -1) {
                        Expr* inv = kf_fn2(SYM_Power, d, expr_new_integer(-1));
                        d = evaluate(inv);
                        expr_free(inv);
                    }
                    expr_free(node);
                    *changed = true;
                    return d;
                }
            }
        }
    }
#endif
    return node;
}

/* Max nesting depth of radicals in e. */
int simp_kummer_depth(const Expr* e);
static int kf_depth(const Expr* e) { return simp_kummer_depth(e); }
int simp_kummer_depth(const Expr* e) {
    if (!e || e->type != EXPR_FUNCTION) return 0;
    int m = 0;
    for (size_t i = 0; i < e->data.function.arg_count; i++) {
        int d = simp_kummer_depth(e->data.function.args[i]);
        if (d > m) m = d;
    }
    if (kf_head(e, SYM_Power) && e->data.function.arg_count == 2
        && kf_head(e->data.function.args[1], SYM_Rational)) m += 1;
    return m;
}

/* Pure radical arithmetic: rationals under Plus, Times and rational powers. */
static bool kf_pure(const Expr* e) {
    if (!e) return false;
    if (e->type == EXPR_INTEGER || e->type == EXPR_BIGINT) return true;
    if (e->type != EXPR_FUNCTION) return false;
    if (kf_head(e, SYM_Rational)) return true;
    if (!(kf_head(e, SYM_Plus) || kf_head(e, SYM_Times) || kf_head(e, SYM_Power))) return false;
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (!kf_pure(e->data.function.args[i])) return false;
    return true;
}

Expr* simp_kummer_denest(const Expr* e) {
    /* only pure radical constants: Simplify returns this form as its answer */
    if (!e || !kf_has_radical(e) || !kf_pure(e)) return NULL;
    bool changed = false;
    Expr* r = kf_walk(e, &changed);
    if (!changed || kf_depth(r) >= kf_depth(e)) { expr_free(r); return NULL; }
    return r;
}

/* ------------------------------------------- Root objects of composed binomials */

/* Coefficients of the polynomial body of Root[body &, k] in Slot[1]. */
static bool kr_poly(const Expr* body, mpq_t* c, int* deg, int cap) {
    for (int i = 0; i < cap; i++) mpq_set_ui(c[i], 0, 1);
    *deg = 0;
    size_t nt = kf_head(body, SYM_Plus) ? body->data.function.arg_count : 1;
    Expr* const* t = kf_head(body, SYM_Plus) ? body->data.function.args : (Expr* const*)&body;
    mpq_t q; mpq_init(q);
    bool ok = true;
    for (size_t i = 0; i < nt && ok; i++) {
        const Expr* ti = t[i];
        mpq_set_ui(q, 1, 1);
        int64_t pw = 0;
        size_t nf = kf_head(ti, SYM_Times) ? ti->data.function.arg_count : 1;
        Expr* const* f = kf_head(ti, SYM_Times) ? ti->data.function.args : (Expr* const*)&ti;
        for (size_t j = 0; j < nf && ok; j++) {
            const Expr* fj = f[j];
            mpq_t r; mpq_init(r);
            if (kf_rational(fj, r)) mpq_mul(q, q, r);
            else if (kf_head(fj, SYM_Slot)) pw += 1;
            else if (kf_head(fj, SYM_Power) && kf_head(fj->data.function.args[0], SYM_Slot)
                     && fj->data.function.args[1]->type == EXPR_INTEGER)
                pw += fj->data.function.args[1]->data.integer;
            else ok = false;
            mpq_clear(r);
        }
        if (ok && (pw < 0 || pw >= cap)) ok = false;
        if (ok) { mpq_add(c[pw], c[pw], q); if (pw > *deg) *deg = (int)pw; }
    }
    mpq_clear(q);
    return ok;
}

static double kr_value(const Expr* e) {
    Expr* n = numericalize_machine_unplanned(e);
    double v = NAN;
    if (n && n->type == EXPR_REAL) v = n->data.real;
    else if (n && n->type == EXPR_INTEGER) v = (double)n->data.integer;
    if (n) expr_free(n);
    return v;
}

Expr* simp_root_to_radicals(const Expr* e) {
    if (!kf_head(e, SYM_Root) || e->data.function.arg_count != 2) return NULL;
    const Expr* fn = e->data.function.args[0];
    if (!kf_head(fn, SYM_Function) || fn->data.function.arg_count != 1) return NULL;
    enum { CAP = 64 };
    mpq_t c[CAP];
    for (int i = 0; i < CAP; i++) mpq_init(c[i]);
    int deg;
    Expr* out = NULL;
    if (!kr_poly(fn->data.function.args[0], c, &deg, CAP) || deg < 2) goto done;
    {
        /* m = gcd of the exponents present; q(y) = sum c_{i m} y^i */
        int m = 0;
        for (int i = 1; i <= deg; i++) if (mpq_sgn(c[i])) m = m ? (int)kf_gcd(m, i) : i;
        int d = deg / m;
        if (d < 1 || d > 3 || (d == 1 && m == 1)) goto done;
        double target = kr_value(e);
        if (target != target) goto done;
        /* q coefficients a0..a3 */
        mpq_t A[4];
        for (int i = 0; i < 4; i++) { mpq_init(A[i]); if (i <= d) mpq_set(A[i], c[i * m]); }
        /* y-roots as Exprs (real only) */
        Expr* ys[3]; int ny = 0;
        if (d == 1) {
            mpq_t y; mpq_init(y); mpq_div(y, A[0], A[1]); mpq_neg(y, y);
            ys[ny++] = kf_q_expr(y); mpq_clear(y);
        } else if (d == 2) {
            /* y = (-b +- Sqrt[b^2 - 4 a c])/(2 a) */
            mpq_t disc, t; mpq_init(disc); mpq_init(t);
            mpq_mul(disc, A[1], A[1]); mpq_mul(t, A[2], A[0]); mpq_mul_2exp(t, t, 2); mpq_sub(disc, disc, t);
            if (mpq_sgn(disc) >= 0)
                for (int s = -1; s <= 1; s += 2) {
                    Expr* y = kf_fn2(SYM_Times, kf_fn2(SYM_Power, kf_fn2(SYM_Times, expr_new_integer(2), kf_q_expr(A[2])), expr_new_integer(-1)),
                                     kf_fn2(SYM_Plus, kf_fn2(SYM_Times, expr_new_integer(-1), kf_q_expr(A[1])),
                                            kf_fn2(SYM_Times, expr_new_integer(s), kf_fn2(SYM_Power, kf_q_expr(disc), make_rational(1, 2)))));
                    ys[ny++] = evaluate(y); expr_free(y);
                }
            mpq_clear(disc); mpq_clear(t);
        } else {
            /* cubic: shift y = z - a2/(3 a3); binomial iff the z coefficient vanishes:
             * a1 - a2^2/(3 a3) == 0; then z^3 = -(const)/a3 with
             * const = a0 - a1 a2/(3 a3) + 2 a2^3/(27 a3^2). */
            mpq_t s, t, u, k0; mpq_inits(s, t, u, k0, (mpq_ptr)0);
            mpq_mul(t, A[2], A[2]); mpq_set_ui(u, 3, 1); mpq_mul(u, u, A[3]); mpq_div(t, t, u);
            mpq_sub(s, A[1], t);                                     /* z coefficient */
            if (!mpq_sgn(s)) {
                mpq_t sh; mpq_init(sh); mpq_div(sh, A[2], u);        /* a2/(3 a3) */
                /* k0 = q(-sh) / a3 = value of q at y = -sh, over a3; z^3 = -k0 */
                mpq_t y0, acc; mpq_init(y0); mpq_init(acc); mpq_neg(y0, sh);
                mpq_set(acc, A[3]);
                for (int i = 2; i >= 0; i--) { mpq_mul(acc, acc, y0); mpq_add(acc, acc, A[i]); }
                mpq_div(k0, acc, A[3]); mpq_neg(k0, k0);             /* z^3 = k0 */
                /* real z = sign(k0) |k0|^(1/3) */
                mpq_t ak; mpq_init(ak); mpq_abs(ak, k0);
                Expr* z = kf_fn2(SYM_Power, kf_q_expr(ak), make_rational(1, 3));
                if (mpq_sgn(k0) < 0) z = kf_fn2(SYM_Times, expr_new_integer(-1), z);
                Expr* y = kf_fn2(SYM_Plus, z, kf_fn2(SYM_Times, expr_new_integer(-1), kf_q_expr(sh)));
                ys[ny++] = evaluate(y); expr_free(y);
                mpq_clear(ak); mpq_clear(y0); mpq_clear(acc); mpq_clear(sh);
            }
            mpq_clears(s, t, u, k0, (mpq_ptr)0);
        }
        /* x: the real m-th roots of each y; the one equal to the Root */
        for (int i = 0; i < ny && !out; i++) {
            double yv = kr_value(ys[i]);
            if (yv != yv) continue;
            Expr* cand[2]; int nc = 0;
            if (m == 1) cand[nc++] = expr_copy(ys[i]);
            else if (yv > 0) {
                cand[nc++] = kf_fn2(SYM_Power, expr_copy(ys[i]), make_rational(1, m));
                if (m % 2 == 0)
                    cand[nc++] = kf_fn2(SYM_Times, expr_new_integer(-1),
                                        kf_fn2(SYM_Power, expr_copy(ys[i]), make_rational(1, m)));
            } else if (m % 2 == 1) {
                cand[nc++] = kf_fn2(SYM_Times, expr_new_integer(-1),
                    kf_fn2(SYM_Power, kf_fn2(SYM_Times, expr_new_integer(-1), expr_copy(ys[i])),
                           make_rational(1, m)));
            }
            for (int j = 0; j < nc; j++) {
                Expr* v = evaluate(cand[j]);
                expr_free(cand[j]);
                double xv = kr_value(v);
                if (!out && fabs(xv - target) <= 1e-9 * (1 + fabs(target))) out = v;
                else expr_free(v);
            }
        }
        for (int i = 0; i < ny; i++) expr_free(ys[i]);
        for (int i = 0; i < 4; i++) mpq_clear(A[i]);
    }
done:
    for (int i = 0; i < CAP; i++) mpq_clear(c[i]);
    if (out) {
        Expr* dn = simp_kummer_denest(out);
        if (dn) { expr_free(out); out = dn; }
    }
    return out;
}
