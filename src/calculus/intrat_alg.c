/* intrat_alg.c — native (FLINT) implementations of the algorithms of
 * Bronstein, Symbolic Integration I, chapter 2, for integrands A/D with
 * A, D in Q[x]:
 *
 *   §2.2  HermiteReduce, quadratic version          hermite_quadratic
 *   §2.2  HermiteReduce, Mack's linear version      hermite_mack
 *   §2.3  Horowitz-Ostrogradsky                     hermite_horowitz
 *   §2.4  Rothstein-Trager log part                 logpart_rt
 *   §2.6  Czichowski log part (lex Groebner basis)  logpart_czichowski
 *   §2.8  LogToReal / LogToAtan (Rioboo)            intrat_alg_log_to_real
 *
 * (The Lazard-Rioboo-Trager log part stays in intrat.c; its subresultant chain
 * is native in intrat_subres.c.)
 *
 * All arithmetic is FLINT: fmpq_poly for Q[x] and Q[t], arrays of fmpq_poly
 * for Q[t][x], fmpq_mat for the Horowitz-Ostrogradsky system, the in-tree
 * FLINT Groebner engine (gb_flint.c) for Czichowski, and qqbar (exact
 * algebraic numbers) for the roots of LogToReal.  The evaluator is touched
 * only to convert the inputs and the final answer.
 *
 * Algebraic-number coefficients: the number-field engine src/poly/nfrat.c
 * holds K[x] = Q(theta)[x] as fmpq_mpoly reduced modulo theta's minimal
 * polynomial, with a certified gcd, but exports only Together / Cancel.  Here
 * a non-rational coefficient declines to the Expr path of intrat.c, whose
 * Together / Cancel steps rat.c already routes to nfrat.  Every routine of this
 * file generalises to K once nfrat exports its K[x] layer (from-Expr
 * conversion, divrem, gcd, xgcd): replace fmpq_poly by that type, Q[t]/q by
 * K[t]/q, and qqbar roots of K-polynomials by roots of their norms.
 */
#include "intrat_alg.h"
#include "expand.h"
#include "eval.h"
#include "sym_names.h"
#include "core.h"

#include <stdlib.h>
#include <string.h>
#include <stdint.h>

#ifdef USE_FLINT
#include <gmp.h>
#include <flint/flint.h>
#include <flint/fmpz.h>
#include <flint/fmpz_vec.h>
#include <flint/fmpq.h>
#include <flint/fmpz_poly.h>
#include <flint/fmpz_poly_factor.h>
#include <flint/fmpq_poly.h>
#include <flint/fmpq_mat.h>
#include <flint/qqbar.h>
#include "groebner.h"
#include "gb_flint.h"

/* Size caps: past them the generic path is no worse off, and a pathological
 * input cannot run away (qqbar arithmetic in particular). */
#define IA_MAX_DEG        400   /* degree in x or t of any input           */
#define IA_MAX_LTR_DEG      6   /* LogToReal with Root objects: degree of one residue factor */
#define IA_MAX_LTR_DEG_RAD  4   /* LogToReal in radicals only (Automatic)                   */
#define IA_MAX_LTR_BITS   256   /* LogToReal: coefficient size of a residue factor          */
#define IA_MAX_RT_DEG       8   /* Rothstein-Trager: degree of one irreducible factor; past
                                 * it the gcd over Q(a) costs more than LRT's subresultants */
#define IA_MAX_QQBAR_DEG   32   /* LogToReal: algebraic degree of a coefficient */

/* ================================================================== */
/* Expr <-> Q[x], Q[t][x]                                              */
/* ================================================================== */

static bool ia_head_is(const Expr* e, const char* h) {
    return e && e->type == EXPR_FUNCTION && e->data.function.head->type == EXPR_SYMBOL
        && e->data.function.head->data.symbol.name == h;
}

static bool ia_sym_is(const Expr* e, const Expr* s) {
    return s && e->type == EXPR_SYMBOL && s->type == EXPR_SYMBOL
        && (e->data.symbol.name == s->data.symbol.name
            || strcmp(e->data.symbol.name, s->data.symbol.name) == 0);
}

static bool ia_int_leaf(const Expr* e, fmpz_t z) {
    if (e->type == EXPR_INTEGER) { fmpz_set_si(z, (slong)e->data.integer); return true; }
    if (e->type == EXPR_BIGINT) { fmpz_set_mpz(z, e->data.bigint); return true; }
    return false;
}

static bool ia_mul_rational(fmpq_t q, const Expr* e) {
    fmpz_t n, d;
    fmpz_init(n); fmpz_init(d);
    bool ok = false;
    if (ia_int_leaf(e, n)) {
        fmpq_mul_fmpz(q, q, n);
        ok = true;
    } else if (ia_head_is(e, SYM_Rational) && e->data.function.arg_count == 2
               && ia_int_leaf(e->data.function.args[0], n)
               && ia_int_leaf(e->data.function.args[1], d) && !fmpz_is_zero(d)) {
        fmpq_t r;
        fmpq_init(r);
        fmpq_set_fmpz_frac(r, n, d);
        fmpq_mul(q, q, r);
        fmpq_clear(r);
        ok = true;
    }
    fmpz_clear(n); fmpz_clear(d);
    return ok;
}

/* One factor of a monomial c x^i t^j (t may be NULL). */
static bool ia_factor(const Expr* f, const Expr* x, const Expr* t, fmpq_t c, slong* i, slong* j) {
    if (ia_sym_is(f, x)) { (*i)++; return true; }
    if (t && ia_sym_is(f, t)) { (*j)++; return true; }
    if (ia_head_is(f, SYM_Power) && f->data.function.arg_count == 2
        && f->data.function.args[1]->type == EXPR_INTEGER) {
        int64_t k = f->data.function.args[1]->data.integer;
        const Expr* b = f->data.function.args[0];
        if (k < 1 || k > IA_MAX_DEG) return false;
        if (ia_sym_is(b, x)) { *i += (slong)k; return true; }
        if (t && ia_sym_is(b, t)) { *j += (slong)k; return true; }
        return false;
    }
    return ia_mul_rational(c, f);
}

/* Dense Q[t][x]: c[k] is the coefficient (in Q[t]) of x^k; len = 0 for 0. */
typedef struct { slong len; slong cap; fmpq_poly_struct* c; } QTX;

static void qtx_init(QTX* p, slong cap) {
    if (cap < 1) cap = 1;
    p->len = 0; p->cap = cap;
    p->c = (fmpq_poly_struct*)flint_malloc(sizeof(fmpq_poly_struct) * (size_t)cap);
    for (slong k = 0; k < cap; k++) fmpq_poly_init(p->c + k);
}

static void qtx_clear(QTX* p) {
    for (slong k = 0; k < p->cap; k++) fmpq_poly_clear(p->c + k);
    flint_free(p->c);
    p->c = NULL; p->len = p->cap = 0;
}

static void qtx_normalise(QTX* p) {
    while (p->len > 0 && fmpq_poly_is_zero(p->c + p->len - 1)) p->len--;
}

/* p (any form) -> Q[t][x]; t == NULL parses a univariate p in x (into c[k]
 * constant polynomials). */
static bool qtx_from_expr(const Expr* p, const Expr* x, const Expr* t, QTX* out) {
    Expr* e = expr_expand((Expr*)p);
    if (!e) return false;
    bool plus = ia_head_is(e, SYM_Plus);
    size_t n = plus ? e->data.function.arg_count : 1;
    slong* is = (slong*)malloc(sizeof(slong) * n);
    slong* js = (slong*)malloc(sizeof(slong) * n);
    fmpq* cs = (fmpq*)flint_malloc(sizeof(fmpq) * n);
    for (size_t k = 0; k < n; k++) fmpq_init(cs + k);
    bool ok = true;
    slong maxi = -1;
    for (size_t k = 0; k < n && ok; k++) {
        const Expr* term = plus ? e->data.function.args[k] : e;
        fmpq_one(cs + k);
        is[k] = js[k] = 0;
        if (ia_head_is(term, SYM_Times)) {
            for (size_t m = 0; m < term->data.function.arg_count && ok; m++)
                ok = ia_factor(term->data.function.args[m], x, t, cs + k, is + k, js + k);
        } else {
            ok = ia_factor(term, x, t, cs + k, is + k, js + k);
        }
        if (ok && (is[k] > IA_MAX_DEG || js[k] > IA_MAX_DEG)) ok = false;
        if (ok && !fmpq_is_zero(cs + k) && is[k] > maxi) maxi = is[k];
    }
    if (ok) {
        qtx_init(out, maxi + 1);
        fmpq_t cur;
        fmpq_init(cur);
        for (size_t k = 0; k < n; k++) {
            if (fmpq_is_zero(cs + k)) continue;
            fmpq_poly_get_coeff_fmpq(cur, out->c + is[k], js[k]);
            fmpq_add(cur, cur, cs + k);
            fmpq_poly_set_coeff_fmpq(out->c + is[k], js[k], cur);
        }
        fmpq_clear(cur);
        out->len = maxi + 1;
        qtx_normalise(out);
    }
    for (size_t k = 0; k < n; k++) fmpq_clear(cs + k);
    flint_free(cs); free(is); free(js);
    expr_free(e);
    return ok;
}

/* p (any form) -> Q[x]; false unless a polynomial in x over Q. */
static bool qx_from_expr(const Expr* p, const Expr* x, fmpq_poly_t out) {
    QTX q;
    if (!qtx_from_expr(p, x, NULL, &q)) return false;
    fmpq_poly_zero(out);
    fmpq_t c;
    fmpq_init(c);
    for (slong k = 0; k < q.len; k++) {
        fmpq_poly_get_coeff_fmpq(c, q.c + k, 0);
        fmpq_poly_set_coeff_fmpq(out, k, c);
    }
    fmpq_clear(c);
    qtx_clear(&q);
    return true;
}

static Expr* ia_expr_fmpz(const fmpz_t z) {
    if (fmpz_fits_si(z)) return expr_new_integer((int64_t)fmpz_get_si(z));
    mpz_t m;
    mpz_init(m);
    fmpz_get_mpz(m, z);
    return expr_new_bigint_take(m);
}

static Expr* ia_expr_fmpq(const fmpq_t q) {
    if (fmpz_is_one(fmpq_denref(q))) return ia_expr_fmpz(fmpq_numref(q));
    return expr_new_function(expr_new_symbol(SYM_Rational),
        (Expr*[]){ ia_expr_fmpz(fmpq_numref(q)), ia_expr_fmpz(fmpq_denref(q)) }, 2);
}

static Expr* ia_pow(const Expr* v, slong k) {
    if (k == 1) return expr_copy((Expr*)v);
    return expr_new_function(expr_new_symbol(SYM_Power),
        (Expr*[]){ expr_copy((Expr*)v), expr_new_integer((int64_t)k) }, 2);
}

/* Unevaluated Plus of c * t^j * x^k over the nonzero coefficients. */
static Expr* ia_sum(Expr** terms, size_t n) {
    if (n == 0) { free(terms); return expr_new_integer(0); }
    Expr* r = n == 1 ? terms[0] : expr_new_function(expr_new_symbol(SYM_Plus), terms, n);
    free(terms);
    return r;
}

static Expr* qtx_to_expr(const QTX* p, const Expr* x, const Expr* t) {
    size_t n = 0;
    for (slong k = 0; k < p->len; k++) n += (size_t)fmpq_poly_length(p->c + k);
    Expr** terms = (Expr**)malloc(sizeof(Expr*) * (n ? n : 1));
    size_t m = 0;
    fmpq_t c;
    fmpq_init(c);
    for (slong k = 0; k < p->len; k++) {
        for (slong j = 0; j < fmpq_poly_length(p->c + k); j++) {
            fmpq_poly_get_coeff_fmpq(c, p->c + k, j);
            if (fmpq_is_zero(c)) continue;
            Expr* fs[3];
            size_t nf = 0;
            fs[nf++] = ia_expr_fmpq(c);
            if (j > 0) fs[nf++] = ia_pow(t, j);
            if (k > 0) fs[nf++] = ia_pow(x, k);
            terms[m++] = nf == 1 ? fs[0] : expr_new_function(expr_new_symbol(SYM_Times), fs, nf);
        }
    }
    fmpq_clear(c);
    return eval_and_free(ia_sum(terms, m));
}

static Expr* qx_to_expr(const fmpq_poly_t p, const Expr* x) {
    slong n = fmpq_poly_length(p);
    Expr** terms = (Expr**)malloc(sizeof(Expr*) * (size_t)(n ? n : 1));
    size_t m = 0;
    fmpq_t c;
    fmpq_init(c);
    for (slong k = 0; k < n; k++) {
        fmpq_poly_get_coeff_fmpq(c, p, k);
        if (fmpq_is_zero(c)) continue;
        Expr* cc = ia_expr_fmpq(c);
        terms[m++] = k == 0 ? cc
            : expr_new_function(expr_new_symbol(SYM_Times), (Expr*[]){ cc, ia_pow(x, k) }, 2);
    }
    fmpq_clear(c);
    return eval_and_free(ia_sum(terms, m));
}

/* n/d as an Expr, after cancelling their gcd and making d monic. */
static Expr* qx_fraction_to_expr(const fmpq_poly_t n, const fmpq_poly_t d, const Expr* x) {
    if (fmpq_poly_is_zero(n)) return expr_new_integer(0);
    fmpq_poly_t g, nn, dd, r;
    fmpq_poly_init(g); fmpq_poly_init(nn); fmpq_poly_init(dd); fmpq_poly_init(r);
    fmpq_poly_gcd(g, n, d);
    fmpq_poly_divrem(nn, r, n, g);
    fmpq_poly_divrem(dd, r, d, g);
    /* scale so the denominator is primitive over Z with a positive lc */
    fmpq_t c;
    fmpq_init(c);
    fmpq_poly_content(c, dd);
    {
        fmpq_t l;
        fmpq_init(l);
        fmpq_poly_get_coeff_fmpq(l, dd, fmpq_poly_degree(dd));
        if (fmpq_sgn(l) < 0) fmpq_neg(c, c);
        fmpq_clear(l);
    }
    fmpq_poly_scalar_div_fmpq(dd, dd, c);
    fmpq_poly_scalar_div_fmpq(nn, nn, c);
    Expr* en = qx_to_expr(nn, x);
    Expr* ed = qx_to_expr(dd, x);
    fmpq_clear(c);
    fmpq_poly_clear(g); fmpq_poly_clear(nn); fmpq_poly_clear(dd); fmpq_poly_clear(r);
    if (ed->type == EXPR_INTEGER && ed->data.integer == 1) { expr_free(ed); return en; }
    return eval_and_free(expr_new_function(expr_new_symbol(SYM_Times), (Expr*[]){ en,
        expr_new_function(expr_new_symbol(SYM_Power), (Expr*[]){ ed, expr_new_integer(-1) }, 2) }, 2));
}

/* ================================================================== */
/* Q[x] helpers                                                        */
/* ================================================================== */

/* q := a / b, exact; false (q untouched) if b does not divide a. */
static bool qx_divexact(fmpq_poly_t q, const fmpq_poly_t a, const fmpq_poly_t b) {
    fmpq_poly_t qq, r;
    fmpq_poly_init(qq); fmpq_poly_init(r);
    fmpq_poly_divrem(qq, r, a, b);
    bool ok = fmpq_poly_is_zero(r);
    if (ok) fmpq_poly_swap(q, qq);
    fmpq_poly_clear(qq); fmpq_poly_clear(r);
    return ok;
}

/* Bronstein's ExtendedEuclidean(a, b, c): s a + t b = c with s = 0 or
 * deg s < deg b.  False when gcd(a, b) does not divide c. */
static bool qx_ext_euclid(fmpq_poly_t s, fmpq_poly_t t, const fmpq_poly_t a,
                          const fmpq_poly_t b, const fmpq_poly_t c) {
    fmpq_poly_t G, S, T, q, r;
    fmpq_poly_init(G); fmpq_poly_init(S); fmpq_poly_init(T);
    fmpq_poly_init(q); fmpq_poly_init(r);
    bool ok = false;
    fmpq_poly_xgcd(G, S, T, a, b);       /* G = S a + T b, G monic */
    fmpq_poly_divrem(q, r, c, G);
    if (fmpq_poly_is_zero(r)) {
        fmpq_poly_mul(S, S, q);
        fmpq_poly_rem(s, S, b);          /* s = q S mod b */
        fmpq_poly_mul(T, s, a);
        fmpq_poly_sub(T, c, T);          /* c - s a */
        ok = qx_divexact(t, T, b);
    }
    fmpq_poly_clear(G); fmpq_poly_clear(S); fmpq_poly_clear(T);
    fmpq_poly_clear(q); fmpq_poly_clear(r);
    return ok;
}

/* gn/gd += b/v (unreduced). */
static void qx_frac_add(fmpq_poly_t gn, fmpq_poly_t gd, const fmpq_poly_t b, const fmpq_poly_t v) {
    fmpq_poly_t t1;
    fmpq_poly_init(t1);
    fmpq_poly_mul(gn, gn, v);
    fmpq_poly_mul(t1, b, gd);
    fmpq_poly_add(gn, gn, t1);
    fmpq_poly_mul(gd, gd, v);
    /* keep it reduced so the products stay small */
    fmpq_poly_gcd(t1, gn, gd);
    if (fmpq_poly_degree(t1) > 0) {
        qx_divexact(gn, gn, t1);
        qx_divexact(gd, gd, t1);
    }
    fmpq_poly_clear(t1);
}

/* The squarefree factorisation D = c * prod_k D_k^k with each D_k monic
 * (D_k = 1 when absent).  Returns m = the largest k; *fac (length m + 1,
 * index 0 unused) is allocated and initialised; c in *lc. */
static slong qx_squarefree(fmpq_poly_struct** fac, fmpq_t lc, const fmpq_poly_t D) {
    fmpz_poly_t Dz;
    fmpz_poly_init(Dz);
    fmpq_poly_get_numerator(Dz, D);
    fmpz_poly_factor_t f;
    fmpz_poly_factor_init(f);
    fmpz_poly_factor_squarefree(f, Dz);
    slong m = 0;
    for (slong i = 0; i < f->num; i++) if (f->exp[i] > m) m = f->exp[i];
    fmpq_poly_struct* F = (fmpq_poly_struct*)flint_malloc(sizeof(fmpq_poly_struct) * (size_t)(m + 1));
    for (slong k = 0; k <= m; k++) { fmpq_poly_init(F + k); fmpq_poly_one(F + k); }
    fmpq_poly_t pq;
    fmpq_poly_init(pq);
    for (slong i = 0; i < f->num; i++) {
        fmpq_poly_set_fmpz_poly(pq, f->p + i);
        fmpq_poly_make_monic(pq, pq);
        fmpq_poly_mul(F + f->exp[i], F + f->exp[i], pq);
    }
    /* lc = D / prod D_k^k (a constant) */
    fmpq_poly_t P, t1;
    fmpq_poly_init(P); fmpq_poly_init(t1);
    fmpq_poly_one(P);
    for (slong k = 1; k <= m; k++) {
        fmpq_poly_pow(t1, F + k, (ulong)k);
        fmpq_poly_mul(P, P, t1);
    }
    fmpq_t a, b;
    fmpq_init(a); fmpq_init(b);
    fmpq_poly_get_coeff_fmpq(a, D, fmpq_poly_degree(D));
    fmpq_poly_get_coeff_fmpq(b, P, fmpq_poly_degree(P));
    fmpq_div(lc, a, b);
    fmpq_clear(a); fmpq_clear(b);
    fmpq_poly_clear(P); fmpq_poly_clear(t1); fmpq_poly_clear(pq);
    fmpz_poly_factor_clear(f);
    fmpz_poly_clear(Dz);
    *fac = F;
    return m;
}

static void qx_free_factors(fmpq_poly_struct* F, slong m) {
    for (slong k = 0; k <= m; k++) fmpq_poly_clear(F + k);
    flint_free(F);
}

/* ================================================================== */
/* Hermite reduction                                                    */
/* ================================================================== */

/* Mack's linear version (Bronstein §2.2, HermiteReduce, linear). */
static bool hermite_mack(fmpq_poly_t gn, fmpq_poly_t gd, fmpq_poly_t hn, fmpq_poly_t hd,
                         const fmpq_poly_t A0, const fmpq_poly_t D) {
    fmpq_poly_t A, Dp, Dbar, Dstar, Dbp, Dbar2, Dbstar, dtil, B, C, t1;
    fmpq_poly_init(A); fmpq_poly_init(Dp); fmpq_poly_init(Dbar); fmpq_poly_init(Dstar);
    fmpq_poly_init(Dbp); fmpq_poly_init(Dbar2); fmpq_poly_init(Dbstar); fmpq_poly_init(dtil);
    fmpq_poly_init(B); fmpq_poly_init(C); fmpq_poly_init(t1);
    bool ok = true;
    fmpq_poly_set(A, A0);
    fmpq_poly_zero(gn); fmpq_poly_one(gd);
    fmpq_poly_derivative(Dp, D);
    fmpq_poly_gcd(Dbar, D, Dp);
    ok = qx_divexact(Dstar, D, Dbar);
    while (ok && fmpq_poly_degree(Dbar) > 0) {
        fmpq_poly_derivative(Dbp, Dbar);
        fmpq_poly_gcd(Dbar2, Dbar, Dbp);
        ok = qx_divexact(Dbstar, Dbar, Dbar2);
        if (!ok) break;
        fmpq_poly_mul(t1, Dstar, Dbp);
        ok = qx_divexact(dtil, t1, Dbar);
        if (!ok) break;
        fmpq_poly_neg(dtil, dtil);
        ok = qx_ext_euclid(B, C, dtil, Dbstar, A);   /* B dtil + C Dbstar = A */
        if (!ok) break;
        fmpq_poly_derivative(t1, B);
        fmpq_poly_mul(t1, t1, Dstar);
        ok = qx_divexact(t1, t1, Dbstar);
        if (!ok) break;
        fmpq_poly_sub(A, C, t1);
        qx_frac_add(gn, gd, B, Dbar);
        fmpq_poly_swap(Dbar, Dbar2);
    }
    if (ok) { fmpq_poly_set(hn, A); fmpq_poly_set(hd, Dstar); }
    fmpq_poly_clear(A); fmpq_poly_clear(Dp); fmpq_poly_clear(Dbar); fmpq_poly_clear(Dstar);
    fmpq_poly_clear(Dbp); fmpq_poly_clear(Dbar2); fmpq_poly_clear(Dbstar); fmpq_poly_clear(dtil);
    fmpq_poly_clear(B); fmpq_poly_clear(C); fmpq_poly_clear(t1);
    return ok;
}

/* The quadratic version (Bronstein §2.2): partial fractions over the
 * squarefree factorisation, then one reduction per power. */
static bool hermite_quadratic(fmpq_poly_t gn, fmpq_poly_t gd, fmpq_poly_t hn, fmpq_poly_t hd,
                              const fmpq_poly_t A0, const fmpq_poly_t D) {
    fmpq_poly_struct* F;
    fmpq_t lc;
    fmpq_init(lc);
    slong m = qx_squarefree(&F, lc, D);
    fmpq_poly_t A, E, Rest, inv, G, T, Bk, V, Vp, B, C, t1, rhs;
    fmpq_poly_init(A); fmpq_poly_init(E); fmpq_poly_init(Rest); fmpq_poly_init(inv);
    fmpq_poly_init(G); fmpq_poly_init(T); fmpq_poly_init(Bk); fmpq_poly_init(V);
    fmpq_poly_init(Vp); fmpq_poly_init(B); fmpq_poly_init(C); fmpq_poly_init(t1);
    fmpq_poly_init(rhs);
    bool ok = true;
    fmpq_poly_scalar_div_fmpq(A, A0, lc);          /* A/D = A1 / prod D_k^k */
    fmpq_poly_zero(gn); fmpq_poly_one(gd);
    fmpq_poly_zero(hn); fmpq_poly_one(hd);
    for (slong k = 1; k <= m && ok; k++) {
        if (fmpq_poly_degree(F + k) <= 0) continue;
        /* B_k = A1 * (prod_{j != k} D_j^j)^(-1) mod D_k^k */
        fmpq_poly_pow(E, F + k, (ulong)k);
        fmpq_poly_one(Rest);
        for (slong j = 1; j <= m; j++) {
            if (j == k || fmpq_poly_degree(F + j) <= 0) continue;
            fmpq_poly_pow(t1, F + j, (ulong)j);
            fmpq_poly_mul(Rest, Rest, t1);
        }
        fmpq_poly_xgcd(G, inv, T, Rest, E);
        if (fmpq_poly_degree(G) != 0) { ok = false; break; }
        fmpq_poly_mul(Bk, A, inv);
        fmpq_poly_rem(Bk, Bk, E);
        /* reduce B_k / V^k, V = D_k */
        fmpq_poly_set(V, F + k);
        fmpq_poly_derivative(Vp, V);
        for (slong j = k - 1; j >= 1 && ok; j--) {
            fmpq_poly_scalar_div_si(rhs, Bk, -(slong)j);    /* -B_k/j */
            ok = qx_ext_euclid(B, C, Vp, V, rhs);         /* B V' + C V = -B_k/j */
            if (!ok) break;
            fmpq_poly_pow(t1, V, (ulong)j);
            qx_frac_add(gn, gd, B, t1);                    /* g += B/V^j */
            fmpq_poly_derivative(t1, B);
            fmpq_poly_scalar_mul_si(Bk, C, -(slong)j);
            fmpq_poly_sub(Bk, Bk, t1);                     /* B_k = -j C - B' */
        }
        if (ok) qx_frac_add(hn, hd, Bk, V);                /* h += B_k / V */
    }
    fmpq_poly_clear(A); fmpq_poly_clear(E); fmpq_poly_clear(Rest); fmpq_poly_clear(inv);
    fmpq_poly_clear(G); fmpq_poly_clear(T); fmpq_poly_clear(Bk); fmpq_poly_clear(V);
    fmpq_poly_clear(Vp); fmpq_poly_clear(B); fmpq_poly_clear(C); fmpq_poly_clear(t1);
    fmpq_poly_clear(rhs);
    fmpq_clear(lc);
    qx_free_factors(F, m);
    return ok;
}

/* Horowitz-Ostrogradsky (Bronstein §2.3): with D^- = gcd(D, D'), D* = D/D^-,
 * solve A = B' D* - B (D* D^-'/D^-) + C D^- for deg B < deg D^-,
 * deg C < deg D*, a square linear system over Q; g = B/D^-, h = C/D*. */
static bool hermite_horowitz(fmpq_poly_t gn, fmpq_poly_t gd, fmpq_poly_t hn, fmpq_poly_t hd,
                             const fmpq_poly_t A, const fmpq_poly_t D) {
    fmpq_poly_t Dp, Dm, Ds, Dmp, H, col, t1;
    fmpq_poly_init(Dp); fmpq_poly_init(Dm); fmpq_poly_init(Ds); fmpq_poly_init(Dmp);
    fmpq_poly_init(H); fmpq_poly_init(col); fmpq_poly_init(t1);
    bool ok = true;
    fmpq_poly_derivative(Dp, D);
    fmpq_poly_gcd(Dm, D, Dp);
    ok = qx_divexact(Ds, D, Dm);
    slong n = fmpq_poly_degree(Dm), m = fmpq_poly_degree(Ds), N = n + m;
    if (ok) {
        fmpq_poly_derivative(Dmp, Dm);
        fmpq_poly_mul(t1, Ds, Dmp);
        ok = qx_divexact(H, t1, Dm);
    }
    if (ok && n == 0) {                    /* D squarefree: nothing to reduce */
        fmpq_poly_zero(gn); fmpq_poly_one(gd);
        fmpq_poly_set(hn, A); fmpq_poly_set(hd, D);
    } else if (ok) {
        fmpq_mat_t M, X, Bv;
        fmpq_mat_init(M, N, N); fmpq_mat_init(X, N, 1); fmpq_mat_init(Bv, N, 1);
        fmpq_t c;
        fmpq_init(c);
        for (slong i = 0; i < n; i++) {                  /* b_i: (x^i)' Ds - x^i H */
            fmpq_poly_zero(col);
            if (i > 0) {
                fmpq_poly_zero(t1);
                fmpq_poly_set_coeff_si(t1, i - 1, i);
                fmpq_poly_mul(col, t1, Ds);
            }
            fmpq_poly_zero(t1);
            fmpq_poly_set_coeff_si(t1, i, 1);
            fmpq_poly_mul(t1, t1, H);
            fmpq_poly_sub(col, col, t1);
            for (slong r = 0; r < N; r++) {
                fmpq_poly_get_coeff_fmpq(c, col, r);
                fmpq_set(fmpq_mat_entry(M, r, i), c);
            }
        }
        for (slong j = 0; j < m; j++) {                  /* c_j: x^j Dm */
            fmpq_poly_zero(t1);
            fmpq_poly_set_coeff_si(t1, j, 1);
            fmpq_poly_mul(col, t1, Dm);
            for (slong r = 0; r < N; r++) {
                fmpq_poly_get_coeff_fmpq(c, col, r);
                fmpq_set(fmpq_mat_entry(M, r, n + j), c);
            }
        }
        for (slong r = 0; r < N; r++) {
            fmpq_poly_get_coeff_fmpq(c, A, r);
            fmpq_set(fmpq_mat_entry(Bv, r, 0), c);
        }
        ok = fmpq_mat_solve(X, M, Bv) != 0;
        if (ok) {
            fmpq_poly_zero(gn);
            for (slong i = 0; i < n; i++) fmpq_poly_set_coeff_fmpq(gn, i, fmpq_mat_entry(X, i, 0));
            fmpq_poly_set(gd, Dm);
            fmpq_poly_zero(hn);
            for (slong j = 0; j < m; j++) fmpq_poly_set_coeff_fmpq(hn, j, fmpq_mat_entry(X, n + j, 0));
            fmpq_poly_set(hd, Ds);
        }
        fmpq_clear(c);
        fmpq_mat_clear(M); fmpq_mat_clear(X); fmpq_mat_clear(Bv);
    }
    fmpq_poly_clear(Dp); fmpq_poly_clear(Dm); fmpq_poly_clear(Ds); fmpq_poly_clear(Dmp);
    fmpq_poly_clear(H); fmpq_poly_clear(col); fmpq_poly_clear(t1);
    return ok;
}

Expr* intrat_alg_hermite(const Expr* Ae, const Expr* De, const Expr* x, int variant) {
    if (!Ae || !De || !x || x->type != EXPR_SYMBOL) return NULL;
    fmpq_poly_t A, D, Q, R, gn, gd, hn, hd;
    fmpq_poly_init(A); fmpq_poly_init(D); fmpq_poly_init(Q); fmpq_poly_init(R);
    fmpq_poly_init(gn); fmpq_poly_init(gd); fmpq_poly_init(hn); fmpq_poly_init(hd);
    Expr* out = NULL;
    if (qx_from_expr(Ae, x, A) && qx_from_expr(De, x, D) && fmpq_poly_degree(D) >= 1) {
        /* the polynomial part rides along in h */
        fmpq_poly_divrem(Q, R, A, D);
        bool ok;
        switch (variant) {
        case INTRAT_HERMITE_QUADRATIC: ok = hermite_quadratic(gn, gd, hn, hd, R, D); break;
        case INTRAT_HERMITE_HOROWITZ:  ok = hermite_horowitz(gn, gd, hn, hd, R, D); break;
        default:                       ok = hermite_mack(gn, gd, hn, hd, R, D); break;
        }
        if (ok) {
            fmpq_poly_mul(Q, Q, hd);
            fmpq_poly_add(hn, hn, Q);
            Expr* g = qx_fraction_to_expr(gn, gd, x);
            Expr* h = qx_fraction_to_expr(hn, hd, x);
            out = expr_new_function(expr_new_symbol(SYM_List), (Expr*[]){ g, h }, 2);
        }
    }
    fmpq_poly_clear(A); fmpq_poly_clear(D); fmpq_poly_clear(Q); fmpq_poly_clear(R);
    fmpq_poly_clear(gn); fmpq_poly_clear(gd); fmpq_poly_clear(hn); fmpq_poly_clear(hd);
    return out;
}

/* ================================================================== */
/* Log part: Rothstein-Trager and Czichowski                           */
/* ================================================================== */

/* Integral A, D with the same ratio (both scaled by the lcm of all
 * denominators). */
static void ia_integral_pair(fmpz_poly_t Az, fmpz_poly_t Dz, const fmpq_poly_t A, const fmpq_poly_t D) {
    fmpz_t L;
    fmpz_init(L);
    fmpz_lcm(L, fmpq_poly_denref(A), fmpq_poly_denref(D));
    fmpq_poly_t t1;
    fmpq_poly_init(t1);
    fmpq_poly_scalar_mul_fmpz(t1, A, L);
    fmpq_poly_get_numerator(Az, t1);
    fmpq_poly_scalar_mul_fmpz(t1, D, L);
    fmpq_poly_get_numerator(Dz, t1);
    fmpq_poly_clear(t1);
    fmpz_clear(L);
}

/* R(t) = Res_x(D, A - t D') by evaluation at deg D + 1 integer points and
 * interpolation (R is in Z[t] for integral A, D).  The points skip the one
 * value of t at which A - t D' drops degree, so every value is the formal
 * Sylvester resultant. */
static void ia_rt_resultant(fmpz_poly_t R, const fmpz_poly_t Az, const fmpz_poly_t Dz) {
    slong n = fmpz_poly_degree(Dz);
    fmpz_poly_t Dp, P;
    fmpz_poly_init(Dp); fmpz_poly_init(P);
    fmpz_poly_derivative(Dp, Dz);
    /* t* = a_{n-1} / (n d_n) */
    fmpq_t tstar;
    fmpq_init(tstar);
    {
        fmpz_t a, b;
        fmpz_init(a); fmpz_init(b);
        fmpz_poly_get_coeff_fmpz(a, Az, n - 1);
        fmpz_poly_get_coeff_fmpz(b, Dp, n - 1);
        fmpq_set_fmpz_frac(tstar, a, b);
        fmpz_clear(a); fmpz_clear(b);
    }
    fmpz* xs = _fmpz_vec_init(n + 1);
    fmpz* ys = _fmpz_vec_init(n + 1);
    slong k = 0, v = 0;
    fmpz_t tk;
    fmpz_init(tk);
    while (k <= n) {
        fmpz_set_si(tk, v);
        v = v >= 0 ? -v - 1 : -v;            /* 0, -1, 1, -2, 2, ... */
        if (fmpz_is_one(fmpq_denref(tstar)) && fmpz_equal(tk, fmpq_numref(tstar))) continue;
        fmpz_poly_scalar_mul_fmpz(P, Dp, tk);
        fmpz_poly_sub(P, Az, P);
        fmpz_set(xs + k, tk);
        fmpz_poly_resultant(ys + k, Dz, P);
        k++;
    }
    fmpz_poly_interpolate_fmpz_vec(R, xs, ys, n + 1);
    fmpz_clear(tk);
    _fmpz_vec_clear(xs, n + 1); _fmpz_vec_clear(ys, n + 1);
    fmpq_clear(tstar);
    fmpz_poly_clear(Dp); fmpz_poly_clear(P);
}

/* In (Q[t]/q)[x]: p := p mod q coefficient-wise. */
static void qtx_reduce_mod(QTX* p, const fmpq_poly_t q) {
    for (slong k = 0; k < p->len; k++) fmpq_poly_rem(p->c + k, p->c + k, q);
    qtx_normalise(p);
}

/* a := a mod b in (Q[t]/q)[x], q irreducible (a field); false if lc(b) is
 * not invertible (cannot happen for irreducible q and b != 0 mod q). */
static bool qtx_rem_mod(QTX* a, const QTX* b, const fmpq_poly_t q) {
    fmpq_poly_t inv, G, T, c, t1;
    fmpq_poly_init(inv); fmpq_poly_init(G); fmpq_poly_init(T); fmpq_poly_init(c); fmpq_poly_init(t1);
    fmpq_poly_xgcd(G, inv, T, b->c + b->len - 1, q);
    bool ok = fmpq_poly_degree(G) == 0;
    while (ok && a->len >= b->len) {
        slong sh = a->len - b->len;
        fmpq_poly_mul(c, a->c + a->len - 1, inv);
        fmpq_poly_rem(c, c, q);
        for (slong k = 0; k < b->len; k++) {
            fmpq_poly_mul(t1, c, b->c + k);
            fmpq_poly_sub(a->c + k + sh, a->c + k + sh, t1);
            fmpq_poly_rem(a->c + k + sh, a->c + k + sh, q);
        }
        fmpq_poly_zero(a->c + a->len - 1);
        qtx_normalise(a);
    }
    fmpq_poly_clear(inv); fmpq_poly_clear(G); fmpq_poly_clear(T); fmpq_poly_clear(c); fmpq_poly_clear(t1);
    return ok;
}

/* gcd(D, A - t D') over Q(alpha) = Q[t]/q, monic in x, as a QTX. */
static bool ia_gcd_mod(QTX* out, const fmpq_poly_t A, const fmpq_poly_t D, const fmpq_poly_t q) {
    fmpq_poly_t Dp;
    fmpq_poly_init(Dp);
    fmpq_poly_derivative(Dp, D);
    QTX a, b;
    slong n = fmpq_poly_length(D);
    qtx_init(&a, n); qtx_init(&b, n);
    fmpq_t c;
    fmpq_init(c);
    for (slong k = 0; k < n; k++) {
        fmpq_poly_get_coeff_fmpq(c, D, k);
        fmpq_poly_set_fmpq(a.c + k, c);
        /* (A - t D')_k = a_k - d'_k t */
        fmpq_poly_get_coeff_fmpq(c, A, k);
        fmpq_poly_set_fmpq(b.c + k, c);
        fmpq_poly_get_coeff_fmpq(c, Dp, k);
        fmpq_neg(c, c);
        fmpq_poly_set_coeff_fmpq(b.c + k, 1, c);
    }
    a.len = b.len = n;
    qtx_reduce_mod(&a, q);
    qtx_reduce_mod(&b, q);
    bool ok = true;
    while (ok && b.len > 0) {
        ok = qtx_rem_mod(&a, &b, q);
        QTX s = a; a = b; b = s;
    }
    if (ok && a.len > 0) {
        /* monic in x */
        fmpq_poly_t inv, G, T;
        fmpq_poly_init(inv); fmpq_poly_init(G); fmpq_poly_init(T);
        fmpq_poly_xgcd(G, inv, T, a.c + a.len - 1, q);
        for (slong k = 0; k < a.len; k++) {
            fmpq_poly_mul(a.c + k, a.c + k, inv);
            fmpq_poly_rem(a.c + k, a.c + k, q);
        }
        fmpq_poly_clear(inv); fmpq_poly_clear(G); fmpq_poly_clear(T);
        *out = a;
        qtx_clear(&b);
    } else {
        ok = false;
        qtx_clear(&a); qtx_clear(&b);
    }
    fmpq_clear(c);
    fmpq_poly_clear(Dp);
    return ok;
}

static Expr* ia_pair(const fmpq_poly_t Q, const QTX* S, const Expr* x, const Expr* t) {
    Expr* qe;
    {
        QTX tmp;
        qtx_init(&tmp, 1);
        fmpq_poly_set(tmp.c, Q);
        tmp.len = 1;
        /* Q as a polynomial in t: store it as the x^0 coefficient */
        qe = qtx_to_expr(&tmp, x, t);
        qtx_clear(&tmp);
    }
    Expr* se = qtx_to_expr(S, x, t);
    return expr_new_function(expr_new_symbol(SYM_List), (Expr*[]){ qe, se }, 2);
}

/* Rothstein-Trager (Bronstein §2.4): R = Res_x(D, A - t D'), factored over
 * Q; for each irreducible factor q, S = gcd(D, A - a D') over Q(a). */
static Expr* logpart_rt(const fmpq_poly_t A, const fmpq_poly_t D, const Expr* x, const Expr* t) {
    fmpz_poly_t Az, Dz, R;
    fmpz_poly_init(Az); fmpz_poly_init(Dz); fmpz_poly_init(R);
    ia_integral_pair(Az, Dz, A, D);
    ia_rt_resultant(R, Az, Dz);
    Expr* out = NULL;
    if (!fmpz_poly_is_zero(R)) {
        fmpz_poly_factor_t f;
        fmpz_poly_factor_init(f);
        fmpz_poly_factor(f, R);
        Expr** pairs = (Expr**)malloc(sizeof(Expr*) * (size_t)(f->num ? f->num : 1));
        size_t np = 0;
        bool ok = true;
        fmpq_poly_t q;
        fmpq_poly_init(q);
        for (slong i = 0; i < f->num && ok; i++)
            if (fmpz_poly_degree(f->p + i) > IA_MAX_RT_DEG) ok = false;
        for (slong i = 0; i < f->num && ok; i++) {
            if (fmpz_poly_degree(f->p + i) < 1) continue;
            fmpq_poly_set_fmpz_poly(q, f->p + i);
            QTX S;
            ok = ia_gcd_mod(&S, A, D, q);
            if (!ok) break;
            pairs[np++] = ia_pair(q, &S, x, t);
            qtx_clear(&S);
        }
        if (ok) out = expr_new_function(expr_new_symbol(SYM_List), pairs, np);
        else for (size_t k = 0; k < np; k++) expr_free(pairs[k]);
        free(pairs);
        fmpq_poly_clear(q);
        fmpz_poly_factor_clear(f);
    }
    fmpz_poly_clear(Az); fmpz_poly_clear(Dz); fmpz_poly_clear(R);
    return out;
}

/* A GBPoly in (x, t) -> QTX by x-degree.  Variable 0 is x, 1 is t. */
static void ia_gb_to_qtx(const GBPoly* p, QTX* out) {
    slong maxi = -1;
    for (size_t k = 0; k < p->n_terms; k++)
        if (p->exps[k * 2] > maxi) maxi = p->exps[k * 2];
    qtx_init(out, maxi + 1);
    fmpq_t c, cur;
    fmpq_init(c); fmpq_init(cur);
    for (size_t k = 0; k < p->n_terms; k++) {
        fmpq_set_mpq(c, p->coefs[k]);
        slong i = p->exps[k * 2], j = p->exps[k * 2 + 1];
        fmpq_poly_get_coeff_fmpq(cur, out->c + i, j);
        fmpq_add(cur, cur, c);
        fmpq_poly_set_coeff_fmpq(out->c + i, j, cur);
    }
    out->len = maxi + 1;
    qtx_normalise(out);
    fmpq_clear(c); fmpq_clear(cur);
}

/* content_x (monic gcd of the x-coefficients) of p. */
static void qtx_content(fmpq_poly_t g, const QTX* p) {
    fmpq_poly_zero(g);
    for (slong k = 0; k < p->len; k++) fmpq_poly_gcd(g, g, p->c + k);
}

static GBPoly* ia_gb_from(const fmpq_poly_t P, const fmpq_poly_t Pt) {
    /* P(x) - t Pt(x) as a GBPoly in (x, t), lex x > t */
    GBPoly* g = gb_poly_new(2, GB_ORDER_LEX, 0);
    mpq_t c;
    mpq_init(c);
    fmpq_t f;
    fmpq_init(f);
    int ex[2];
    for (slong k = 0; k < fmpq_poly_length(P); k++) {
        fmpq_poly_get_coeff_fmpq(f, P, k);
        if (fmpq_is_zero(f)) continue;
        fmpq_get_mpq(c, f);
        ex[0] = (int)k; ex[1] = 0;
        gb_poly_push_term(g, ex, c);
    }
    if (Pt) {
        for (slong k = 0; k < fmpq_poly_length(Pt); k++) {
            fmpq_poly_get_coeff_fmpq(f, Pt, k);
            if (fmpq_is_zero(f)) continue;
            fmpq_neg(f, f);
            fmpq_get_mpq(c, f);
            ex[0] = (int)k; ex[1] = 1;
            gb_poly_push_term(g, ex, c);
        }
    }
    gb_poly_normalize(g);
    fmpq_clear(f);
    mpq_clear(c);
    return g;
}

/* Czichowski (Bronstein §2.6, Theorem 2.6.1): the reduced lex Groebner
 * basis (P_1, ..., P_m) of (D, A - t D') with x > t, sorted by increasing
 * leading term; Q_i = content_x(P_i)/content_x(P_{i+1}), S_i = pp_x(P_{i+1}). */
static Expr* logpart_czichowski(const fmpq_poly_t A, const fmpq_poly_t D, const Expr* x, const Expr* t) {
    fmpq_poly_t Dp;
    fmpq_poly_init(Dp);
    fmpq_poly_derivative(Dp, D);
    GBPoly* F[2];
    F[0] = ia_gb_from(D, NULL);
    F[1] = ia_gb_from(A, Dp);
    GBPoly** G = NULL;
    size_t ng = 0;
    Expr* out = NULL;
    if (gb_flint_buchberger((GBPoly* const*)F, 2, &G, &ng) && ng >= 2) {
        QTX* P = (QTX*)malloc(sizeof(QTX) * ng);
        for (size_t i = 0; i < ng; i++) ia_gb_to_qtx(G[i], P + i);
        /* the basis comes sorted ascending by leading monomial; P_1 must be
         * free of x (the resultant's squarefree part) */
        bool ok = P[0].len == 1;
        for (size_t i = 1; i < ng && ok; i++) ok = P[i].len > P[i - 1].len;
        if (ok) {
            Expr** pairs = (Expr**)malloc(sizeof(Expr*) * ng);
            size_t np = 0;
            fmpq_poly_t ci, cn, Q;
            fmpq_poly_init(ci); fmpq_poly_init(cn); fmpq_poly_init(Q);
            qtx_content(ci, P + 0);
            for (size_t i = 0; i + 1 < ng && ok; i++) {
                qtx_content(cn, P + i + 1);
                ok = qx_divexact(Q, ci, cn);
                if (!ok) break;
                if (fmpq_poly_degree(Q) >= 1) {
                    QTX S;
                    qtx_init(&S, P[i + 1].len);
                    for (slong k = 0; k < P[i + 1].len; k++)
                        qx_divexact(S.c + k, P[i + 1].c + k, cn);
                    S.len = P[i + 1].len;
                    fmpq_poly_make_monic(Q, Q);
                    pairs[np++] = ia_pair(Q, &S, x, t);
                    qtx_clear(&S);
                }
                fmpq_poly_swap(ci, cn);
            }
            if (ok) out = expr_new_function(expr_new_symbol(SYM_List), pairs, np);
            else for (size_t k = 0; k < np; k++) expr_free(pairs[k]);
            free(pairs);
            fmpq_poly_clear(ci); fmpq_poly_clear(cn); fmpq_poly_clear(Q);
        }
        for (size_t i = 0; i < ng; i++) qtx_clear(P + i);
        free(P);
    }
    if (G) {
        for (size_t i = 0; i < ng; i++) gb_poly_free(G[i]);
        free(G);
    }
    gb_poly_free(F[0]); gb_poly_free(F[1]);
    fmpq_poly_clear(Dp);
    return out;
}

Expr* intrat_alg_logpart(const Expr* Ae, const Expr* De, const Expr* x, const Expr* t, int method) {
    if (!Ae || !De || !x || !t || x->type != EXPR_SYMBOL || t->type != EXPR_SYMBOL) return NULL;
    fmpq_poly_t A, D, g;
    fmpq_poly_init(A); fmpq_poly_init(D); fmpq_poly_init(g);
    Expr* out = NULL;
    if (qx_from_expr(Ae, x, A) && qx_from_expr(De, x, D)
        && fmpq_poly_degree(D) >= 1 && fmpq_poly_degree(A) < fmpq_poly_degree(D)
        && !fmpq_poly_is_zero(A)) {
        /* preconditions: D squarefree, gcd(A, D) = 1 */
        fmpq_poly_t Dp;
        fmpq_poly_init(Dp);
        fmpq_poly_derivative(Dp, D);
        fmpq_poly_gcd(g, D, Dp);
        bool sqfree = fmpq_poly_degree(g) == 0;
        fmpq_poly_gcd(g, A, D);
        fmpq_poly_clear(Dp);
        if (sqfree && fmpq_poly_degree(g) == 0) {
            if (method == INTRAT_LOGPART_RT) out = logpart_rt(A, D, x, t);
            else if (method == INTRAT_LOGPART_CZICHOWSKI) out = logpart_czichowski(A, D, x, t);
        }
    }
    fmpq_poly_clear(A); fmpq_poly_clear(D); fmpq_poly_clear(g);
    return out;
}

/* ================================================================== */
/* LogToReal over the real algebraic numbers (qqbar)                   */
/* ================================================================== */

/* Dense polynomial in x with qqbar coefficients. */
typedef struct { slong len; slong cap; qqbar_struct* c; } QBP;

static void qbp_init(QBP* p, slong cap) {
    if (cap < 1) cap = 1;
    p->len = 0; p->cap = cap;
    p->c = (qqbar_struct*)flint_malloc(sizeof(qqbar_struct) * (size_t)cap);
    for (slong k = 0; k < cap; k++) qqbar_init(p->c + k);
}
static void qbp_clear(QBP* p) {
    for (slong k = 0; k < p->cap; k++) qqbar_clear(p->c + k);
    flint_free(p->c);
    p->c = NULL; p->len = p->cap = 0;
}
static void qbp_fit(QBP* p, slong len) {
    if (len <= p->cap) return;
    p->c = (qqbar_struct*)flint_realloc(p->c, sizeof(qqbar_struct) * (size_t)len);
    for (slong k = p->cap; k < len; k++) qqbar_init(p->c + k);
    p->cap = len;
}
static void qbp_norm(QBP* p) {
    while (p->len > 0 && qqbar_is_zero(p->c + p->len - 1)) p->len--;
}
static void qbp_set(QBP* r, const QBP* a) {
    qbp_fit(r, a->len);
    for (slong k = 0; k < a->len; k++) qqbar_set(r->c + k, a->c + k);
    r->len = a->len;
}
static void qbp_neg(QBP* r, const QBP* a) {
    qbp_set(r, a);
    for (slong k = 0; k < r->len; k++) qqbar_neg(r->c + k, r->c + k);
}
/* r = a + s * b (s = +1 or -1) */
static void qbp_addsub(QBP* r, const QBP* a, const QBP* b, int s) {
    slong n = a->len > b->len ? a->len : b->len;
    QBP t;
    qbp_init(&t, n);
    for (slong k = 0; k < n; k++) {
        if (k < a->len) qqbar_set(t.c + k, a->c + k); else qqbar_zero(t.c + k);
        if (k < b->len) {
            if (s > 0) qqbar_add(t.c + k, t.c + k, b->c + k);
            else       qqbar_sub(t.c + k, t.c + k, b->c + k);
        }
    }
    t.len = n;
    qbp_norm(&t);
    qbp_clear(r);
    *r = t;
}
static void qbp_mul(QBP* r, const QBP* a, const QBP* b) {
    QBP t;
    if (a->len == 0 || b->len == 0) { qbp_init(&t, 1); qbp_clear(r); *r = t; return; }
    slong n = a->len + b->len - 1;
    qbp_init(&t, n);
    qqbar_t p;
    qqbar_init(p);
    for (slong k = 0; k < n; k++) qqbar_zero(t.c + k);
    for (slong i = 0; i < a->len; i++) {
        if (qqbar_is_zero(a->c + i)) continue;
        for (slong j = 0; j < b->len; j++) {
            qqbar_mul(p, a->c + i, b->c + j);
            qqbar_add(t.c + i + j, t.c + i + j, p);
        }
    }
    qqbar_clear(p);
    t.len = n;
    qbp_norm(&t);
    qbp_clear(r);
    *r = t;
}
/* q, r with a = q b + r, deg r < deg b (b != 0). */
static void qbp_divrem(QBP* q, QBP* r, const QBP* a, const QBP* b) {
    QBP rr, qq;
    qbp_init(&rr, a->len); qbp_set(&rr, a);
    slong ql = a->len - b->len + 1;
    qbp_init(&qq, ql > 0 ? ql : 1);
    for (slong k = 0; k < qq.cap; k++) qqbar_zero(qq.c + k);
    qq.len = ql > 0 ? ql : 0;
    qqbar_t c, p;
    qqbar_init(c); qqbar_init(p);
    while (rr.len >= b->len && rr.len > 0) {
        slong sh = rr.len - b->len;
        qqbar_div(c, rr.c + rr.len - 1, b->c + b->len - 1);
        qqbar_set(qq.c + sh, c);
        for (slong k = 0; k < b->len; k++) {
            qqbar_mul(p, c, b->c + k);
            qqbar_sub(rr.c + k + sh, rr.c + k + sh, p);
        }
        qqbar_zero(rr.c + rr.len - 1);
        qbp_norm(&rr);
    }
    qqbar_clear(c); qqbar_clear(p);
    qbp_norm(&qq);
    if (q) { qbp_clear(q); *q = qq; } else qbp_clear(&qq);
    if (r) { qbp_clear(r); *r = rr; } else qbp_clear(&rr);
}

/* g = s a + t b (extended Euclid). */
static void qbp_xgcd(QBP* g, QBP* s, QBP* t, const QBP* a, const QBP* b) {
    QBP r0, r1, s0, s1, t0, t1, q, tmp, tmp2;
    qbp_init(&r0, 1); qbp_set(&r0, a);
    qbp_init(&r1, 1); qbp_set(&r1, b);
    qbp_init(&s0, 1); qqbar_one(s0.c); s0.len = 1;
    qbp_init(&s1, 1);
    qbp_init(&t0, 1);
    qbp_init(&t1, 1); qqbar_one(t1.c); t1.len = 1;
    qbp_init(&q, 1); qbp_init(&tmp, 1); qbp_init(&tmp2, 1);
    while (r1.len > 0) {
        qbp_divrem(&q, &tmp, &r0, &r1);
        qbp_clear(&r0); r0 = r1; r1 = tmp; qbp_init(&tmp, 1);
        qbp_mul(&tmp2, &q, &s1);
        qbp_addsub(&tmp2, &s0, &tmp2, -1);
        qbp_clear(&s0); s0 = s1; s1 = tmp2; qbp_init(&tmp2, 1);
        qbp_mul(&tmp2, &q, &t1);
        qbp_addsub(&tmp2, &t0, &tmp2, -1);
        qbp_clear(&t0); t0 = t1; t1 = tmp2; qbp_init(&tmp2, 1);
    }
    qbp_clear(g); *g = r0;
    qbp_clear(s); *s = s0;
    qbp_clear(t); *t = t0;
    qbp_clear(&r1); qbp_clear(&s1); qbp_clear(&t1);
    qbp_clear(&q); qbp_clear(&tmp); qbp_clear(&tmp2);
}

/* An algebraic number as an Expr: a rational, (a + b Sqrt[c])/q, or (when
 * allowed) Root[minpoly, k] with k found by comparing enclosures with the
 * evaluator's Root objects.  NULL to decline. */
static Expr* ia_qqbar_expr(const qqbar_t a, bool allow_root) {
    if (qqbar_is_rational(a)) {
        fmpq_t q;
        fmpq_init(q);
        qqbar_get_fmpq(q, a);
        Expr* r = ia_expr_fmpq(q);
        fmpq_clear(q);
        return eval_and_free(r);
    }
    if (qqbar_degree(a) == 2) {
        fmpz_t A, B, C, Q;
        fmpz_init(A); fmpz_init(B); fmpz_init(C); fmpz_init(Q);
        qqbar_get_quadratic(A, B, C, Q, a, 0);       /* (A + B Sqrt[C]) / Q */
        Expr* sq = expr_new_function(expr_new_symbol(SYM_Sqrt), (Expr*[]){ ia_expr_fmpz(C) }, 1);
        Expr* num = expr_new_function(expr_new_symbol(SYM_Plus), (Expr*[]){ ia_expr_fmpz(A),
            expr_new_function(expr_new_symbol(SYM_Times), (Expr*[]){ ia_expr_fmpz(B), sq }, 2) }, 2);
        Expr* r = expr_new_function(expr_new_symbol(SYM_Times), (Expr*[]){ num,
            expr_new_function(expr_new_symbol(SYM_Power), (Expr*[]){ ia_expr_fmpz(Q), expr_new_integer(-1) }, 2) }, 2);
        fmpz_clear(A); fmpz_clear(B); fmpz_clear(C); fmpz_clear(Q);
        return eval_and_free(r);
    }
    if (!allow_root || !qqbar_is_real(a)) return NULL;
    /* Root[p, k]: the real roots of p are indexed in increasing order */
    slong d = qqbar_degree(a);
    qqbar_ptr rs = _qqbar_vec_init(d);
    qqbar_roots_fmpz_poly(rs, QQBAR_POLY(a), QQBAR_ROOTS_IRREDUCIBLE);
    slong below = 0;
    for (slong i = 0; i < d; i++)
        if (qqbar_is_real(rs + i) && qqbar_cmp_re(rs + i, a) < 0) below++;
    _qqbar_vec_clear(rs, d);
    QTX tmp;
    qtx_init(&tmp, 1);
    fmpq_poly_set_fmpz_poly(tmp.c, QQBAR_POLY(a));
    tmp.len = 1;
    Expr* slot = expr_new_function(expr_new_symbol(SYM_Slot), (Expr*[]){ expr_new_integer(1) }, 1);
    Expr* body = qtx_to_expr(&tmp, slot, slot);   /* univariate in the slot */
    expr_free(slot);
    qtx_clear(&tmp);
    Expr* fn = expr_new_function(expr_new_symbol(SYM_Function), (Expr*[]){ body }, 1);
    return eval_and_free(expr_new_function(expr_new_symbol(SYM_Root),
        (Expr*[]){ fn, expr_new_integer((int64_t)(below + 1)) }, 2));
}

static Expr* qbp_to_expr(const QBP* p, const Expr* x, bool allow_root) {
    Expr** terms = (Expr**)malloc(sizeof(Expr*) * (size_t)(p->len ? p->len : 1));
    size_t m = 0;
    for (slong k = 0; k < p->len; k++) {
        if (qqbar_is_zero(p->c + k)) continue;
        Expr* c = ia_qqbar_expr(p->c + k, allow_root);
        if (!c) {
            for (size_t i = 0; i < m; i++) expr_free(terms[i]);
            free(terms);
            return NULL;
        }
        terms[m++] = k == 0 ? c
            : expr_new_function(expr_new_symbol(SYM_Times), (Expr*[]){ c, ia_pow(x, k) }, 2);
    }
    return eval_and_free(ia_sum(terms, m));
}

/* An ArcTan argument with its common denominator pulled out: for rational
 * and quadratic-surd coefficients (a + b Sqrt[c])/q the lcm L of the q's is
 * factored, x/2 - 3x^3/2 + x^5/2 -> (x - 3x^3 + x^5)/2, the reference's
 * spelling; otherwise as qbp_to_expr. */
static Expr* qbp_to_atan_arg(const QBP* p, const Expr* x, bool allow_root) {
    fmpz_t L, A, B, C, Q;
    fmpz_init(L); fmpz_init(A); fmpz_init(B); fmpz_init(C); fmpz_init(Q);
    fmpz_one(L);
    bool small = true;
    for (slong k = 0; k < p->len && small; k++) {
        if (qqbar_is_zero(p->c + k)) continue;
        slong d = qqbar_degree(p->c + k);
        if (d > 2) { small = false; break; }
        if (d == 1) {
            fmpq_t q;
            fmpq_init(q);
            qqbar_get_fmpq(q, p->c + k);
            fmpz_lcm(L, L, fmpq_denref(q));
            fmpq_clear(q);
        } else {
            qqbar_get_quadratic(A, B, C, Q, p->c + k, 0);
            fmpz_lcm(L, L, Q);
        }
    }
    Expr* out;
    if (!small || fmpz_is_one(L)) {
        out = qbp_to_expr(p, x, allow_root);
    } else {
        QBP s;
        qbp_init(&s, p->len);
        for (slong k = 0; k < p->len; k++) qqbar_mul_fmpz(s.c + k, p->c + k, L);
        s.len = p->len;
        Expr* num = qbp_to_expr(&s, x, allow_root);
        qbp_clear(&s);
        out = num ? eval_and_free(expr_new_function(expr_new_symbol(SYM_Times), (Expr*[]){ num,
            expr_new_function(expr_new_symbol(SYM_Power), (Expr*[]){ ia_expr_fmpz(L), expr_new_integer(-1) }, 2) }, 2))
            : NULL;
    }
    fmpz_clear(L); fmpz_clear(A); fmpz_clear(B); fmpz_clear(C); fmpz_clear(Q);
    return out;
}

static bool qbp_degree_ok(const QBP* p) {
    for (slong k = 0; k < p->len; k++)
        if (qqbar_degree(p->c + k) > IA_MAX_QQBAR_DEG) return false;
    return true;
}

/* Rioboo's LogToAtan (Bronstein §2.8): a sum of 2 ArcTan[polynomial] whose
 * derivative is that of I Log[(A + I B)/(A - I B)].  Appends to *acc. */
static bool ia_log_to_atan(const QBP* A, const QBP* B, const Expr* x, bool allow_root,
                           Expr*** acc, size_t* n, size_t* cap, int depth) {
    if (B->len == 0) return true;
    if (depth > 64 || !qbp_degree_ok(A) || !qbp_degree_ok(B)) return false;
    QBP q, r;
    qbp_init(&q, 1); qbp_init(&r, 1);
    qbp_divrem(&q, &r, A, B);
    bool ok = true;
    if (r.len == 0) {                              /* B | A: 2 ArcTan[A/B] */
        Expr* pe = qbp_to_atan_arg(&q, x, allow_root);
        ok = pe != NULL;
        if (ok) {
            if (*n == *cap) { *cap = *cap * 2 + 4; *acc = (Expr**)realloc(*acc, sizeof(Expr*) * *cap); }
            (*acc)[(*n)++] = expr_new_function(expr_new_symbol(SYM_Times), (Expr*[]){ expr_new_integer(2),
                expr_new_function(expr_new_symbol(SYM_ArcTan), (Expr*[]){ pe }, 1) }, 2);
        }
    } else if (A->len < B->len) {                  /* LogToAtan(-B, A) */
        QBP nb;
        qbp_init(&nb, 1);
        qbp_neg(&nb, B);
        ok = ia_log_to_atan(&nb, A, x, allow_root, acc, n, cap, depth + 1);
        qbp_clear(&nb);
    } else {
        /* (D, C, G) <- ExtendedEuclidean(B, -A):  B D - A C = G */
        QBP na, G, Dd, C, t1, t2, P, rem;
        qbp_init(&na, 1); qbp_init(&G, 1); qbp_init(&Dd, 1); qbp_init(&C, 1);
        qbp_init(&t1, 1); qbp_init(&t2, 1); qbp_init(&P, 1); qbp_init(&rem, 1);
        qbp_neg(&na, A);
        qbp_xgcd(&G, &Dd, &C, B, &na);
        qbp_mul(&t1, A, &Dd);
        qbp_mul(&t2, B, &C);
        qbp_addsub(&t1, &t1, &t2, +1);            /* A D + B C */
        qbp_divrem(&P, &rem, &t1, &G);
        ok = rem.len == 0;
        if (ok) {
            Expr* pe = qbp_to_atan_arg(&P, x, allow_root);
            ok = pe != NULL;
            if (ok) {
                if (*n == *cap) { *cap = *cap * 2 + 4; *acc = (Expr**)realloc(*acc, sizeof(Expr*) * *cap); }
                (*acc)[(*n)++] = expr_new_function(expr_new_symbol(SYM_Times), (Expr*[]){ expr_new_integer(2),
                    expr_new_function(expr_new_symbol(SYM_ArcTan), (Expr*[]){ pe }, 1) }, 2);
                ok = ia_log_to_atan(&Dd, &C, x, allow_root, acc, n, cap, depth + 1);
            }
        }
        qbp_clear(&na); qbp_clear(&G); qbp_clear(&Dd); qbp_clear(&C);
        qbp_clear(&t1); qbp_clear(&t2); qbp_clear(&P); qbp_clear(&rem);
    }
    qbp_clear(&q); qbp_clear(&r);
    return ok;
}

static void ia_push(Expr*** acc, size_t* n, size_t* cap, Expr* e) {
    if (*n == *cap) { *cap = *cap * 2 + 4; *acc = (Expr**)realloc(*acc, sizeof(Expr*) * *cap); }
    (*acc)[(*n)++] = e;
}

/* One pair (Q, S) over Q: the real form of Sum_{Q(a)=0} a Log[S(a, x)]. */
static bool ia_ltr_pair(const fmpq_poly_t Q, const QTX* S, const Expr* x, bool allow_root,
                        Expr*** acc, size_t* n, size_t* cap) {
    slong d = fmpq_poly_degree(Q);
    if (d < 1 || d > (allow_root ? IA_MAX_LTR_DEG : IA_MAX_LTR_DEG_RAD)) return false;
    fmpz_poly_t Qz;
    fmpz_poly_init(Qz);
    fmpq_poly_get_numerator(Qz, Q);
    /* the cost of qqbar arithmetic (Re and Im of the roots, then LogToAtan
     * over their field) grows with the coefficient size too, and cannot be
     * interrupted: decline before it starts */
    if (FLINT_ABS(fmpz_poly_max_bits(Qz)) > IA_MAX_LTR_BITS) { fmpz_poly_clear(Qz); return false; }
    qqbar_ptr roots = _qqbar_vec_init(d);
    qqbar_roots_fmpz_poly(roots, Qz, 0);
    bool ok = true;
    qqbar_t u, v, sv;
    qqbar_init(u); qqbar_init(v); qqbar_init(sv);
    for (slong i = 0; i < d && ok; i++) {
        const qqbar_struct* a = roots + i;
        if (qqbar_is_real(a)) {
            /* a Log[S(a, x)] */
            QBP s;
            qbp_init(&s, S->len);
            for (slong k = 0; k < S->len; k++) qqbar_evaluate_fmpq_poly(s.c + k, S->c + k, a);
            s.len = S->len;
            qbp_norm(&s);
            ok = qbp_degree_ok(&s);
            Expr* ae = ok ? ia_qqbar_expr(a, allow_root) : NULL;
            Expr* se = ae ? qbp_to_expr(&s, x, allow_root) : NULL;
            ok = ae && se;
            if (ok) ia_push(acc, n, cap, expr_new_function(expr_new_symbol(SYM_Times), (Expr*[]){ ae,
                    expr_new_function(expr_new_symbol(SYM_Log), (Expr*[]){ se }, 1) }, 2));
            else { if (ae) expr_free(ae); if (se) expr_free(se); }
            qbp_clear(&s);
            continue;
        }
        if (qqbar_sgn_im(a) < 0) continue;         /* the conjugate's partner */
        qqbar_re(u, a);
        qqbar_im(v, a);
        QBP Ap, Bp, A2, B2;
        qbp_init(&Ap, S->len); qbp_init(&Bp, S->len); qbp_init(&A2, 1); qbp_init(&B2, 1);
        for (slong k = 0; k < S->len; k++) {
            qqbar_evaluate_fmpq_poly(sv, S->c + k, a);
            qqbar_re(Ap.c + k, sv);
            qqbar_im(Bp.c + k, sv);
        }
        Ap.len = Bp.len = S->len;
        qbp_norm(&Ap); qbp_norm(&Bp);
        ok = qbp_degree_ok(&Ap) && qbp_degree_ok(&Bp);
        if (ok && !qqbar_is_zero(u)) {
            /* u Log[A^2 + B^2] */
            qbp_mul(&A2, &Ap, &Ap);
            qbp_mul(&B2, &Bp, &Bp);
            qbp_addsub(&A2, &A2, &B2, +1);
            Expr* ue = ia_qqbar_expr(u, allow_root);
            Expr* me = ue ? qbp_to_expr(&A2, x, allow_root) : NULL;
            ok = ue && me;
            if (ok) ia_push(acc, n, cap, expr_new_function(expr_new_symbol(SYM_Times), (Expr*[]){ ue,
                    expr_new_function(expr_new_symbol(SYM_Log), (Expr*[]){ me }, 1) }, 2));
            else { if (ue) expr_free(ue); if (me) expr_free(me); }
        }
        if (ok) {
            /* v LogToAtan(A, B) */
            Expr** at = NULL;
            size_t na = 0, ca = 0;
            ok = ia_log_to_atan(&Ap, &Bp, x, allow_root, &at, &na, &ca, 0);
            Expr* ve = ok ? ia_qqbar_expr(v, allow_root) : NULL;
            ok = ok && ve;
            if (ok) {
                Expr* sum = na == 0 ? expr_new_integer(0)
                          : na == 1 ? at[0]
                          : expr_new_function(expr_new_symbol(SYM_Plus), at, na);
                if (na > 1) { /* adopted */ }
                ia_push(acc, n, cap, expr_new_function(expr_new_symbol(SYM_Times), (Expr*[]){ ve, sum }, 2));
            } else {
                if (ve) expr_free(ve);
                for (size_t k = 0; k < na; k++) expr_free(at[k]);
            }
            free(at);
        }
        qbp_clear(&Ap); qbp_clear(&Bp); qbp_clear(&A2); qbp_clear(&B2);
    }
    qqbar_clear(u); qqbar_clear(v); qqbar_clear(sv);
    _qqbar_vec_clear(roots, d);
    fmpz_poly_clear(Qz);
    return ok;
}

Expr* intrat_alg_log_to_real(const Expr* pairs, const Expr* x, const Expr* t, bool allow_root) {
    if (!pairs || !ia_head_is(pairs, SYM_List) || !x || !t) return NULL;
    Expr** acc = NULL;
    size_t n = 0, cap = 0;
    bool ok = pairs->data.function.arg_count > 0;
    for (size_t i = 0; i < pairs->data.function.arg_count && ok; i++) {
        const Expr* pr = pairs->data.function.args[i];
        if (!ia_head_is(pr, SYM_List) || pr->data.function.arg_count != 2) { ok = false; break; }
        QTX Qt, S;
        /* Q(t): parse as a polynomial in t (t in the x slot) */
        if (!qtx_from_expr(pr->data.function.args[0], t, NULL, &Qt)) { ok = false; break; }
        if (!qtx_from_expr(pr->data.function.args[1], x, t, &S)) { qtx_clear(&Qt); ok = false; break; }
        fmpq_poly_t Q;
        fmpq_poly_init(Q);
        fmpq_t c;
        fmpq_init(c);
        for (slong k = 0; k < Qt.len; k++) {
            fmpq_poly_get_coeff_fmpq(c, Qt.c + k, 0);
            fmpq_poly_set_coeff_fmpq(Q, k, c);
        }
        fmpq_clear(c);
        ok = ia_ltr_pair(Q, &S, x, allow_root, &acc, &n, &cap);
        fmpq_poly_clear(Q);
        qtx_clear(&Qt); qtx_clear(&S);
    }
    if (!ok) {
        for (size_t k = 0; k < n; k++) expr_free(acc[k]);
        free(acc);
        return NULL;
    }
    Expr* sum = n == 0 ? expr_new_integer(0) : n == 1 ? acc[0]
              : expr_new_function(expr_new_symbol(SYM_Plus), acc, n);
    if (n > 1) { /* adopted */ }
    free(acc);
    return eval_and_free(sum);
}

#else  /* !USE_FLINT */

Expr* intrat_alg_hermite(const Expr* A, const Expr* D, const Expr* x, int variant) {
    (void)A; (void)D; (void)x; (void)variant; return NULL;
}
Expr* intrat_alg_logpart(const Expr* A, const Expr* D, const Expr* x, const Expr* t, int method) {
    (void)A; (void)D; (void)x; (void)t; (void)method; return NULL;
}
Expr* intrat_alg_log_to_real(const Expr* pairs, const Expr* x, const Expr* t, bool allow_root) {
    (void)pairs; (void)x; (void)t; (void)allow_root; return NULL;
}

#endif
