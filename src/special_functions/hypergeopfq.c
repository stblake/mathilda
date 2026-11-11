/*
 * hypergeopfq.c -- HypergeometricPFQ for Mathilda.
 *
 * HypergeometricPFQ[{a1,...,ap}, {b1,...,bq}, z] is the generalized
 * hypergeometric function pFq(a;b;z) with series
 *
 *     pFq(a;b;z) = Sum_{k>=0}  ( prod_i (a_i)_k / prod_j (b_j)_k ) z^k / k!
 *
 * where (x)_k is the Pochhammer (rising-factorial) symbol.  The evaluation
 * order in builtin_hypergeometric_pfq is:
 *
 *   1. validate (arity 3; first two args are Lists);
 *   2. thread over a List third argument;
 *   3. cancel parameters common to the upper and lower lists (generic case);
 *   4. z == 0 -> 1;
 *   5. terminate to an explicit polynomial when an upper parameter is a
 *      non-positive integer (exact, symbolic-z capable);
 *   6. numeric evaluation by direct series summation (machine / MPFR /
 *      complex) in the convergent regime;
 *   7. otherwise return NULL (stay unevaluated).
 *
 * The convergence regime is: p<=q (entire, all finite z); p==q+1 (|z|<1).
 * For p==q+1 with |z|>=1 and for p>q+1 (asymptotic only) the series is not
 * summed directly -- the call stays unevaluated.  Analytic continuation is a
 * deliberate follow-on (see docs/spec).
 *
 * Memory contract: builtin takes ownership of res but must not free it.
 */

#include "hypergeopfq.h"
#include "ndarray.h"
#include "sf_machine.h"    /* NDArray buffer path for the elementwise argument */
#include "sf_machine_hyper.h" /* complex twins of the machine kernels */
#include "sf_scalar.h"     /* sf_machine_scalar_args, sf_machine_leaf */
#include "symtab.h"
#include "attr.h"
#include "eval.h"
#include "expr.h"
#include "sym_names.h"
#include "arithmetic.h"   /* is_rational, is_complex, make_complex */
#include "numeric.h"      /* numeric_combined_bits, get_approx_mpfr, ... */
#include "flint_num_bridge.h" /* flint_num_hyp2f1: the arbitrary-precision 2F1 */
#include "parse.h"        /* parse_expression: closed-form templates */
#include "simp_sf_engine.h" /* sf_contains_head */
#ifdef USE_MPFR
#include "numeric_complex.h"
#include <mpfr.h>
#endif
#include <complex.h>
#include <math.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>

/* Hard caps on the number of series terms summed before giving up. */
#define HGPFQ_MACHINE_MAX_TERMS 200000
#define HGPFQ_MPFR_MAX_TERMS    1000000

/* ------------------------------------------------------------------ */
/*  Small structural helpers                                           */
/* ------------------------------------------------------------------ */

static bool is_list(const Expr* e) {
    return e && e->type == EXPR_FUNCTION
        && e->data.function.head->type == EXPR_SYMBOL
        && e->data.function.head->data.symbol.name == SYM_List;
}

/* True if e is a non-positive integer (..., -2, -1, 0); fills *neg with
 * its magnitude n>=0 (param == -n) when true. */
static bool is_nonpos_int(const Expr* e, long* neg) {
    if (e && e->type == EXPR_INTEGER && e->data.integer <= 0) {
        if (neg) *neg = (long)(-e->data.integer);
        return true;
    }
    return false;
}

/* Build HypergeometricPFQ[a, b, z] from owned arguments and evaluate. */
static Expr* rebuild_eval(Expr* a, Expr* b, Expr* z) {
    Expr* args[3] = { a, b, z };
    Expr* call = expr_new_function(expr_new_symbol(SYM_HypergeometricPFQ), args, 3);
    Expr* r = evaluate(call);
    expr_free(call);
    return r;
}

/* Pochhammer product (x)_k = x (x+1) ... (x+k-1) as a symbolic Times,
 * or Integer 1 when k == 0.  `x` is copied. */
static Expr* poch_product(Expr* x, long k) {
    if (k <= 0) return expr_new_integer(1);
    Expr** fac = malloc(sizeof(Expr*) * (size_t)k);
    for (long i = 0; i < k; i++) {
        if (i == 0) {
            fac[i] = expr_copy(x);
        } else {
            Expr* pa[2] = { expr_copy(x), expr_new_integer(i) };
            fac[i] = expr_new_function(expr_new_symbol(SYM_Plus), pa, 2);
        }
    }
    Expr* t = expr_new_function(expr_new_symbol(SYM_Times), fac, (size_t)k);
    free(fac);
    return t;
}

/* ------------------------------------------------------------------ */
/*  Parameter cancellation (generic): drop one element common to the   */
/*  upper and lower lists, unless it is a non-positive integer (which   */
/*  makes the series terminate and must be kept).                       */
/* ------------------------------------------------------------------ */

static Expr* try_cancel(Expr* a, Expr* b, Expr* z) {
    size_t p = a->data.function.arg_count;
    size_t q = b->data.function.arg_count;
    for (size_t i = 0; i < p; i++) {
        Expr* ai = a->data.function.args[i];
        long dummy;
        if (is_nonpos_int(ai, &dummy)) continue;   /* keep terminating params */
        for (size_t j = 0; j < q; j++) {
            if (!expr_eq(ai, b->data.function.args[j])) continue;
            /* Build new lists with element i (upper) / j (lower) removed. */
            Expr** na = malloc(sizeof(Expr*) * (p - 1));
            Expr** nb = malloc(sizeof(Expr*) * (q - 1));
            size_t w = 0;
            for (size_t k = 0; k < p; k++)
                if (k != i) na[w++] = expr_copy(a->data.function.args[k]);
            w = 0;
            for (size_t k = 0; k < q; k++)
                if (k != j) nb[w++] = expr_copy(b->data.function.args[k]);
            Expr* nal = expr_new_function(expr_new_symbol(SYM_List), na, p - 1);
            Expr* nbl = expr_new_function(expr_new_symbol(SYM_List), nb, q - 1);
            free(na); free(nb);
            return rebuild_eval(nal, nbl, expr_copy(z));
        }
    }
    return NULL;
}

/* ------------------------------------------------------------------ */
/*  Native builders for the closed-form stages below.  Each hg_* call  */
/*  builder takes ownership of its Expr arguments.                     */
/* ------------------------------------------------------------------ */

static Expr* hg_fn(const char* head, Expr** args, size_t n) {
    return expr_new_function(expr_new_symbol(head), args, n);
}
static Expr* hg_plus(Expr* x, Expr* y)  { return hg_fn(SYM_Plus,  (Expr*[]){ x, y }, 2); }
static Expr* hg_times(Expr* x, Expr* y) { return hg_fn(SYM_Times, (Expr*[]){ x, y }, 2); }
static Expr* hg_neg(Expr* x)            { return hg_times(expr_new_integer(-1), x); }
static Expr* hg_half(Expr* x)           { return hg_times(make_rational(1, 2), x); }
static Expr* hg_inv(Expr* x) {
    return hg_fn(SYM_Power, (Expr*[]){ x, expr_new_integer(-1) }, 2);
}
static Expr* hg_gamma(Expr* x) { return hg_fn(SYM_Gamma, (Expr*[]){ x }, 1); }
#define HG_CP(x) expr_copy((Expr*)(x))

/* Evaluate e and free it. */
static Expr* hg_eval(Expr* e) {
    Expr* r = evaluate(e);
    expr_free(e);
    return r;
}

/* u - v, expanded: parameters arrive in either distributed (1 + a/2) or
 * factored (1/2 (2 + a)) form, and Plus alone does not distribute. */
static Expr* hg_diff(const Expr* u, const Expr* v) {
    return hg_eval(hg_fn(SYM_Expand, (Expr*[]){ hg_plus(HG_CP(u), hg_neg(HG_CP(v))) }, 1));
}

/* True iff u - v is exactly 0.  Classical parameter relations are linear in
 * the parameters (1 + a - b, a/2, ...), so the expanded difference decides
 * them; no Simplify is needed. */
static bool hg_eq(const Expr* u, const Expr* v) {
    if (expr_eq(u, v)) return true;
    Expr* d = hg_diff(u, v);
    bool z = d && d->type == EXPR_INTEGER && d->data.integer == 0;
    if (d) expr_free(d);
    return z;
}

/* The sign of Re(e): 1 if e is a numeric quantity with Re(e) > 0, 0 if it is
 * numeric with Re(e) <= 0, -1 if e is symbolic (undecided: the theorem is
 * then the formal value, as for symbolic parameters). */
static int hg_re_positive(const Expr* e) {
    if (!expr_is_numeric_quantity(e)) return -1;
#ifndef USE_MPFR
    /* No-MPFR build: decide on the machine value's real part. */
    Expr* mv = numericalize_machine_unplanned(e);
    const Expr* r = mv;
    if (r && r->type == EXPR_FUNCTION && r->data.function.arg_count == 2
        && r->data.function.head->type == EXPR_SYMBOL
        && r->data.function.head->data.symbol.name == SYM_Complex)
        r = r->data.function.args[0];
    int ms = -1;
    if (r && r->type == EXPR_REAL)    ms = r->data.real > 0.0 ? 1 : 0;
    if (r && r->type == EXPR_INTEGER) ms = r->data.integer > 0 ? 1 : 0;
    if (mv) expr_free(mv);
    return ms;
#else
    mpfr_t re, im;
    mpfr_init2(re, 128); mpfr_init2(im, 128);
    bool inexact = false;
    Expr* nv = NULL;
    bool ok = get_approx_mpfr(e, re, im, &inexact);
    if (!ok) {
        nv = numericalize_machine_unplanned(e);
        ok = nv && get_approx_mpfr(nv, re, im, &inexact);
    }
    int s = ok ? (mpfr_sgn(re) > 0 ? 1 : 0) : -1;
    mpfr_clear(re); mpfr_clear(im);
    if (nv) expr_free(nv);
    return s;
#endif
}

/* Gate a closed form on the convergence condition Re(s) > 0, s owned. */
static bool hg_converges(Expr* s) {
    Expr* v = hg_eval(s);
    bool ok = v && hg_re_positive(v) != 0;
    if (v) expr_free(v);
    return ok;
}

/* L with the element at index skip removed, as a new List. */
static Expr* hg_delete(const Expr* l, size_t skip) {
    size_t n = l->data.function.arg_count;
    Expr** v = malloc(sizeof(Expr*) * (n ? n : 1));
    size_t w = 0;
    for (size_t k = 0; k < n; k++)
        if (k != skip) v[w++] = expr_copy(l->data.function.args[k]);
    Expr* r = hg_fn(SYM_List, v, w);
    free(v);
    return r;
}

/* ------------------------------------------------------------------ */
/*  Integer-shifted parameter pair: an upper a_i = b_j + n (n a        */
/*  positive integer) against a lower b_j.                             */
/*                                                                     */
/*  (b+n)_k/(b)_k = P(k) = prod_{m<n} (k+b+m)/(b+m) is a polynomial    */
/*  in k of degree n, and k z^k = theta z^k with theta = z d/dz, so    */
/*      pFq[U + {b+n}, L + {b}, z] = P(theta) pFq[U, L, z].            */
/*  Applied only when the reduced series has a closed form (it is      */
/*  evaluated at a fresh symbol t first): 2F2[{2,2},{1,1},z] ->        */
/*  e^z (1 + 3z + z^2), the Dobinski/Touchard sums Sum k^m z^k/k!,     */
/*  2F1[{a, c+1},{c},z] -> (1-z)^(-a-1) (1 - z + a z/c).  A pair whose */
/*  lower b is a non-positive integer is a pole, never reduced.        */
/*                                                                     */
/*  P's coefficients are the elementary symmetric functions of the     */
/*  b+m, built by the recurrence e_d <- e_{d-1} + (b+m) e_d.           */
/* ------------------------------------------------------------------ */
#define HG_SHIFT_MAX 12

static Expr* hg_shift_apply(const Expr* a, const Expr* b, const Expr* z,
                            size_t i, size_t j, int64_t n) {
    Expr* t = expr_new_symbol("Hyp`t");
    Expr* F = rebuild_eval(hg_delete(a, i), hg_delete(b, j), HG_CP(t));
    if (!F || sf_contains_head(F, SYM_HypergeometricPFQ)) {
        if (F) expr_free(F);
        expr_free(t);
        return NULL;
    }
    const Expr* bj = b->data.function.args[j];
    Expr* e[HG_SHIFT_MAX + 1];            /* e[d] = coefficient of k^d */
    e[0] = expr_new_integer(1);
    for (int64_t d = 1; d <= n; d++) e[d] = expr_new_integer(0);
    Expr* norm = expr_new_integer(1);     /* prod (b+m) */
    for (int64_t m = 0; m < n; m++) {
        Expr* cm = hg_eval(hg_plus(HG_CP(bj), expr_new_integer(m)));
        for (int64_t d = m + 1; d >= 1; d--)
            e[d] = hg_eval(hg_plus(e[d - 1] ? HG_CP(e[d - 1]) : expr_new_integer(0),
                                   hg_times(HG_CP(cm), e[d])));
        e[0] = hg_eval(hg_times(HG_CP(cm), e[0]));
        norm = hg_eval(hg_times(norm, cm));
    }
    /* sum_d e[d] theta^d F / norm */
    Expr** terms = malloc(sizeof(Expr*) * (size_t)(n + 1));
    Expr* G = F;
    for (int64_t d = 0; d <= n; d++) {
        terms[d] = hg_times(e[d], HG_CP(G));
        if (d < n) {
            Expr* dG = hg_eval(hg_fn(SYM_D, (Expr*[]){ HG_CP(G), HG_CP(t) }, 2));
            expr_free(G);
            G = hg_eval(hg_times(HG_CP(t), dG));
        }
    }
    expr_free(G);
    Expr* sum = hg_times(hg_fn(SYM_Plus, terms, (size_t)(n + 1)), hg_inv(norm));
    free(terms);
    Expr* fac = hg_eval(hg_fn(SYM_Factor, (Expr*[]){ sum }, 1));
    Expr* rule = hg_fn(SYM_Rule, (Expr*[]){ t, HG_CP(z) }, 2);
    return hg_eval(hg_fn(SYM_ReplaceAll, (Expr*[]){ fac, rule }, 2));
}

static Expr* try_shift_pair(Expr* a, Expr* b, Expr* z) {
    size_t p = a->data.function.arg_count, q = b->data.function.arg_count;
    for (size_t i = 0; i < p; i++) {
        for (size_t j = 0; j < q; j++) {
            long dummy;
            if (is_nonpos_int(b->data.function.args[j], &dummy)) continue;
            Expr* d = hg_diff(a->data.function.args[i], b->data.function.args[j]);
            bool ok = d && d->type == EXPR_INTEGER && d->data.integer >= 1
                   && d->data.integer <= HG_SHIFT_MAX;
            int64_t n = ok ? d->data.integer : 0;
            if (d) expr_free(d);
            if (!ok) continue;
            Expr* r = hg_shift_apply(a, b, z, i, j, n);
            if (r) return r;
        }
    }
    return NULL;
}

/* ------------------------------------------------------------------ */
/*  Classical summation theorems at z = 1 and z = -1.                  */
/*                                                                     */
/*  Gauss     2F1(a, b; c; 1)            = G(c)G(c-a-b)/(G(c-a)G(c-b))  */
/*            (Re(c-a-b) > 0)                                           */
/*  Kummer    2F1(a, b; 1+a-b; -1)       = G(1+a-b)G(1+a/2)/            */
/*            (Re b < 1)                   (G(1+a)G(1+a/2-b))           */
/*  Pfaff-Saalschutz (terminating, balanced)                            */
/*            3F2(a, b, -n; c, 1+a+b-c-n; 1) = (c-a)_n (c-b)_n/         */
/*                                         ((c)_n (c-a-b)_n)            */
/*  Dixon     3F2(a, b, c; 1+a-b, 1+a-c; 1)   (Re(1+a/2-b-c) > 0)       */
/*  Dougall   5F4(a, 1+a/2, b, c, d; a/2, 1+a-b, 1+a-c, 1+a-d; 1)       */
/*            (very-well-poised; Ramanujan's quarter series;            */
/*             Re(1+a-b-c-d) > 0)                                       */
/*                                                                     */
/*  Parameters are matched up to permutation by exact comparison of     */
/*  the linear relations (hg_eq): an upper parameter is chosen for each */
/*  role and the lower list is checked against the relations, so the    */
/*  work is O(p^2 q) differences, not O(p! q!) Simplify calls.  A       */
/*  numeric parameter set outside the region of convergence is          */
/*  declined; symbolic parameters give the formal value.                */
/* ------------------------------------------------------------------ */

/* 1 + a - x, unevaluated. */
static Expr* hg_1pa_m(const Expr* a, const Expr* x) {
    return hg_fn(SYM_Plus, (Expr*[]){ expr_new_integer(1), HG_CP(a), hg_neg(HG_CP(x)) }, 3);
}
/* Gamma-ratio  prod Gamma(num) / prod Gamma(den), all owned, evaluated. */
static Expr* hg_gamma_ratio(Expr** num, size_t nn, Expr** den, size_t nd) {
    Expr* f[8];
    size_t w = 0;
    for (size_t k = 0; k < nn; k++) f[w++] = hg_gamma(num[k]);
    for (size_t k = 0; k < nd; k++) f[w++] = hg_inv(hg_gamma(den[k]));
    return hg_eval(hg_fn(SYM_Times, f, w));
}

static Expr* hg_gauss(const Expr* a, const Expr* b, const Expr* c) {
    Expr* s = hg_fn(SYM_Plus, (Expr*[]){ HG_CP(c), hg_neg(HG_CP(a)), hg_neg(HG_CP(b)) }, 3);
    if (!hg_converges(expr_copy(s))) { expr_free(s); return NULL; }
    return hg_gamma_ratio((Expr*[]){ HG_CP(c), s }, 2,
        (Expr*[]){ hg_plus(HG_CP(c), hg_neg(HG_CP(a))),
                   hg_plus(HG_CP(c), hg_neg(HG_CP(b))) }, 2);
}

static Expr* hg_kummer(const Expr* A, const Expr* c) {
    for (size_t i = 0; i < 2; i++) {
        const Expr* a = A->data.function.args[i];
        const Expr* b = A->data.function.args[1 - i];
        Expr* rel = hg_1pa_m(a, b);
        bool m = hg_eq(c, rel);
        expr_free(rel);
        if (!m) continue;
        if (!hg_converges(hg_plus(expr_new_integer(1), hg_neg(HG_CP(b))))) continue;
        Expr* ha = hg_plus(expr_new_integer(1), hg_half(HG_CP(a)));
        return hg_gamma_ratio(
            (Expr*[]){ hg_1pa_m(a, b), expr_copy(ha) }, 2,
            (Expr*[]){ hg_plus(expr_new_integer(1), HG_CP(a)),
                       hg_plus(ha, hg_neg(HG_CP(b))) }, 2);
    }
    return NULL;
}

/* -u is a non-negative integer or a bare symbol: the series terminates at
 * the (possibly symbolic) degree n = -u. */
static bool hg_terminating(const Expr* u) {
    if (u->type == EXPR_INTEGER) return u->data.integer <= 0;
    return u->type == EXPR_FUNCTION && u->data.function.head->type == EXPR_SYMBOL
        && u->data.function.head->data.symbol.name == SYM_Times
        && u->data.function.arg_count == 2
        && u->data.function.args[0]->type == EXPR_INTEGER
        && u->data.function.args[0]->data.integer == -1
        && u->data.function.args[1]->type == EXPR_SYMBOL;
}

static Expr* hg_poch(Expr* x, Expr* n) { return hg_fn(SYM_Pochhammer, (Expr*[]){ x, n }, 2); }

static Expr* hg_saalschutz(const Expr* A, const Expr* B) {
    for (size_t k = 0; k < 3; k++) {
        const Expr* u = A->data.function.args[k];
        if (!hg_terminating(u)) continue;
        const Expr* a = A->data.function.args[(k + 1) % 3];
        const Expr* b = A->data.function.args[(k + 2) % 3];
        for (size_t j = 0; j < 2; j++) {
            const Expr* c = B->data.function.args[j];
            Expr* rel = hg_fn(SYM_Plus, (Expr*[]){ expr_new_integer(1), HG_CP(a), HG_CP(b),
                                                   hg_neg(HG_CP(c)), HG_CP(u) }, 5);
            bool m = hg_eq(B->data.function.args[1 - j], rel);
            expr_free(rel);
            if (!m) continue;
            Expr* n = hg_eval(hg_neg(HG_CP(u)));
            Expr* ca = hg_plus(HG_CP(c), hg_neg(HG_CP(a)));
            Expr* cb = hg_plus(HG_CP(c), hg_neg(HG_CP(b)));
            Expr* cab = hg_fn(SYM_Plus, (Expr*[]){ HG_CP(c), hg_neg(HG_CP(a)), hg_neg(HG_CP(b)) }, 3);
            Expr* r = hg_eval(hg_fn(SYM_Times, (Expr*[]){
                hg_poch(ca, HG_CP(n)), hg_poch(cb, HG_CP(n)),
                hg_inv(hg_poch(HG_CP(c), HG_CP(n))), hg_inv(hg_poch(cab, HG_CP(n))) }, 4));
            expr_free(n);
            return r;
        }
    }
    return NULL;
}

static Expr* hg_dixon(const Expr* A, const Expr* B) {
    for (size_t k = 0; k < 3; k++) {
        const Expr* a = A->data.function.args[k];
        for (size_t s = 0; s < 2; s++) {
            const Expr* b = A->data.function.args[(k + 1 + s) % 3];
            const Expr* c = A->data.function.args[(k + 2 - s) % 3];
            Expr* r1 = hg_1pa_m(a, b), *r2 = hg_1pa_m(a, c);
            bool m = hg_eq(B->data.function.args[0], r1) && hg_eq(B->data.function.args[1], r2);
            expr_free(r1); expr_free(r2);
            if (!m) continue;
            Expr* ha = hg_plus(expr_new_integer(1), hg_half(HG_CP(a)));
            Expr* hbc = hg_fn(SYM_Plus, (Expr*[]){ expr_copy(ha), hg_neg(HG_CP(b)), hg_neg(HG_CP(c)) }, 3);
            if (!hg_converges(expr_copy(hbc))) { expr_free(ha); expr_free(hbc); return NULL; }
            return hg_gamma_ratio(
                (Expr*[]){ expr_copy(ha), hg_1pa_m(a, b), hg_1pa_m(a, c), hbc }, 4,
                (Expr*[]){ hg_plus(expr_new_integer(1), HG_CP(a)),
                           hg_plus(expr_copy(ha), hg_neg(HG_CP(b))),
                           hg_plus(ha, hg_neg(HG_CP(c))),
                           hg_fn(SYM_Plus, (Expr*[]){ expr_new_integer(1), HG_CP(a),
                                 hg_neg(HG_CP(b)), hg_neg(HG_CP(c)) }, 4) }, 4);
        }
    }
    return NULL;
}

static Expr* hg_dougall(const Expr* A, const Expr* B) {
    for (size_t ia = 0; ia < 5; ia++) {
        const Expr* a = A->data.function.args[ia];
        if (a->type == EXPR_INTEGER && a->data.integer == 0) continue;
        Expr* ha = hg_plus(expr_new_integer(1), hg_half(HG_CP(a)));
        Expr* a2 = hg_half(HG_CP(a));
        for (size_t ih = 0; ih < 5; ih++) {
            if (ih == ia || !hg_eq(A->data.function.args[ih], ha)) continue;
            size_t ub[3], w = 0;
            for (size_t k = 0; k < 5; k++) if (k != ia && k != ih) ub[w++] = k;
            for (size_t jl = 0; jl < 4; jl++) {
                if (!hg_eq(B->data.function.args[jl], a2)) continue;
                bool used[3] = { false, false, false }, all = true;
                for (size_t l = 0; l < 4 && all; l++) {
                    if (l == jl) continue;
                    bool found = false;
                    for (size_t u = 0; u < 3 && !found; u++) {
                        if (used[u]) continue;
                        Expr* rel = hg_1pa_m(a, A->data.function.args[ub[u]]);
                        if (hg_eq(B->data.function.args[l], rel)) used[u] = found = true;
                        expr_free(rel);
                    }
                    all = found;
                }
                if (!all) continue;
                const Expr* b = A->data.function.args[ub[0]];
                const Expr* c = A->data.function.args[ub[1]];
                const Expr* d = A->data.function.args[ub[2]];
                expr_free(ha); expr_free(a2);
                Expr* s = hg_fn(SYM_Plus, (Expr*[]){ expr_new_integer(1), HG_CP(a),
                          hg_neg(HG_CP(b)), hg_neg(HG_CP(c)), hg_neg(HG_CP(d)) }, 5);
                if (!hg_converges(expr_copy(s))) { expr_free(s); return NULL; }
                #define HG_1A(x, y) hg_fn(SYM_Plus, (Expr*[]){ expr_new_integer(1), HG_CP(a), \
                                        hg_neg(HG_CP(x)), hg_neg(HG_CP(y)) }, 4)
                Expr* r = hg_gamma_ratio(
                    (Expr*[]){ hg_1pa_m(a, b), hg_1pa_m(a, c), hg_1pa_m(a, d), s }, 4,
                    (Expr*[]){ hg_plus(expr_new_integer(1), HG_CP(a)),
                               HG_1A(c, d), HG_1A(b, d), HG_1A(b, c) }, 4);
                #undef HG_1A
                return r;
            }
        }
        expr_free(ha); expr_free(a2);
    }
    return NULL;
}

static Expr* try_classical_sum(Expr* a, Expr* b, Expr* z) {
    bool z1 = z->type == EXPR_INTEGER && z->data.integer == 1;
    bool zm1 = z->type == EXPR_INTEGER && z->data.integer == -1;
    if (!z1 && !zm1) return NULL;
    size_t p = a->data.function.arg_count, q = b->data.function.arg_count;
    const Expr* const* A = (const Expr* const*)a->data.function.args;
    if (p == 2 && q == 1)
        return z1 ? hg_gauss(A[0], A[1], b->data.function.args[0])
                  : hg_kummer(a, b->data.function.args[0]);
    if (!z1) return NULL;
    if (p == 3 && q == 2) {
        Expr* r = hg_saalschutz(a, b);
        return r ? r : hg_dixon(a, b);
    }
    if (p == 5 && q == 4) return hg_dougall(a, b);
    return NULL;
}

/* ------------------------------------------------------------------ */
/*  Termination to a polynomial                                        */
/* ------------------------------------------------------------------ */

/* If any upper parameter is a non-positive integer -n, the series truncates
 * at k = n (the smallest such n).  Returns the explicit Sum_{k=0}^{n} c_k z^k
 * evaluated, or NULL when termination does not apply or a lower non-positive
 * integer would force a 0 denominator inside the truncation range. */
static Expr* try_terminate(Expr* a, Expr* b, Expr* z) {
    size_t p = a->data.function.arg_count;
    size_t q = b->data.function.arg_count;

    long n = -1;
    for (size_t i = 0; i < p; i++) {
        long m;
        if (is_nonpos_int(a->data.function.args[i], &m)) {
            if (n < 0 || m < n) n = m;
        }
    }
    if (n < 0) return NULL;   /* no terminating upper parameter */

    /* A lower non-positive integer -m gives (b_j)_k == 0 for k > m; if m < n
     * the denominator vanishes inside the range -> undefined.  Bail. */
    for (size_t j = 0; j < q; j++) {
        long m;
        if (is_nonpos_int(b->data.function.args[j], &m) && m < n) return NULL;
    }

    Expr** terms = malloc(sizeof(Expr*) * (size_t)(n + 1));
    for (long k = 0; k <= n; k++) {
        /* numerator: prod_i (a_i)_k */
        size_t nfac = p + 2;                 /* uppers + denom-power + z-power */
        Expr** mul = malloc(sizeof(Expr*) * nfac);
        size_t w = 0;
        for (size_t i = 0; i < p; i++)
            mul[w++] = poch_product(a->data.function.args[i], k);

        /* denominator: prod_j (b_j)_k * k!  -> Power[Times[...], -1] */
        size_t nden = q + 1;
        Expr** den = malloc(sizeof(Expr*) * nden);
        size_t dw = 0;
        for (size_t j = 0; j < q; j++)
            den[dw++] = poch_product(b->data.function.args[j], k);
        den[dw++] = expr_new_function(expr_new_symbol(SYM_Factorial),
                        (Expr*[]){ expr_new_integer(k) }, 1);
        Expr* denom = expr_new_function(expr_new_symbol(SYM_Times), den, nden);
        free(den);
        mul[w++] = expr_new_function(expr_new_symbol(SYM_Power),
                        (Expr*[]){ denom, expr_new_integer(-1) }, 2);

        /* z^k */
        mul[w++] = expr_new_function(expr_new_symbol(SYM_Power),
                        (Expr*[]){ expr_copy(z), expr_new_integer(k) }, 2);

        terms[k] = expr_new_function(expr_new_symbol(SYM_Times), mul, w);
        free(mul);
    }
    Expr* sum = expr_new_function(expr_new_symbol(SYM_Plus), terms, (size_t)(n + 1));
    free(terms);
    Expr* r = evaluate(sum);
    expr_free(sum);
    return r;
}

/* ------------------------------------------------------------------ */
/*  Symbolic reductions to elementary functions                        */
/* ------------------------------------------------------------------ */

static bool is_int_val(const Expr* e, int64_t v) {
    return e && e->type == EXPR_INTEGER && e->data.integer == v;
}
/* True if e == odd/2 (half-odd-integer); fills *twice with the numerator. */
static bool is_half_int(const Expr* e, int64_t* twice) {
    int64_t n, d;
    if (is_rational((Expr*)e, &n, &d) && d == 2) { if (twice) *twice = n; return true; }
    return false;
}

/* 1 - z as Plus[1, Times[-1, z]]. */
static Expr* one_minus(Expr* z) {
    Expr* neg = expr_new_function(expr_new_symbol(SYM_Times),
                    (Expr*[]){ expr_new_integer(-1), expr_copy(z) }, 2);
    return expr_new_function(expr_new_symbol(SYM_Plus),
                    (Expr*[]){ expr_new_integer(1), neg }, 2);
}
/* Sqrt[z] as Power[z, 1/2]. */
static Expr* sqrt_of(Expr* z) {
    return expr_new_function(expr_new_symbol(SYM_Power),
               (Expr*[]){ expr_copy(z), make_rational(1, 2) }, 2);
}

/* Recognise a handful of closed forms that need no new function heads.
 * Returns the evaluated elementary expression, or NULL if no rule matches.
 * Valid as analytic continuations for all z, so safe to run before the
 * numeric series. */
/* Does e equal the (reduced) rational num/den? */
static bool hg_is_rat(Expr* e, int64_t num, int64_t den) {
    Expr* r = make_rational(num, den);
    bool ok = expr_eq(e, r);
    expr_free(r);
    return ok;
}

/* `fn` -- the source of a pure function of one argument, e.g.
 * "ArcTan[Sqrt[#]]/Sqrt[#] &" -- applied to a copy of z, UNevaluated (the
 * caller of try_reduce evaluates its result). Templates keep the closed forms
 * below readable; they are parsed only when a shape has already matched. */
static Expr* hg_apply_tmpl(const char* fn, const Expr* z) {
    Expr* f = parse_expression(fn);
    if (!f) return NULL;
    return expr_new_function(f, (Expr*[]){ expr_copy((Expr*)z) }, 1);
}

/* Is z syntactically a negated quantity (a negative real number, or a
 * product with a negative numeric coefficient)? Picks between the two
 * equally valid spellings of an even-in-Sqrt[z] closed form -- Erf[Sqrt[-z]]
 * versus Erfi[Sqrt[z]], ArcTan versus ArcTanh -- so the printed form is the
 * real-looking one: 1F1(1/2;3/2;-x) -> Erf[Sqrt[x]], not Erfi[Sqrt[-x]]. */
static bool hg_negated(const Expr* z) {
    if (z->type == EXPR_INTEGER) return z->data.integer < 0;
    if (z->type == EXPR_REAL) return z->data.real < 0;
    int64_t n, d;
    if (is_rational((Expr*)z, &n, &d)) return n < 0;
    if (z->type == EXPR_FUNCTION && z->data.function.head->type == EXPR_SYMBOL &&
        z->data.function.head->data.symbol.name == SYM_Times &&
        z->data.function.arg_count >= 1)
        return hg_negated(z->data.function.args[0]);
    return false;
}

/* l - u evaluated equals the integer `v` exactly. */
static bool hg_diff_is(const Expr* l, const Expr* u, int64_t v) {
    Expr* neg = expr_new_function(expr_new_symbol(SYM_Times),
                    (Expr*[]){ expr_new_integer(-1), expr_copy((Expr*)u) }, 2);
    Expr* d = eval_and_free(expr_new_function(expr_new_symbol(SYM_Plus),
                    (Expr*[]){ expr_copy((Expr*)l), neg }, 2));
    bool ok = is_int_val(d, v);
    expr_free(d);
    return ok;
}

/* 1F1(a; a+1; z) for a half-integer or a positive integer, a = n/d.
 *
 * F_a(z) = 1F1(a; a+1; z) = a z^-a Int_0^z t^(a-1) e^t dt, and integrating
 * by parts gives the two-term recurrence
 *     F_(a+1) = (a+1) (E^z - F_a) / z,      F_(a-1) = E^z - z F_a / a,
 * so every half-integer a reaches F_(1/2) = Sqrt[Pi] Erfi[Sqrt[z]]/(2 Sqrt[z])
 * (written Erf[Sqrt[-z]] when z is a negated quantity) and every positive
 * integer reaches F_1 = (E^z - 1)/z. These match the reference system's
 * forms (1F1(-1/2;1/2;z) = E^z - Sqrt[Pi] Sqrt[z] Erfi[Sqrt[z]]).
 *
 * Other a are declined: the general closed form a (-z)^-a Gamma[a, 0, -z]
 * is correct but Gamma[a, w] has no numerics on its cut w < 0, so it would
 * turn a summable series at real z > 0 into an expression N cannot evaluate.
 * Evaluated and owned, or NULL. */
#define HG_1F1_MAX_STEPS 24
static Expr* hg_1f1_a_ap1(const Expr* a, const Expr* z) {
    int64_t n, d;
    if (!is_rational((Expr*)a, &n, &d)) return NULL;
    Expr* F;
    Expr* cur;                      /* the parameter F currently belongs to */
    int64_t steps;
    if (d == 2) {
        F = hg_negated(z)
            ? hg_apply_tmpl("(Sqrt[Pi] Erf[Sqrt[-#]]/(2 Sqrt[-#]))&", z)
            : hg_apply_tmpl("(Sqrt[Pi] Erfi[Sqrt[#]]/(2 Sqrt[#]))&", z);
        cur = make_rational(1, 2);
        steps = (n - 1) / 2;        /* n odd: a = 1/2 + steps */
    } else if (d == 1 && n >= 1) {
        F = hg_apply_tmpl("((E^# - 1)/#)&", z);
        cur = expr_new_integer(1);
        steps = n - 1;
    } else {
        return NULL;
    }
    if (!F || steps > HG_1F1_MAX_STEPS || steps < -HG_1F1_MAX_STEPS) {
        expr_free(F); expr_free(cur);
        return NULL;
    }
    F = eval_and_free(F);
    for (int64_t k = 0; k != steps; k += (steps > 0 ? 1 : -1)) {
        const char* rec = steps > 0 ? "((#1 + 1) (E^#3 - #2)/#3)&"      /* up   */
                                    : "(E^#3 - #3 #2/#1)&";              /* down */
        Expr* f = parse_expression(rec);
        Expr* call = expr_new_function(f, (Expr*[]){ expr_copy(cur), F, expr_copy((Expr*)z) }, 3);
        F = eval_and_free(call);
        Expr* nc = expr_new_function(expr_new_symbol(SYM_Plus),
                       (Expr*[]){ cur, expr_new_integer(steps > 0 ? 1 : -1) }, 2);
        cur = eval_and_free(nc);
    }
    expr_free(cur);
    return F;
}

static Expr* try_reduce(Expr* a, Expr* b, Expr* z) {
    size_t p = a->data.function.arg_count;
    size_t q = b->data.function.arg_count;
    Expr** au = a->data.function.args;
    Expr** bl = b->data.function.args;
    Expr* out = NULL;

    if (p == 0 && q == 0) {
        /* 0F0(;;z) = E^z */
        out = expr_new_function(expr_new_symbol(SYM_Exp), (Expr*[]){ expr_copy(z) }, 1);
    } else if (p == 1 && q == 0) {
        /* 1F0(a;;z) = (1 - z)^(-a) */
        Expr* nega = expr_new_function(expr_new_symbol(SYM_Times),
                         (Expr*[]){ expr_new_integer(-1), expr_copy(au[0]) }, 2);
        out = expr_new_function(expr_new_symbol(SYM_Power),
                  (Expr*[]){ one_minus(z), nega }, 2);
    } else if (p == 0 && q == 1) {
        int64_t tw;
        if (is_half_int(bl[0], &tw) && tw == 1) {
            /* 0F1(;1/2;z) = Cosh[2 Sqrt[z]] */
            Expr* arg = expr_new_function(expr_new_symbol(SYM_Times),
                            (Expr*[]){ expr_new_integer(2), sqrt_of(z) }, 2);
            out = expr_new_function(expr_new_symbol(SYM_Cosh), (Expr*[]){ arg }, 1);
        } else if (is_half_int(bl[0], &tw) && tw == 3) {
            /* 0F1(;3/2;z) = Sinh[2 Sqrt[z]] / (2 Sqrt[z]) */
            Expr* twosq = expr_new_function(expr_new_symbol(SYM_Times),
                              (Expr*[]){ expr_new_integer(2), sqrt_of(z) }, 2);
            Expr* sinh = expr_new_function(expr_new_symbol(SYM_Sinh),
                             (Expr*[]){ expr_copy(twosq) }, 1);
            Expr* inv = expr_new_function(expr_new_symbol(SYM_Power),
                            (Expr*[]){ twosq, expr_new_integer(-1) }, 2);
            out = expr_new_function(expr_new_symbol(SYM_Times),
                      (Expr*[]){ sinh, inv }, 2);
        }
    } else if (p == 1 && q == 1) {
        /* 1F1(1;2;z) = (E^z - 1) / z */
        if (is_int_val(au[0], 1) && is_int_val(bl[0], 2)) {
            Expr* em1 = expr_new_function(expr_new_symbol(SYM_Plus),
                            (Expr*[]){ expr_new_function(expr_new_symbol(SYM_Exp),
                                           (Expr*[]){ expr_copy(z) }, 1),
                                       expr_new_integer(-1) }, 2);
            Expr* invz = expr_new_function(expr_new_symbol(SYM_Power),
                             (Expr*[]){ expr_copy(z), expr_new_integer(-1) }, 2);
            out = expr_new_function(expr_new_symbol(SYM_Times),
                      (Expr*[]){ em1, invz }, 2);
        } else if (hg_diff_is(bl[0], au[0], 1) && !is_nonpos_int(au[0], NULL)) {
            /* 1F1(a; a+1; z): incomplete gamma / error function. */
            out = hg_1f1_a_ap1(au[0], z);
        } else if (is_int_val(au[0], 1) && !is_nonpos_int(bl[0], NULL)) {
            /* Kummer: 1F1(1; c; z) = E^z 1F1(c-1; c; -z), the family above
             * with a = c - 1 (c = 2 was handled first; c - 1 a non-positive
             * integer would be a terminating series on the right but not on
             * the left, and is excluded with c itself). */
            Expr* cm1 = eval_and_free(expr_new_function(expr_new_symbol(SYM_Plus),
                            (Expr*[]){ expr_copy(bl[0]), expr_new_integer(-1) }, 2));
            Expr* mz = expr_new_function(expr_new_symbol(SYM_Times),
                            (Expr*[]){ expr_new_integer(-1), expr_copy(z) }, 2);
            Expr* mze = eval_and_free(mz);
            if (!is_nonpos_int(cm1, NULL)) {
                Expr* inner = hg_1f1_a_ap1(cm1, mze);
                if (inner)
                    out = expr_new_function(expr_new_symbol(SYM_Times),
                              (Expr*[]){ expr_new_function(expr_new_symbol(SYM_Exp),
                                             (Expr*[]){ expr_copy(z) }, 1), inner }, 2);
            }
            expr_free(cm1); expr_free(mze);
        }
    } else if (p == 2 && q == 1) {
        int64_t tw;
        bool u11 = is_int_val(au[0], 1) && is_int_val(au[1], 1);
        bool u12 = (is_int_val(au[0], 1) && is_int_val(au[1], 2)) ||
                   (is_int_val(au[0], 2) && is_int_val(au[1], 1));
        int64_t gn, gd;
        bool rat3 = is_int_val(z, 1) && is_rational(au[0], &gn, &gd) &&
                    is_rational(au[1], &gn, &gd) && is_rational(bl[0], &gn, &gd) &&
                    !is_nonpos_int(bl[0], NULL);
        Expr* excess = NULL;
        if (rat3) {
            /* Gauss: 2F1(a,b;c;1) = G(c) G(c-a-b) / (G(c-a) G(c-b)) for
             * c - a - b > 0 (rational parameters, as the reference system
             * sums it; symbolic ones stay). */
            Expr* t = parse_expression("(#3 - #1 - #2)&");
            excess = eval_and_free(expr_new_function(t,
                         (Expr*[]){ expr_copy(au[0]), expr_copy(au[1]), expr_copy(bl[0]) }, 3));
            int64_t en, ed;
            if (!(is_rational(excess, &en, &ed) && en > 0)) rat3 = false;
            expr_free(excess);
        }
        if (rat3) {
            Expr* f = parse_expression(
                "(Gamma[#3] Gamma[#3 - #1 - #2] / (Gamma[#3 - #1] Gamma[#3 - #2]))&");
            out = expr_new_function(f, (Expr*[]){ expr_copy(au[0]), expr_copy(au[1]),
                                                  expr_copy(bl[0]) }, 3);
        } else if (u11 && is_int_val(bl[0], 2)) {
            /* 2F1(1,1;2;z) = -Log[1 - z] / z */
            Expr* lg = expr_new_function(expr_new_symbol(SYM_Log),
                           (Expr*[]){ one_minus(z) }, 1);
            Expr* invz = expr_new_function(expr_new_symbol(SYM_Power),
                             (Expr*[]){ expr_copy(z), expr_new_integer(-1) }, 2);
            out = expr_new_function(expr_new_symbol(SYM_Times),
                      (Expr*[]){ expr_new_integer(-1), lg, invz }, 3);
        } else if (u11 && is_half_int(bl[0], &tw) && tw == 1) {
            /* 2F1(1,1;1/2;z) = 1/(1-z) + Sqrt[z] ArcSin[Sqrt[z]] / (1-z)^(3/2).
             * (Central-binomial family, e.g. Sum[2^k/Binomial[2k,k],{k,0,Inf}].) */
            Expr* t1 = expr_new_function(expr_new_symbol(SYM_Power),
                           (Expr*[]){ one_minus(z), expr_new_integer(-1) }, 2);
            Expr* asin = expr_new_function(expr_new_symbol(SYM_ArcSin),
                             (Expr*[]){ sqrt_of(z) }, 1);
            Expr* omz32 = expr_new_function(expr_new_symbol(SYM_Power),
                              (Expr*[]){ one_minus(z), make_rational(-3, 2) }, 2);
            Expr* t2 = expr_new_function(expr_new_symbol(SYM_Times),
                           (Expr*[]){ sqrt_of(z), asin, omz32 }, 3);
            out = expr_new_function(expr_new_symbol(SYM_Plus), (Expr*[]){ t1, t2 }, 2);
        } else if (u11 && is_half_int(bl[0], &tw) && tw == 3) {
            /* 2F1(1,1;3/2;z) = ArcSin[Sqrt[z]] / (Sqrt[z] Sqrt[1-z]). */
            Expr* asin = expr_new_function(expr_new_symbol(SYM_ArcSin),
                             (Expr*[]){ sqrt_of(z) }, 1);
            Expr* invsq = expr_new_function(expr_new_symbol(SYM_Power),
                              (Expr*[]){ sqrt_of(z), expr_new_integer(-1) }, 2);
            Expr* invomz = expr_new_function(expr_new_symbol(SYM_Power),
                               (Expr*[]){ one_minus(z), make_rational(-1, 2) }, 2);
            out = expr_new_function(expr_new_symbol(SYM_Times),
                      (Expr*[]){ asin, invsq, invomz }, 3);
        } else if (u12 && is_half_int(bl[0], &tw) && tw == 3) {
            /* 2F1(2,1;3/2;z) = 1/(2(1-z)) + ArcSin[Sqrt[z]] / (2 Sqrt[z] (1-z)^(3/2)).
             * (e.g. Sum[2^k/Binomial[2k,k],{k,1,Inf}] = 1 + Pi/2.) */
            Expr* t1 = expr_new_function(expr_new_symbol(SYM_Times),
                           (Expr*[]){ make_rational(1, 2),
                                      expr_new_function(expr_new_symbol(SYM_Power),
                                          (Expr*[]){ one_minus(z), expr_new_integer(-1) }, 2) }, 2);
            Expr* asin = expr_new_function(expr_new_symbol(SYM_ArcSin),
                             (Expr*[]){ sqrt_of(z) }, 1);
            Expr* invsq = expr_new_function(expr_new_symbol(SYM_Power),
                              (Expr*[]){ sqrt_of(z), expr_new_integer(-1) }, 2);
            Expr* omz32 = expr_new_function(expr_new_symbol(SYM_Power),
                              (Expr*[]){ one_minus(z), make_rational(-3, 2) }, 2);
            Expr* t2 = expr_new_function(expr_new_symbol(SYM_Times),
                           (Expr*[]){ make_rational(1, 2), asin, invsq, omz32 }, 4);
            out = expr_new_function(expr_new_symbol(SYM_Plus), (Expr*[]){ t1, t2 }, 2);
        } else if (is_half_int(bl[0], &tw) && tw == 3 &&
                   ((hg_is_rat(au[0], 1, 2) && is_int_val(au[1], 1)) ||
                    (is_int_val(au[0], 1) && hg_is_rat(au[1], 1, 2)))) {
            /* 2F1(1/2, 1; 3/2; z) = ArcTanh[Sqrt[z]]/Sqrt[z]
             *                     = ArcTan[Sqrt[-z]]/Sqrt[-z]  (even in Sqrt). */
            out = hg_negated(z) ? hg_apply_tmpl("(ArcTan[Sqrt[-#]]/Sqrt[-#])&", z)
                                : hg_apply_tmpl("(ArcTanh[Sqrt[#]]/Sqrt[#])&", z);
        } else if (is_half_int(bl[0], &tw) && tw == 3 &&
                   hg_is_rat(au[0], 1, 2) && hg_is_rat(au[1], 1, 2)) {
            /* 2F1(1/2, 1/2; 3/2; z) = ArcSin[Sqrt[z]]/Sqrt[z]
             *                       = ArcSinh[Sqrt[-z]]/Sqrt[-z]. */
            out = hg_negated(z) ? hg_apply_tmpl("(ArcSinh[Sqrt[-#]]/Sqrt[-#])&", z)
                                : hg_apply_tmpl("(ArcSin[Sqrt[#]]/Sqrt[#])&", z);
        }
    } else if (p == 4 && q == 3) {
        /* Very-well-poised Ramanujan 1/Pi series (table-backed; general 1/Pi
         * summation is an open problem).  The classic central-binomial-cubed
         * family Sum[(4k+1) Binomial[2k,k]^3/(-64)^k] maps to
         *     4F3({1/2,1/2,1/2,5/4}, {1/4,1,1}, -1) = 2/Pi.
         * Matched as a multiset (parameter order is not canonical). */
        if (is_int_val(z, -1)) {
            int c_half = 0, c_54 = 0, other_u = 0;
            for (size_t i = 0; i < 4; i++) {
                if (hg_is_rat(au[i], 1, 2))      c_half++;
                else if (hg_is_rat(au[i], 5, 4)) c_54++;
                else                             other_u++;
            }
            int c_14 = 0, c_1 = 0, other_l = 0;
            for (size_t j = 0; j < 3; j++) {
                if (hg_is_rat(bl[j], 1, 4))   c_14++;
                else if (is_int_val(bl[j], 1)) c_1++;
                else                           other_l++;
            }
            if (c_half == 3 && c_54 == 1 && other_u == 0
                && c_14 == 1 && c_1 == 2 && other_l == 0) {
                out = expr_new_function(expr_new_symbol(SYM_Times),
                          (Expr*[]){ expr_new_integer(2),
                                     expr_new_function(expr_new_symbol(SYM_Power),
                                         (Expr*[]){ expr_new_symbol(SYM_Pi),
                                                    expr_new_integer(-1) }, 2) }, 2);
            }
        }
    }

    if (!out) return NULL;
    Expr* r = evaluate(out);
    expr_free(out);
    return r;
}

/* ------------------------------------------------------------------ */
/*  Numeric: machine double-complex                                    */
/* ------------------------------------------------------------------ */

/* Machine-precision numeric approximation of e into *out; sets *inexact if
 * any contributing leaf was inexact (Real / MPFR).  Returns false if e is
 * not a concrete number. */
static bool approx_machine(const Expr* e, double complex* out, bool* inexact) {
    if (e->type == EXPR_INTEGER) { *out = (double)e->data.integer; return true; }
    if (e->type == EXPR_REAL)    { *out = e->data.real; if (inexact) *inexact = true; return true; }
#ifdef USE_MPFR
    if (e->type == EXPR_MPFR)    { *out = mpfr_get_d(e->data.mpfr, MPFR_RNDN);
                                   if (inexact) *inexact = true;
                                   return true; }
#endif
    int64_t nn, dd;
    if (is_rational((Expr*)e, &nn, &dd)) { *out = (double)nn / (double)dd; return true; }
    Expr *re, *im;
    if (is_complex((Expr*)e, &re, &im)) {
        double complex r, i;
        if (approx_machine(re, &r, inexact) && approx_machine(im, &i, inexact)) {
            *out = creal(r) + creal(i) * I;
            return true;
        }
    }
    return false;
}

/* True if z is (near) a non-positive integer -- a lower-parameter pole. */
static bool near_nonpos_int(double complex z) {
    if (fabs(cimag(z)) > 1e-9) return false;
    double r = creal(z);
    return r <= 1e-9 && fabs(r - round(r)) < 1e-9;
}

static Expr* machine_sum(double complex zc,
                         const double complex* ac, size_t p,
                         const double complex* bc, size_t q,
                         bool all_real, double* lost_bits) {
    double complex term = 1.0, sum = 1.0;
    double maxmag = 1.0;
    int settled = 0;
    *lost_bits = 0.0;
    for (long k = 0; k < HGPFQ_MACHINE_MAX_TERMS; k++) {
        double complex fac = zc / (double)(k + 1);
        for (size_t i = 0; i < p; i++) fac *= (ac[i] + (double)k);
        for (size_t j = 0; j < q; j++) {
            double complex d = bc[j] + (double)k;
            if (cabs(d) < 1e-300) return NULL;   /* denominator pole */
            fac /= d;
        }
        term *= fac;
        double tm = cabs(term);
        if (tm > maxmag) maxmag = tm;
        sum += term;
        if (cabs(term) <= cabs(sum) * 1e-16) {
            if (++settled >= 2) break;
        } else {
            settled = 0;
        }
        if (!isfinite(creal(sum)) || !isfinite(cimag(sum))) return NULL;
    }
    /* How much of the double's precision this summation has already destroyed.
     * The series is alternating for negative real z, so its terms grow far past
     * the value they sum to (|term|_max ~ e^|z| for the confluent cases); every
     * bit of that excess magnitude is a bit of the answer that cancelled away
     * and cannot be recovered by summing more carefully in the same precision.
     * The caller re-sums in MPFR when this is large -- reporting it here rather
     * than guessing from |z| keeps the estimate exact for every p, q and every
     * parameter set, including the ones where no cancellation happens at all. */
    double as = cabs(sum);
    if (as > 0.0 && maxmag > as) *lost_bits = log2(maxmag / as);

    if (all_real && cimag(sum) == 0.0) return expr_new_real(creal(sum));
    return make_complex(expr_new_real(creal(sum)), expr_new_real(cimag(sum)));
}

/* ------------------------------------------------------------------ */
/*  Numeric: MPFR (real & complex via (re,im) pairs)                   */
/* ------------------------------------------------------------------ */

#ifdef USE_MPFR
typedef struct { mpfr_t re, im; } cpx_t;

static void cpx_init(cpx_t* c, mpfr_prec_t bits) { mpfr_init2(c->re, bits); mpfr_init2(c->im, bits); }
static void cpx_clear(cpx_t* c) { mpfr_clear(c->re); mpfr_clear(c->im); }
static void cpx_set(cpx_t* d, const cpx_t* s) {
    mpfr_set(d->re, s->re, MPFR_RNDN); mpfr_set(d->im, s->im, MPFR_RNDN);
}
static void cpx_add(cpx_t* o, const cpx_t* x, const cpx_t* y) {
    mpfr_add(o->re, x->re, y->re, MPFR_RNDN);
    mpfr_add(o->im, x->im, y->im, MPFR_RNDN);
}
/* o = x * y  (o must not alias x or y) */
static void cpx_mul(cpx_t* o, const cpx_t* x, const cpx_t* y, mpfr_prec_t bits) {
    mpfr_t t1, t2;
    mpfr_init2(t1, bits); mpfr_init2(t2, bits);
    mpfr_mul(t1, x->re, y->re, MPFR_RNDN);
    mpfr_mul(t2, x->im, y->im, MPFR_RNDN);
    mpfr_sub(t1, t1, t2, MPFR_RNDN);          /* re = xr*yr - xi*yi */
    mpfr_mul(o->im, x->re, y->im, MPFR_RNDN);
    mpfr_mul(t2, x->im, y->re, MPFR_RNDN);
    mpfr_add(o->im, o->im, t2, MPFR_RNDN);    /* im = xr*yi + xi*yr */
    mpfr_set(o->re, t1, MPFR_RNDN);
    mpfr_clear(t1); mpfr_clear(t2);
}
/* o = x / y */
static void cpx_div(cpx_t* o, const cpx_t* x, const cpx_t* y) {
    mpfr_complex_div(o->re, o->im, x->re, x->im, y->re, y->im);
}
/* magnitude proxy max(|re|,|im|) into m */
static void cpx_absmax(mpfr_t m, const cpx_t* x) {
    mpfr_abs(m, x->re, MPFR_RNDN);
    mpfr_t t; mpfr_init2(t, mpfr_get_prec(m));
    mpfr_abs(t, x->im, MPFR_RNDN);
    if (mpfr_cmp(t, m) > 0) mpfr_set(m, t, MPFR_RNDN);
    mpfr_clear(t);
}

/* Extract e into a cpx pair at `bits` precision; false if not numeric. */
static bool cpx_from_expr(const Expr* e, cpx_t* c, mpfr_prec_t bits) {
    cpx_init(c, bits);
    bool inexact = false;
    if (!get_approx_mpfr(e, c->re, c->im, &inexact)) { cpx_clear(c); return false; }
    return true;
}

static Expr* hgpfq_mpfr_sum(const Expr* z, Expr* a, Expr* b,
                            size_t p, size_t q, long target_bits) {
    mpfr_prec_t work = (mpfr_prec_t)(target_bits + 64);

    cpx_t zc, *ac = NULL, *bc = NULL;
    if (!cpx_from_expr(z, &zc, work)) return NULL;
    if (p) ac = malloc(sizeof(cpx_t) * p);
    if (q) bc = malloc(sizeof(cpx_t) * q);
    bool ok = true;
    /* Track how many entries were successfully initialised; a failing
     * cpx_from_expr clears its own slot, so on the error path we must clear
     * exactly the prior survivors (else their mpfr_t limbs leak). */
    size_t na = 0, nb = 0;
    for (size_t i = 0; i < p && ok; i++) {
        ok = cpx_from_expr(a->data.function.args[i], &ac[i], work);
        if (ok) na++;
    }
    for (size_t j = 0; j < q && ok; j++) {
        ok = cpx_from_expr(b->data.function.args[j], &bc[j], work);
        if (ok) nb++;
    }
    if (!ok) {
        cpx_clear(&zc);
        for (size_t i = 0; i < na; i++) cpx_clear(&ac[i]);
        for (size_t j = 0; j < nb; j++) cpx_clear(&bc[j]);
        free(ac); free(bc);
        return NULL;
    }

    cpx_t term, sum, fac, tmp, num;
    cpx_init(&term, work); cpx_init(&sum, work);
    cpx_init(&fac, work);  cpx_init(&tmp, work); cpx_init(&num, work);
    mpfr_set_ui(term.re, 1, MPFR_RNDN); mpfr_set_ui(term.im, 0, MPFR_RNDN);
    mpfr_set_ui(sum.re, 1, MPFR_RNDN);  mpfr_set_ui(sum.im, 0, MPFR_RNDN);

    mpfr_t tmag, smag, thr;
    mpfr_init2(tmag, work); mpfr_init2(smag, work); mpfr_init2(thr, work);

    bool failed = false;
    int settled = 0;
    long k = 0;
    for (; k < HGPFQ_MPFR_MAX_TERMS; k++) {
        /* fac = z / (k+1) * prod(a_i + k) / prod(b_j + k) */
        cpx_set(&fac, &zc);
        mpfr_div_ui(fac.re, fac.re, (unsigned long)(k + 1), MPFR_RNDN);
        mpfr_div_ui(fac.im, fac.im, (unsigned long)(k + 1), MPFR_RNDN);
        for (size_t i = 0; i < p; i++) {
            cpx_set(&tmp, &ac[i]);
            mpfr_add_ui(tmp.re, tmp.re, (unsigned long)k, MPFR_RNDN);
            cpx_mul(&num, &fac, &tmp, work);
            cpx_set(&fac, &num);
        }
        for (size_t j = 0; j < q; j++) {
            cpx_set(&tmp, &bc[j]);
            mpfr_add_ui(tmp.re, tmp.re, (unsigned long)k, MPFR_RNDN);
            cpx_div(&num, &fac, &tmp);
            cpx_set(&fac, &num);
            if (mpfr_nan_p(fac.re) || mpfr_inf_p(fac.re) ||
                mpfr_nan_p(fac.im) || mpfr_inf_p(fac.im)) { failed = true; break; }
        }
        if (failed) break;
        cpx_mul(&num, &term, &fac, work);   /* term *= fac */
        cpx_set(&term, &num);
        cpx_add(&sum, &sum, &term);

        cpx_absmax(tmag, &term);
        cpx_absmax(smag, &sum);
        mpfr_mul_2si(thr, smag, -(long)work, MPFR_RNDN);   /* |sum| * 2^-work */
        if (mpfr_cmp(tmag, thr) <= 0) { if (++settled >= 2) break; }
        else settled = 0;
    }

    Expr* result = NULL;
    if (!failed && k < HGPFQ_MPFR_MAX_TERMS) {
        mpfr_t orr, oii;
        mpfr_init2(orr, (mpfr_prec_t)target_bits);
        mpfr_init2(oii, (mpfr_prec_t)target_bits);
        mpfr_set(orr, sum.re, MPFR_RNDN);
        mpfr_set(oii, sum.im, MPFR_RNDN);
        result = numeric_mpfr_make_complex(orr, oii);
        mpfr_clear(orr); mpfr_clear(oii);
    }

    mpfr_clear(tmag); mpfr_clear(smag); mpfr_clear(thr);
    cpx_clear(&term); cpx_clear(&sum); cpx_clear(&fac); cpx_clear(&tmp); cpx_clear(&num);
    cpx_clear(&zc);
    for (size_t i = 0; i < p; i++) cpx_clear(&ac[i]);
    for (size_t j = 0; j < q; j++) cpx_clear(&bc[j]);
    free(ac); free(bc);
    return result;
}
/* Round a high-precision sum back down to machine precision, preserving the
 * shape (bare real, or Complex[re, im]).  Returns NULL if `e` is not one of
 * those forms, in which case the caller keeps the plain machine result. */
static double hgpfq_num_to_d(const Expr* e) {
    switch (e->type) {
        case EXPR_REAL:    return e->data.real;
        case EXPR_INTEGER: return (double)e->data.integer;
        case EXPR_MPFR:    return mpfr_get_d(e->data.mpfr, MPFR_RNDN);
        default:           return 0.0;
    }
}
static bool hgpfq_is_num(const Expr* e) {
    return e->type == EXPR_REAL || e->type == EXPR_INTEGER || e->type == EXPR_MPFR;
}
static Expr* hgpfq_demote_to_machine(Expr* e) {
    if (hgpfq_is_num(e)) return expr_new_real(hgpfq_num_to_d(e));
    Expr *re = NULL, *im = NULL;
    if (is_complex(e, &re, &im) && hgpfq_is_num(re) && hgpfq_is_num(im))
        return make_complex(expr_new_real(hgpfq_num_to_d(re)),
                            expr_new_real(hgpfq_num_to_d(im)));
    return NULL;
}
#endif /* USE_MPFR */

/* ------------------------------------------------------------------ */
/*  Numeric dispatch                                                   */
/* ------------------------------------------------------------------ */

static Expr* try_numeric(Expr* a, Expr* b, Expr* z) {
    size_t p = a->data.function.arg_count;
    size_t q = b->data.function.arg_count;

    /* Machine extraction of z and all parameters: gating + machine path. */
    double complex zc; bool inexact = false;
    if (!approx_machine(z, &zc, &inexact)) return NULL;

    double complex* ac = p ? malloc(sizeof(double complex) * p) : NULL;
    double complex* bc = q ? malloc(sizeof(double complex) * q) : NULL;
    bool all_real = (cimag(zc) == 0.0);
    bool ok = true;
    for (size_t i = 0; i < p && ok; i++) {
        ok = approx_machine(a->data.function.args[i], &ac[i], &inexact);
        if (ok && cimag(ac[i]) != 0.0) all_real = false;
    }
    for (size_t j = 0; j < q && ok; j++) {
        ok = approx_machine(b->data.function.args[j], &bc[j], &inexact);
        if (ok && cimag(bc[j]) != 0.0) all_real = false;
        if (ok && near_nonpos_int(bc[j])) ok = false;   /* lower-parameter pole */
    }
    if (!ok || !inexact) { free(ac); free(bc); return NULL; }

#ifdef USE_MPFR
    /* 2F1 at arbitrary precision OUTSIDE the unit disc: arb's kernel is the
     * principal-branch analytic continuation, so it is the value wherever that
     * is single-valued -- everywhere but the cut [1, Infinity), where the side
     * is a convention left to the caller. The machine twin already covers
     * |z| >= 1 (hg_machine_pfq); without this, N[2F1(...; -1), 20] stayed
     * unevaluated where the machine N gave a number. */
    if (p == 2 && q == 1 && cabs(zc) >= 1.0 - 1e-12 &&
        !(cimag(zc) == 0.0 && creal(zc) >= 1.0 - 1e-12)) {
        long gov = numeric_min_inexact_bits(z);
        for (size_t i = 0; i < p; i++) {
            long bi = numeric_min_inexact_bits(a->data.function.args[i]);
            if (bi > 0 && (gov == 0 || bi < gov)) gov = bi;
        }
        long bj = numeric_min_inexact_bits(b->data.function.args[0]);
        if (bj > 0 && (gov == 0 || bj < gov)) gov = bj;
        Expr* v = flint_num_hyp2f1(a->data.function.args[0], a->data.function.args[1],
                                   b->data.function.args[0], z);
        /* A machine-precision call (the machine kernel declined it, e.g. a
         * complex parameter) still answers a machine number. */
        if (v && gov <= 53) {
            Expr* dn = hgpfq_demote_to_machine(v);
            expr_free(v);
            v = dn;
        }
        if (v) { free(ac); free(bc); return v; }
    }
#endif

    /* Convergence gate (radius of convergence by p vs q). */
    double az = cabs(zc);
    if (p > q + 1) { free(ac); free(bc); return NULL; }        /* divergent */
    if (p == q + 1 && az >= 1.0 - 1e-12) { free(ac); free(bc); return NULL; }

    Expr* result = NULL;
#ifdef USE_MPFR
    /* Precision contagion: the governing precision is the MINIMUM among the
     * inexact (Real/MPFR) leaves of z and every parameter. A machine-precision
     * argument forces the machine path even when others are high-precision
     * MPFR (matching Mathematica). Only when every inexact argument exceeds
     * machine precision do we take the arbitrary-precision MPFR sum. */
    {
        long gov = numeric_min_inexact_bits(z);
        for (size_t i = 0; i < p; i++) {
            long bi = numeric_min_inexact_bits(a->data.function.args[i]);
            if (bi > 0 && (gov == 0 || bi < gov)) gov = bi;
        }
        for (size_t j = 0; j < q; j++) {
            long bj = numeric_min_inexact_bits(b->data.function.args[j]);
            if (bj > 0 && (gov == 0 || bj < gov)) gov = bj;
        }
        if (gov > 53) {
            /* 2F1 inside the unit disc: arb's transformation-based kernel,
             * not the direct series, whose term count blows up as |z| -> 1
             * (the convergence gate above has already excluded |z| >= 1). */
            if (p == 2 && q == 1) {
                result = flint_num_hyp2f1(a->data.function.args[0], a->data.function.args[1],
                                          b->data.function.args[0], z);
                if (result) { free(ac); free(bc); return result; }
            }
            result = hgpfq_mpfr_sum(z, a, b, p, q, gov);
            free(ac); free(bc);
            return result;
        }
    }
#endif
    double lost = 0.0;
    result = machine_sum(zc, ac, p, bc, q, all_real, &lost);
#ifdef USE_MPFR
    /* Machine-precision arguments must still give a machine-precision ANSWER,
     * and a doubles-only sum of a badly cancelling series does not: 1F1(1;2;-40)
     * came back with 5e-2 relative error because its largest term is ~2e17 times
     * the result.  Re-sum with enough working precision to absorb the loss and
     * round back, so the argument precision -- not the conditioning of the
     * series -- decides the output type. */
    if (result && lost > 4.0 && lost < 4096.0) {
        Expr* hi = hgpfq_mpfr_sum(z, a, b, p, q, 53 + (long)lost + 16);
        if (hi) {
            Expr* dn = hgpfq_demote_to_machine(hi);
            expr_free(hi);
            if (dn) { expr_free(result); result = dn; }
        }
    }
#endif
    free(ac); free(bc);
    return result;
}

/* ------------------------------------------------------------------ */
/*  Machine-precision fast path (scalar and buffer)                    */
/* ------------------------------------------------------------------ */

#define HG_MAXPARAMS 16

/* One argument as a machine (re, im) pair: a machine-scalar leaf (Real,
 * Integer, BigInt, Rational) or Complex[] of two such. *inexact is OR-ed. */
static bool hg_leaf(const Expr* e, double* re, double* im, bool* inexact, bool* cplx) {
    bool ix = false;
    if (sf_machine_leaf(e, re, &ix)) { *im = 0.0; *inexact |= ix; return true; }
    if (e->type == EXPR_FUNCTION && e->data.function.head->type == EXPR_SYMBOL &&
        e->data.function.head->data.symbol.name == SYM_Complex &&
        e->data.function.arg_count == 2) {
        bool ia = false, ib = false;
        if (!sf_machine_leaf(e->data.function.args[0], re, &ia) ||
            !sf_machine_leaf(e->data.function.args[1], im, &ib)) return false;
        *inexact |= (ia || ib);
        *cplx = true;
        return true;
    }
    return false;
}

typedef struct {
    size_t p, q;
    double ar[HG_MAXPARAMS], ai[HG_MAXPARAMS], br[HG_MAXPARAMS], bi[HG_MAXPARAMS];
    bool inexact, cplx;
} HgParams;

static bool hg_params(const Expr* a, const Expr* b, HgParams* hp) {
    hp->p = a->data.function.arg_count;
    hp->q = b->data.function.arg_count;
    if (hp->p > HG_MAXPARAMS || hp->q > HG_MAXPARAMS) return false;
    hp->inexact = false; hp->cplx = false;
    for (size_t i = 0; i < hp->p; i++)
        if (!hg_leaf(a->data.function.args[i], &hp->ar[i], &hp->ai[i],
                     &hp->inexact, &hp->cplx)) return false;
    for (size_t j = 0; j < hp->q; j++)
        if (!hg_leaf(b->data.function.args[j], &hp->br[j], &hp->bi[j],
                     &hp->inexact, &hp->cplx)) return false;
    return true;
}

/* Evaluate one point. Real parameters and a real z use the real kernel; for
 * 2F1 on its cut (real z > 1) the value is complex and comes from the complex
 * kernel as the limit from below (zi = -0), the Mathematica / mpmath branch.
 * *is_c reports whether the result is complex. */
static bool hg_point(const HgParams* hp, double zr, double zi, bool zc,
                     double* rr, double* ri, bool* is_c) {
    if (!hp->cplx && !zc) {
        *is_c = false; *ri = 0.0;
        if (sf_machine_pfq(hp->ar, hp->p, hp->br, hp->q, zr, rr)) return true;
        if (hp->p == 2 && hp->q == 1 && zr > 1.0) {
            *is_c = true;
            return sf_machine_pfq_c(hp->ar, hp->ai, 2, hp->br, hp->bi, 1, zr, -0.0, rr, ri);
        }
        return false;
    }
    *is_c = true;
    return sf_machine_pfq_c(hp->ar, hp->ai, hp->p, hp->br, hp->bi, hp->q, zr, zi, rr, ri);
}

/* The scalar machine gate: every parameter and z a machine-scalar leaf, at
 * least one of them a machine Real. NULL when the call is not machine
 * precision or the kernel declined (the caller's exact / MPFR route answers). */
static Expr* hg_machine_pfq(const Expr* a, const Expr* b, const Expr* z) {
    HgParams hp;
    if (!hg_params(a, b, &hp)) return NULL;
    double zr, zi; bool zc = false;
    if (!hg_leaf(z, &zr, &zi, &hp.inexact, &zc)) return NULL;
    if (!hp.inexact) return NULL;                 /* exact call stays exact */
    double rr, ri; bool is_c;
    if (!hg_point(&hp, zr, zi, zc, &rr, &ri, &is_c)) return NULL;
    return is_c ? sf_machine_complex_result(rr, ri) : expr_new_real(rr);
}

/* The buffer path: z an NDArray / packed List, parameters as above. Returns
 * NULL when the buffer or a parameter is not a machine case, or any element
 * declines. The output dtype is the input's, widened to complex when the
 * parameters are complex or a 2F1 element lies on the cut. */
static Expr* hg_pfq_ndarray(const Expr* a, const Expr* b, const Expr* z) {
    NDType zdt = z->data.ndarray.dtype;
    if (zdt == NDT_INT64) return NULL;            /* exact data: element-wise */
    HgParams hp;
    if (!hg_params(a, b, &hp)) return NULL;
    bool zc = ndt_is_complex(zdt);
    size_t sz = ndarray_size(z);
    /* The common case -- real parameters over a float64 buffer -- goes straight
     * through the vector driver, which hoists the parameter work (2F1
     * coefficient tables) out of the loop and is bit-identical to the scalar
     * kernel. A partial run (some element declined, or a 2F1 element on the
     * cut, whose value is complex) falls to the general loop below. */
    if (!hp.cplx && zdt == NDT_FLOAT64) {
        double* out = malloc(sizeof(double) * (sz ? sz : 1));
        if (!out) return NULL;
        size_t done = sf_machine_pfq_vec(hp.ar, hp.p, hp.br, hp.q,
                                         (const double*)z->data.ndarray.data, out, sz);
        if (done == sz)
            return expr_new_ndarray_like(z, z->data.ndarray.rank, z->data.ndarray.dims,
                                         out, NDT_FLOAT64);
        free(out);
    }
    /* Compute into a complex scratch buffer; narrow to real at the end when
     * nothing came back complex. */
    double* buf = malloc(sizeof(double) * 2 * (sz ? sz : 1));
    if (!buf) return NULL;
    bool any_c = zc || hp.cplx;
    for (size_t i = 0; i < sz; i++) {
        double xr, xi, rr, ri; bool is_c;
        ndt_get(z->data.ndarray.data, i, zdt, &xr, &xi);
        if (!hg_point(&hp, xr, xi, zc, &rr, &ri, &is_c)) { free(buf); return NULL; }
        buf[2 * i] = rr; buf[2 * i + 1] = ri;
        any_c |= is_c;
    }
    NDType odt = any_c ? (zdt == NDT_FLOAT32 || zdt == NDT_COMPLEX32 ? NDT_COMPLEX32
                                                                      : NDT_COMPLEX64)
                       : zdt;
    void* out = malloc(ndt_elem_size(odt) * (sz ? sz : 1));
    if (!out) { free(buf); return NULL; }
    for (size_t i = 0; i < sz; i++) ndt_set(out, i, odt, buf[2 * i], buf[2 * i + 1]);
    free(buf);
    return expr_new_ndarray_like(z, z->data.ndarray.rank, z->data.ndarray.dims, out, odt);
}

/* ------------------------------------------------------------------ */
/*  Entry point                                                        */
/* ------------------------------------------------------------------ */

Expr* builtin_hypergeometric_pfq(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 3) return NULL;
    Expr* a = res->data.function.args[0];
    Expr* b = res->data.function.args[1];
    Expr* z = res->data.function.args[2];
    if (!is_list(a) || !is_list(b)) return NULL;

    /* (1b) A machine buffer as the argument: map the machine kernel over it.
     *
     * This is what Hypergeometric0F1/1F1/2F1[..., arr] turn into -- they
     * rewrite to PFQ when the argument is a buffer -- so without it a visible
     * NDArray came back unevaluated and a packed List was materialised into one
     * Expr per element, each going round the evaluator (1.2-2.6 s per 10^6).
     *
     * Parameters may be exact (Integer / Rational) or machine Real, real or
     * complex: a machine-real buffer element makes every element call a
     * machine-precision call, and the scalar builtin answers such a call from
     * the SAME kernel (hg_machine_pfq, ahead of every exact rewrite), so
     * the two surfaces agree bit-for-bit. An int64 buffer is exact data and an
     * MPFR parameter makes the call arbitrary-precision: both take the
     * element-wise route, as does any element the kernel declines. */
    if (is_ndarray(z)) {
        Expr* r = hg_pfq_ndarray(a, b, z);
        if (r) return r;
        /* Declined: the buffer must not reach the scalar steps below, which
         * test z for Integer/Real and would leave the call unevaluated. */
        return ndarray_delist_and_reeval(res);
    }

    /* (2) Thread over a List third argument. */
    if (is_list(z)) {
        size_t n = z->data.function.arg_count;
        Expr** out = malloc(sizeof(Expr*) * (n ? n : 1));
        for (size_t i = 0; i < n; i++)
            out[i] = rebuild_eval(expr_copy(a), expr_copy(b),
                                  expr_copy(z->data.function.args[i]));
        Expr* lst = expr_new_function(expr_new_symbol(SYM_List), out, n);
        free(out);
        return lst;
    }

    /* (2b) Machine-precision call: the raw-double kernel, before any of the
     * exact rewrites below (which would re-evaluate a whole expression tree
     * per call, and would not match the buffer/Compile surfaces bit-for-bit). */
    {
        Expr* m = hg_machine_pfq(a, b, z);
        if (m) return m;
    }

    /* (3) Generic parameter cancellation. */
    Expr* c = try_cancel(a, b, z);
    if (c) return c;

    /* (3b) An upper parameter exceeding a lower one by a positive integer. */
    Expr* sp = try_shift_pair(a, b, z);
    if (sp) return sp;

    /* (4) z == 0 -> 1. */
    if (z->type == EXPR_INTEGER && z->data.integer == 0) return expr_new_integer(1);

    /* (5) Termination to a polynomial. */
    Expr* t = try_terminate(a, b, z);
    if (t) return t;

    /* (5b) Classical summation theorems at z = +-1. */
    Expr* cs = try_classical_sum(a, b, z);
    if (cs) return cs;

    /* (6) Reduction to elementary closed forms (analytic continuations). */
    Expr* red = try_reduce(a, b, z);
    if (red) return red;

    /* (7) Numeric evaluation. */
    Expr* n = try_numeric(a, b, z);
    if (n) return n;

    /* (7b) The machine series bails out (returns NULL) once a partial sum
     * overflows a double -- Hypergeometric1F1[1., 2., 900.] ~ 8*10^387. Promote
     * the whole call to MPFR (parameters inside the lists included) and retry at
     * the extended exponent range. */
    { Expr* p = numeric_promote_machine_call(res); if (p) return p; }

    /* (8) Stay unevaluated. */
    return NULL;
}

/* ------------------------------------------------------------------ */
/*  Convenience heads: rewrite to HypergeometricPFQ                    */
/* ------------------------------------------------------------------ */

/* Wrap a scalar parameter into a one-element List (copying). */
static Expr* list1(Expr* x) {
    return expr_new_function(expr_new_symbol(SYM_List), (Expr*[]){ expr_copy(x) }, 1);
}

/* The machine gate for the fixed-arity heads: all arguments machine-scalar
 * leaves (sf_machine_scalar_args), at least one a machine Real. The same
 * kernels as HypergeometricPFQ's gate (sf_machine_pfq dispatches these shapes
 * to them), so every spelling agrees bit-for-bit with the buffer and Compile
 * surfaces. */
static Expr* hg_machine_head(Expr* res, size_t p, size_t q) {
    SFMachineArgs ma;
    if (!sf_machine_scalar_args(res, &ma)) return NULL;
    HgParams hp;
    hp.p = p; hp.q = q; hp.inexact = true; hp.cplx = false;
    size_t k = 0;
    for (size_t i = 0; i < p; i++, k++) {
        hp.ar[i] = ma.re[k]; hp.ai[i] = ma.im[k]; hp.cplx |= ma.is_complex[k];
    }
    for (size_t j = 0; j < q; j++, k++) {
        hp.br[j] = ma.re[k]; hp.bi[j] = ma.im[k]; hp.cplx |= ma.is_complex[k];
    }
    double rr, ri; bool is_c;
    if (!hg_point(&hp, ma.re[k], ma.im[k], ma.is_complex[k], &rr, &ri, &is_c))
        return NULL;
    return is_c ? sf_machine_complex_result(rr, ri) : expr_new_real(rr);
}

Expr* builtin_hypergeometric_0f1(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 2) return NULL;
    Expr* m = hg_machine_head(res, 0, 1);
    if (m) return m;
    Expr* b = res->data.function.args[0];
    Expr* z = res->data.function.args[1];
    Expr* empty = expr_new_function(expr_new_symbol(SYM_List), NULL, 0);
    return rebuild_eval(empty, list1(b), expr_copy(z));
}

Expr* builtin_hypergeometric_1f1(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 3) return NULL;
    Expr* m = hg_machine_head(res, 1, 1);
    if (m) return m;
    Expr* a = res->data.function.args[0];
    Expr* b = res->data.function.args[1];
    Expr* z = res->data.function.args[2];
    return rebuild_eval(list1(a), list1(b), expr_copy(z));
}

Expr* builtin_hypergeometric_2f1(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 4) return NULL;
    Expr* m = hg_machine_head(res, 2, 1);
    if (m) return m;
    Expr* a = res->data.function.args[0];
    Expr* b = res->data.function.args[1];
    Expr* cc = res->data.function.args[2];
    Expr* z = res->data.function.args[3];
    Expr* up = expr_new_function(expr_new_symbol(SYM_List),
                   (Expr*[]){ expr_copy(a), expr_copy(b) }, 2);
    return rebuild_eval(up, list1(cc), expr_copy(z));
}

/* ------------------------------------------------------------------ */
/*  Registration                                                       */
/* ------------------------------------------------------------------ */

void hypergeopfq_init(void) {
    symtab_add_builtin("HypergeometricPFQ", builtin_hypergeometric_pfq);
    symtab_get_def("HypergeometricPFQ")->attributes |=
        ATTR_NUMERICFUNCTION | ATTR_PROTECTED;

    /* The fixed-arity heads are Listable, as in Mathematica (HypergeometricPFQ
     * is not: its first two arguments ARE lists). */
    symtab_add_builtin("Hypergeometric0F1", builtin_hypergeometric_0f1);
    symtab_get_def("Hypergeometric0F1")->attributes |=
        ATTR_LISTABLE | ATTR_NUMERICFUNCTION | ATTR_PROTECTED;

    symtab_add_builtin("Hypergeometric1F1", builtin_hypergeometric_1f1);
    symtab_get_def("Hypergeometric1F1")->attributes |=
        ATTR_LISTABLE | ATTR_NUMERICFUNCTION | ATTR_PROTECTED;

    symtab_add_builtin("Hypergeometric2F1", builtin_hypergeometric_2f1);
    symtab_get_def("Hypergeometric2F1")->attributes |=
        ATTR_LISTABLE | ATTR_NUMERICFUNCTION | ATTR_PROTECTED;
}
