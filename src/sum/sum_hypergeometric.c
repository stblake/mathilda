/*
 * sum_hypergeometric.c -- Sum`Hypergeometric: infinite hypergeometric series.
 *
 * For a summand t(k) whose term ratio t(k+1)/t(k) is a rational function of k,
 * the infinite sum from a concrete lower bound is a generalized hypergeometric
 * function:
 *
 *     Sum_{k=imin}^{Infinity} t(k) = t(imin) * pFq({a}, {b}, z),
 *
 * where the upper/lower parameters and the argument z are read off by factoring
 * the (reindexed) term ratio into monic linear factors over Q:
 *
 *     t(imin+m+1)/t(imin+m) = z * prod_i (m + a_i) / ( prod_j (m + b_j) (m+1) ).
 *
 * The (m+1) is the canonical factorial factor.  Numerator roots give the upper
 * parameters, the remaining denominator roots the lower parameters, and the
 * ratio of leading coefficients is z.  The result is handed to HypergeometricPFQ,
 * which reduces the common cases to elementary closed forms (e.g. 0F0 -> E^z,
 * the geometric 1F0(1) -> 1/(1-z)) and otherwise returns the pFq itself.
 *
 *   Sum`Hypergeometric[f, k, imin, Infinity]  -> t(imin) pFq({a},{b},z)
 *
 * Only the infinite case with a concrete integer lower bound is handled; every
 * other shape (and a non-rational term ratio, an irrational/complex root, or a
 * divergent p > q+1 series) falls through unevaluated so the cascade can try the
 * next stage or leave Sum[...] held.
 */

#include "sum_internal.h"
#include "symtab.h"
#include "eval.h"
#include "expr.h"
#include "attr.h"
#include "poly.h"
#include "arithmetic.h"   /* is_infinity_sym */
#include "sym_names.h"
#include <string.h>
#include <stdlib.h>
#include <stdint.h>

#define SUM_HG_MAX_DEGREE 256   /* cap on the linear-factor deflation loop */

/* True iff e still contains a HypergeometricPFQ head anywhere (i.e. the pFq did
 * not reduce to a closed form). */
static bool hg_contains_pfq(Expr* e) {
    if (!e || e->type != EXPR_FUNCTION) return false;
    if (e->data.function.head->type == EXPR_SYMBOL
        && e->data.function.head->data.symbol.name == SYM_HypergeometricPFQ) return true;
    if (hg_contains_pfq(e->data.function.head)) return true;
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (hg_contains_pfq(e->data.function.args[i])) return true;
    return false;
}

/* Rewrite Binomial[a, b] -> a! / (b! (a-b)!) throughout e (returns an owned
 * copy; e is not consumed).  Mathilda's Simplify reduces Factorial ratios but
 * not Binomial/Gamma/Pochhammer ratios, so this exposes the rational term ratio
 * of a binomial hypergeometric term (e.g. 2^k/Binomial[2k,k]) to the detector. */
static Expr* hg_expand_binomial(const Expr* e) {
    if (e->type != EXPR_FUNCTION) return expr_copy((Expr*)e);
    const Expr* h = e->data.function.head;
    size_t n = e->data.function.arg_count;
    if (h->type == EXPR_SYMBOL && h->data.symbol.name == SYM_Binomial && n == 2) {
        Expr* a = hg_expand_binomial(e->data.function.args[0]);
        Expr* b = hg_expand_binomial(e->data.function.args[1]);
        Expr* amb = expr_new_function(expr_new_symbol(SYM_Plus),
            (Expr*[]){ expr_copy(a),
                       expr_new_function(expr_new_symbol(SYM_Times),
                           (Expr*[]){ expr_new_integer(-1), expr_copy(b) }, 2) }, 2);
        Expr* fa   = expr_new_function(expr_new_symbol(SYM_Factorial), (Expr*[]){ a }, 1);
        Expr* fb   = expr_new_function(expr_new_symbol(SYM_Factorial), (Expr*[]){ b }, 1);
        Expr* famb = expr_new_function(expr_new_symbol(SYM_Factorial), (Expr*[]){ amb }, 1);
        return expr_new_function(expr_new_symbol(SYM_Times),
            (Expr*[]){ fa,
                       expr_new_function(expr_new_symbol(SYM_Power),
                           (Expr*[]){ fb, expr_new_integer(-1) }, 2),
                       expr_new_function(expr_new_symbol(SYM_Power),
                           (Expr*[]){ famb, expr_new_integer(-1) }, 2) }, 3);
    }
    /* Pochhammer[a, k] = Gamma[a + k]/Gamma[a]: the term ratio of a
     * Pochhammer then simplifies to the rational a + k. */
    if (h->type == EXPR_SYMBOL && h->data.symbol.name == SYM_Pochhammer && n == 2) {
        Expr* a = hg_expand_binomial(e->data.function.args[0]);
        Expr* b = hg_expand_binomial(e->data.function.args[1]);
        Expr* apb = expr_new_function(expr_new_symbol(SYM_Plus), (Expr*[]){ expr_copy(a), b }, 2);
        return expr_new_function(expr_new_symbol(SYM_Times),
            (Expr*[]){ expr_new_function(expr_new_symbol(SYM_Gamma), (Expr*[]){ apb }, 1),
                       expr_new_function(expr_new_symbol(SYM_Power),
                           (Expr*[]){ expr_new_function(expr_new_symbol(SYM_Gamma), (Expr*[]){ a }, 1),
                                      expr_new_integer(-1) }, 2) }, 2);
    }
    Expr** nargs = malloc(sizeof(Expr*) * (n ? n : 1));
    for (size_t i = 0; i < n; i++) nargs[i] = hg_expand_binomial(e->data.function.args[i]);
    Expr* out = expr_new_function(expr_copy((Expr*)h), nargs, n);
    free(nargs);
    return out;
}

/* e is a polynomial in var? */
static bool hg_poly_in(Expr* e, Expr* var) {
    Expr* vars[1] = { var };
    return is_polynomial(e, vars, 1);
}

/* degree of e in var (expanding first). */
static int hg_pdeg(Expr* e, Expr* var) {
    Expr* ea[1] = { expr_copy(e) };
    Expr* ex = sum_eval("Expand", ea, 1);
    int d = get_degree_poly(ex, var);
    expr_free(ex);
    return d;
}

/* Pull the first root r out of a Solve result List[List[Rule[var, r]], ...].
 * Returns an owned copy of r, or NULL if the shape is anything else (e.g. an
 * unsolved Solve[...] or an empty solution set). */
static Expr* hg_first_root(Expr* sol) {
    if (!sol || sol->type != EXPR_FUNCTION) return NULL;
    if (sol->data.function.head->type != EXPR_SYMBOL
        || sol->data.function.head->data.symbol.name != SYM_List
        || sol->data.function.arg_count < 1) return NULL;
    Expr* s0 = sol->data.function.args[0];
    if (s0->type != EXPR_FUNCTION
        || s0->data.function.head->type != EXPR_SYMBOL
        || s0->data.function.head->data.symbol.name != SYM_List
        || s0->data.function.arg_count < 1) return NULL;
    Expr* rule = s0->data.function.args[0];
    if (rule->type != EXPR_FUNCTION
        || rule->data.function.head->type != EXPR_SYMBOL
        || rule->data.function.head->data.symbol.name != SYM_Rule
        || rule->data.function.arg_count != 2) return NULL;
    return expr_copy(rule->data.function.args[1]);
}

/* Factor polynomial P (in var) as lc * prod_i (var + a_i) with each a_i free
 * of var.  On success returns true and fills *params (owned array, *np entries,
 * each owned) and *lc (owned, the remaining constant).  Returns false if P does
 * not split into rational/parametric linear factors (e.g. an irreducible
 * quadratic). */
static bool hg_linear_params(Expr* P, Expr* var,
                             Expr*** params, size_t* np, Expr** lc) {
    size_t cap = 4, n = 0;
    Expr** ps = malloc(sizeof(Expr*) * cap);
    Expr* work = expr_copy(P);

    for (int guard = 0; guard < SUM_HG_MAX_DEGREE; guard++) {
        int d = hg_pdeg(work, var);
        if (d <= 0) { *params = ps; *np = n; *lc = work; return true; }

        Expr* eq = expr_new_function(expr_new_symbol(SYM_Equal),
                       (Expr*[]){ expr_copy(work), sum_int(0) }, 2);
        Expr* sol = sum_eval("Solve", (Expr*[]){ eq, expr_copy(var) }, 2);
        Expr* r = hg_first_root(sol);
        expr_free(sol);
        if (!r || !sum_free_of(r, var)) {
            if (r) expr_free(r);
            goto fail;
        }

        /* param = -r */
        Expr* param = sum_eval("Times", (Expr*[]){ sum_int(-1), expr_copy(r) }, 2);
        if (n == cap) { cap *= 2; ps = realloc(ps, sizeof(Expr*) * cap); }
        ps[n++] = param;

        /* deflate: work = PolynomialQuotient[work, var - r, var] */
        Expr* lin = expr_new_function(expr_new_symbol(SYM_Plus),
                        (Expr*[]){ expr_copy(var),
                                   expr_new_function(expr_new_symbol(SYM_Times),
                                       (Expr*[]){ sum_int(-1), r }, 2) }, 2);
        Expr* q = sum_eval("PolynomialQuotient",
                      (Expr*[]){ expr_copy(work), lin, expr_copy(var) }, 3);
        expr_free(work);
        work = q;
        if (hg_pdeg(work, var) >= d) goto fail;   /* degree must drop */
    }

fail:
    for (size_t i = 0; i < n; i++) expr_free(ps[i]);
    free(ps);
    expr_free(work);
    return false;
}

/* List[elems...] from an owned array (adopts the elements; copies pointers). */
static Expr* hg_make_list(Expr** elems, size_t n) {
    Expr* l = expr_new_function(expr_new_symbol(SYM_List), elems, n);
    return l;
}

Expr* builtin_sum_hypergeometric(Expr* res);

/* The series from the (first nonzero) index imin; see builtin_sum_hypergeometric. */
static Expr* hg_sum_from(Expr* f, Expr* var, Expr* imin, Expr* nterm) {

    /* Term ratio t(k+1)/t(k), reduced to a rational function num/den in var.
     * Binomials are first expanded to factorials so the ratio simplifies. */
    Expr* fx = hg_expand_binomial(f);
    Expr* nv = expr_new_function(expr_new_symbol(SYM_Plus),
                   (Expr*[]){ expr_copy(var), sum_int(1) }, 2);
    Expr* tshift = sum_subst(fx, var, nv);
    expr_free(nv);
    Expr* ratio_raw = expr_new_function(expr_new_symbol(SYM_Times),
                          (Expr*[]){ tshift,
                                     expr_new_function(expr_new_symbol(SYM_Power),
                                         (Expr*[]){ fx, sum_int(-1) }, 2) }, 2);
    Expr* simp = sum_eval("Simplify", (Expr*[]){ ratio_raw }, 1);
    Expr* ratio = sum_eval("Together", (Expr*[]){ simp }, 1);
    Expr* num = sum_eval("Numerator", (Expr*[]){ expr_copy(ratio) }, 1);
    Expr* den = sum_eval("Denominator", (Expr*[]){ ratio }, 1);

    if (!hg_poly_in(num, var) || !hg_poly_in(den, var)) {
        expr_free(num); expr_free(den);
        return NULL;
    }

    /* Reindex to m >= 0: substitute var -> var + imin (no-op when imin == 0). */
    Expr* shift = expr_new_function(expr_new_symbol(SYM_Plus),
                      (Expr*[]){ expr_copy(var), expr_copy(imin) }, 2);
    Expr* num_m = sum_subst(num, var, shift);
    Expr* den_m = sum_subst(den, var, shift);
    expr_free(shift);
    expr_free(num); expr_free(den);

    Expr **A = NULL, **B = NULL, *lcN = NULL, *lcD = NULL;
    size_t p = 0, q = 0;
    if (!hg_linear_params(num_m, var, &A, &p, &lcN)) { expr_free(num_m); expr_free(den_m); return NULL; }
    if (!hg_linear_params(den_m, var, &B, &q, &lcD)) {
        for (size_t i = 0; i < p; i++) expr_free(A[i]);
        free(A); expr_free(lcN); expr_free(num_m); expr_free(den_m);
        return NULL;
    }
    expr_free(num_m); expr_free(den_m);

    /* Upper parameters start as the numerator roots. */
    size_t up_cap = p + 1, up_n = 0;
    Expr** up = malloc(sizeof(Expr*) * (up_cap ? up_cap : 1));
    for (size_t i = 0; i < p; i++) up[up_n++] = A[i];
    free(A);

    /* Lower parameters: the denominator roots minus one canonical (m+1)
     * factorial factor (param == 1).  If none is present, introduce the
     * factorial by adding an upper parameter 1 (the (m+1) in the denominator). */
    size_t lo_cap = q ? q : 1, lo_n = 0;
    Expr** lo = malloc(sizeof(Expr*) * lo_cap);
    Expr* one = sum_int(1);
    bool factorial_removed = false;
    for (size_t j = 0; j < q; j++) {
        if (!factorial_removed && expr_eq(B[j], one)) {
            factorial_removed = true;
            expr_free(B[j]);
        } else {
            lo[lo_n++] = B[j];
        }
    }
    free(B);
    expr_free(one);
    if (!factorial_removed) up[up_n++] = sum_int(1);

    /* z = lcN / lcD. */
    Expr* z = sum_eval("Times",
                  (Expr*[]){ lcN,
                             expr_new_function(expr_new_symbol(SYM_Power),
                                 (Expr*[]){ lcD, sum_int(-1) }, 2) }, 2);

    /* Convergence gate.  p <= q is entire (sums for all z).  p == q+1 converges
     * only for |z| < 1 -- outside that, the elementary closed form is the
     * analytic continuation, NOT the value of the (divergent) series, so we must
     * NOT return it (e.g. Sum[2^k, {k,0,Infinity}] must stay Infinity/unevaluated,
     * not 1/(1-2) = -1).  p > q+1 always diverges.  When z is symbolic the gate
     * passes and the closed form is returned as the conventional CAS answer. */
    /* A finite symbolic upper bound: the sum to imax is the whole series
     * exactly when the series terminates there, i.e. an upper parameter is
     * -(imax - imin) (Pfaff-Saalschutz's -n, Chu-Vandermonde's -r). */
    if (nterm) {
        bool term = false;
        for (size_t i = 0; i < up_n && !term; i++) {
            Expr* d = sum_eval("Simplify", (Expr*[]){ expr_new_function(expr_new_symbol(SYM_Plus),
                                  (Expr*[]){ expr_copy(up[i]), expr_copy(nterm) }, 2) }, 1);
            term = d && d->type == EXPR_INTEGER && d->data.integer == 0;
            if (d) expr_free(d);
        }
        if (!term) {
            for (size_t i = 0; i < up_n; i++) expr_free(up[i]);
            for (size_t j = 0; j < lo_n; j++) expr_free(lo[j]);
            free(up); free(lo); expr_free(z);
            return NULL;
        }
    }
    bool diverges = false, boundary = false;
    if (nterm) {
        /* terminating: a polynomial, no convergence question */
    } else if (up_n > lo_n + 1) {
        diverges = true;
    } else if (up_n == lo_n + 1) {
        Expr* absz = sum_eval("Abs", (Expr*[]){ expr_copy(z) }, 1);
        Expr* nz = sum_eval("N", (Expr*[]){ absz }, 1);
        if (nz->type == EXPR_REAL) {
            if (nz->data.real > 1.0 + 1e-9) {
                diverges = true;                                    /* |z| > 1 */
            } else if (nz->data.real > 1.0 - 1e-9) {
                /* |z| == 1: the series converges only for z != 1 (Dirichlet:
                 * partial sums of z^k stay bounded) AND when the term magnitude
                 * decays, i.e. the growth exponent sigma = Sum(up) - Sum(lo) - 1
                 * is < 0.  Otherwise it diverges and the pFq closed form would be
                 * a spurious analytic continuation. */
                Expr** st = malloc(sizeof(Expr*) * (up_n + lo_n + 1));
                size_t sn = 0;
                for (size_t i = 0; i < up_n; i++) st[sn++] = expr_copy(up[i]);
                for (size_t j = 0; j < lo_n; j++)
                    st[sn++] = expr_new_function(expr_new_symbol(SYM_Times),
                                   (Expr*[]){ sum_int(-1), expr_copy(lo[j]) }, 2);
                st[sn++] = sum_int(-1);
                Expr* tot = expr_new_function(expr_new_symbol(SYM_Plus), st, sn);
                free(st);
                /* With s = Sum(lo) - Sum(up) = -1 - sigma, the series converges
                 * at z == 1 iff Re s > 0 (sigma < -1: Gauss, Dixon, Dougall at 1)
                 * and elsewhere on |z| == 1 iff Re s > -1 (sigma < 0: Kummer at
                 * -1).  A symbolic sigma is left to the assumptions: the boundary
                 * path below still returns only a genuine closed form. */
                Expr* nsig = sum_eval("N", (Expr*[]){ tot }, 1);
                bool z_is_one = (z->type == EXPR_INTEGER && z->data.integer == 1);
                bool numeric = nsig->type == EXPR_REAL;
                double lim = z_is_one ? -1.0 : 0.0;
                bool converges = numeric ? (nsig->data.real < lim - 1e-9) : true;
                expr_free(nsig);
                if (converges) boundary = true;
                else diverges = true;
            }
        }
        expr_free(nz);
    }
    if (diverges) {
        for (size_t i = 0; i < up_n; i++) expr_free(up[i]);
        for (size_t j = 0; j < lo_n; j++) expr_free(lo[j]);
        free(up); free(lo); expr_free(z);
        return NULL;
    }

    /* prefactor t(imin). */
    Expr* prefactor = sum_subst(f, var, imin);

    Expr* upList = hg_make_list(up, up_n); free(up);
    Expr* loList = hg_make_list(lo, lo_n); free(lo);
    Expr* pfq = expr_new_function(expr_new_symbol(SYM_HypergeometricPFQ),
                    (Expr*[]){ upList, loList, z }, 3);
    Expr* result = sum_eval("Times", (Expr*[]){ prefactor, pfq }, 2);

    /* On the |z| == 1 boundary the pFq series is only conditionally convergent
     * (or divergent).  The closed form is safe to return ONLY when HypergeometricPFQ
     * actually reduced it to a genuine closed form (e.g. a very-well-poised
     * Ramanujan 4F3 -> 2/Pi); if it stayed a HypergeometricPFQ, returning it would
     * fabricate the divergent analytic continuation, so we leave Sum[...] held. */
    if (boundary && hg_contains_pfq(result)) {
        expr_free(result);
        return NULL;
    }
    return result;
}

static bool hgr_head(const Expr* e, const char* h) {
    return e && e->type == EXPR_FUNCTION && e->data.function.head->type == EXPR_SYMBOL
        && e->data.function.head->data.symbol.name == h;
}

/* g == c + m v with machine-integer c, m (g evaluated, borrowed). */
static bool hgr_linear(const Expr* g, const Expr* v, int64_t* c, int64_t* m) {
    *c = 0; *m = 0;
    size_t nt = hgr_head(g, SYM_Plus) ? g->data.function.arg_count : 1;
    Expr* const* t = hgr_head(g, SYM_Plus) ? g->data.function.args : (Expr* const*)&g;
    for (size_t i = 0; i < nt; i++) {
        const Expr* ti = t[i];
        if (ti->type == EXPR_INTEGER) *c += ti->data.integer;
        else if (expr_eq(ti, v)) *m += 1;
        else if (hgr_head(ti, SYM_Times) && ti->data.function.arg_count == 2
                 && ti->data.function.args[0]->type == EXPR_INTEGER
                 && expr_eq(ti->data.function.args[1], v))
            *m += ti->data.function.args[0]->data.integer;
        else return false;
    }
    return true;
}

/* A factor Gamma[c + m v]^e (e a machine integer). */
static bool hgr_gamma_factor(const Expr* f, const Expr* v, int64_t* c, int64_t* m, int64_t* e) {
    *e = 1;
    if (hgr_head(f, SYM_Power) && f->data.function.arg_count == 2
        && f->data.function.args[1]->type == EXPR_INTEGER) {
        *e = f->data.function.args[1]->data.integer;
        f = f->data.function.args[0];
    }
    return hgr_head(f, SYM_Gamma) && f->data.function.arg_count == 1
        && hgr_linear(f->data.function.args[0], v, c, m);
}

/* Gamma[c + x] (x = m v, m > 0) as Gamma[1 + x] times a rational factor. */
static Expr* hgr_shift_gamma(int64_t c, int64_t m, const Expr* v) {
    Expr* x = expr_new_function(expr_new_symbol(SYM_Times),
                                (Expr*[]){ expr_new_integer(m), expr_copy((Expr*)v) }, 2);
    int64_t lo = c >= 1 ? 1 : c, hi = c >= 1 ? c - 1 : 0;
    size_t nf = (size_t)(hi >= lo ? hi - lo + 1 : 0) + 1, w = 0;
    Expr** f = malloc(sizeof(Expr*) * nf);
    f[w++] = expr_new_function(expr_new_symbol(SYM_Gamma), (Expr*[]){
        expr_new_function(expr_new_symbol(SYM_Plus), (Expr*[]){ expr_new_integer(1), expr_copy(x) }, 2) }, 1);
    for (int64_t j = lo; j <= hi; j++) {
        Expr* xj = expr_new_function(expr_new_symbol(SYM_Plus), (Expr*[]){ expr_copy(x), expr_new_integer(j) }, 2);
        f[w++] = c >= 1 ? xj : expr_new_function(expr_new_symbol(SYM_Power),
                                   (Expr*[]){ xj, expr_new_integer(-1) }, 2);
    }
    expr_free(x);
    Expr* r = expr_new_function(expr_new_symbol(SYM_Times), f, w);
    free(f);
    return r;
}

/* Pole-pair regularisation for a closed form in the summation bounds' symbols.
 * The bounds of an integer range are integers, so Gamma[c + m n] with integer
 * c and negative integer m is a pole at every admissible n -- Dixon's
 * Sum_{k=-n}^{n} (-1)^k Binomial[2n, k+n]^3 comes out of the 3F2 theorem as
 * Gamma[1-n] Gamma[1+3n]/((-1)^n Gamma[1-2n] Gamma[1+n]^2 Gamma[1+2n]),
 * finite only as a ratio.  Reflection writes each pole factor as
 * Pi/(Sin[Pi (c + m n)] Gamma[1 - c - m n]).  Near an integer n,
 * Sin[Pi (c + m n)] ~ (-1)^(c + m n) m Pi eps, so when the Sin factors balance
 * (as many above as below the fraction bar) the eps and Pi cancel and the ratio
 * is prod (m (-1)^(c + m n))^(-e).  The signs are collected into one exponent
 * C + M n, reduced mod 2 (n is an integer), and every remaining Gamma[c' + m' n]
 * (m' > 0) is shifted to Gamma[1 + m' n]: (3n)!/(n!)^3 for Dixon.
 * Unbalanced poles, and any non-product form, are left as they were.
 * The factors are classified structurally in C; r is consumed. */
static Expr* hg_regularise_poles(Expr* r, Expr* vars) {
    if (!hgr_head(vars, SYM_List) || vars->data.function.arg_count != 1) return r;
    const Expr* v = vars->data.function.args[0];
    size_t nf = hgr_head(r, SYM_Times) ? r->data.function.arg_count : 1;
    Expr* const* f = hgr_head(r, SYM_Times) ? r->data.function.args : &r;
    int64_t balance = 0, signC = 0, signM = 0, npole = 0;
    for (size_t i = 0; i < nf; i++) {
        int64_t c, m, e;
        if (hgr_gamma_factor(f[i], v, &c, &m, &e) && m < 0) {
            balance += e; npole++;
        }
    }
    if (npole == 0 || balance != 0) return r;
    Expr** out = malloc(sizeof(Expr*) * (nf + 2));
    size_t w = 0;
    for (size_t i = 0; i < nf; i++) {
        int64_t c, m, e;
        if (!hgr_gamma_factor(f[i], v, &c, &m, &e)) { out[w++] = expr_copy(f[i]); continue; }
        Expr* g;
        if (m < 0) {
            /* (Pi/(m (-1)^(c+mn) Gamma[1-c-mn]))^e; Pi^e cancels over the
             * balanced set and the sign goes into (signC, signM). */
            signC += e * c; signM += e * m;
            g = expr_new_function(expr_new_symbol(SYM_Times), (Expr*[]){
                expr_new_integer(m),
                hgr_shift_gamma(1 - c, -m, v) }, 2);
            e = -e;
        } else if (m > 0) {
            g = hgr_shift_gamma(c, m, v);
        } else {
            out[w++] = expr_copy(f[i]);
            continue;
        }
        out[w++] = expr_new_function(expr_new_symbol(SYM_Power), (Expr*[]){ g, expr_new_integer(e) }, 2);
    }
    if (signC % 2) out[w++] = expr_new_integer(-1);
    if (signM % 2)
        out[w++] = expr_new_function(expr_new_symbol(SYM_Power),
                       (Expr*[]){ expr_new_integer(-1), expr_copy((Expr*)v) }, 2);
    Expr* prod = expr_new_function(expr_new_symbol(SYM_Times), out, w);
    free(out);
    Expr* res = sum_eval("Together", (Expr*[]){ prod }, 1);
    if (!res) return r;
    expr_free(r);
    return res;
}

Expr* builtin_sum_hypergeometric(Expr* res) {
    Expr *f, *var, *imin, *imax;
    bool definite;
    if (!sum_stage_args(res, &f, &var, &imin, &imax, &definite)) return NULL;
    if (!definite) return NULL;

    /* A symbolic lower bound is shifted to 0 (Dixon's Sum_{k=-n}^{n}):
     * Sum_{k=a}^{b} f(k) = Sum_{j=0}^{b-a} f(j + a). */
    if (imin->type != EXPR_INTEGER) {
        if (is_infinity_sym(imin)) return NULL;
        Expr* sh = expr_new_function(expr_new_symbol(SYM_Plus),
                       (Expr*[]){ expr_copy(var), expr_copy(imin) }, 2);
        Expr* f2 = sum_subst(f, var, sh);
        expr_free(sh);
        Expr* hi = is_infinity_sym(imax) ? expr_copy(imax)
                 : sum_eval("Simplify", (Expr*[]){ expr_new_function(expr_new_symbol(SYM_Plus),
                       (Expr*[]){ expr_copy(imax), expr_new_function(expr_new_symbol(SYM_Times),
                           (Expr*[]){ sum_int(-1), expr_copy(imin) }, 2) }, 2) }, 1);
        Expr* zero = sum_int(0);
        Expr* r = hg_sum_from(f2, var, zero, is_infinity_sym(hi) ? NULL : hi);
        expr_free(zero); expr_free(f2); expr_free(hi);
        if (r) {
            Expr* vs = sum_eval("Variables", (Expr*[]){ expr_copy(imin) }, 1);
            r = hg_regularise_poles(r, vs);
            expr_free(vs);
        }
        return r;
    }
    /* A finite upper bound must be symbolic (a concrete one is enumerated
     * by the finite-sum stage); the series must terminate there. */
    Expr* nterm = NULL;
    if (!is_infinity_sym(imax)) {
        if (imax->type == EXPR_INTEGER) return NULL;
        nterm = sum_eval("Simplify", (Expr*[]){ expr_new_function(expr_new_symbol(SYM_Plus),
                   (Expr*[]){ expr_copy(imax), expr_new_function(expr_new_symbol(SYM_Times),
                       (Expr*[]){ sum_int(-1), expr_copy(imin) }, 2) }, 2) }, 1);
    }

    /* The series is read off from its first term t(imin): the closed form is
     * t(imin) pFq(..).  When the leading terms VANISH -- k^m x^k/k! at k = 0,
     * the Dobinski/Touchard sums -- that product is identically 0 (the bug was
     * Sum[k^10/k!, {k, 0, Infinity}] == 0, true value 115975 E).  Start instead
     * at the first index whose term is nonzero; the vanishing terms contribute
     * nothing.  A term that is zero for every index scanned is left alone. */
    Expr* imin_eff = NULL;
    for (int64_t j = imin->data.integer; j < imin->data.integer + 64; j++) {
        Expr* jj = sum_int(j);
        Expr* tj = sum_subst(f, var, jj);
        bool zero = tj && tj->type == EXPR_INTEGER && tj->data.integer == 0;
        if (tj) expr_free(tj);
        if (!zero) {
            if (j != imin->data.integer) imin_eff = jj; else expr_free(jj);
            break;
        }
        expr_free(jj);
    }
    if (imin_eff && nterm)
        nterm = sum_eval("Plus", (Expr*[]){ nterm, sum_int(imin->data.integer - imin_eff->data.integer) }, 2);
    Expr* r = hg_sum_from(f, var, imin_eff ? imin_eff : imin, nterm);
    if (imin_eff) expr_free(imin_eff);
    if (nterm) expr_free(nterm);
    return r;
}

void sum_hypergeometric_init(void) {
    symtab_add_builtin("Sum`Hypergeometric", builtin_sum_hypergeometric);
}
