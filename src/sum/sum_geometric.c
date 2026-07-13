/*
 * sum_geometric.c -- Sum`Geometric: summation of p(i) * r^i.
 *
 * Handles f = p(i) r^i where p is a polynomial in i and r is free of i (the
 * base of one or more factors r^i).  The antidifference has the same shape,
 * F(i) = q(i) r^i with deg q = deg p, because
 *
 *     Delta(q(i) r^i) = r^i ( r q(i+1) - q(i) ).
 *
 * So q is the polynomial solving  r q(i+1) - q(i) = p(i).  Writing
 * q = sum_k c_k i^k, this is a triangular linear system in the c_k (the i^m
 * coefficient fixes c_m, the next fixes c_{m-1}, ...), solved here by reading
 * off each coefficient with CoefficientList and back-substituting -- no
 * Bernoulli numbers and no general linear solver required.
 *
 *   Sum`Geometric[f, i]              -> F(i)                 (indefinite)
 *   Sum`Geometric[f, i, imin, imax]  -> F(imax+1) - F(imin)  (definite, finite)
 *
 * Pure geometric (p constant) gives r^i/(r-1), matching Sum[a^i,i] = a^i/(-1+a).
 * Infinite bounds (convergence-dependent) are deferred and fall through.
 */

#include "sum_internal.h"
#include "symtab.h"
#include "eval.h"
#include "expr.h"
#include "attr.h"
#include "poly.h"
#include "sym_names.h"
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>

static int geom_counter = 0;   /* unique coefficient-symbol suffix */

/* Return the per-step geometric ratio of a single factor e(var), or NULL if e
 * is not a (constant times) geometric factor in var.
 *
 * A factor g(var) is geometric with ratio rho iff  g(var+1)/g(var)  is free of
 * var -- then g(var) = C rho^var for some var-free C.  This ratio test is
 * representation-independent, so it catches every surface form of the same
 * factor alike:  r^i,  r^(-i),  (r^i)^(-1),  (1/r)^i,  r^(a i + b),  2^(-k)
 * (held as (2^k)^(-1)), ...  Polynomial factors fail the test -- e.g.
 * (k+1)^2/k^2 still depends on k -- and are left in the polynomial part by the
 * caller.  The var-free C stays in the polynomial part too, since divide_off
 * removes exactly rho^var.
 *
 * The plain case e = Power[base, var] (base free of var) is short-circuited to
 * avoid the substitute/evaluate round-trip on the hot path. */
static Expr* var_power_base(Expr* e, Expr* var) {
    /* Constant factor: not geometric. */
    if (sum_free_of(e, var)) return NULL;

    /* Fast path: Power[base, var] with base free of var and not 0/1. */
    if (e->type == EXPR_FUNCTION && e->data.function.head->type == EXPR_SYMBOL
        && e->data.function.head->data.symbol.name == SYM_Power
        && e->data.function.arg_count == 2
        && expr_eq(e->data.function.args[1], var)
        && sum_free_of(e->data.function.args[0], var)) {
        Expr* base = e->data.function.args[0];
        if (!(base->type == EXPR_INTEGER &&
              (base->data.integer == 0 || base->data.integer == 1)))
            return expr_copy(base);
    }

    /* General ratio test: rho = g(var+1) / g(var). */
    Expr* vp1 = expr_new_function(expr_new_symbol(SYM_Plus),
                    (Expr*[]){ expr_copy(var), sum_int(1) }, 2);
    Expr* eshift = sum_subst(e, var, vp1);            /* g(var+1) */
    expr_free(vp1);
    Expr* einv = expr_new_function(expr_new_symbol(SYM_Power),
                     (Expr*[]){ expr_copy(e), sum_int(-1) }, 2);
    Expr* ratio_raw = expr_new_function(expr_new_symbol(SYM_Times),
                          (Expr*[]){ eshift, einv }, 2);   /* consumes both */
    Expr* ratio = evaluate(ratio_raw);
    expr_free(ratio_raw);
    if (!sum_free_of(ratio, var) ||
        (ratio->type == EXPR_INTEGER &&
         (ratio->data.integer == 0 || ratio->data.integer == 1))) {
        expr_free(ratio); return NULL;
    }
    return ratio;
}

/* Detect f = p(i) * r^i.  Sets *r_out (product of the ratios of all geometric
 * factors) and *p_out (the polynomial part, which must be a polynomial in var).
 * Returns false if f has no geometric factor.
 *
 * Each factor g that is geometric (var_power_base returns its ratio rho) equals
 * C * rho^var for a var-free constant C = g|_{var->0}; we accumulate rho into r
 * and multiply C into the polynomial part.  Multiplying (rather than dividing
 * the whole f by a reconstructed rho^var) sidesteps base-representation
 * mismatches -- e.g. the factor 2^(-k) has ratio 1/2, and 2^(-k)/(1/2)^k would
 * not cancel structurally though it is 1.  q1^i q2^i still works: each power is
 * its own geometric factor. */
static bool find_geometric(Expr* f, Expr* var, Expr** r_out, Expr** p_out) {
    bool is_times = (f->type == EXPR_FUNCTION
                     && f->data.function.head->type == EXPR_SYMBOL
                     && f->data.function.head->data.symbol.name == SYM_Times);
    size_t n = is_times ? f->data.function.arg_count : 1;

    Expr* r = NULL;                 /* product of ratios */
    Expr* p = sum_int(1);           /* accumulated polynomial part */
    for (size_t k = 0; k < n; k++) {
        Expr* g = is_times ? f->data.function.args[k] : f;
        Expr* rho = var_power_base(g, var);
        Expr* factor;               /* what to fold into the polynomial part */
        if (rho) {
            if (!r) { r = rho; }
            else {
                Expr* tr = expr_new_function(expr_new_symbol(SYM_Times),
                               (Expr*[]){ r, rho }, 2);
                r = evaluate(tr);
                expr_free(tr);
            }
            /* C = g |_{var -> 0}, the var-free constant with g = C rho^var. */
            Expr* zero = sum_int(0);
            factor = sum_subst(g, var, zero);
            expr_free(zero);
        } else {
            factor = expr_copy(g);
        }
        Expr* tp = expr_new_function(expr_new_symbol(SYM_Times),
                       (Expr*[]){ p, factor }, 2);
        p = evaluate(tp);
        expr_free(tp);
    }
    if (!r) { expr_free(p); return false; }
    *r_out = r;
    *p_out = p;
    return true;
}

/* Solve r q(i+1) - q(i) = p for the polynomial q (degree = deg p); return q,
 * or NULL on failure (r == 1, or p not a polynomial). */
static Expr* solve_q(Expr* p, Expr* r, Expr* var) {
    Expr* pe[1] = { expr_copy(p) };
    Expr* pexp = sum_eval("Expand", pe, 1);
    Expr* vars[1] = { var };
    if (!is_polynomial(pexp, vars, 1)) { expr_free(pexp); return NULL; }
    int m = get_degree_poly(pexp, var);
    if (m < 0) m = 0;

    /* q = sum_k c_k var^k with fresh coefficient symbols. */
    int id = ++geom_counter;
    char name[64];
    Expr** csym = malloc(sizeof(Expr*) * (m + 1));
    Expr** qterms = malloc(sizeof(Expr*) * (m + 1));
    for (int k = 0; k <= m; k++) {
        snprintf(name, sizeof name, "Sum`g`%d`%d", id, k);
        csym[k] = expr_new_symbol(name);
        Expr* pw = expr_new_function(expr_new_symbol(SYM_Power),
                       (Expr*[]){ expr_copy(var), sum_int(k) }, 2);
        qterms[k] = expr_new_function(expr_new_symbol(SYM_Times),
                       (Expr*[]){ expr_copy(csym[k]), pw }, 2);
    }
    /* q = Plus[qterms...] (takes ownership of qterms elements). */
    Expr* q = expr_new_function(expr_new_symbol(SYM_Plus), qterms, m + 1);
    free(qterms);

    /* E = r * q(var+1) - q(var) - p */
    Expr* vp1 = expr_new_function(expr_new_symbol(SYM_Plus),
                    (Expr*[]){ expr_copy(var), sum_int(1) }, 2);
    Expr* qshift = sum_subst(q, var, vp1);
    expr_free(vp1);
    Expr* rqshift = expr_new_function(expr_new_symbol(SYM_Times),
                        (Expr*[]){ expr_copy(r), qshift }, 2);
    Expr* Eraw = expr_new_function(expr_new_symbol(SYM_Plus),
                     (Expr*[]){ rqshift,
                                expr_new_function(expr_new_symbol(SYM_Times),
                                    (Expr*[]){ sum_int(-1), expr_copy(q) }, 2),
                                expr_new_function(expr_new_symbol(SYM_Times),
                                    (Expr*[]){ sum_int(-1), expr_copy(pexp) }, 2) }, 3);
    expr_free(pexp);
    Expr* Eargs[1] = { Eraw };
    Expr* E = sum_eval("Expand", Eargs, 1);

    /* coeffs = CoefficientList[E, var] */
    Expr* clargs[2] = { expr_copy(E), expr_copy(var) };
    Expr* cl = sum_eval("CoefficientList", clargs, 2);
    expr_free(E);

    bool ok = true;
    Expr* rules = expr_new_function(expr_new_symbol(SYM_List), (Expr*[]){ NULL }, 0);

    size_t ncoeff = (cl->type == EXPR_FUNCTION) ? cl->data.function.arg_count : 0;
    /* Walk coefficients high degree -> low; each is linear in its own c_k. */
    for (int k = m; k >= 0 && ok; k--) {
        Expr* ck = csym[k];
        Expr* coeff_raw = ((size_t)k < ncoeff)
            ? expr_copy(cl->data.function.args[k])
            : sum_int(0);
        /* resolve already-known higher c's */
        Expr* rr[2] = { coeff_raw, expr_copy(rules) };
        Expr* coeff = sum_eval("ReplaceRepeated", rr, 2);

        Expr* zero = sum_int(0);
        Expr* one  = sum_int(1);
        Expr* B = sum_subst(coeff, ck, zero);          /* coeff at c_k = 0 */
        Expr* C1 = sum_subst(coeff, ck, one);          /* coeff at c_k = 1 */
        expr_free(zero); expr_free(one);
        expr_free(coeff);
        Expr* A = sum_sub(C1, B);                      /* A = dcoeff/dc_k    */
        expr_free(C1);

        bool A_zero = (A->type == EXPR_INTEGER && A->data.integer == 0);
        if (A_zero) {
            /* c_k absent: equation must already be satisfiable; treat c_k = 0. */
            expr_free(A); expr_free(B);
            continue;
        }
        /* c_k = -B / A */
        Expr* negB = expr_new_function(expr_new_symbol(SYM_Times),
                         (Expr*[]){ sum_int(-1), B }, 2);
        Expr* invA = expr_new_function(expr_new_symbol(SYM_Power),
                         (Expr*[]){ A, sum_int(-1) }, 2);
        Expr* valraw = expr_new_function(expr_new_symbol(SYM_Times),
                           (Expr*[]){ negB, invA }, 2);
        Expr* vargs[1] = { valraw };
        Expr* val = sum_eval("Cancel", vargs, 1);

        Expr* rule = expr_new_function(expr_new_symbol(SYM_Rule),
                         (Expr*[]){ expr_copy(ck), val }, 2);
        Expr* aargs[2] = { rules, rule };
        rules = sum_eval("Append", aargs, 2);
    }
    expr_free(cl);

    Expr* result = NULL;
    if (ok) {
        Expr* rr[2] = { expr_copy(q), expr_copy(rules) };
        result = sum_eval("ReplaceRepeated", rr, 2);
    }
    expr_free(rules);
    expr_free(q);
    for (int k = 0; k <= m; k++) expr_free(csym[k]);
    free(csym);
    return result;
}

Expr* builtin_sum_geometric(Expr* res) {
    Expr *f, *var, *imin, *imax;
    bool definite;
    if (!sum_stage_args(res, &f, &var, &imin, &imax, &definite)) return NULL;

    bool infinite = (definite && imax->type == EXPR_SYMBOL
                     && imax->data.symbol.name == SYM_Infinity);

    /* Canonicalise the (held) summand so exponential factors are in base^exp
     * form: the parser leaves 1/2^k as (2^k)^(-1), whose base 2^k contains the
     * variable, and divide_off then cannot cancel it back out of the polynomial
     * part.  evaluate rewrites it to 2^(-k) (a var-free base), which cancels
     * cleanly.  This does not expand or run Together, so it is safe on the
     * symbolic r^var forms the memory warns about. */
    Expr* fc = expr_copy(f);
    Expr* fn = evaluate(fc);
    expr_free(fc);   /* evaluate() borrows; free the copy we made */

    Expr *r = NULL, *p = NULL;
    if (!find_geometric(fn, var, &r, &p)) { expr_free(fn); return NULL; }
    expr_free(fn);

    Expr* vars[1] = { var };
    if (!is_polynomial(p, vars, 1)) { expr_free(p); expr_free(r); return NULL; }

    /* Infinite range converges iff |r| < 1: the summand p(k) r^k is a polynomial
     * times a decaying exponential, so the boundary term q(n+1) r^(n+1) -> 0 and
     * the sum collapses to -F(imin).  An undecidable ratio (symbolic r) or
     * |r| >= 1 (including the divergent r = -1 and r = 1 edges) is deferred:
     * we return NULL so the sum is left unevaluated rather than fabricated. */
    if (infinite) {
        Expr* lt = expr_new_function(expr_new_symbol(SYM_Less),
                       (Expr*[]){ expr_new_function(expr_new_symbol(SYM_Abs),
                                      (Expr*[]){ expr_copy(r) }, 1),
                                  sum_int(1) }, 2);
        Expr* conv = evaluate(lt);
        expr_free(lt);
        bool converges = (conv->type == EXPR_SYMBOL && conv->data.symbol.name == SYM_True);
        expr_free(conv);
        if (!converges) { expr_free(p); expr_free(r); return NULL; }
    }

    Expr* q = solve_q(p, r, var);
    expr_free(p);
    if (!q) { expr_free(r); return NULL; }

    /* Combine the rational coefficient over a common denominator.  We must NOT
     * run Together/Factor on the full F = q * r^var: both loop on the symbolic
     * exponential r^var.  q itself is rational in (i, r) and safe to Together. */
    Expr* qt_arg[1] = { q };
    Expr* qt = sum_eval("Together", qt_arg, 1);

    /* F(var) = qt(var) * r^var */
    Expr* rv = expr_new_function(expr_new_symbol(SYM_Power),
                   (Expr*[]){ expr_copy(r), expr_copy(var) }, 2);
    Expr* Fraw = expr_new_function(expr_new_symbol(SYM_Times),
                     (Expr*[]){ expr_copy(qt), rv }, 2);
    Expr* F = evaluate(Fraw);
    expr_free(Fraw);
    expr_free(r);

    Expr* out;
    if (!definite) {
        out = F;
        F = NULL;
    } else if (infinite) {
        /* |r| < 1 was verified above, so F(n+1) = q(n+1) r^(n+1) -> 0 and the
         * sum from imin to infinity is  lim_n [F(n+1) - F(imin)] = -F(imin). */
        Expr* Flo = sum_subst(F, var, imin);
        Expr* neg = expr_new_function(expr_new_symbol(SYM_Times),
                        (Expr*[]){ sum_int(-1), Flo }, 2);   /* Flo consumed */
        out = evaluate(neg);
        expr_free(neg);
    } else {
        Expr* up = expr_new_function(expr_new_symbol(SYM_Plus),
                       (Expr*[]){ expr_copy(imax), sum_int(1) }, 2);
        Expr* Fhi = sum_subst(F, var, up);
        Expr* Flo = sum_subst(F, var, imin);
        expr_free(up);
        out = sum_sub(Fhi, Flo);   /* difference; no Together (r^pow present) */
        expr_free(Fhi); expr_free(Flo);
    }
    expr_free(qt);
    if (F) expr_free(F);
    return out;
}

void sum_geometric_init(void) {
    symtab_add_builtin("Sum`Geometric", builtin_sum_geometric);
    symtab_get_def("Sum`Geometric")->attributes |= ATTR_PROTECTED | ATTR_READPROTECTED;
    symtab_set_docstring("Sum`Geometric",
        "Sum`Geometric[f, i] gives the indefinite sum of f = p(i) r^i where p is "
        "a polynomial in i and r is free of i; Sum`Geometric[f, i, imin, imax] "
        "gives the finite definite sum. Returns unevaluated otherwise.");
}
