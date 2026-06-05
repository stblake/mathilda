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

/* base^var with exponent structurally equal to var, base free of var and not
 * the trivial 0/1?  If so return base (copied), else NULL. */
static Expr* var_power_base(Expr* e, Expr* var) {
    if (e->type != EXPR_FUNCTION) return NULL;
    if (e->data.function.head->type != EXPR_SYMBOL) return NULL;
    if (e->data.function.head->data.symbol != SYM_Power) return NULL;
    if (e->data.function.arg_count != 2) return NULL;
    Expr* base = e->data.function.args[0];
    Expr* exp  = e->data.function.args[1];
    if (!expr_eq(exp, var)) return NULL;
    if (!sum_free_of(base, var)) return NULL;
    if (base->type == EXPR_INTEGER &&
        (base->data.integer == 0 || base->data.integer == 1)) return NULL;
    return expr_copy(base);
}

/* Multiply acc (owned) by base^(-var) and evaluate, cancelling base^var. */
static Expr* divide_off(Expr* acc, Expr* base, Expr* var) {
    Expr* inv = expr_new_function(expr_new_symbol("Power"),
                    (Expr*[]){ expr_copy(base),
                               expr_new_function(expr_new_symbol("Times"),
                                   (Expr*[]){ sum_int(-1), expr_copy(var) }, 2) }, 2);
    Expr* prod = expr_new_function(expr_new_symbol("Times"),
                     (Expr*[]){ acc, inv }, 2);
    Expr* r = evaluate(prod);
    expr_free(prod);
    return r;
}

/* Detect f = p(i) * r^i.  Sets *r_out (product of the bases of all r^i factors)
 * and *p_out (f with those factors divided off, which must be a polynomial in
 * var -- each base is cancelled individually so q1^i q2^i works).  Returns
 * false if f has no exponential factor or the remainder is not polynomial. */
static bool find_geometric(Expr* f, Expr* var, Expr** r_out, Expr** p_out) {
    Expr* direct = var_power_base(f, var);
    if (direct) {
        *r_out = direct;
        *p_out = divide_off(expr_copy(f), direct, var);
        return true;
    }
    if (!(f->type == EXPR_FUNCTION && f->data.function.head->type == EXPR_SYMBOL
          && f->data.function.head->data.symbol == SYM_Times))
        return false;

    Expr* r = NULL;
    Expr* p = expr_copy(f);
    for (size_t k = 0; k < f->data.function.arg_count; k++) {
        Expr* b = var_power_base(f->data.function.args[k], var);
        if (!b) continue;
        p = divide_off(p, b, var);                 /* cancel b^i out of p */
        if (!r) { r = b; }
        else {
            Expr* tr = expr_new_function(expr_new_symbol("Times"),
                           (Expr*[]){ r, b }, 2);
            r = evaluate(tr);
            expr_free(tr);
        }
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
        Expr* pw = expr_new_function(expr_new_symbol("Power"),
                       (Expr*[]){ expr_copy(var), sum_int(k) }, 2);
        qterms[k] = expr_new_function(expr_new_symbol("Times"),
                       (Expr*[]){ expr_copy(csym[k]), pw }, 2);
    }
    /* q = Plus[qterms...] (takes ownership of qterms elements). */
    Expr* q = expr_new_function(expr_new_symbol("Plus"), qterms, m + 1);
    free(qterms);

    /* E = r * q(var+1) - q(var) - p */
    Expr* vp1 = expr_new_function(expr_new_symbol("Plus"),
                    (Expr*[]){ expr_copy(var), sum_int(1) }, 2);
    Expr* qshift = sum_subst(q, var, vp1);
    expr_free(vp1);
    Expr* rqshift = expr_new_function(expr_new_symbol("Times"),
                        (Expr*[]){ expr_copy(r), qshift }, 2);
    Expr* Eraw = expr_new_function(expr_new_symbol("Plus"),
                     (Expr*[]){ rqshift,
                                expr_new_function(expr_new_symbol("Times"),
                                    (Expr*[]){ sum_int(-1), expr_copy(q) }, 2),
                                expr_new_function(expr_new_symbol("Times"),
                                    (Expr*[]){ sum_int(-1), expr_copy(pexp) }, 2) }, 3);
    expr_free(pexp);
    Expr* Eargs[1] = { Eraw };
    Expr* E = sum_eval("Expand", Eargs, 1);

    /* coeffs = CoefficientList[E, var] */
    Expr* clargs[2] = { expr_copy(E), expr_copy(var) };
    Expr* cl = sum_eval("CoefficientList", clargs, 2);
    expr_free(E);

    bool ok = true;
    Expr* rules = expr_new_function(expr_new_symbol("List"), (Expr*[]){ NULL }, 0);

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
        Expr* negB = expr_new_function(expr_new_symbol("Times"),
                         (Expr*[]){ sum_int(-1), B }, 2);
        Expr* invA = expr_new_function(expr_new_symbol("Power"),
                         (Expr*[]){ A, sum_int(-1) }, 2);
        Expr* valraw = expr_new_function(expr_new_symbol("Times"),
                           (Expr*[]){ negB, invA }, 2);
        Expr* vargs[1] = { valraw };
        Expr* val = sum_eval("Cancel", vargs, 1);

        Expr* rule = expr_new_function(expr_new_symbol("Rule"),
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

    /* Infinite (convergence-dependent) sums are deferred. */
    if (definite && imax->type == EXPR_SYMBOL && imax->data.symbol == SYM_Infinity)
        return NULL;

    Expr *r = NULL, *p = NULL;
    if (!find_geometric(f, var, &r, &p)) return NULL;

    Expr* vars[1] = { var };
    if (!is_polynomial(p, vars, 1)) { expr_free(p); expr_free(r); return NULL; }

    Expr* q = solve_q(p, r, var);
    expr_free(p);
    if (!q) { expr_free(r); return NULL; }

    /* Combine the rational coefficient over a common denominator.  We must NOT
     * run Together/Factor on the full F = q * r^var: both loop on the symbolic
     * exponential r^var.  q itself is rational in (i, r) and safe to Together. */
    Expr* qt_arg[1] = { q };
    Expr* qt = sum_eval("Together", qt_arg, 1);

    /* F(var) = qt(var) * r^var */
    Expr* rv = expr_new_function(expr_new_symbol("Power"),
                   (Expr*[]){ expr_copy(r), expr_copy(var) }, 2);
    Expr* Fraw = expr_new_function(expr_new_symbol("Times"),
                     (Expr*[]){ expr_copy(qt), rv }, 2);
    Expr* F = evaluate(Fraw);
    expr_free(Fraw);
    expr_free(r);

    Expr* out;
    if (!definite) {
        out = F;
        F = NULL;
    } else {
        Expr* up = expr_new_function(expr_new_symbol("Plus"),
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
