/*
 * sum_zetaseries.c -- Sum`ZetaSeries: sums of functions of Zeta[a k + b].
 *
 * Closes  Sum[c(k) (Zeta[a k + b] - 1), {k, k0, Infinity}]  (and the (-1)^k
 * variant) by interchanging the order of summation.  Since
 *
 *     Zeta[a k + b] - 1 = Sum_{n>=2} n^{-(a k + b)},
 *
 * we have
 *
 *     Sum_{k>=k0} c(k) (Zeta[a k+b] - 1)
 *        = Sum_{n>=2} n^{-b} ( Sum_{k>=k0} c(k) (n^{-a})^k ),
 *
 * where the inner geometric series (|n^{-a}| < 1 for n >= 2) closes to a
 * rational function R(n), and the outer Sum over n closes via Sum`Rational.
 * E.g.  Sum[Zeta[k]-1, {k,2,Inf}] = Sum_{n>=2} 1/(n(n-1)) = 1.
 *
 * Scope: the subtracted part must be exactly the constant -1 (removing the
 * n = 1 term, so the interchanged series starts at n0 = 2 and converges);
 * c(k) must be a k-free constant or that times (-1)^k.  Anything else falls
 * through (returns NULL) rather than risk a divergent interchange.
 *
 * This unblocks Product[Exp[Zeta[k]-1], {k,2,Inf}] = E via Product`LogSum.
 *
 * Memory contract: takes ownership of res but must not free it; returns an
 * owned closed form or NULL to fall through.
 */

#include "sum_internal.h"
#include "symtab.h"
#include "eval.h"
#include "expr.h"
#include "attr.h"
#include "arithmetic.h"   /* is_infinity_sym */
#include "sym_names.h"
#include <string.h>
#include <stdlib.h>

/* e is the one-argument call Zeta[arg] with arg linear (a var + b), a a
 * positive integer and b free of var?  On success a_out and b_out are owned. */
static bool zs_zeta_linear(Expr* e, Expr* var, Expr** a_out, Expr** b_out) {
    if (e->type != EXPR_FUNCTION) return false;
    Expr* h = e->data.function.head;
    if (h->type != EXPR_SYMBOL || h->data.symbol.name != SYM_Zeta) return false;
    if (e->data.function.arg_count != 1) return false;
    Expr* arg = e->data.function.args[0];
    if (sum_free_of(arg, var)) return false;

    Expr* zero = sum_int(0), *one = sum_int(1);
    Expr* b  = sum_subst(arg, var, zero);          /* arg /. var->0 */
    Expr* a1 = sum_subst(arg, var, one);           /* arg /. var->1 */
    expr_free(zero); expr_free(one);
    Expr* a = sum_sub(a1, b);                       /* a = arg(1) - arg(0) */
    expr_free(a1);

    if (!(a->type == EXPR_INTEGER && a->data.integer > 0) || !sum_free_of(b, var)) {
        expr_free(a); expr_free(b); return false;
    }
    /* verify arg == a var + b. */
    Expr* lin = sum_eval("Plus", (Expr*[]){ expr_copy(b),
                    sum_eval("Times", (Expr*[]){ expr_copy(a), expr_copy(var) }, 2) }, 2);
    Expr* diff = sum_sub(arg, lin); expr_free(lin);
    Expr* ds = sum_eval("Simplify", (Expr*[]){ diff }, 1);
    bool linear = (ds->type == EXPR_INTEGER && ds->data.integer == 0);
    expr_free(ds);
    if (!linear) { expr_free(a); expr_free(b); return false; }

    *a_out = a; *b_out = b;
    return true;
}

/* Within a Plus, find the (single) Zeta[linear] child; return its index and
 * a,b, or -1. */
static int zs_find_zeta_in_plus(Expr* p, Expr* var, Expr** a_out, Expr** b_out) {
    for (size_t i = 0; i < p->data.function.arg_count; i++) {
        if (zs_zeta_linear(p->data.function.args[i], var, a_out, b_out))
            return (int)i;
    }
    return -1;
}

static bool zs_is_plus(Expr* e) {
    return e->type == EXPR_FUNCTION && e->data.function.head->type == EXPR_SYMBOL
        && e->data.function.head->data.symbol.name == SYM_Plus;
}

Expr* builtin_sum_zetaseries(Expr* res);

Expr* builtin_sum_zetaseries(Expr* res) {
    Expr *f, *var, *imin, *imax;
    bool definite;
    if (!sum_stage_args(res, &f, &var, &imin, &imax, &definite)) return NULL;
    if (!definite) return NULL;
    if (!is_infinity_sym(imax)) return NULL;
    if (imin->type != EXPR_INTEGER) return NULL;
    int64_t k0 = imin->data.integer;

    Expr* fc = expr_copy(f);
    Expr* fn = evaluate(fc);
    expr_free(fc);   /* evaluate() borrows; free the copy we made */

    /* Split fn = C * G, where G is the additive group holding Zeta[linear]. */
    Expr* G = NULL;        /* alias into fn */
    Expr* C = NULL;        /* owned coefficient, or NULL == 1 */
    Expr* a = NULL, *b = NULL;   /* owned */
    Expr* rest0 = NULL;    /* alias: the single subtracted term (must be -1) */

    if (zs_is_plus(fn)) {
        int zi = zs_find_zeta_in_plus(fn, var, &a, &b);
        if (zi >= 0) {
            G = fn;
            /* rest = all Plus args except zi; require exactly one, == -1. */
            if (fn->data.function.arg_count == 2)
                rest0 = fn->data.function.args[(zi == 0) ? 1 : 0];
        }
    } else if (fn->type == EXPR_FUNCTION
               && fn->data.function.head->type == EXPR_SYMBOL
               && fn->data.function.head->data.symbol.name == SYM_Times) {
        /* find the factor that is a Plus[..Zeta[linear]..]; rest = other factors. */
        int gfi = -1;
        for (size_t i = 0; i < fn->data.function.arg_count; i++) {
            Expr* fac = fn->data.function.args[i];
            if (zs_is_plus(fac)) {
                int zi = zs_find_zeta_in_plus(fac, var, &a, &b);
                if (zi >= 0) {
                    gfi = (int)i; G = fac;
                    if (fac->data.function.arg_count == 2)
                        rest0 = fac->data.function.args[(zi == 0) ? 1 : 0];
                    break;
                }
            }
        }
        if (gfi >= 0) {
            /* C = product of the other factors. */
            size_t m = fn->data.function.arg_count;
            if (m == 2) {
                C = expr_copy(fn->data.function.args[(gfi == 0) ? 1 : 0]);
            } else {
                Expr** parts = malloc(sizeof(Expr*) * (m - 1));
                size_t j = 0;
                for (size_t i = 0; i < m; i++)
                    if ((int)i != gfi) parts[j++] = expr_copy(fn->data.function.args[i]);
                C = sum_eval("Times", parts, m - 1);
                free(parts);
            }
        }
    }

    if (!G || !a || !rest0) {
        expr_free(fn); if (a) expr_free(a); if (b) expr_free(b); if (C) expr_free(C);
        return NULL;
    }

    /* Subtracted part must be exactly -1 (removes n = 1 -> interchange from n0 = 2). */
    Expr* rp1 = sum_eval("Plus", (Expr*[]){ expr_copy(rest0), sum_int(1) }, 1 + 1);
    bool rest_is_m1 = (rp1->type == EXPR_INTEGER && rp1->data.integer == 0);
    expr_free(rp1);
    if (!rest_is_m1) {
        expr_free(fn); expr_free(a); expr_free(b); if (C) expr_free(C);
        return NULL;
    }

    /* Coefficient c(k): k-free constant, optionally times (-1)^k. */
    Expr* coef = NULL; bool alt = false;
    if (!C) {
        coef = sum_int(1);
    } else if (sum_free_of(C, var)) {
        coef = C; C = NULL;
    } else {
        /* try C = coef * (-1)^k. */
        Expr* invsign = expr_new_function(expr_new_symbol(SYM_Power),
                            (Expr*[]){ expr_new_function(expr_new_symbol(SYM_Power),
                                          (Expr*[]){ sum_int(-1), expr_copy(var) }, 2),
                                       sum_int(-1) }, 2);   /* ((-1)^k)^-1 */
        Expr* cred = sum_eval("Simplify",
                        (Expr*[]){ expr_new_function(expr_new_symbol(SYM_Times),
                                       (Expr*[]){ expr_copy(C), invsign }, 2) }, 1);
        if (sum_free_of(cred, var)) { coef = cred; alt = true; }
        else { expr_free(cred); }
        expr_free(C); C = NULL;
    }
    expr_free(fn);   /* G, rest0 aliased into fn; done with them */

    if (!coef) { expr_free(a); expr_free(b); return NULL; }

    int64_t aval = a->data.integer;

    /* Build R(n) = coef * n^{-b} * ( ratio^{k0} / (1 - ratio) ),  ratio = (+/-) n^{-a}. */
    Expr* n = expr_new_symbol("Sum`ZetaSeries`n");
    Expr* base = expr_new_function(expr_new_symbol(SYM_Power),
                     (Expr*[]){ expr_copy(n), sum_int(-aval) }, 2);   /* n^{-a} */
    Expr* ratio = alt
        ? expr_new_function(expr_new_symbol(SYM_Times), (Expr*[]){ sum_int(-1), base }, 2)
        : base;                                                       /* adopts base */
    Expr* num = expr_new_function(expr_new_symbol(SYM_Power),
                    (Expr*[]){ expr_copy(ratio), sum_int(k0) }, 2);   /* ratio^k0 */
    Expr* den = expr_new_function(expr_new_symbol(SYM_Plus),
                    (Expr*[]){ sum_int(1),
                               expr_new_function(expr_new_symbol(SYM_Times),
                                   (Expr*[]){ sum_int(-1), ratio }, 2) }, 2); /* 1 - ratio */
    Expr* inner = expr_new_function(expr_new_symbol(SYM_Times),
                      (Expr*[]){ num,
                                 expr_new_function(expr_new_symbol(SYM_Power),
                                     (Expr*[]){ den, sum_int(-1) }, 2) }, 2);
    Expr* nmb = expr_new_function(expr_new_symbol(SYM_Power),
                    (Expr*[]){ expr_copy(n),
                               expr_new_function(expr_new_symbol(SYM_Times),
                                   (Expr*[]){ sum_int(-1), expr_copy(b) }, 2) }, 2); /* n^{-b} */
    Expr* Rraw = expr_new_function(expr_new_symbol(SYM_Times),
                     (Expr*[]){ coef, nmb, inner }, 3);
    expr_free(a); expr_free(b);

    Expr* R = sum_eval("Together", (Expr*[]){ Rraw }, 1);   /* adopts Rraw */

    /* Outer sum: Sum[R, {n, 2, Infinity}]  (n0 = 2). */
    Expr* spec = expr_new_function(expr_new_symbol(SYM_List),
                     (Expr*[]){ expr_copy(n), sum_int(2), expr_new_symbol(SYM_Infinity) }, 3);
    expr_free(n);
    Expr* out = sum_eval("Sum", (Expr*[]){ R, spec }, 2);   /* adopts R, spec */

    if (out->type == EXPR_FUNCTION && out->data.function.head->type == EXPR_SYMBOL
        && out->data.function.head->data.symbol.name == SYM_Sum) {
        expr_free(out); return NULL;   /* outer sum did not close */
    }
    return out;
}

void sum_zetaseries_init(void) {
    symtab_add_builtin("Sum`ZetaSeries", builtin_sum_zetaseries);
    symtab_get_def("Sum`ZetaSeries")->attributes |= ATTR_PROTECTED | ATTR_READPROTECTED;
    symtab_set_docstring("Sum`ZetaSeries",
        "Sum`ZetaSeries[f, i, k0, Infinity] evaluates sums of c[i](Zeta[a i+b]-1) "
        "by interchanging the order of summation into a rational sum over the "
        "zeta argument.");
}
