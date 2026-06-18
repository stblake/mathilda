/* Mathilda -- QPochhammer, the q-Pochhammer symbol (q-shifted factorial).
 *
 *   QPochhammer[a, q, n] = prod_{k=0}^{n-1} (1 - a q^k)
 *   QPochhammer[a, q]     = prod_{k=0}^{Infinity} (1 - a q^k)   ((a;q)_inf)
 *
 * Finite form (3 args): for a non-negative integer n the product is built and
 * handed to the evaluator, which reduces it exactly for exact a, q and at
 * machine / MPFR precision for inexact a, q (so N works through it).  A
 * symbolic / non-integer n is left unevaluated -- which is exactly what Product
 * relies on to emit QPochhammer[a, q, n] as a closed form.
 *
 * Infinite form (2 args): evaluated for machine-real a, q with |q| < 1 by
 * accumulating factors until they fall below machine epsilon.  Symbolic or
 * |q| >= 1 inputs stay unevaluated.
 *
 * Memory: honours the builtin ownership contract (never frees res).
 */

#include "qpochhammer.h"
#include "eval.h"
#include "symtab.h"
#include "attr.h"
#include "arithmetic.h"   /* is_rational */
#include "sym_names.h"
#include <gmp.h>
#include <math.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdlib.h>

/* (1 - a q^k) as an owned, unevaluated tree. */
static Expr* qfactor(Expr* a, Expr* q, int64_t k) {
    Expr* qk = expr_new_function(expr_new_symbol(SYM_Power),
                   (Expr*[]){ expr_copy(q), expr_new_integer(k) }, 2);
    Expr* aqk = expr_new_function(expr_new_symbol(SYM_Times),
                    (Expr*[]){ expr_new_integer(-1), expr_copy(a), qk }, 3);
    return expr_new_function(expr_new_symbol(SYM_Plus),
               (Expr*[]){ expr_new_integer(1), aqk }, 2);
}

/* Best-effort machine double for an exact/inexact real Expr; returns false for
 * non-real / symbolic inputs. */
static bool to_double(const Expr* e, double* out) {
    int64_t n, d;
    if (e->type == EXPR_REAL) { *out = e->data.real; return true; }
#ifdef USE_MPFR
    if (e->type == EXPR_MPFR) { *out = mpfr_get_d(e->data.mpfr, MPFR_RNDN); return true; }
#endif
    if (is_rational((Expr*)e, &n, &d) && d != 0) { *out = (double)n / (double)d; return true; }
    return false;
}

Expr* builtin_qpochhammer(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    Expr** a = res->data.function.args;

    if (argc == 3) {
        Expr *av = a[0], *qv = a[1], *nv = a[2];
        int64_t n, d;
        if (!(is_rational(nv, &n, &d) && d == 1)) return NULL;  /* symbolic n stays */
        if (n < 0) return NULL;
        if (n == 0) return expr_new_integer(1);
        Expr** facs = malloc(sizeof(Expr*) * (size_t)n);
        for (int64_t k = 0; k < n; k++) facs[k] = qfactor(av, qv, k);
        Expr* prod = expr_new_function(expr_new_symbol(SYM_Times), facs, (size_t)n);
        free(facs);
        Expr* out = evaluate(prod);
        expr_free(prod);
        return out;
    }

    if (argc == 2) {
        double av, qv;
        if (!to_double(a[0], &av) || !to_double(a[1], &qv)) return NULL;
        if (!(fabs(qv) < 1.0)) return NULL;                     /* diverges / undefined */
        /* Only commit to a machine value when at least one input is inexact;
         * an all-exact (a;q)_inf is left symbolic (no exact closed form). */
        if (a[0]->type != EXPR_REAL && a[1]->type != EXPR_REAL
#ifdef USE_MPFR
            && a[0]->type != EXPR_MPFR && a[1]->type != EXPR_MPFR
#endif
           ) return NULL;
        double prod = 1.0, qk = 1.0;
        for (int k = 0; k < 100000; k++) {
            double factor = 1.0 - av * qk;
            prod *= factor;
            qk *= qv;
            if (fabs(av * qk) < 1e-18) break;
        }
        return expr_new_real(prod);
    }

    return NULL;
}

void qpochhammer_init(void) {
    symtab_add_builtin("QPochhammer", builtin_qpochhammer);
    symtab_get_def("QPochhammer")->attributes |=
        (ATTR_LISTABLE | ATTR_NUMERICFUNCTION | ATTR_PROTECTED);
}
