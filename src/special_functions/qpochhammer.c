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
 * Infinite form (2 args): evaluated for machine-real a, q with |q| < 1 by the
 * shared kernel sf_machine_qpochhammer (also the NDArray / Compile one), and
 * for complex or arbitrary-precision a, q through arb at the call's precision.
 * Symbolic or |q| >= 1 inputs stay unevaluated.
 *
 * A negative integer n gives 1 / prod_{k=1}^{|n|} (1 - a q^-k); a non-integer
 * n, with an inexact argument, (a; q)_inf / (a q^n; q)_inf.
 *
 * Memory: honours the builtin ownership contract (never frees res).
 */

#include <float.h>
#include <math.h>
#include "qpochhammer.h"
#include "eval.h"
#include "symtab.h"
#include "attr.h"
#include "arithmetic.h"   /* is_rational */
#include "sym_names.h"
#include "flint_num_bridge.h"   /* complex / MPFR and non-integer n */
#include "sf_machine.h"         /* sf_machine_qpochhammer: the shared kernel */
#include "sf_scalar.h"          /* sf_machine_scalar_args */
#include <gmp.h>
#include <math.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdlib.h>

/* The most factors the finite form expands into an explicit product, for
 * either sign of n.  Past it an exact or symbolic call stays in closed form
 * (QPochhammer[a, q, 10^7] used to build and evaluate a 10^7-factor product,
 * minutes and gigabytes -- audit QPochhammer#6); an inexact one goes through
 * arb's (a; q)_inf / (a q^n; q)_inf, as a non-integer n does. */
#define QP_EXPAND_MAX 100000

/* (1 - a q^k) as an owned, unevaluated tree. */
static Expr* qfactor(Expr* a, Expr* q, int64_t k) {
    Expr* qk = expr_new_function(expr_new_symbol(SYM_Power),
                   (Expr*[]){ expr_copy(q), expr_new_integer(k) }, 2);
    Expr* aqk = expr_new_function(expr_new_symbol(SYM_Times),
                    (Expr*[]){ expr_new_integer(-1), expr_copy(a), qk }, 3);
    return expr_new_function(expr_new_symbol(SYM_Plus),
               (Expr*[]){ expr_new_integer(1), aqk }, 2);
}

Expr* builtin_qpochhammer(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    Expr** a = res->data.function.args;

    if (argc == 3) {
        Expr *av = a[0], *qv = a[1], *nv = a[2];
        int64_t n, d;
        if (!(is_rational(nv, &n, &d) && d == 1)) {
            /* A non-integer n: (a; q)_inf / (a q^n; q)_inf, for an inexact call. */
            if (flint_num_call_inexact(res))
                return flint_num_for_call(flint_num_qpochhammer(av, qv, nv), res);
            return NULL;                                 /* symbolic n stays */
        }
        if (n > QP_EXPAND_MAX || n < -QP_EXPAND_MAX) {
            if (flint_num_call_inexact(res))
                return flint_num_for_call(flint_num_qpochhammer(av, qv, nv), res);
            return NULL;                                 /* closed form stays */
        }
        if (n < 0) {
            /* (a; q)_-m = 1 / prod_{k=1}^{m} (1 - a q^-k). */
            int64_t m = -n;
            Expr** facs = malloc(sizeof(Expr*) * (size_t)m);
            for (int64_t k = 0; k < m; k++) facs[k] = qfactor(av, qv, -(k + 1));
            Expr* prod = expr_new_function(expr_new_symbol(SYM_Times), facs, (size_t)m);
            free(facs);
            return eval_and_free(expr_new_function(expr_new_symbol(SYM_Power),
                       (Expr*[]){ prod, expr_new_integer(-1) }, 2));
        }
        if (n == 0) return expr_new_integer(1);
        Expr** facs = malloc(sizeof(Expr*) * (size_t)n);
        for (int64_t k = 0; k < n; k++) facs[k] = qfactor(av, qv, k);
        Expr* prod = expr_new_function(expr_new_symbol(SYM_Times), facs, (size_t)n);
        free(facs);
        Expr* out = evaluate(prod);
        expr_free(prod);
        return out;
    }

    if (argc == 1) {
        /* QPochhammer[q] is the Euler function (q; q)_inf = QPochhammer[q, q]
         * (QPochhammer#1 of the 2026-11 audit): an inexact q takes the
         * two-argument kernels; an exact q stays symbolic, as the
         * two-argument (a; q)_inf does. */
        if (!flint_num_call_inexact(res)) return NULL;
        return eval_and_free(expr_new_function(expr_new_symbol("QPochhammer"),
                   (Expr*[]){ expr_copy(a[0]), expr_copy(a[0]) }, 2));
    }

    if (argc == 2) {
        /* (0; q)_inf = 1: every factor is 1 (QPochhammer#3). */
        if (a[0]->type == EXPR_INTEGER && a[0]->data.integer == 0) return expr_new_integer(1);
        /* Machine-real a, q: the kernel the NDArray map and Compile[] share.
         * An all-exact (a;q)_inf stays symbolic (no exact closed form).  A
         * product that underflows the double range (a subnormal or a zero
         * the factors do not explain) goes to arb, whose 53-bit result keeps
         * the exponent -- the numeric tower's policy for a machine value
         * outside [DBL_MIN, DBL_MAX] (QPochhammer#5). */
        SFMachineArgs ma;
        if (sf_machine_scalar_args(res, &ma) && !ma.any_complex) {
            double v;
            if (sf_machine_qpochhammer(ma.re[0], ma.re[1], &v) &&
                (fabs(v) >= DBL_MIN || ma.re[0] == 1.0))
                return expr_new_real(v);
        }
        /* Complex or arbitrary-precision a, q: the same product through arb,
         * at the call's precision (|q| >= 1 diverges and stays unevaluated). */
        if (flint_num_call_inexact(res))
            return flint_num_for_call(flint_num_qpochhammer(a[0], a[1], NULL), res);
        return NULL;
    }

    return NULL;
}

void qpochhammer_init(void) {
    symtab_add_builtin("QPochhammer", builtin_qpochhammer);
    symtab_get_def("QPochhammer")->attributes |=
        (ATTR_LISTABLE | ATTR_NUMERICFUNCTION | ATTR_PROTECTED);
}
