/*
 * sum_logzeta.c -- Sum`LogZeta: log-weighted zeta series.
 *
 *   Sum[c * Log[k] / k^s, {k, 1, Infinity}]  =  -c * Zeta'[s].
 *
 * The system has no symbolic Zeta' (D[Zeta[s],s] is an inert Derivative), so a
 * usable elementary closed form is emitted only for s == 2, via the Glaisher
 * bridge
 *
 *   -Zeta'[2] = (Pi^2/6) * (12 Log[Glaisher] - EulerGamma - Log[2 Pi]).
 *
 * This unblocks Product[k^(1/k^2), {k,1,Inf}] = Exp[-Zeta'[2]] through the
 * Product`LogSum bridge.  For any other s the stage returns NULL (leaves the
 * Sum held) rather than fabricate a Zeta' symbol.
 *
 * Memory contract: takes ownership of res but must not free it; returns an
 * owned closed form or NULL to fall through.
 */

#include "sum_internal.h"
#include "symtab.h"
#include "eval.h"
#include "expr.h"
#include "attr.h"
#include "sym_names.h"
#include <string.h>
#include <stdlib.h>

static bool is_infinity_sym(const Expr* e) {
    return e->type == EXPR_SYMBOL && e->data.symbol == SYM_Infinity;
}

Expr* builtin_sum_logzeta(Expr* res) {
    Expr *f, *var, *imin, *imax;
    bool definite;
    if (!sum_stage_args(res, &f, &var, &imin, &imax, &definite)) return NULL;
    if (!definite || !is_infinity_sym(imax)) return NULL;
    if (!(imin->type == EXPR_INTEGER && imin->data.integer == 1)) return NULL;

    /* Body must carry a Log[var] factor. */
    Expr* logk = expr_new_function(expr_new_symbol("Log"),
                     (Expr*[]){ expr_copy(var) }, 1);
    if (sum_free_of(f, logk)) { expr_free(logk); return NULL; }

    /* h = Simplify[f / Log[var]]; must be free of Log[var] (a pure power part). */
    Expr* invlog = expr_new_function(expr_new_symbol(SYM_Power),
                       (Expr*[]){ expr_copy(logk), expr_new_integer(-1) }, 2);
    Expr* prod = expr_new_function(expr_new_symbol(SYM_Times),
                     (Expr*[]){ expr_copy(f), invlog }, 2);        /* adopts invlog */
    Expr* h = sum_eval("Simplify", (Expr*[]){ prod }, 1);          /* adopts prod */
    if (!sum_free_of(h, logk)) { expr_free(h); expr_free(logk); return NULL; }
    expr_free(logk);

    /* s = -var * D[Log[h], var]  (for h = c*var^(-s) this is exactly s). */
    Expr* logh = expr_new_function(expr_new_symbol("Log"),
                     (Expr*[]){ expr_copy(h) }, 1);
    Expr* dlogh = sum_eval("D", (Expr*[]){ logh, expr_copy(var) }, 2);  /* adopts logh */
    Expr* smul = expr_new_function(expr_new_symbol(SYM_Times),
                     (Expr*[]){ expr_new_integer(-1), expr_copy(var), dlogh }, 3); /* adopts dlogh */
    Expr* s = sum_eval("Simplify", (Expr*[]){ smul }, 1);          /* adopts smul */

    /* Accept any integer s >= 2 (convergent).  s == 2 gets an elementary
     * Glaisher closed form; higher s emit the inert -c Zeta'[s]. */
    bool s_ok = (s->type == EXPR_INTEGER && s->data.integer >= 2 &&
                 sum_free_of(s, var));
    if (!s_ok) { expr_free(s); expr_free(h); return NULL; }
    int64_t sval = s->data.integer;

    /* c = Simplify[h * var^s]  (the constant coefficient, must be free of var). */
    Expr* vs = expr_new_function(expr_new_symbol(SYM_Power),
                   (Expr*[]){ expr_copy(var), expr_copy(s) }, 2);
    Expr* hc = expr_new_function(expr_new_symbol(SYM_Times),
                   (Expr*[]){ expr_copy(h), vs }, 2);              /* adopts vs */
    Expr* c = sum_eval("Simplify", (Expr*[]){ hc }, 1);            /* adopts hc */
    expr_free(h);
    if (!sum_free_of(c, var)) { expr_free(c); expr_free(s); return NULL; }

    Expr* negZs;
    if (sval == 2) {
        /* -Zeta'[2] = (Pi^2/6)(12 Log[Glaisher] - EulerGamma - Log[2 Pi]). */
        Expr* pi2 = expr_new_function(expr_new_symbol(SYM_Power),
                        (Expr*[]){ expr_new_symbol("Pi"), expr_new_integer(2) }, 2);
        Expr* pi2over6 = expr_new_function(expr_new_symbol(SYM_Times),
                             (Expr*[]){ expr_new_function(expr_new_symbol(SYM_Rational),
                                            (Expr*[]){ expr_new_integer(1),
                                                       expr_new_integer(6) }, 2),
                                        pi2 }, 2);                     /* adopts pi2 */
        Expr* logG = expr_new_function(expr_new_symbol("Log"),
                         (Expr*[]){ expr_new_symbol("Glaisher") }, 1);
        Expr* t1 = expr_new_function(expr_new_symbol(SYM_Times),
                       (Expr*[]){ expr_new_integer(12), logG }, 2);    /* adopts logG */
        Expr* t2 = expr_new_function(expr_new_symbol(SYM_Times),
                       (Expr*[]){ expr_new_integer(-1), expr_new_symbol("EulerGamma") }, 2);
        Expr* twopi = expr_new_function(expr_new_symbol(SYM_Times),
                          (Expr*[]){ expr_new_integer(2), expr_new_symbol("Pi") }, 2);
        Expr* log2pi = expr_new_function(expr_new_symbol("Log"),
                           (Expr*[]){ twopi }, 1);                     /* adopts twopi */
        Expr* t3 = expr_new_function(expr_new_symbol(SYM_Times),
                       (Expr*[]){ expr_new_integer(-1), log2pi }, 2);  /* adopts log2pi */
        Expr* bracket = expr_new_function(expr_new_symbol(SYM_Plus),
                            (Expr*[]){ t1, t2, t3 }, 3);               /* adopts t1,t2,t3 */
        negZs = expr_new_function(expr_new_symbol(SYM_Times),
                    (Expr*[]){ pi2over6, bracket }, 2);               /* adopts both */
    } else {
        /* General s: -Zeta'[s] = -Derivative[1][Zeta][s] (inert first zeta
         * derivative, matching the object D[Zeta[s],s] produces). */
        Expr* op = expr_new_function(expr_new_symbol("Derivative"),
                       (Expr*[]){ expr_new_integer(1) }, 1);          /* Derivative[1] */
        Expr* op_z = expr_new_function(op,
                         (Expr*[]){ expr_new_symbol(SYM_Zeta) }, 1);  /* Derivative[1][Zeta] */
        Expr* zprime = expr_new_function(op_z,
                           (Expr*[]){ expr_copy(s) }, 1);             /* Derivative[1][Zeta][s] */
        negZs = expr_new_function(expr_new_symbol(SYM_Times),
                    (Expr*[]){ expr_new_integer(-1), zprime }, 2);    /* -Zeta'[s] */
    }
    expr_free(s);

    /* result = c * (-Zeta'[s]). */
    Expr* result_expr = expr_new_function(expr_new_symbol(SYM_Times),
                            (Expr*[]){ c, negZs }, 2);             /* adopts c, negZs */
    Expr* out = evaluate(result_expr);
    expr_free(result_expr);
    return out;
}

void sum_logzeta_init(void) {
    symtab_add_builtin("Sum`LogZeta", builtin_sum_logzeta);
    symtab_get_def("Sum`LogZeta")->attributes |= ATTR_PROTECTED | ATTR_READPROTECTED;
    symtab_set_docstring("Sum`LogZeta",
        "Sum`LogZeta[f, i, 1, Infinity] evaluates Sum[c Log[i]/i^s] = -c Zeta'[s] "
        "for integer s>=2. s==2 yields an elementary Glaisher-constant closed form; "
        "higher s yield the inert -c Derivative[1][Zeta][s].");
}
