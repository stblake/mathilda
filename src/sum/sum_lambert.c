/* sum_lambert.c — Sum`Lambert: the Lambert series Sum_{n>=1} n^s/(E^(2 Pi n) - 1)
 * for s = 1 (mod 4), the values the Eisenstein series give at tau = i.
 *
 * With q = e^{2 pi i tau}, the weight-2m Eisenstein series is
 *     E_2m(tau) = 1 - (4m/B_2m) Sum_{n>=1} n^{2m-1} q^n/(1 - q^n),
 * and Sum n^{2m-1} q^n/(1-q^n) = Sum n^{2m-1}/(q^{-n} - 1).  At tau = i,
 * q^{-n} = E^(2 Pi n), and the modular relation E_2m(-1/tau) = tau^2m E_2m(tau)
 * forces E_2m(i) = 0 whenever i^2m = -1, i.e. 2m = 2 (mod 4), m >= 3:
 *
 *     Sum_{n>=1} n^{4k+1}/(E^(2 Pi n) - 1) = B_{4k+2}/(8k + 4),   k >= 1
 *
 * (Ramanujan; k = 3 gives 1/24 for n^13).  Weight 2 is quasi-modular,
 * E_2(-1/tau) = tau^2 E_2(tau) + 12 tau/(2 pi i), so E_2(i) = 3/pi and
 *
 *     Sum_{n>=1} n/(E^(2 Pi n) - 1) = 1/24 - 1/(8 Pi)        (Schlomilch).
 *
 * The term must be exactly c n^s/(E^(2 Pi n) - 1) (c free of n) summed from 1
 * to Infinity; anything else declines.  The factors are classified
 * structurally; the denominator is confirmed by an exact difference.
 */

#include "sum_internal.h"
#include "arithmetic.h"   /* is_infinity_sym */
#include "symtab.h"
#include "eval.h"
#include "expr.h"
#include "sym_names.h"
#include <stdlib.h>

Expr* builtin_sum_lambert(Expr* res);

/* True iff d evaluates to E^(2 Pi n) - 1 (borrowed d, n). */
static bool lam_is_den(const Expr* d, const Expr* n) {
    Expr* ref = expr_new_function(expr_new_symbol(SYM_Plus), (Expr*[]){
        expr_new_integer(-1),
        expr_new_function(expr_new_symbol(SYM_Power), (Expr*[]){ expr_new_symbol(SYM_E),
            expr_new_function(expr_new_symbol(SYM_Times), (Expr*[]){ expr_new_integer(2),
                expr_new_symbol(SYM_Pi), expr_copy((Expr*)n) }, 3) }, 2) }, 2);
    Expr* diff = expr_new_function(expr_new_symbol(SYM_Plus), (Expr*[]){ expr_copy((Expr*)d),
        expr_new_function(expr_new_symbol(SYM_Times), (Expr*[]){ expr_new_integer(-1), ref }, 2) }, 2);
    Expr* r = evaluate(diff);
    expr_free(diff);
    bool z = r && r->type == EXPR_INTEGER && r->data.integer == 0;
    if (r) expr_free(r);
    return z;
}

Expr* builtin_sum_lambert(Expr* res) {
    Expr *f, *var, *imin, *imax;
    bool definite;
    if (!sum_stage_args(res, &f, &var, &imin, &imax, &definite)) return NULL;
    if (!definite || !is_infinity_sym(imax) || var->type != EXPR_SYMBOL) return NULL;
    if (!(imin->type == EXPR_INTEGER && imin->data.integer == 1)) return NULL;
    /* Classify the factors of the term: c * n^s * (E^(2 Pi n) - 1)^-1. */
    bool times = f->type == EXPR_FUNCTION && f->data.function.head->type == EXPR_SYMBOL
              && f->data.function.head->data.symbol.name == SYM_Times;
    size_t nf = times ? f->data.function.arg_count : 1;
    Expr* const* fac = times ? f->data.function.args : &f;
    Expr** cf = malloc(sizeof(Expr*) * nf);
    size_t nc = 0;
    int64_t s = 0;
    int nden = 0, npow = 0;
    bool ok = true;
    for (size_t i = 0; i < nf && ok; i++) {
        const Expr* g = fac[i];
        if (sum_free_of((Expr*)g, var)) { cf[nc++] = expr_copy((Expr*)g); continue; }
        if (expr_eq(g, var)) { s = 1; npow++; continue; }
        bool pw = g->type == EXPR_FUNCTION && g->data.function.head->type == EXPR_SYMBOL
               && g->data.function.head->data.symbol.name == SYM_Power
               && g->data.function.arg_count == 2;
        if (pw && expr_eq(g->data.function.args[0], var)
            && g->data.function.args[1]->type == EXPR_INTEGER) {
            s = g->data.function.args[1]->data.integer; npow++; continue;
        }
        if (pw && g->data.function.args[1]->type == EXPR_INTEGER
            && g->data.function.args[1]->data.integer == -1
            && lam_is_den(g->data.function.args[0], var)) { nden++; continue; }
        ok = false;
    }
    Expr* out = NULL;
    if (ok && nden == 1 && npow == 1 && s >= 1 && s % 4 == 1) {
        Expr* val;
        if (s == 1)   /* 1/24 - 1/(8 Pi) */
            val = expr_new_function(expr_new_symbol(SYM_Plus), (Expr*[]){ make_rational(1, 24),
                expr_new_function(expr_new_symbol(SYM_Times), (Expr*[]){ make_rational(-1, 8),
                    expr_new_function(expr_new_symbol(SYM_Power), (Expr*[]){
                        expr_new_symbol(SYM_Pi), expr_new_integer(-1) }, 2) }, 2) }, 2);
        else          /* B_{s+1} / (2 (s+1)) */
            val = expr_new_function(expr_new_symbol(SYM_Times), (Expr*[]){
                make_rational(1, 2 * (s + 1)),
                expr_new_function(expr_new_symbol(SYM_BernoulliB),
                    (Expr*[]){ expr_new_integer(s + 1) }, 1) }, 2);
        cf[nc++] = val;
        Expr* prod = expr_new_function(expr_new_symbol(SYM_Times), cf, nc);
        out = evaluate(prod);
        expr_free(prod);
    } else {
        for (size_t i = 0; i < nc; i++) expr_free(cf[i]);
    }
    free(cf);
    return out;
}

void sum_lambert_init(void) {
    symtab_add_builtin("Sum`Lambert", builtin_sum_lambert);
}
