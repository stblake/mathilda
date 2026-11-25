/* Mathilda -- sub-machine precision for special functions that compute at
 * no fewer than 53 bits.
 *
 * Precision contagion says f[1.5`3] is a 3-digit number, as Sin[1.5`3] is.
 * Several special functions floor their working precision at 53 bits (their
 * kernels are machine or arb evaluators), and then build the result at that
 * floor: CosIntegral[1.5`3] came back as a MachinePrecision Real (CosIntegral#1,
 * CoshIntegral#1, Erfi#1, HurwitzZeta#3, LogIntegral#1,
 * MathieuCharacteristicB#3, SinhIntegral#1 of the 2026-11 audit).
 *
 * Rather than re-plumb each kernel, the heads listed below are wrapped once,
 * after registration: when every inexact argument is arbitrary precision and
 * the least of them carries p < 53 bits (a machine Real anywhere makes the
 * result machine, so no wrap applies), the original builtin runs unchanged and
 * its numeric result -- a Real, an MPFR, or a Complex of them -- is re-tagged
 * at p bits.  Any other call goes straight to the original builtin; the cost
 * of the wrapper is one table scan and an argument scan.
 */

#include "expr.h"
#include "symtab.h"
#include "numeric.h"
#include "arithmetic.h"        /* is_complex, make_complex */
#include "sym_intern.h"
#include <stddef.h>

void sf_lowprec_init(void);

#ifdef USE_MPFR
#include <mpfr.h>

/* The heads whose kernels floor at 53 bits and build the result there. */
static const char* const LOWPREC_HEADS[] = {
    "CosIntegral", "CoshIntegral", "SinhIntegral", "Erfi",
    "HurwitzZeta", "LogIntegral", "MathieuCharacteristicB",
};
#define N_LOWPREC (sizeof LOWPREC_HEADS / sizeof LOWPREC_HEADS[0])

static const char* g_names[N_LOWPREC];     /* interned head names */
static BuiltinFunc g_orig[N_LOWPREC];

/* The sub-machine precision p (< 53 bits) that `res`'s arguments impose, or 0
 * when none does (no inexact argument, a machine one, or all at >= 53). */
static long lowprec_bits(const Expr* res) {
    long lo = 0;
    for (size_t i = 0; i < res->data.function.arg_count; i++) {
        const Expr* a = res->data.function.args[i];
        long b = numeric_min_inexact_bits(a);
        if (b == 0) continue;                      /* exact or symbolic */
        if (a->type == EXPR_REAL) return 0;
        if (a->type == EXPR_FUNCTION) {            /* Complex: a machine part wins */
            Expr *re, *im;
            if (is_complex((Expr*)a, &re, &im) &&
                (re->type == EXPR_REAL || im->type == EXPR_REAL))
                return 0;
        }
        if (lo == 0 || b < lo) lo = b;
    }
    return (lo > 0 && lo < 53) ? lo : 0;
}

/* `e` at p bits when it is an inexact real or complex number; otherwise `e`
 * unchanged.  Consumes e. */
static Expr* lowprec_retag(Expr* e, long p) {
    if (!e) return NULL;
    if (e->type == EXPR_REAL) {
        Expr* r = expr_new_mpfr_from_d(e->data.real, (mpfr_prec_t)p);
        expr_free(e);
        return r;
    }
    if (e->type == EXPR_MPFR) {
        if ((long)mpfr_get_prec(e->data.mpfr) <= p) return e;
        Expr* r = expr_new_mpfr_bits((mpfr_prec_t)p);
        mpfr_set(r->data.mpfr, e->data.mpfr, MPFR_RNDN);
        expr_free(e);
        return r;
    }
    Expr *re, *im;
    if (e->type == EXPR_FUNCTION && is_complex(e, &re, &im) &&
        (re->type == EXPR_REAL || re->type == EXPR_MPFR) &&
        (im->type == EXPR_REAL || im->type == EXPR_MPFR)) {
        Expr* r = make_complex(lowprec_retag(expr_copy(re), p), lowprec_retag(expr_copy(im), p));
        expr_free(e);
        return r;
    }
    return e;
}

static Expr* lowprec_dispatch(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.head->type != EXPR_SYMBOL) return NULL;
    const char* h = res->data.function.head->data.symbol.name;
    BuiltinFunc orig = NULL;
    for (size_t k = 0; k < N_LOWPREC; k++)
        if (g_names[k] == h) { orig = g_orig[k]; break; }
    if (!orig) return NULL;
    long p = lowprec_bits(res);
    Expr* r = orig(res);
    return p ? lowprec_retag(r, p) : r;
}

void sf_lowprec_init(void) {
    for (size_t k = 0; k < N_LOWPREC; k++) {
        SymbolDef* d = symtab_get_def(LOWPREC_HEADS[k]);
        if (!d->builtin_func || d->builtin_func == lowprec_dispatch) continue;
        g_names[k] = intern_symbol(LOWPREC_HEADS[k]);
        g_orig[k] = d->builtin_func;
        d->builtin_func = lowprec_dispatch;
    }
}

#else
void sf_lowprec_init(void) {}
#endif
