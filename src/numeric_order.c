/* numeric_order.c -- certified numeric ordering of expressions; see
 * numeric_order.h for the order and its callers. */

#include "numeric_order.h"
#include "comparisons.h"       /* comparisons_numeric_order: the certified comparator */
#include "numeric.h"           /* numericalize, numeric_machine_spec, numeric_machine_value */
#include "list/list_common.h"  /* is_real_numeric, is_infinity, is_minus_infinity */
#include "arithmetic.h"        /* is_complex */
#include "expr.h"              /* expr_compare */
#include "sym_names.h"
#include <float.h>
#include <math.h>
#ifdef USE_MPFR
#include <mpfr.h>
#endif

/* +1 / -1 for the real infinities (Infinity, -Infinity, DirectedInfinity[+-1]),
 * else 0. */
static int numord_inf_sign(Expr* e) {
    if (is_infinity(e)) return 1;
    if (is_minus_infinity(e)) return -1;
    if (e->type == EXPR_FUNCTION && e->data.function.head->type == EXPR_SYMBOL &&
        e->data.function.head->data.symbol.name == SYM_DirectedInfinity &&
        e->data.function.arg_count == 1) {
        Expr* d = e->data.function.args[0];
        if (d->type == EXPR_INTEGER) return (d->data.integer > 0) - (d->data.integer < 0);
    }
    return 0;
}

/* The double of a real number atom (Integer, BigInt, Rational, Real, MPFR, or
 * a Complex with an exact-zero imaginary part). False for anything else. */
static bool numord_atom_value(Expr* e, double* v) {
    Expr *re, *im;
    if (is_complex(e, &re, &im)) e = re;
    if (numeric_machine_value(e, v)) return true;
#ifdef USE_MPFR
    if (e->type == EXPR_MPFR) { *v = mpfr_get_d(e->data.mpfr, MPFR_RNDN); return true; }
#endif
    return false;
}

/* True when no Plus occurs in e: a machine reading of such an expression keeps
 * its relative accuracy (no cancellation), so a clear gap between two such
 * readings decides their order. A sum's reading carries no relative accuracy
 * (PolyLog[2, 1/3] - 2006534954/5479143815 is +6.9*10^-21 and reads
 * -5.6*10^-17). Bounded: a deeper tree is conservatively "may cancel". */
static bool numord_plus_free(const Expr* e, int depth) {
    if (e->type != EXPR_FUNCTION) return true;
    if (depth > 32) return false;
    const Expr* h = e->data.function.head;
    if (h->type == EXPR_SYMBOL && h->data.symbol.name == SYM_Plus) return false;
    if (!numord_plus_free(h, depth + 1)) return false;
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (!numord_plus_free(e->data.function.args[i], depth + 1)) return false;
    return true;
}

void numord_key(Expr* e, NumOrdKey* k) {
    k->v = 0.0; k->cls = NUMORD_OTHER; k->exact = 0; k->sep = 0;
    /* O(1) leaves: the packed / numeric data path. */
    if (e->type == EXPR_INTEGER) {
        int64_t x = e->data.integer;
        k->v = (double)x;
        k->cls = NUMORD_REAL;
        k->exact = (x <= INT64_C(9007199254740992) && x >= -INT64_C(9007199254740992));
        k->sep = 1;
        return;
    }
    if (e->type == EXPR_REAL) {
        if (isnan(e->data.real)) return;
        k->v = e->data.real; k->cls = NUMORD_REAL; k->exact = 1; k->sep = 1;
        return;
    }
    int s = numord_inf_sign(e);
    if (s) { k->v = s > 0 ? HUGE_VAL : -HUGE_VAL; k->cls = NUMORD_REAL; k->exact = 1; k->sep = 1; return; }
    if (is_real_numeric(e)) {                    /* BigInt, Rational, MPFR, real Complex */
        double v;
        if (numord_atom_value(e, &v) && !isnan(v)) {
            k->v = v; k->cls = NUMORD_REAL; k->sep = 1;
        }
        return;
    }
    if (e->type != EXPR_FUNCTION && e->type != EXPR_SYMBOL) return;   /* a string, ... */
    /* A symbolic real (Pi, E, Sqrt[2], -Pi, Log[3], E + 1): numericalize once;
     * a free symbol stays symbolic and ranks after the reals. */
    Expr* nv = numericalize(e, numeric_machine_spec());
    if (!nv) return;
    double v;
    if (is_real_numeric(nv) && numord_atom_value(nv, &v) && !isnan(v)) {
        k->v = v; k->cls = NUMORD_REAL; k->sep = numord_plus_free(e, 0);
    }
    expr_free(nv);
}

int numord_cmp(Expr* a, const NumOrdKey* ka, Expr* b, const NumOrdKey* kb) {
    if (ka->cls != kb->cls) return ka->cls < kb->cls ? -1 : 1;
    if (ka->cls == NUMORD_OTHER) return expr_compare(a, b);
    double va = ka->v, vb = kb->v;
    if (ka->exact && kb->exact) return (va > vb) - (va < vb);
    /* A clear gap between two accurate readings in the normal range decides. */
    if (ka->sep && kb->sep && isfinite(va) && isfinite(vb) &&
        (va == 0.0 || fabs(va) >= DBL_MIN) && (vb == 0.0 || fabs(vb) >= DBL_MIN)) {
        double mx = fabs(va) > fabs(vb) ? fabs(va) : fabs(vb);
        if (fabs(va - vb) > mx * 0x1p-40) return (va > vb) - (va < vb);
    }
    /* Close, cancelling, or out of the double range: the certified comparator
     * (exact GMP for rationals and big integers, MPFR at increasing precision
     * for exact quantities). Undecided means no order to give: a tie. */
    bool decided = false;
    int c = comparisons_numeric_order(a, b, &decided);
    return decided ? c : 0;
}

int numord_compare(Expr* a, Expr* b) {
    NumOrdKey ka, kb;
    numord_key(a, &ka);
    numord_key(b, &kb);
    return numord_cmp(a, &ka, b, &kb);
}
