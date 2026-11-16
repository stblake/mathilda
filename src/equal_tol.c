/*
 * equal_tol.c -- Internal`$EqualTolerance and Internal`$SameQTolerance.
 * See equal_tol.h.
 */
#include "equal_tol.h"
#include "expr.h"
#include "symtab.h"
#include "sym_names.h"
#include "attr.h"
#include "arithmetic.h"
#include <float.h>
#include <stdint.h>

#ifdef USE_MPFR
#include <mpfr.h>
#endif

/* log2 10 and its reciprocal, for digits <-> bits. */
#define EQT_LOG2_10   3.3219280948873623
#define EQT_LOG10_2   0.30102999566398120

#define EQT_DEFAULT_EQUAL_BITS 7.0
#define EQT_DEFAULT_SAMEQ_BITS 1.0

/* The two symbols' definitions, captured once at init: a SymbolDef lives as
 * long as the symbol table, and its own_values are read live, so Set, Unset
 * and Block on the variable take effect at the next comparison. */
static SymbolDef* g_eqt_equal_def = NULL;
static SymbolDef* g_eqt_sameq_def = NULL;

/* The numeric value of a symbol's OwnValue, in decimal digits; false when it
 * has none or it is not a real number. */
static bool eqt_read_digits(const SymbolDef* d, double* out) {
    if (!d || !d->own_values || !d->own_values->replacement) return false;
    const Expr* v = d->own_values->replacement;
    switch (v->type) {
        case EXPR_INTEGER: *out = (double)v->data.integer; return true;
        case EXPR_REAL:    *out = v->data.real; return isfinite(*out);
#ifdef USE_MPFR
        case EXPR_MPFR:    *out = mpfr_get_d(v->data.mpfr, MPFR_RNDN); return isfinite(*out);
#endif
        default: break;
    }
    int64_t n, q;
    if (is_rational((Expr*)v, &n, &q) && q != 0) { *out = (double)n / (double)q; return true; }
    return false;
}

static double eqt_bits(const SymbolDef* d, double dflt) {
    double digits;
    if (!eqt_read_digits(d, &digits)) return dflt;
    double bits = digits * EQT_LOG2_10;
    if (!(bits > 0.0)) return 0.0;
    if (bits > (double)DBL_MANT_DIG) return (double)DBL_MANT_DIG;
    return bits;
}

double mth_equal_tolerance_bits(void) {
    return eqt_bits(g_eqt_equal_def, EQT_DEFAULT_EQUAL_BITS);
}

double mth_sameq_tolerance_bits(void) {
    return eqt_bits(g_eqt_sameq_def, EQT_DEFAULT_SAMEQ_BITS);
}

/* The default 7 bits returns the exact 2^-46 the comparison heads always used
 * (the digit value 2.10721 is 7 log10 2 rounded, and its bits are 7 only to
 * within a few ulps). */
static double eqt_rel(double bits, double dflt_bits) {
    if (fabs(bits - dflt_bits) < 1e-9) return ldexp(1.0, (int)dflt_bits - DBL_MANT_DIG);
    return exp2(bits - (double)DBL_MANT_DIG);
}

double mth_machine_equal_rel(void) {
    return eqt_rel(mth_equal_tolerance_bits(), EQT_DEFAULT_EQUAL_BITS);
}

double mth_machine_sameq_rel(void) {
    return eqt_rel(mth_sameq_tolerance_bits(), EQT_DEFAULT_SAMEQ_BITS);
}

static void eqt_register(const char* name, double digits, const char* doc) {
    Expr* sym = expr_new_symbol(name);
    Expr* val = expr_new_real(digits);
    symtab_add_own_value(name, sym, val);
    expr_free(sym);
    expr_free(val);
    symtab_set_docstring(name, doc);
}

void equal_tol_init(void) {
    eqt_register(SYM_InternalEqualTolerance, EQT_DEFAULT_EQUAL_BITS * EQT_LOG10_2,
        "Internal`$EqualTolerance is the number of trailing decimal digits Equal (and "
        "Unequal, Less, ... and compiled comparisons) ignores when comparing approximate "
        "numbers of at least machine precision. Default Log10[2^7] (7 bits).");
    eqt_register(SYM_InternalSameQTolerance, EQT_DEFAULT_SAMEQ_BITS * EQT_LOG10_2,
        "Internal`$SameQTolerance is the number of trailing decimal digits SameQ and "
        "UnsameQ ignore when comparing approximate numbers. Default Log10[2] (1 bit).");
    g_eqt_equal_def = symtab_get_def(SYM_InternalEqualTolerance);
    g_eqt_sameq_def = symtab_get_def(SYM_InternalSameQTolerance);
}
