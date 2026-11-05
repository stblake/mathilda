/* Mathilda -- shared numeric dispatch for two-argument special functions.
 * See sf_dispatch2.h. */
#include "sf_dispatch2.h"

#include "common.h"
#include "numeric.h"
#include "arithmetic.h"         /* is_complex, make_complex */

#ifdef USE_MPFR
#include <mpfr.h>
#endif

/* A real numeric leaf (machine/MPFR Real, Integer, Rational) as a double. */
static bool sf2_leaf_d(const Expr* e, double* out) {
#ifdef USE_MPFR
    if (e && e->type == EXPR_MPFR) { *out = mpfr_get_d(e->data.mpfr, MPFR_RNDN); return true; }
#endif
    return common_machine_real_value(e, out);
}

/* Demote an arb result (53-bit MPFR leaves, possibly Complex) to machine Reals.
 * Takes ownership of `e`; returns it unchanged if it is not a plain number. */
static Expr* sf2_demote_machine(Expr* e) {
    double d;
    if (sf2_leaf_d(e, &d)) { expr_free(e); return expr_new_real(d); }
    Expr *re, *im;
    double dr, di;
    if (is_complex(e, &re, &im) && sf2_leaf_d(re, &dr) && sf2_leaf_d(im, &di)) {
        expr_free(e);
        return make_complex(expr_new_real(dr), expr_new_real(di));
    }
    return e;
}

/* Minimum precision (bits) of e's inexact components; 0 when e is exact. */
static long sf2_inexact_bits(const Expr* e) {
#ifdef USE_MPFR
    return numeric_min_inexact_bits(e);
#else
    Expr *re, *im;
    if (e && e->type == EXPR_REAL) return 53;
    if (e && is_complex((Expr*)e, &re, &im)
            && (re->type == EXPR_REAL || im->type == EXPR_REAL)) return 53;
    return 0;
#endif
}

Expr* sf2_numeric(const Expr* nu, const Expr* z, int which,
                  sf2_machine_fn machine, sf2_arb_fn arb, bool* handled) {
    long bnu = sf2_inexact_bits(nu), bz = sf2_inexact_bits(z);
    *handled = false;
    if (bnu == 0 && bz == 0) return NULL;
    if (!expr_is_numeric_quantity(nu) || !expr_is_numeric_quantity(z)) return NULL;
    *handled = true;
    long outb = bnu;
    if (bz > 0 && (outb == 0 || bz < outb)) outb = bz;

    if (outb > 53) return arb(nu, z, which);

    double nuv, zv;
    if (common_machine_real_value(nu, &nuv) && common_machine_real_value(z, &zv)) {
        double v;
        if (machine(which, nuv, zv, &v)) return expr_new_real(v);
    }
    Expr* v = arb(nu, z, which);
    return v ? sf2_demote_machine(v) : NULL;
}

bool sf2_exact_integer(const Expr* e) {
    return e && (e->type == EXPR_INTEGER || e->type == EXPR_BIGINT);
}

bool sf2_exact_zero(const Expr* e) {
    return e && e->type == EXPR_INTEGER && e->data.integer == 0;
}
