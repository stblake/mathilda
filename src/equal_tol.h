/*
 * equal_tol.h -- the tolerance of Equal and SameQ on approximate numbers.
 *
 * Two approximate numbers of at least machine precision are Equal when they
 * differ in at most their last k binary digits, k = Internal`$EqualTolerance
 * * log2 10 (default 2.10721 decimal digits = 7 bits); SameQ ignores
 * Internal`$SameQTolerance (default 0.30103 = 1 bit) the same way.  Both
 * variables are ordinary user-settable symbols, so Block[{Internal`$EqualTolerance
 * = 0}, ...] makes Equal exact.
 *
 * This header is the ONE place the machine-real rule is written down: the
 * Equal/Less family (comparisons.c), the compiled VM (compile_vm.c), the
 * numloop predicate fast path (numloop.c) and the sparse Equal kernel
 * (sparse_reduce.c) all call mth_machine_tolerant_cmp, so a packed, sparse,
 * compiled and scalar comparison of the same two doubles cannot disagree.
 *
 * Reading the variables costs two pointer loads and a type test (no cache, so
 * no shared mutable state under threads); hot loops read the relative
 * tolerance once per call and pass it down.
 */
#ifndef MATHILDA_EQUAL_TOL_H
#define MATHILDA_EQUAL_TOL_H

#include <float.h>
#include <math.h>
#include <stdbool.h>

/* Bits ignored by Equal / SameQ: the current variable, read as a real number
 * of decimal digits and converted to bits, clamped to [0, 53]. A value that is
 * not a real number falls back to the default (7 bits / 1 bit). */
double mth_equal_tolerance_bits(void);
double mth_sameq_tolerance_bits(void);

/* The relative tolerance for two machine reals: 2^(k - 53). */
double mth_machine_equal_rel(void);
double mth_machine_sameq_rel(void);

/* True iff two machine reals are equal within relative tolerance `rel`:
 * identical (two infinities of one sign included), or both finite and apart by
 * at most rel * max(|a|, |b|).  False whenever either is NaN. */
static inline bool mth_machine_tol_equal(double a, double b, double rel) {
    if (a == b) return true;
    double ma = fabs(a), mb = fabs(b);
    double mx = ma > mb ? ma : mb;
    if (!(mx <= DBL_MAX)) return false;          /* an infinity or a NaN */
    return fabs(a - b) <= mx * rel;
}

/* Order of two machine reals under relative tolerance `rel`: 0 when equal
 * within it, else the sign of a - b.  A NaN operand is unordered: callers that
 * must treat it specially (IEEE) test for it first; here it yields +1. */
static inline int mth_machine_tolerant_cmp(double a, double b, double rel) {
    if (mth_machine_tol_equal(a, b, rel)) return 0;
    return a < b ? -1 : 1;
}

/* The six comparison heads on two machine reals, IEEE-faithful for NaN (every
 * ordering and Equal False, Unequal True) -- the form the compiled VM uses. */
static inline bool mth_tol_eq(double a, double b, double r) { return mth_machine_tol_equal(a, b, r); }
static inline bool mth_tol_ne(double a, double b, double r) { return !mth_machine_tol_equal(a, b, r); }
static inline bool mth_tol_lt(double a, double b, double r) { return a < b && !mth_machine_tol_equal(a, b, r); }
static inline bool mth_tol_gt(double a, double b, double r) { return a > b && !mth_machine_tol_equal(a, b, r); }
static inline bool mth_tol_le(double a, double b, double r) { return a <= b || mth_machine_tol_equal(a, b, r); }
static inline bool mth_tol_ge(double a, double b, double r) { return a >= b || mth_machine_tol_equal(a, b, r); }

/* Register the two variables with their default values (core init). */
void equal_tol_init(void);

#endif
