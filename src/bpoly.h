/*
 * bpoly.h
 * -------
 * Bivariate polynomials over the integers Z, designed as the working
 * representation for the multivariate factoring pipeline (mvfactor.c).
 *
 * Storage: a BPoly represents
 *
 *      P(x, y) = sum_{i = 0}^{deg_x} cx[i] * x^i
 *
 * where each cx[i] is a ZUPoly in y.  The "main" variable x is dense
 * (cx is a contiguous array); the secondary variable y is dense within
 * each x-coefficient.  This layout is well suited to:
 *
 *   - Hensel lifting in y (we routinely truncate every cx[i] modulo
 *     y^k -- a coefficient-by-coefficient operation).
 *   - Polynomial division viewing P as Z[y][x] (each step reduces the
 *     leading x-coefficient using ZUPoly arithmetic).
 *   - Substitution y -> alpha (collapse each cx[i] to a constant in
 *     Z, yielding a univariate polynomial in x).
 *
 * Memory model:
 *   - Every cx slot up to `cap_x` holds an owning ZUPoly* (or NULL).
 *   - Slots beyond `deg_x` are NULL or zero-polynomials.
 *   - Construction routines transfer ownership of ZUPoly* arguments
 *     where indicated.
 *
 * Invariants:
 *   - deg_x >= -1 (the zero polynomial has deg_x = -1).
 *   - For deg_x >= 0, cx[deg_x] is non-zero.
 *   - cap_x > deg_x (so the leading slot exists).
 */

#ifndef BPOLY_H
#define BPOLY_H

#include <stdbool.h>
#include <stdint.h>
#include <stddef.h>
#include "zupoly.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    int       deg_x;   /* x-degree; -1 for zero */
    int       cap_x;   /* allocated capacity for cx[] */
    ZUPoly**  cx;      /* cx[i] is the y-polynomial coefficient of x^i.
                        * Always allocated for 0 <= i < cap_x.  May be
                        * NULL or a zero-poly for missing terms. */
} BPoly;

/* ---------------------------------------------------------------------- */
/*  Construction / destruction                                            */
/* ---------------------------------------------------------------------- */

/* Allocate a zero polynomial with capacity for at least `cap` x-terms. */
BPoly* bpoly_new(int cap);

/* The zero polynomial. */
BPoly* bpoly_zero(void);

/* Deep copy (each cx[i] is also copied). */
BPoly* bpoly_copy(const BPoly* p);

/* Free.  Safe with NULL. */
void bpoly_free(BPoly* p);

/* ---------------------------------------------------------------------- */
/*  Coefficient access                                                    */
/* ---------------------------------------------------------------------- */

/* Set cx[i] to `c`.  TRANSFERS OWNERSHIP of c -- caller must not free
 * it after.  Pass NULL or zero-poly to clear.  Updates deg_x. */
void bpoly_set_xcoef(BPoly* p, int i, ZUPoly* c);

/* Get cx[i] as a borrowed pointer (or NULL if out of range / missing). */
const ZUPoly* bpoly_get_xcoef(const BPoly* p, int i);

/* Re-establish deg_x invariant after direct mutation of cx[]. */
void bpoly_normalize(BPoly* p);

/* ---------------------------------------------------------------------- */
/*  Predicates and degree                                                 */
/* ---------------------------------------------------------------------- */

bool bpoly_is_zero(const BPoly* p);
bool bpoly_eq(const BPoly* a, const BPoly* b);

/* Degree of P as a polynomial in x: cx index of the leading term. */
int bpoly_deg_x(const BPoly* p);

/* Maximum y-degree across all cx[i] (or -1 for the zero polynomial). */
int bpoly_deg_y(const BPoly* p);

/* Leading coefficient in x: a borrowed pointer to cx[deg_x] (or NULL
 * for the zero polynomial). */
const ZUPoly* bpoly_lc_x(const BPoly* p);

/* ---------------------------------------------------------------------- */
/*  Arithmetic                                                            */
/* ---------------------------------------------------------------------- */

BPoly* bpoly_add(const BPoly* a, const BPoly* b);
BPoly* bpoly_sub(const BPoly* a, const BPoly* b);
BPoly* bpoly_mul(const BPoly* a, const BPoly* b);
BPoly* bpoly_neg(const BPoly* a);

/* Multiply, then truncate every y-coefficient modulo y^max_y_plus_one.
 * Useful in Hensel lifting where we work modulo successive powers of y. */
BPoly* bpoly_mul_truncate_y(const BPoly* a, const BPoly* b,
                            int max_y_plus_one);

/* Multiply each cx[i] (y-poly) by a scalar y-poly s -- equivalent to
 * a * (s viewed as a bivariate constant in x). */
BPoly* bpoly_mul_zupoly(const BPoly* a, const ZUPoly* s);

/* Exact bivariate division viewing as Z[y][x].  Returns NULL if not
 * exact.  The leading-coefficient ZUPoly division must succeed at
 * every step. */
BPoly* bpoly_divexact(const BPoly* a, const BPoly* b);

/* ---------------------------------------------------------------------- */
/*  Truncation, evaluation, shift                                         */
/* ---------------------------------------------------------------------- */

/* Truncate every y-coefficient of P modulo y^k.  Returns a fresh BPoly. */
BPoly* bpoly_truncate_y(const BPoly* p, int k);

/* Substitute y -> alpha.  Returns a univariate ZUPoly in x. */
ZUPoly* bpoly_eval_y_si(const BPoly* p, int64_t alpha);

/* Substitute y -> y + alpha (a translation).  Returns a fresh BPoly. */
BPoly* bpoly_shift_y_si(const BPoly* p, int64_t alpha);

/* ---------------------------------------------------------------------- */
/*  Conversion to / from picocas Expr                                     */
/* ---------------------------------------------------------------------- */

struct Expr;

/* Try to interpret `e` as a polynomial in (x_var, y_var) with integer
 * coefficients in both.  Returns NULL if any coefficient is non-integer
 * or any variable other than x_var/y_var appears with non-zero
 * exponent.  Caller owns the returned BPoly. */
BPoly* expr_to_bpoly(const struct Expr* e,
                     const struct Expr* x_var,
                     const struct Expr* y_var);

/* Build a fresh Expr representing P(x_var, y_var).  Caller owns. */
struct Expr* bpoly_to_expr(const BPoly* p,
                           const struct Expr* x_var,
                           const struct Expr* y_var);

/* ---------------------------------------------------------------------- */
/*  Debug                                                                 */
/* ---------------------------------------------------------------------- */

/* Print a human-readable form to stderr. */
void bpoly_print(const BPoly* p, const char* x_name, const char* y_name);

#ifdef __cplusplus
}
#endif

#endif /* BPOLY_H */
