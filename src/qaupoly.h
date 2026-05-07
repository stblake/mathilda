/* qaupoly.h — Q(α)[x] univariate polynomials (Phase G2, Trager).
 *
 * QAUPoly represents a polynomial in `x` whose coefficients are
 * elements of an algebraic extension Q(α) = Q[α]/(P_α(α)).  This is
 * the substrate Phase G3 (norm via resultant) and G4 (Trager's
 * sqfr_norm + alg_factor) build on.
 *
 * Storage parallels ZUPoly: a dense array `c[]` of `QANum*`, with
 * `c[i]` the coefficient of x^i.  All slots in the underlying capacity
 * are valid (initialised to qa_zero(ext)); slots above `deg` are
 * available for in-place updates.  `deg = -1` denotes the zero
 * polynomial.
 *
 * Q(α) is a field whenever P_α is irreducible, so QAUPoly admits
 * standard Euclidean division and gcd.  The implementation does not
 * verify P_α irreducibility — the caller is responsible for that
 * contract (and Phase G5's picocas-level API only builds extensions
 * where it holds, e.g. Sqrt[c] / c^(1/n)). */

#ifndef PICOCAS_QAUPOLY_H
#define PICOCAS_QAUPOLY_H

#include <stddef.h>
#include <stdbool.h>
#include "qa.h"

typedef struct QAUPoly {
    const QAExt* ext;   /* borrowed: extension definition */
    int     deg;        /* -1 for zero, otherwise degree of x */
    int     cap;        /* allocated slots in c[]: each >= 1 */
    QANum** c;          /* c[i] is owned; non-NULL after init */
} QAUPoly;

/* ============================ Construction ============================ */

/* Allocate a polynomial with at least `cap` coefficient slots
 * (capacity grows automatically on setcoef when needed).  All slots
 * start at qa_zero(ext); deg = -1 (zero polynomial). */
QAUPoly* qaupoly_new(const QAExt* ext, int cap);

/* Convenience: the zero polynomial. */
QAUPoly* qaupoly_zero(const QAExt* ext);

/* Constant polynomial `v`.  Takes a borrowed QANum (deep-copies it). */
QAUPoly* qaupoly_from_qa(const QANum* v);

/* Constant polynomial in Q (no α dependency). */
QAUPoly* qaupoly_from_si(const QAExt* ext, long num, long den);

/* p(x) = x   (the indeterminate). */
QAUPoly* qaupoly_x(const QAExt* ext);

/* Deep copy. */
QAUPoly* qaupoly_copy(const QAUPoly* p);

void qaupoly_free(QAUPoly* p);

/* ====================== Coefficient access ======================== */

/* Replace c[i] with a copy of `v`.  Grows capacity if i >= cap; updates
 * deg if necessary (does not normalize trailing zeros — call
 * qaupoly_normalize after building). */
void qaupoly_setcoef(QAUPoly* p, int i, const QANum* v);

/* Returns a borrowed pointer to c[i].  If i > deg or i out of range,
 * returns the zero element (also borrowed; valid until the polynomial
 * is freed or its capacity grown). */
const QANum* qaupoly_getcoef(const QAUPoly* p, int i);

/* Recompute deg by scanning down from cap-1 for the highest non-zero
 * coefficient. */
void qaupoly_normalize(QAUPoly* p);

/* ============================== Predicates ============================== */

bool qaupoly_is_zero(const QAUPoly* p);
bool qaupoly_eq(const QAUPoly* a, const QAUPoly* b);

/* ============================== Arithmetic ============================== */

QAUPoly* qaupoly_add(const QAUPoly* a, const QAUPoly* b);
QAUPoly* qaupoly_sub(const QAUPoly* a, const QAUPoly* b);
QAUPoly* qaupoly_neg(const QAUPoly* a);
QAUPoly* qaupoly_mul(const QAUPoly* a, const QAUPoly* b);

QAUPoly* qaupoly_scale_qa(const QAUPoly* a, const QANum* k);
QAUPoly* qaupoly_scale_si(const QAUPoly* a, long num, long den);

/* Long division over Q(α)[x].  On exit `*q_out` is the quotient and
 * `*r_out` is the remainder; both are freshly allocated and owned by
 * the caller.  Returns false (with both outputs untouched) if `b` is
 * zero or its leading coefficient is non-invertible in Q(α). */
bool qaupoly_divrem(const QAUPoly* a, const QAUPoly* b,
                    QAUPoly** q_out, QAUPoly** r_out);

/* Euclidean gcd over Q(α)[x].  Result is monic (leading coefficient
 * normalised to 1 in Q(α)).  Returns NULL if the gcd is zero, which
 * only happens when both inputs are zero. */
QAUPoly* qaupoly_gcd(const QAUPoly* a, const QAUPoly* b);

/* Make `p` monic by dividing all coefficients by the leading
 * coefficient.  Returns NULL if `p` is zero. */
QAUPoly* qaupoly_make_monic(const QAUPoly* p);

/* ============================ Substitution ============================ */

/* Evaluate p at x = c (an element of Q(α)).  Horner. */
QANum* qaupoly_eval(const QAUPoly* p, const QANum* c);

/* Compute p(x + c) for c in Q(α).  Used by Trager's sqfr_norm to
 * apply the linear shift x → x − sα. */
QAUPoly* qaupoly_shift(const QAUPoly* p, const QANum* c);

#endif
