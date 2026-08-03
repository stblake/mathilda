#ifndef MATHILDA_NDINTEGER_H
#define MATHILDA_NDINTEGER_H

/* ---------------------------------------------------------------------------
 * ndinteger.c — exact-integer NDArray kernels for the number-theoretic heads.
 *
 * ndkernels.c owns the libc-expressible elementary and special functions: every
 * kernel there computes in `double` and its exactness question is only ever
 * "does this real answer round-trip". These are the opposite kind. GCD, LCM,
 * EulerPhi, MoebiusMu, DivisorSigma, IntegerLength, PowerMod, Prime and
 * IntegerDigits are defined on Z, their answers are exact Integers, and a
 * `double` intermediate is not an approximation of the answer but a different
 * answer. So they are written in int64 end to end, through the ci_*_i64
 * overflow-checked helpers, and ABANDON the whole array the moment a result
 * would not fit — the List path then re-runs the call and GMP answers exactly.
 * Never wrap, never round: the same contract as the narrowing kernels in
 * ndkernels.c and Total's int64 accumulate.
 *
 * WHY THIS EXISTS. `Mod[Range[10^6] 7919, 1000]` is how an integer array is
 * made, and until it kept its buffer every consumer downstream paid. The eighth
 * round fixed Mod and Quotient and left the rest of the integer domain with no
 * kernel at all: `GCD[cv, 1234]` and its seven siblings answered a visible
 * NDArray with the unevaluated call and a packed List with one Expr per
 * element. The natural next batch, and the same reasoning — an integer pipeline
 * is only as packed as its least-aware link.
 * -------------------------------------------------------------------------- */

#include "expr.h"
#include <stdbool.h>
#include <stdint.h>

/* Register the unary/binary integer kernels (GCD, LCM, DivisorSigma, EulerPhi,
 * MoebiusMu, IntegerLength) on their symbols. Called from core_init. */
void ndinteger_init(void);

/* ---- head-level buffer paths ------------------------------------------- *
 * Three heads cannot be expressed as an element kernel and get their own entry
 * point, called from the head's builtin before it does anything else. Each
 * returns NULL for any shape it does not handle, leaving the builtin's ordinary
 * path untouched.
 *
 *   Prime      per-element Prime[n] is a Meissel-style prime count each time.
 *              Over an array the answer is a SIEVE: one pass to the largest
 *              requested index's upper bound, then a gather. That is a
 *              different algorithm, not a different loop, which is why it
 *              cannot be a kernel.
 *   PowerMod   ternary — PowerMod[b, e, m] — and NDUnary/NDBinaryKernel cannot
 *              express three arguments.
 *   IntegerDigits  the result is RAGGED (a list per element, of differing
 *              lengths), so no buffer holds it. The win is on the INPUT side:
 *              read the buffer directly instead of materialising 10^6 Integers
 *              only to read each one back out again.
 */
Expr* ndint_prime(Expr* res);          /* Prime[arr] */
Expr* ndint_powermod(Expr* res);       /* PowerMod[arr, e, m] */
Expr* ndint_integerdigits(Expr* res);  /* IntegerDigits[arr] */

/* ---- elementwise sign predicates ---------------------------------------- *
 * Positive / Negative / NonNegative / NonPositive. Listable, answering with a
 * List of True/False. Both sides now stay on the buffer: the INPUT is read
 * straight off it and the OUTPUT is a packed one-byte-per-element NDT_BOOL array
 * (performance.md §13 gap C.1, now closed), inheriting the input's presentation.
 * No per-element Expr and no evaluator round-trip on either side.
 *
 * `which` selects the comparison; see NDSignPred. Returns NULL unless the sole
 * argument is a rank-1 real or int64 NDArray (a complex or Indeterminate element
 * has no order, so those degrade to the List path element for element). */
typedef enum {
    NDSP_POSITIVE, NDSP_NEGATIVE, NDSP_NONNEGATIVE, NDSP_NONPOSITIVE
} NDSignPred;
Expr* ndint_sign_predicate(Expr* res, NDSignPred which);

#endif /* MATHILDA_NDINTEGER_H */
