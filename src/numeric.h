/* Mathilda — numeric evaluation (N[], SetPrecision, …)
 *
 * This module implements numeric approximation of symbolic expressions.
 * The public entry point is `N[expr]` and `N[expr, prec]`. Phase 1 supports
 * machine-precision doubles; Phase 2 adds MPFR arbitrary precision behind
 * the `USE_MPFR` compile flag.
 *
 * Design overview
 * ---------------
 * `N` does NOT reimplement numeric math. Its job is to replace *leaves*
 * (integers, rationals, named constants like Pi/E/EulerGamma) with their
 * numeric equivalents and then hand the rebuilt expression back to the
 * evaluator. The existing builtins — Plus/Times/Power/Sin/Cos/Exp/Log —
 * already take the numeric fast path whenever they see a numeric argument,
 * so re-evaluation naturally produces a numeric result.
 *
 * This "descend + re-evaluate" strategy keeps `numeric.c` small and places
 * the numeric-math responsibility in the modules that own each function.
 *
 * Extensibility
 * -------------
 * • New constants (e.g. `Glaisher`): append to `kConstants` in numeric.c.
 * • New backends (e.g. GSL for special functions, libquadmath for
 *   __float128): extend `NumericMode` below and thread through
 *   `numericalize`. Each per-function numeric branch (in trig.c, logexp.c,
 *   etc.) is the natural place to hook in backend-specific calls — no
 *   central dispatcher is needed.
 * • New special functions (e.g. `BesselJ`, `Gamma`, `Zeta`): add their
 *   numeric branch to the function's own builtin (as Sin already does for
 *   `csin`). `N` will drive it automatically because the evaluator reaches
 *   that builtin after we rebuild with numericalized arguments.
 */
#ifndef NUMERIC_H
#define NUMERIC_H

#include "expr.h"
#include <float.h>
#include <stdbool.h>

/* Evaluation backends. MACHINE uses IEEE 754 doubles; MPFR uses
 * arbitrary-precision reals (Phase 2). Additional backends (e.g. GSL,
 * __float128) would be appended here. */
typedef enum {
    NUMERIC_MODE_MACHINE = 0
#ifdef USE_MPFR
    , NUMERIC_MODE_MPFR
#endif
} NumericMode;

/* A precision request. For MACHINE, `bits` is ignored. For MPFR, `bits` is
 * the target precision in binary bits (MPFR's native unit). Converters
 * below bridge the Mathematica-style decimal-digit specification. */
typedef struct {
    NumericMode mode;
    long bits;
    /* Bare-N[expr] intent: when true (and mode is MACHINE), already-inexact
     * arbitrary-precision leaves keep their existing precision instead of
     * being down-converted to a machine double. Mathematica's N[expr]
     * numericalizes only *exact* quantities and leaves approximate numbers
     * alone, so N[N[Pi, 100]] must stay 100 digits. Consulted solely in the
     * MACHINE/EXPR_MPFR branch of numericalize(); every MPFR-mode spec leaves
     * it untouched (harmless) and every MACHINE-mode spec is built through
     * numeric_machine_spec(), which clears it. */
    bool preserve_inexact;
} NumericSpec;

/* Decimal-digit equivalent of one full machine double's mantissa,
 * derived from the platform's DBL_MANT_DIG so it tracks the local float
 * representation (53 for IEEE 754, giving the familiar ~15.955). This is
 * the C-side value backing the $MachinePrecision system variable. */
#define NUMERIC_LOG10_2 0.30102999566398119521
#define NUMERIC_MACHINE_PRECISION_DIGITS ((double)DBL_MANT_DIG * NUMERIC_LOG10_2)

/* Digit / bit conversions for MPFR precision. Both round up so the target
 * carries at least the requested amount of information. */
long numeric_digits_to_bits(double digits);
double numeric_bits_to_digits(long bits);

/* Make a default-machine spec. */
static inline NumericSpec numeric_machine_spec(void) {
    NumericSpec s;
    s.mode = NUMERIC_MODE_MACHINE;
    s.bits = 0;
    s.preserve_inexact = false;
    return s;
}

/* Numericalize an expression under `spec`, returning a *newly allocated*
 * Expr tree. The input `e` is read-only and not consumed by this function.
 *
 * Semantics (matching Mathematica):
 *   - Exact numeric leaves (Integer, BigInt, Rational) → numeric value.
 *   - Inexact numeric leaves (Real; MPFR in Phase 2) → re-rounded to spec.
 *   - Known constants (Pi, E, EulerGamma, Catalan, GoldenRatio, Degree)
 *     → their numeric value at the requested precision.
 *   - Unknown symbols → returned unchanged (Mathematica's N[x] = x).
 *   - Hold / HoldForm / Unevaluated → returned unchanged.
 *   - Complex[re, im] → recursively numericalize components.
 *   - Rational[n, d] → direct (double)n/d (or MPFR equivalent).
 *   - Any other function f[a, b, …] → f[N[a], N[b], …] then re-evaluate.
 *
 * Ownership: returned Expr* is owned by the caller. */
Expr* numericalize(const Expr* e, NumericSpec spec);

/* Look up a named constant (Pi, E, …). Returns true and fills *out if the
 * symbol is a recognized numeric constant; false otherwise. Exposed so the
 * Phase-2 MPFR path can implement its own constant fill while sharing the
 * registry. */
bool numeric_constant_machine_value(const char* name, double* out);

#ifdef USE_MPFR
/* The target precision (in bits) for combining two numeric expressions.
 * If either operand is an MPFR value, we use the max of their precisions;
 * otherwise the caller's `default_bits` is used. Exact inputs (Integer,
 * Rational) contribute no precision and accept the ambient precision. */
long numeric_combined_bits(const Expr* a, const Expr* b, long default_bits);

/* Governing precision (bits) for a numeric result built from `e`, following
 * Mathematica's precision contagion: the MINIMUM precision among the inexact
 * (Real/MPFR) leaves of `e`, descending into Complex[...]. A machine `Real`
 * counts as 53 bits, an `MPFR` as its own bit precision, and exact atoms
 * (Integer/Rational/BigInt) impose no constraint. Returns 0 if `e` has no
 * inexact leaf.
 *
 * Unlike `numeric_combined_bits` (which takes the MAX and treats a machine
 * Real as no constraint), this is the rule a multi-argument special function
 * needs: e.g. PolyLog[N[1/2], 1.3429`50] must yield a machine-precision
 * result because one argument is machine precision. Combine several arguments
 * by taking the min of their non-zero results, then floor at 53. */
long numeric_min_inexact_bits(const Expr* e);

/* Binary MPFR-level numeric operators. Each returns a freshly allocated
 * EXPR_MPFR (real-valued) on success, or NULL if either input is not
 * representable as a real MPFR value (e.g. it's a Complex[...] or symbolic).
 * Caller must free the return with expr_free.
 *
 * The output precision is `numeric_combined_bits(a, b, default_bits)`.
 * These helpers are used by Plus, Times, Power, and the transcendental
 * builtins to take the arbitrary-precision numeric path when any operand
 * carries MPFR precision. */
Expr* numeric_mpfr_add (const Expr* a, const Expr* b, long default_bits);
Expr* numeric_mpfr_sub (const Expr* a, const Expr* b, long default_bits);
Expr* numeric_mpfr_mul (const Expr* a, const Expr* b, long default_bits);
Expr* numeric_mpfr_div (const Expr* a, const Expr* b, long default_bits);
Expr* numeric_mpfr_pow (const Expr* a, const Expr* b, long default_bits);

/* Complex power at MPFR precision via polar form, used by `Power` when
 * `numeric_mpfr_pow` can't apply (negative real base with fractional
 * exponent, complex base, or complex exponent). Returns Complex[MPFR, MPFR],
 * a pure MPFR when the imaginary part rounds to zero, or NULL if either
 * operand isn't MPFR-promotable or `base` is zero. */
Expr* numeric_mpfr_complex_pow(const Expr* base, const Expr* exp,
                               long default_bits);

/* Apply an MPFR unary real function to `e`. Returns a fresh EXPR_MPFR,
 * or NULL if `e` isn't a real-valued MPFR-promotable number (e.g. it's
 * complex or symbolic). Used by Sin/Cos/Tan/Sinh/Cosh/Log/Exp/…
 *
 * Output precision is `max(mpfr_prec in e, default_bits)`. */
typedef int (*MpfrUnaryOp)(mpfr_t, const mpfr_t, mpfr_rnd_t);
Expr* numeric_mpfr_apply_unary(const Expr* e, long default_bits, MpfrUnaryOp op);

/* Convenience: true iff `e` should trigger the MPFR numeric path (carries
 * any MPFR subvalue). */
bool numeric_expr_is_mpfr(const Expr* e);

/* True iff any of a or b is an EXPR_MPFR — trigger for the MPFR code
 * paths in plus/times/power. Complex-of-MPFR is also detected so that
 * MPFR propagation happens through the Complex wrapper. */
bool numeric_any_mpfr(const Expr* a, const Expr* b);

/* Fill the recognized constant `name` into `out` at `bits` precision.
 * Returns true if the name is a known numeric constant, false otherwise.
 * `out` must be a valid (initialized) mpfr_t on entry; its precision is
 * reset to `bits` before filling. */
bool numeric_constant_mpfr_value(const char* name, mpfr_t out, long bits);

/* Extract an MPFR approximation of `e` into `(re, im)`.
 *
 * The caller must have already `mpfr_init2`'d both `re` and `im` to the
 * desired output precision. `re` and `im` are written with the real and
 * imaginary components of the value. Returns true if extraction was
 * possible (i.e. `e` is some kind of numeric or Complex-of-numerics), and
 * sets *is_inexact according to whether any input was inexact (EXPR_REAL
 * or EXPR_MPFR — "inexact" in the Mathematica sense).
 *
 * Exact inputs (Integer, BigInt, Rational) are converted at the target
 * precision. Complex[re, im] recurses on components.
 *
 * When `e` is unrecognized (symbolic, or a non-numeric function), returns
 * false without touching `re` or `im`.
 */
bool get_approx_mpfr(const Expr* e,
                     mpfr_t re, mpfr_t im,
                     bool* is_inexact);
#endif

/* Mathematica inexact-contagion pre-pass for Plus / Times / other numeric
 * heads. Scans `args[0..n]` for any inexact numeric leaf (Real, MPFR, or
 * Complex[...] containing one). If found, returns a freshly allocated
 * array of `n` Expr* where each entry is numericalize(args[i]) at the
 * maximum MPFR precision observed among the inexact args (or machine
 * precision if only EXPR_REAL inexactness is present). This is how
 * `1. Pi` collapses to `3.14159` — the Pi is numericalized in-line.
 *
 * Returns NULL if no arg is inexact; the caller should then proceed with
 * the original arg list. The returned array and all its elements are
 * caller-owned (free each with expr_free and the array with free()). */
Expr** numeric_contagion_args(Expr* const* args, size_t n);

/* Builtin entry points */
Expr* builtin_n(Expr* res);

/* MachineNumberQ[expr] -> True/False. True iff expr is a finite machine
 * (IEEE double) real, or a Complex of two finite machine reals. Always
 * returns a True/False symbol when arity is 1; returns NULL only on a
 * structurally invalid call (arity != 1). */
Expr* builtin_machinenumberq(Expr* res);

/* Registers N (and, in Phase 2, Precision/Accuracy/SetPrecision/SetAccuracy
 * and MachinePrecision) with the symbol table. */
void numeric_init(void);

#endif /* NUMERIC_H */
