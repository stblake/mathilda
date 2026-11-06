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
    /* Like NUMERIC_MODE_MPFR, but for the two-argument N[expr, p]: it never
     * *increases* the precision of an already-inexact leaf. A machine Real
     * stays a machine Real (it cannot gain digits it doesn't have); an MPFR
     * value is capped at min(existing, requested) bits. Exact leaves
     * (Integer/Rational/constants) are still produced at the full requested
     * precision. SetPrecision/SetAccuracy and every internal
     * working-precision numericalization keep plain NUMERIC_MODE_MPFR (which
     * pads up), so this capping variant is produced solely by the N builtin.
     * See numeric_spec_is_mpfr() / numeric_spec_caps_inexact() below. */
    , NUMERIC_MODE_MPFR_CAP
#endif
} NumericMode;

/* A precision request. For MACHINE, `bits` is ignored. For MPFR, `bits` is
 * the target precision in binary bits (MPFR's native unit). Converters
 * below bridge the Mathematica-style decimal-digit specification. */
typedef struct {
    NumericMode mode;
    long bits;
    /* Two fields only, and keep it that way: 30+ call sites declare a
     * NumericSpec uninitialised and set .mode / .bits, so any further field
     * is stack garbage there (a `contagion` flag here was garbage-true at
     * -O2 and stopped N threading into List: Fit[..., WorkingPrecision -> 24]
     * returned the exact fit). Per-call state lives in numeric.c instead. */
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

/* The IEEE-754 double nearest the exact rational q, ties to even, with
 * gradual underflow: the same double a correctly rounded decimal literal of
 * q would give.  Exactly one rounding, done in GMP integers.  GMP's
 * mpq_get_d truncates; rounding to a 53-bit MPFR and then to a double rounds
 * TWICE, and below DBL_MIN, where a double holds fewer than 53 bits, the
 * second rounding could land on the wrong neighbour (N[r] for r just under
 * 3/2 2^-1074 gave 2^-1073).  Magnitudes past DBL_MAX give +-HUGE_VAL;
 * magnitudes below half the smallest subnormal give a zero of q's sign.
 * Declared with an opaque pointer so this header need not include gmp.h:
 * pass an mpq_t (it decays to mpq_ptr). */
double numeric_mpq_to_double(const void* q);

/* Make a default-machine spec. */
static inline NumericSpec numeric_machine_spec(void) {
    NumericSpec s;
    s.mode = NUMERIC_MODE_MACHINE;
    s.bits = 0;
    return s;
}

/* True iff `spec` requests an arbitrary-precision (MPFR) result, whether or
 * not it caps inexact leaves. Use in place of a bare
 * `spec.mode == NUMERIC_MODE_MPFR` test so the capping N variant
 * (NUMERIC_MODE_MPFR_CAP) still takes the MPFR production paths for exact
 * leaves. `mode` is always initialized, so this is UB-free at every call
 * site (unlike reading a bool field a partial initializer may leave unset). */
static inline bool numeric_spec_is_mpfr(NumericSpec spec) {
#ifdef USE_MPFR
    return spec.mode == NUMERIC_MODE_MPFR || spec.mode == NUMERIC_MODE_MPFR_CAP;
#else
    (void)spec;
    return false;
#endif
}

/* True iff `spec` must not increase the precision of an already-inexact leaf
 * (the two-argument N[expr, p] semantics). */
static inline bool numeric_spec_caps_inexact(NumericSpec spec) {
#ifdef USE_MPFR
    return spec.mode == NUMERIC_MODE_MPFR_CAP;
#else
    (void)spec;
    return false;
#endif
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
/* Machine precision without the working-precision plan (for screening; see numeric.c). */
Expr* numericalize_machine_unplanned(const Expr* e);
/* Re-evaluations N[expr, p]'s adaptive precision has run since startup (the
 * first evaluation of each N is not counted). Monotone; for tests that pin
 * "an expression that cancels nothing is evaluated once". */
long numeric_adaptive_reevaluations(void);

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
/* Record that an MPFR sum computed outside numeric_mpfr_add (a fused dot
 * product, say) lost `lost` leading bits to cancellation and kept `left`
 * (its result's precision less `lost`; 0 when it cancelled to an exact zero),
 * so N[expr, p]'s adaptive precision sees it as it sees a cancelling Plus:
 * `lost` sizes the next step, `left` must reach p to certify a pass. */
void numeric_note_cancellation_bits(long lost, long left);
Expr* numeric_mpfr_add (const Expr* a, const Expr* b, long default_bits);
Expr* numeric_mpfr_sub (const Expr* a, const Expr* b, long default_bits);
Expr* numeric_mpfr_mul (const Expr* a, const Expr* b, long default_bits);
Expr* numeric_mpfr_div (const Expr* a, const Expr* b, long default_bits);
Expr* numeric_mpfr_pow (const Expr* a, const Expr* b, long default_bits);
/* Real MPFR x^y with mpfr_pow's signature: faithful exp(y log x) for x > 0 and
 * non-integer y, mpfr_pow otherwise (see numeric.c). */
int numeric_mpfr_pow_fast(mpfr_ptr rop, mpfr_srcptr x, mpfr_srcptr y,
                          mpfr_rnd_t rnd);

/* Complex power at MPFR precision via polar form, used by `Power` when
 * `numeric_mpfr_pow` can't apply (negative real base with fractional
 * exponent, complex base, or complex exponent). Returns Complex[MPFR, MPFR],
 * a pure MPFR when the imaginary part rounds to zero, or NULL if either
 * operand isn't MPFR-promotable or `base` is zero. */
Expr* numeric_mpfr_complex_pow(const Expr* base, const Expr* exp,
                               long default_bits);

/* f(I y) = I g(y) on the imaginary axis, for an odd f whose partner g is
 * real there (Erf/Erfi, Sin/Sinh, Tan/Tanh, both ways). When `arg` is
 * Complex[0, y] with an exact zero real part and an MPFR y, returns
 * Complex[0, partner[y]] (the zero part stays exact, as it is in the
 * argument); NULL otherwise -- including a machine y, which keeps the machine
 * rule. */
Expr* numeric_mpfr_imag_axis(const Expr* arg, const char* partner);

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
 * Complex[...] containing one). If found, writes into the caller-provided
 * `out[0..n]` (which must have room for `n` entries) each numericalize(args[i])
 * at the maximum MPFR precision observed among the inexact args (or machine
 * precision if only EXPR_REAL inexactness is present) and returns true. This
 * is how `1. Pi` collapses to `3.14159` — the Pi is numericalized in-line.
 *
 * Returns false (leaving `out` untouched) if no arg is inexact, or if every
 * arg is already a bare machine Real (nothing to do); the caller should then
 * proceed with the original arg list. When true, each written element is
 * caller-owned (free with expr_free); the `out` buffer itself is the caller's
 * (may be stack-allocated). */
bool numeric_contagion_args(Expr* const* args, size_t n, Expr** out);

/* True iff `e` is (or contains, inside Complex[...]) an inexact numeric leaf —
 * an EXPR_REAL or, under USE_MPFR, an EXPR_MPFR value. This is exactly the
 * predicate numeric_contagion_args scans for, exported so a caller can cheaply
 * pre-detect whether the inexact-contagion pass could do anything before paying
 * for it (e.g. Plus's fused special-case guard). */
bool arg_is_inexact(const Expr* e);

/* A numeric builtin given inexact and exact numeric arguments together
 * (Zeta[2., Sqrt[2]], BesselJ[Sqrt[2], 2.], HypergeometricPFQ[{1., 2}, {3},
 * Sqrt[2]/3]) numericalizes the exact ones at the LOWEST precision among the
 * inexact ones (machine precision is the floor), as Plus / Times do. Returns
 * that rewritten call, or NULL when the call has no such mixture, an argument is
 * not a number / numeric quantity / List of them, or an exact argument does not
 * become a number under N. Positions held by NHoldFirst / NHoldRest in `attrs`
 * are neither scanned nor converted. `call` is read-only; the result is the
 * caller's. The evaluator calls it once, when a NumericFunction builtin has
 * declined the original call. */
Expr* numeric_mixed_args_call(const Expr* call, uint32_t attrs);

/* True iff `e` is a numeric quantity in the NumericQ sense: a number, a named
 * numeric constant (Pi, E, I, EulerGamma, ...), or a numeric function (a head
 * with the NumericFunction attribute) applied to numeric-quantity arguments.
 * This is the C predicate behind the NumericQ builtin, exported so comparison
 * heads can gate an exact zero-test on closed-form constant operands (where a
 * free symbol like `x` must NOT be treated as numeric). Read-only. */
bool expr_is_numeric_quantity(const Expr* e);

/* Builtin entry points */
Expr* builtin_n(Expr* res);

/* MachineNumberQ[expr] -> True/False. True iff expr is a finite machine
 * (IEEE double) real, or a Complex of two finite machine reals. Always
 * returns a True/False symbol when arity is 1; returns NULL only on a
 * structurally invalid call (arity != 1). */
Expr* builtin_machinenumberq(Expr* res);

/* inf/nan are libm artifacts, never Mathilda values. The evaluator routes
 * every NumericFunction result through this backstop so a leaked non-finite
 * machine real (a double that overflowed to +/-Inf, or a NaN) is rewritten:
 *   +Inf -> Overflow[]   -Inf -> -Overflow[]   NaN -> Indeterminate
 * numeric_result_has_nonfinite() is the cheap gate (a type check on the
 * common finite path); numeric_sanitize_nonfinite() takes ownership of `e`
 * and returns the rewrite (freeing `e`), or `e` unchanged. Both always
 * compiled — the invariant holds with or without MPFR. */
bool  numeric_result_has_nonfinite(const Expr* e);
Expr* numeric_sanitize_nonfinite(Expr* e);

/* L2 correctness: a head whose machine-real fast path produced a degenerate
 * double (Inf/NaN from overflow, or 0 from underflow) calls this from its own
 * degeneracy branch, with `call` = the builtin's own `res` (head[args...], args
 * still live). It rebuilds the call with every machine-real argument (anywhere,
 * including inside List parameters) promoted to MPFR, evaluates it (reaching the
 * head's MPFR path — the same one N[f[exact]] uses), and returns the finite
 * result at machine precision, or Overflow[]/Indeterminate when the magnitude is
 * past even MPFR's exponent range. Returns NULL when there is nothing to
 * promote, the value is a genuine zero, or the MPFR path itself declines — the
 * caller then keeps its existing behaviour (and the evaluator's inf/nan net
 * remains the final backstop). Without MPFR it is a no-op returning NULL. */
Expr* numeric_promote_machine_call(const Expr* call);

/* Convenience for the shared `get_approx` + libm pure-real branch: wrap a
 * machine result `v` of the unary `call` (input `xin`), promoting a degenerate
 * result via numeric_promote_machine_call. Always compiled; without MPFR it is
 * expr_new_real(v). */
Expr* numeric_machine_real_or_promote(double v, double xin, const Expr* call);

/* Wrap the result `r` of a numeric `call`: if `r` is a degenerate machine real
 * (Inf/NaN, or a flush-to-zero while the call has a non-zero real argument),
 * recompute in MPFR (numeric_promote_machine_call) and return that (freeing
 * `r`); else return `r`. For heads whose numeric core only sees the argument
 * (Erfc/Erfi/ExpIntegralEi/...). Always compiled; a no-op without MPFR. */
Expr* numeric_promote_result_if_degenerate(Expr* r, const Expr* call);

/* Lower a (consumed) 53-bit MPFR / Complex[MPFR, MPFR] value to machine Reals
 * when every part fits a normal double; otherwise return it unchanged (the
 * extended-exponent form). For machine calls a head evaluates through arb.
 * Without MPFR it returns `r`. */
Expr* numeric_lower_mpfr_to_machine(Expr* r);

/* Registers N (and, in Phase 2, Precision/Accuracy/SetPrecision/SetAccuracy
 * and MachinePrecision) with the symbol table. */
void numeric_init(void);

#endif /* NUMERIC_H */
