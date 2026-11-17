/* compile_underflow.h -- detect a machine underflow inside one compiled run.
 *
 * The machine/MPFR boundary policy (docs/design/numeric_tower_review.md): a
 * result that underflows is kept by the interpreter at machine precision in
 * MPFR (1.*^-200 * 1.*^-200 is 1.*10^-400), never flushed to 0. or left
 * subnormal. A compiled run cannot do that, and its final value alone cannot
 * show the underflow (0. may be the right answer), so the IEEE underflow flag
 * is read instead: clear it before the run, test it after, and treat a raised
 * flag -- or a subnormal result -- as a declined run, which falls back to the
 * interpreter.
 *
 * Used only where the compiled value IS the answer the user sees (a
 * CompiledFunction call, Table's elements). Inside NIntegrate / Plot /
 * FindMinimum / NDSolve an underflowed term only feeds a sum or a sample, and
 * a per-sample fallback there would cost speed for no change in the result.
 *
 * <fenv.h> is C99; a platform without FE_UNDERFLOW keeps the old behaviour.
 * The flag is per thread, and these runs execute on the calling thread. */
#ifndef MTH_COMPILE_UNDERFLOW_H
#define MTH_COMPILE_UNDERFLOW_H

#include <fenv.h>
#include <float.h>
#include <math.h>
#include <stdbool.h>

/* Clear the flag only when it is set: IEEE flags are sticky, so in steady
 * state this is one cheap status read, and the (slow, on x86 a full
 * environment rewrite) feclearexcept runs only after an underflow somewhere.
 * Clearing unconditionally per call cost 60 ns, +45% on a compiled call. */
static inline void cvm_underflow_clear(void) {
#ifdef FE_UNDERFLOW
    if (fetestexcept(FE_UNDERFLOW)) feclearexcept(FE_UNDERFLOW);
#endif
}

static inline bool cvm_underflowed(void) {
#ifdef FE_UNDERFLOW
    return fetestexcept(FE_UNDERFLOW) != 0;
#else
    return false;
#endif
}

/* A nonzero double below DBL_MIN: not a machine number under the policy. */
static inline bool cvm_subnormal(double v) {
    return v != 0.0 && fabs(v) < DBL_MIN;
}

#endif /* MTH_COMPILE_UNDERFLOW_H */
