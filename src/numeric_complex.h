/* Mathilda — MPFR-complex transcendental helpers.
 *
 * Mathilda represents complex numbers structurally as `Complex[re, im]`
 * (an EXPR_FUNCTION whose head is the symbol "Complex"). When the
 * components are MPFR values, every transcendental builtin should be
 * able to compute its result at MPFR precision rather than silently
 * coercing the pair to double + libc's complex math (`csin`, `cexp`,
 * `cpow`, etc.).
 *
 * This module owns the per-function MPFR-complex routines and a small
 * dispatch helper so per-builtin call sites stay one line. The
 * existing real-only `numeric_mpfr_apply_unary` (numeric.h) returns
 * NULL for any Complex[..] input; the complex-aware analogue lives
 * here so the responsibility split stays clean.
 *
 * Identities used internally:
 *   exp(a + bi) = exp(a)(cos(b) + i sin(b))
 *   log(a + bi) = log(hypot(a, b)) + i atan2(b, a)
 *   sin(a + bi) = sin(a) cosh(b) + i cos(a) sinh(b)
 *   cos(a + bi) = cos(a) cosh(b) - i sin(a) sinh(b)
 *   tan(a + bi) = (sin(2a) + i sinh(2b)) / (cos(2a) + cosh(2b))
 *                 (numerically stable; avoids cosh overflow for large |b|)
 *   sinh(z)     = -i sin(i z)         (via the swap a <-> -b, b <-> a)
 *   cosh(z)     =       cos(i z)
 *   tanh(z)     = -i tan(i z)
 *   asin(z)     = -i log(i z + sqrt(1 - z^2))
 *   acos(z)     = pi/2 - asin(z)
 *   atan(z)     = (i/2)(log(1 - i z) - log(1 + i z))
 *   asinh(z)    = log(z + sqrt(z^2 + 1))
 *   acosh(z)    = log(z + sqrt((z - 1)(z + 1)))   (branch-stable form)
 *   atanh(z)    = (1/2) log((1 + z)/(1 - z))
 *
 * All ops respect MPFR_RNDN throughout. Output precision is taken from
 * the caller's already-init2'd `out_re`/`out_im` MPFRs; ops never resize
 * the output. Internal scratch MPFRs are init2'd at the same precision
 * as `out_re`.
 */
#ifndef NUMERIC_COMPLEX_H
#define NUMERIC_COMPLEX_H

#ifdef USE_MPFR

#include "expr.h"
#include <mpfr.h>
#include <stdbool.h>

/* Construct a result expression from an (out_re, out_im) MPFR pair. If
 * `out_im` rounds to exactly zero, returns a pure EXPR_MPFR carrying
 * `out_re`; otherwise returns Complex[EXPR_MPFR, EXPR_MPFR]. Always
 * COPIES the inputs — caller still owns and must clear `out_re` /
 * `out_im` after the call. Precision of the result matches the inputs. */
Expr* numeric_mpfr_make_complex(const mpfr_t out_re, const mpfr_t out_im);

/* Signature shared by every MPFR-complex unary op below.
 *
 *   void op(mpfr_t out_re, mpfr_t out_im,
 *           const mpfr_t in_re, const mpfr_t in_im);
 *
 * The op writes into `out_re`/`out_im` (already init2'd by the caller at
 * the working precision). It must not alias `in_re`/`in_im` with the
 * outputs in a way that breaks under MPFR aliasing rules — the
 * implementations below use their own scratch MPFRs for any value
 * that needs to outlive an mpfr_* call. */
typedef void (*MpfrComplexUnaryOp)(mpfr_t, mpfr_t,
                                   const mpfr_t, const mpfr_t);

/* Apply an MPFR-complex op to expression `e`. Returns a freshly
 * allocated Expr (EXPR_MPFR if the result's imaginary part rounds to
 * zero, otherwise Complex[EXPR_MPFR, EXPR_MPFR]) on success, NULL if
 * `e` is not MPFR-promotable (e.g. it carries a free symbol).
 *
 * Working precision is `numeric_combined_bits(e, NULL, default_bits)`
 * — the same convention used by numeric_mpfr_apply_unary, so MPFR
 * precision is inherited from the input when default_bits == 0.
 *
 * Caller owns the returned Expr*. */
Expr* numeric_mpfr_apply_complex_unary(const Expr* e, long default_bits,
                                       MpfrComplexUnaryOp op);

/* MPFR-complex division: (out_re, out_im) = (a_re + a_im i) / (b_re + b_im i)
 * using the Smith algorithm for numerical stability. Outputs are
 * already init2'd. Defined here because the inverse trig / reciprocal
 * trig (Cot/Sec/Csc and the corresponding Arc* variants) compose
 * naturally via "1 / forward". */
void mpfr_complex_div(mpfr_t out_re, mpfr_t out_im,
                      const mpfr_t a_re, const mpfr_t a_im,
                      const mpfr_t b_re, const mpfr_t b_im);

/* Concrete unary ops — one per transcendental. Each implements the
 * identity listed at the top of this file. */
void mpfr_complex_exp  (mpfr_t out_re, mpfr_t out_im, const mpfr_t in_re, const mpfr_t in_im);
void mpfr_complex_log  (mpfr_t out_re, mpfr_t out_im, const mpfr_t in_re, const mpfr_t in_im);
void mpfr_complex_sin  (mpfr_t out_re, mpfr_t out_im, const mpfr_t in_re, const mpfr_t in_im);
void mpfr_complex_cos  (mpfr_t out_re, mpfr_t out_im, const mpfr_t in_re, const mpfr_t in_im);
void mpfr_complex_tan  (mpfr_t out_re, mpfr_t out_im, const mpfr_t in_re, const mpfr_t in_im);
void mpfr_complex_sinh (mpfr_t out_re, mpfr_t out_im, const mpfr_t in_re, const mpfr_t in_im);
void mpfr_complex_cosh (mpfr_t out_re, mpfr_t out_im, const mpfr_t in_re, const mpfr_t in_im);
void mpfr_complex_tanh (mpfr_t out_re, mpfr_t out_im, const mpfr_t in_re, const mpfr_t in_im);
void mpfr_complex_asin (mpfr_t out_re, mpfr_t out_im, const mpfr_t in_re, const mpfr_t in_im);
void mpfr_complex_acos (mpfr_t out_re, mpfr_t out_im, const mpfr_t in_re, const mpfr_t in_im);
void mpfr_complex_atan (mpfr_t out_re, mpfr_t out_im, const mpfr_t in_re, const mpfr_t in_im);
void mpfr_complex_asinh(mpfr_t out_re, mpfr_t out_im, const mpfr_t in_re, const mpfr_t in_im);
void mpfr_complex_acosh(mpfr_t out_re, mpfr_t out_im, const mpfr_t in_re, const mpfr_t in_im);
void mpfr_complex_atanh(mpfr_t out_re, mpfr_t out_im, const mpfr_t in_re, const mpfr_t in_im);

/* Reciprocal trig / hyperbolic composites — implemented as 1 / forward.
 * Provided as named ops so callers stay one-liners. */
void mpfr_complex_cot  (mpfr_t out_re, mpfr_t out_im, const mpfr_t in_re, const mpfr_t in_im);
void mpfr_complex_sec  (mpfr_t out_re, mpfr_t out_im, const mpfr_t in_re, const mpfr_t in_im);
void mpfr_complex_csc  (mpfr_t out_re, mpfr_t out_im, const mpfr_t in_re, const mpfr_t in_im);
void mpfr_complex_coth (mpfr_t out_re, mpfr_t out_im, const mpfr_t in_re, const mpfr_t in_im);
void mpfr_complex_sech (mpfr_t out_re, mpfr_t out_im, const mpfr_t in_re, const mpfr_t in_im);
void mpfr_complex_csch (mpfr_t out_re, mpfr_t out_im, const mpfr_t in_re, const mpfr_t in_im);
void mpfr_complex_acot (mpfr_t out_re, mpfr_t out_im, const mpfr_t in_re, const mpfr_t in_im);
void mpfr_complex_asec (mpfr_t out_re, mpfr_t out_im, const mpfr_t in_re, const mpfr_t in_im);
void mpfr_complex_acsc (mpfr_t out_re, mpfr_t out_im, const mpfr_t in_re, const mpfr_t in_im);
void mpfr_complex_acoth(mpfr_t out_re, mpfr_t out_im, const mpfr_t in_re, const mpfr_t in_im);
void mpfr_complex_asech(mpfr_t out_re, mpfr_t out_im, const mpfr_t in_re, const mpfr_t in_im);
void mpfr_complex_acsch(mpfr_t out_re, mpfr_t out_im, const mpfr_t in_re, const mpfr_t in_im);

/* --------------------------------------------------------------------
 *  ncpx — a raw paired-mpfr_t complex value with explicit-working-precision
 *  arithmetic.
 *
 *  This is the tight-inner-loop counterpart to the Expr-level
 *  `mpfr_complex_*` ops above. Special-function kernels that sum series or
 *  asymptotic expansions (BesselJ today; BesselY/BesselI/BesselK and friends
 *  in future) build directly on `ncpx` instead of re-rolling a file-local
 *  complex struct — the role played by the older `acx`/`ecx`/`gcx` toolkits
 *  in airyai.c / erf.c / gamma.c, hoisted here so siblings share one copy.
 *
 *  Each value is `mpfr_init2`'d at a precision (`ncpx_init`); binary ops that
 *  need scratch take an explicit working precision `wp`. Every op is
 *  alias-safe (the output may alias an input). Rounding is MPFR_RNDN. The
 *  memory contract matches the rest of this file: every ncpx_init pairs with
 *  an ncpx_clear; no allocation survives a call. */
typedef struct { mpfr_t re, im; } ncpx;

void ncpx_init (ncpx* z, mpfr_prec_t p);
void ncpx_clear(ncpx* z);
void ncpx_set   (ncpx* d, const ncpx* s);            /* d = s             */
void ncpx_set_d (ncpx* d, double re, double im);     /* d = re + i*im     */
void ncpx_set_ui(ncpx* d, unsigned long re);         /* d = re (real)     */

void ncpx_add  (ncpx* out, const ncpx* a, const ncpx* b);
void ncpx_sub  (ncpx* out, const ncpx* a, const ncpx* b);
void ncpx_neg  (ncpx* out, const ncpx* a);
void ncpx_mul  (ncpx* out, const ncpx* a, const ncpx* b, mpfr_prec_t wp);
void ncpx_div  (ncpx* out, const ncpx* a, const ncpx* b, mpfr_prec_t wp);
void ncpx_scale(ncpx* out, const ncpx* z, const mpfr_t s);   /* out = z*s, real s */

void ncpx_abs(mpfr_t mag, const ncpx* z);   /* mag = |z|              */
void ncpx_arg(mpfr_t out, const ncpx* z);   /* out = arg z, (-pi, pi] */

void ncpx_exp  (ncpx* out, const ncpx* z, mpfr_prec_t wp);
void ncpx_log  (ncpx* out, const ncpx* z, mpfr_prec_t wp);
void ncpx_sin  (ncpx* out, const ncpx* z, mpfr_prec_t wp);
void ncpx_cos  (ncpx* out, const ncpx* z, mpfr_prec_t wp);
void ncpx_sqrt (ncpx* out, const ncpx* z, mpfr_prec_t wp);
void ncpx_pow_d(ncpx* out, const ncpx* z, double e, mpfr_prec_t wp);          /* z^e, real e        */
void ncpx_pow  (ncpx* out, const ncpx* z, const ncpx* w, mpfr_prec_t wp);     /* z^w = exp(w log z) */

#endif /* USE_MPFR */

#endif /* NUMERIC_COMPLEX_H */
