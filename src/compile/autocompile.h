/* Mathilda — autocompile: a thin adapter that lets numeric builtins (Plot,
 * NIntegrate, FindRoot, Table, ...) transparently compile a held body once and
 * evaluate it at many machine-number sample points via the bytecode VM instead
 * of running the full symbolic evaluator per point.
 *
 * Contract: `autocompile_new` returns NULL when the body is outside the
 * compilable subset — the caller then keeps using its interpreter path.  When
 * it succeeds, the compiled result agrees with the interpreter to machine
 * precision at every point where the compiled program returns a finite value.
 * A compiled REAL program returns `false` (non-finite) exactly where the
 * interpreter would instead produce a complex or singular value; callers that
 * need that contribution (integrals, roots) fall back to the interpreter at that
 * point, while callers that exclude non-real points (Plot) simply drop it.
 *
 * All inputs are real doubles.  A body that is genuinely complex-valued for real
 * inputs (e.g. contains an explicit I) compiles to a complex program and is read
 * via autocompiled_eval_complex.
 */
#ifndef MATHILDA_AUTOCOMPILE_H
#define MATHILDA_AUTOCOMPILE_H

#include <stdbool.h>
#include <stddef.h>
#include "../expr.h"
/* Deliberately does NOT include <complex.h> or compile.h: `double _Complex` is a
 * built-in C99 type needing no header, so including this file does not leak the
 * <complex.h> `I` macro into callers (which would clash with any identifier
 * named `I`, e.g. in plot3d.c). */

typedef struct AutoCompiled AutoCompiled;

/* Master switch, mirroring the user-visible $AutoCompilation. Defaults on; off
 * when the environment sets MATHILDA_NO_AUTOCOMPILE (read once, lazily) or the
 * user assigns $AutoCompilation = False.
 *
 * Auto-compilation is invisible by construction -- Plot, Table, NIntegrate and
 * friends compile behind the caller's back purely as an optimisation, and the
 * compiled answer is contracted to agree with the interpreter's. So the switch
 * is not there to change behaviour but to make the two paths comparable: a
 * differential run flips it and diffs, and a user who suspects the compiled path
 * has it wrong can turn it off and confirm in one line. */
bool autocompile_enabled(void);
void autocompile_set_enabled(bool on);

/* Compile `body` (borrowed) as a function of `nvars` variable symbols (borrowed
 * EXPR_SYMBOL nodes — e.g. the plot / integration / iteration variable).  All
 * inputs are typed CT_REAL.  Returns NULL if the body does not compile or any
 * var is not a symbol. */
AutoCompiled* autocompile_new(const Expr* body, const Expr* const* vars, size_t nvars);

/* As autocompile_new, but the variables are declared CT_COMPLEX — the program
 * takes complex INPUTS, not merely a complex result.  For samplers whose
 * variable ranges over a region of the plane rather than an interval
 * (ComplexPlot), where binding a real and hoping for a complex answer would be
 * the wrong function.
 *
 * Compiles a strictly smaller set of bodies than the real version: a head needs
 * a complex kernel, not just a real one, and many special functions have only
 * the latter (see NUMERIC_FUNCTION_MISSING.md).  Read it back with
 * autocompiled_eval_z. */
AutoCompiled* autocompile_new_z(const Expr* body, const Expr* const* vars, size_t nvars);

/* As autocompile_new, but compiles the body in MPFR at `prec_bits` bits (real
 * inputs), for a sampler running at a high WorkingPrecision.  The compilable
 * subset is the arbitrary-precision one (straight-line arithmetic + elementary
 * functions; see docs/design/compile-arbitrary-precision.md).  Returns NULL for
 * prec_bits <= 0 (use autocompile_new for the machine path) or a non-compilable
 * body.  Read back with autocompiled_eval_mpfr. */
AutoCompiled* autocompile_new_prec(const Expr* body, const Expr* const* vars,
                                   size_t nvars, long prec_bits);

/* As autocompile_new_prec, but the variables are declared COMPLEX (MPFR): the
 * program takes complex INPUTS at `prec_bits` bits.  For a sampler that evaluates
 * the body on a contour in the plane at high precision (e.g. NSum's
 * Euler–Maclaurin derivative-by-Cauchy step).  Compiles a strictly smaller set of
 * bodies than the real version (a head needs a complex ncpx kernel).  Read back
 * with autocompiled_eval_mpfr, passing Complex[MPFR,MPFR] / real sample points. */
AutoCompiled* autocompile_new_prec_z(const Expr* body, const Expr* const* vars,
                                     size_t nvars, long prec_bits);

/* Evaluate a program built by autocompile_new_prec at the numeric sample points
 * `xs[nvars]` (each an EXPR_MPFR / Integer / Real / Rational, borrowed).  Returns
 * a freshly allocated result Expr (EXPR_MPFR, or Complex[MPFR,MPFR]) the caller
 * owns, or NULL when the caller should interpret this point (non-finite result,
 * or a non-managed program). */
Expr* autocompiled_eval_mpfr(const AutoCompiled* ac, const Expr* const* xs);

size_t autocompiled_num_vars(const AutoCompiled* ac);

/* True when the program has an all-real SIGNATURE and a CT_REAL result type, so
 * every value it produces is a machine real and autocompiled_eval_boxed can only
 * ever return expr_new_real(...) or NULL.
 *
 * The gate a collection builtin needs before filling a float64 buffer directly:
 * without it a body like If[p, 1, 2] or a complex-valued one would have its
 * element HEADS changed by the buffer (Integer or Complex becoming Real), which
 * is user-visible. A `false` result still means "fall back for this element",
 * so a producer using this must handle abandonment, not assume success. */
bool autocompiled_result_is_real(const AutoCompiled* ac);

/* Evaluate at real inputs xs[nvars].  Writes a real *out and returns true, or
 * returns false when the result is non-finite OR not real-valued (a complex
 * result with nonzero imaginary part).  Callers that only handle reals treat
 * false as "no usable value at this point". */
bool autocompiled_eval_real(const AutoCompiled* ac, const double* xs, double* out);

/* Evaluate at real inputs xs[nvars] → complex *out (a real/integer result is
 * returned with zero imaginary part).  Returns false only on a non-finite
 * result. */
bool autocompiled_eval_complex(const AutoCompiled* ac, const double* xs, double _Complex* out);

/* Evaluate at real inputs xs[nvars] → a freshly allocated Expr of the result's
 * own type, or NULL when the caller should fall back to the interpreter for this
 * point.  Unlike autocompiled_eval_real this does NOT flatten an integer result
 * to a double: a body like If[p, 1, 2] yields Integer 1, matching the
 * interpreter exactly.  For collection builtins (Table) where the element type
 * is user-visible, prefer this over autocompiled_eval_real. */
Expr* autocompiled_eval_boxed(const AutoCompiled* ac, const double* xs);

/* Evaluate a program built by autocompile_new_z at complex inputs zs[nvars].
 * Returns false on a non-finite result, or if `ac` is a real-argument program.
 * (Contrast autocompiled_eval_complex, whose INPUTS are real and only whose
 * result may be complex.) */
bool autocompiled_eval_z(const AutoCompiled* ac, const double _Complex* zs,
                         double _Complex* out);

void autocompiled_free(AutoCompiled* ac);

#endif /* MATHILDA_AUTOCOMPILE_H */
