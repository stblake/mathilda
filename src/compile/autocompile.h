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

/* Compile `body` (borrowed) as a function of `nvars` variable symbols (borrowed
 * EXPR_SYMBOL nodes — e.g. the plot / integration / iteration variable).  All
 * inputs are typed CT_REAL.  Returns NULL if the body does not compile or any
 * var is not a symbol. */
AutoCompiled* autocompile_new(const Expr* body, const Expr* const* vars, size_t nvars);

size_t autocompiled_num_vars(const AutoCompiled* ac);

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

void autocompiled_free(AutoCompiled* ac);

#endif /* MATHILDA_AUTOCOMPILE_H */
