/* Mathilda — CompiledFunction: the user-facing wrapper around the numeric
 * compiler (compile.{c,h}).
 *
 * `Compile[argspec, body]` builds a CompiledFunction object (an EXPR_COMPILED
 * atom).  Applying it to numeric arguments runs the compiled bytecode with no
 * Expr allocation; symbolic / non-numeric arguments and bodies outside the
 * compilable subset transparently fall back to the interpreter, so a
 * CompiledFunction always behaves like the function it was built from.
 *
 * The payload is reference-counted and immutable after construction, so
 * Expr-level sharing (expr_copy / expr_unshare) just bumps the refcount.
 */
#ifndef MATHILDA_COMPILED_FUNCTION_H
#define MATHILDA_COMPILED_FUNCTION_H

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>
#include "compile.h"
#include "../expr.h"

typedef struct CompiledFunction CompiledFunction;

/* Build from `Compile[argspec, body]` (both borrowed; deep-copied as needed).
 * `argspec` is a List of either bare symbols (defaulting to Real) or
 * `{sym, typespec}` pairs where typespec is `_Real`/`_Integer`/`_Complex`
 * (Blank[Real] etc.) or the bare type symbol.  Never fails on an uncompilable
 * body — it keeps the body for the interpreter fallback.  Returns NULL only on a
 * malformed argspec (wrong shape / duplicate or non-symbol parameter). */
CompiledFunction* compiled_function_new(const Expr* argspec, const Expr* body);

CompiledFunction* compiled_function_ref(CompiledFunction* cf);
void              compiled_function_free(CompiledFunction* cf);

/* Apply to `nargs` already-evaluated arguments.  Returns a fresh boxed Expr on
 * success, or NULL only when the object cannot handle the call at all (arity
 * mismatch), so the evaluator leaves the application unevaluated.  Numeric args
 * run the bytecode; symbolic args or an uncompilable body use the interpreter
 * fallback (which still returns a value). */
Expr* compiled_function_apply(const CompiledFunction* cf, Expr* const* args, size_t nargs);

/* Accessors for the printer. */
size_t             compiled_function_num_args(const CompiledFunction* cf);
bool               compiled_function_is_compiled(const CompiledFunction* cf);
const Expr*        compiled_function_body(const CompiledFunction* cf);
const char* const* compiled_function_arg_names(const CompiledFunction* cf);
const CompileType* compiled_function_arg_types(const CompiledFunction* cf);

/* Registers the `Compile` builtin (attributes + docstring). Called from
 * core_init. */
void compiled_function_init(void);

#endif /* MATHILDA_COMPILED_FUNCTION_H */
