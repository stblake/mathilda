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
 * malformed argspec (wrong shape / duplicate or non-symbol parameter).
 *
 * `runtime_attrs` carries the `RuntimeAttributes` option as an attribute
 * bitmask; the only setting the object acts on is ATTR_LISTABLE (pass ATTR_NONE
 * for the default).  A Listable object threads over List arguments the way a
 * Listable symbol does — see compiled_function_apply.
 *
 * `compile_flags` are compile-engine flags (compile.h) from the `RuntimeOptions`
 * option (COMPILE_WRAP_INT) and the `"BigIntegers"` option (COMPILE_BIGINT).
 * `prec_bits` is the MPFR working precision from `WorkingPrecision -> n` (0 for
 * the default machine precision).  Pass 0 flags / 0 prec_bits for the defaults. */
CompiledFunction* compiled_function_new(const Expr* argspec, const Expr* body,
                                        uint32_t runtime_attrs, unsigned compile_flags,
                                        long prec_bits);

CompiledFunction* compiled_function_ref(CompiledFunction* cf);
void              compiled_function_free(CompiledFunction* cf);

/* Apply to `nargs` already-evaluated arguments (borrowed).  Returns a fresh
 * boxed Expr on success, or NULL only when the object cannot handle the call at
 * all (arity mismatch, or Listable threading over unequal-length lists), so the
 * evaluator leaves the application unevaluated.  Numeric args run the bytecode;
 * symbolic args or an uncompilable body use the interpreter fallback (which
 * still returns a value).
 *
 * A `RuntimeAttributes -> Listable` object first threads over any argument that
 * is a List deeper than that parameter's declared rank, so `cf[{1., 2.}]`
 * gives `{cf[1.], cf[2.]}`. */
Expr* compiled_function_apply(const CompiledFunction* cf, Expr* const* args, size_t nargs);

/* Accessors for the printer. */
size_t             compiled_function_num_args(const CompiledFunction* cf);
bool               compiled_function_is_compiled(const CompiledFunction* cf);
/* The RuntimeAttributes bitmask the object was built with (ATTR_LISTABLE or
 * ATTR_NONE). */
uint32_t           compiled_function_attributes(const CompiledFunction* cf);
const Expr*        compiled_function_body(const CompiledFunction* cf);
const char* const* compiled_function_arg_names(const CompiledFunction* cf);
const CompileType* compiled_function_arg_types(const CompiledFunction* cf);

/* The compiled program, or NULL if the body was outside the compilable subset.
 * Exposed so a compiled CALLER can invoke this callee directly through the VM
 * instead of routing every call back through the evaluator. */
const CompiledProgram* compiled_function_program(const CompiledFunction* cf);

/* The object's signature followed by compiled_disassemble of its program — what
 * CompilePrint[] prints.  For an object whose body was outside the compilable
 * subset there is no program to show, so it reports the bail reason and the
 * subexpression that caused it instead: "why does this run interpreted" is the
 * same question, asked of an object rather than of a body.
 *
 * Caller owns the returned string and must free() it; NULL on allocation
 * failure.  Implemented in disasm.c. */
char* compiled_function_disassemble(const CompiledFunction* cf);

/* Registers the `Compile` builtin (attributes + docstring). Called from
 * core_init. */
void compiled_function_init(void);

#endif /* MATHILDA_COMPILED_FUNCTION_H */
