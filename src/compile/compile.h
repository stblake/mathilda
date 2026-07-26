/* Mathilda — Compile[]: numeric compiler engine (scalar core).
 *
 * Compiles a numeric Expr body, as a function of declared-typed argument
 * symbols, into a typed register-machine program that evaluates over raw machine
 * numbers with NO Expr allocation and NO runtime type dispatch (the opcode
 * carries the type; a Real add is one instruction).  Any construct outside the
 * compilable subset makes compilation bail (returns NULL); the caller keeps
 * using the symbolic interpreter.  This is the reusable substrate for NDSolve,
 * Plot, NIntegrate, a user-facing Compile[], etc.  See docs/design/compile.md.
 *
 * Scope of this module (M0): scalar lattice Bool/Int/Real/Complex, arithmetic,
 * comparisons, boolean logic, and elementary functions.  Control flow, arrays,
 * and the generic special-function kernel path are later milestones.
 */
#ifndef MATHILDA_COMPILE_H
#define MATHILDA_COMPILE_H

#include <stdbool.h>
#include <stddef.h>
#include <complex.h>
#include "../expr.h"

/* Scalar type lattice, ordered by widening (Bool is not numeric). */
typedef enum { CT_BOOL = 0, CT_INT = 1, CT_REAL = 2, CT_COMPLEX = 3 } CompileType;

/* A boxed scalar value (compile-time-known type). */
typedef struct {
    CompileType type;
    union { long long i; double r; double _Complex z; unsigned char b; } v;
} CompileValue;

typedef struct CompiledProgram CompiledProgram;

/* Compile `body` (borrowed) as a function of `nargs` argument symbols with the
 * given interned names and declared types.  Returns NULL if any construct is not
 * compilable (caller falls back to the interpreter).  The result is independent
 * of `body` and must be freed with compiled_free. */
CompiledProgram* compile_expr(const Expr* body,
                              const char* const* arg_names,
                              const CompileType* arg_types, size_t nargs);

CompileType compiled_result_type(const CompiledProgram* p);
size_t      compiled_num_args(const CompiledProgram* p);

/* Evaluate with `nargs` boxed argument values (coerced to the declared arg
 * types).  Writes *out.  Returns false if a numeric result is non-finite. */
bool compiled_eval(const CompiledProgram* p, const CompileValue* args, CompileValue* out);

/* Fast path for an all-Real signature (every arg and the result CT_REAL):
 * args/out are plain doubles, no boxing.  Returns false if the signature is not
 * all-Real or the result is non-finite. */
bool compiled_eval_real(const CompiledProgram* p, const double* args, double* out);

/* Evaluate `nprogs` all-Real programs that share the SAME argument layout
 * (compiled with identical arg_names/arg_types), writing out[i].  The shared
 * arguments are loaded once into a single scratch frame (the largest program's),
 * so a system of N components costs O(nargs + total instructions), not O(N·nargs)
 * — the fast path for a vector RHS (e.g. NDSolve).  Returns false on the first
 * non-real signature or non-finite result. */
bool compiled_eval_real_batch(const CompiledProgram* const* progs, size_t nprogs,
                              const double* args, size_t nargs, double* out);

/* Which argument indices the program actually reads (sorted).  Writes up to
 * `cap` entries into `deps`, returns the count.  Used by clients (e.g. NDSolve
 * colored-FD Jacobian) that need the sparsity of the compiled function. */
size_t compiled_arg_deps(const CompiledProgram* p, int* deps, size_t cap);

void compiled_free(CompiledProgram* p);

#endif /* MATHILDA_COMPILE_H */
