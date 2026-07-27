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
 * Scope of this module: scalar lattice Bool/Int/Real/Complex with arithmetic,
 * comparisons, boolean logic, elementary + special-function kernels, control
 * flow and procedural constructs (M0–M2), plus rank-1 machine arrays delegating
 * to the NDArray subsystem (M3a).  Rank >= 2, Dot/Part, and array locals are
 * later milestones.
 */
#ifndef MATHILDA_COMPILE_H
#define MATHILDA_COMPILE_H

#include <stdbool.h>
#include <stddef.h>
#include <complex.h>
#include "../expr.h"

/* Type lattice.
 *
 * Scalars occupy 0..3, ordered by widening (Bool is not numeric).  Array types
 * (M3) are packed into the same integer above CT_ARR as
 * `CT_ARR + 4*(rank-1) + elem`, so every field that already carries a
 * CompileType — infer_type's result, a register's static type, the declared
 * argument types — carries array types with no parallel plumbing, and the
 * scalar comparisons `t == CT_REAL` / `t < CT_REAL` keep their meanings.
 *
 * CT_ERR is the "no common type" sentinel.  Keeping it in the enum forces a
 * signed underlying type, so both the `(int)t < 0` error checks and the
 * `(int)t >= CT_ARR` array test are well defined. */
typedef enum {
    CT_ERR     = -1,
    CT_BOOL    = 0,
    CT_INT     = 1,
    CT_REAL    = 2,
    CT_COMPLEX = 3,
    CT_ARR     = 4       /* == CT_ARRAY(CT_BOOL, 1); first array encoding */
} CompileType;

/* Highest rank the packed encoding supports (M3a implements rank 1 only). */
#define CT_MAX_RANK 8
#define CT_ARRAY(elem, rank) ((CompileType)((int)CT_ARR + 4 * ((rank) - 1) + (int)(elem)))
#define CT_IS_ARRAY(t)       ((int)(t) >= (int)CT_ARR)
#define CT_ELEM(t)           ((CompileType)(((int)(t) - (int)CT_ARR) & 3))
#define CT_RANK(t)           ((((int)(t) - (int)CT_ARR) >> 2) + 1)

/* A boxed value (compile-time-known type).  For an array type `a` holds an
 * EXPR_NDARRAY node: BORROWED when passed in as an argument (the program never
 * frees an argument array), OWNED by the caller when returned as a result. */
typedef struct {
    CompileType type;
    union { long long i; double r; double _Complex z; unsigned char b; Expr* a; } v;
} CompileValue;

typedef struct CompiledProgram CompiledProgram;

/* Compile `body` (borrowed) as a function of `nargs` argument symbols with the
 * given interned names and declared types.  Returns NULL if any construct is not
 * compilable (caller falls back to the interpreter).  The result is independent
 * of `body` and must be freed with compiled_free.  An argument may be declared
 * `CT_ARRAY(elem, 1)` — a rank-1 machine vector supplied as an EXPR_NDARRAY. */
CompiledProgram* compile_expr(const Expr* body,
                              const char* const* arg_names,
                              const CompileType* arg_types, size_t nargs);

/* Compile-time folding of non-argument symbols that currently hold a machine
 * number as their OwnValue (e.g. the outer iteration variable of a nested
 * Table, which the interpreter binds through the symbol table).  Without it such
 * a symbol is simply uncompilable and the whole body bails.
 *
 * ONLY for programs whose lifetime ends before the symbol can be reassigned —
 * i.e. the throwaway programs built by autocompile.{c,h} inside a single builtin
 * call.  A user `Compile[]` object outlives its defining scope, so it must NOT
 * fold: it keeps the bail-to-interpreter behaviour, which stays correct however
 * the global is later redefined. */
#define COMPILE_FOLD_GLOBALS 0x1u

/* Skip the bytecode optimiser (constant folding, CSE, copy propagation, dead-code
 * elimination, loop-invariant code motion).  For A/B testing only: the optimiser
 * is required to be result-preserving, so any body that answers differently with
 * and without this flag is a bug in a pass. */
#define COMPILE_NO_OPT       0x2u

/* Fuse an elementwise array chain into ONE pass over the buffers, instead of
 * delegating each operation to the NDArray layer (which makes one full-length
 * pass and allocates one temporary buffer per operation).
 *
 * OPT-IN, and off by default, because it is not yet a win: the fused loop runs
 * the scalar VM once per element, and at ~4 ns per bytecode instruction that
 * costs more than the temporary buffers it saves — measured 70 ns/element fused
 * against 61 ns/element delegated for `Total[Sin[v] Exp[-v] + Sqrt[v]]` at
 * length 65536.  The lowering itself is correct, rank-general and parity-tested;
 * what it is waiting for is BLOCK strip-mining (each opcode processing a tile of
 * ~64 elements in a vectorisable C loop, so dispatch is amortised 64x and the
 * temporaries stay in L1).  That is where the order of magnitude is.
 *
 * Rank > 1 does NOT depend on this flag: the delegated NDArray path is already
 * rank-general. */
#define COMPILE_FUSE         0x4u

CompiledProgram* compile_expr_ex(const Expr* body,
                                 const char* const* arg_names,
                                 const CompileType* arg_types, size_t nargs,
                                 unsigned flags);

CompileType compiled_result_type(const CompiledProgram* p);
size_t      compiled_num_args(const CompiledProgram* p);

/* Number of bytecode instructions in the finished program.  Exposed so tests and
 * benchmarks can measure what the optimiser removed; not needed to run a program. */
size_t      compiled_num_instructions(const CompiledProgram* p);

/* True when every argument and the result are CT_REAL with no array temporaries,
 * i.e. compiled_eval_real applies.  Callers that hold a generic program (a user
 * CompiledFunction, say) use this to take the unboxed entry point instead of
 * paying for CompileValue boxing on every call. */
bool        compiled_program_all_real(const CompiledProgram* p);

/* Evaluate with `nargs` boxed argument values (coerced to the declared arg
 * types).  Writes *out.  Returns false if the call could not produce a usable
 * value: a non-finite numeric result, an argument that does not match its
 * declared type, or an array operation that left the promised element type
 * (a real-typed program whose buffer went complex, exactly mirroring the scalar
 * "returns non-finite where the interpreter would go complex" rule).  On false
 * the caller falls back to the interpreter.
 *
 * When the result type is an array, `out->v.a` is a NEW EXPR_NDARRAY the caller
 * owns and must expr_free; every array temporary the program allocated is
 * released before returning, on both the success and failure paths. */
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
