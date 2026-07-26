# Compile[] engine — M1b: user-facing Compile[] / CompiledFunction object

Goal: expose the numeric compiler to the REPL. `Compile[argspec, body]` returns a
`CompiledFunction` object; applying it to numeric args runs the bytecode (Expr-free);
symbolic args / uncompilable bodies fall back to the interpreter.

## Representation
- New `EXPR_COMPILED` atom; union member `struct CompiledFunction* compiled`
  (a pointer — does NOT grow the union, low ABI risk). Own refcount on the
  payload (immutable after build) → expr_copy shares the node, expr_unshare refs
  the payload, expr_free dec-refs.
- `CompiledFunction` owns: `CompiledProgram* prog` (NULL if body uncompilable),
  interned `arg_names`, `arg_types`, `nargs`, and `Expr* body` (fallback + print).

## Tasks
- [ ] expr.h: EXPR_COMPILED enum, forward-decl, union member, `expr_new_compiled`.
- [ ] expr.c: constructor; expr_free / expr_unshare / expr_eq / expr_hash /
      expr_compare cases (identity by payload pointer).
- [ ] src/compile/compiled_function.{h,c}: CompiledFunction struct + new/ref/free/
      apply + accessors; argspec parser ({x}|{{x,_Real}}...); numeric box/unbox;
      interpreter fallback (replace_bindings + eval_and_free).
- [ ] Compile builtin + `compiled_function_init()` (HoldAll, Protected, docstring);
      call from core_init.
- [ ] eval.c: EXPR_COMPILED atomic self-eval; application dispatch (head is
      EXPR_COMPILED → compiled_function_apply).
- [ ] print.c: EXPR_COMPILED → `CompiledFunction[{args}, body]` (or -compiled-).
- [ ] tests/CMakeLists.txt COMMON_SRC += compiled_function.c; tests.
- [ ] docs/spec/builtins + changelog + memory.

## Deferred (later milestones)
- Array/NDArray arg types ({x,_Real,rank}) — M3.
- Compile options (RuntimeAttributes, parallelization) — later.

## Review — DONE (2026-07-26)

Shipped `Compile[]` / `CompiledFunction` (M1b). New `EXPR_COMPILED` atom with a
reference-counted, immutable `CompiledFunction` payload (bare pointer in the
union), wired through expr copy/unshare/free/eq/hash/compare, print (+TeX), and
the evaluator (atomic self-eval + `object[args]` dispatch). `Compile[argspec,
body]` (HoldAll|Protected) parses `_Real`/`_Integer`/`_Complex` argspecs, compiles
the raw body, and keeps it for the interpreter fallback. Application runs the
bytecode for numeric args and boxes the result; symbolic args / uncompilable
bodies (e.g. `Zeta`) fall back to substitution+evaluate; wrong arity stays
unevaluated. `src/compile/compiled_function.{c,h}`;
`tests/test_compiledfunction.c` (all pass); silenced 6 new `-Wswitch` sites for
the enum; `leaks`-clean; `simplify_tests`' 1 failure is pre-existing (verified
against HEAD). Docs + changelog updated.

Next: M3 arrays/NDArray, or wire Plot/NIntegrate/FindRoot/Table to auto-compile.
