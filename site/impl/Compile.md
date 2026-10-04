---
source: src/compile/compiled_function.c
---
**Algorithm.** `Compile` is `HoldAll, Protected`, so the body is compiled in its
raw unevaluated form and the argument symbols stay local (a global value for one
does not leak in). `builtin_compile` seeds `RuntimeAttributes`/`RuntimeOptions`
from the registered defaults (so `SetOptions[Compile, …]` takes effect), then
strips trailing options right-to-left — `RuntimeAttributes`, `RuntimeOptions`,
`WorkingPrecision`, `"BigIntegers"` — the rightmost setting winning; an option it
has no meaning for leaves `Compile[…]` unevaluated rather than being silently
ignored. The remaining two arguments (argspec, body) go to
`compiled_function_new`, which parses the spec and type-infers and lowers the body
to typed bytecode; a malformed argspec makes it return `NULL`, again leaving
`Compile[…]` unevaluated. Success wraps the program in an `EXPR_COMPILED` atom.

**Data structures.** The result is an `EXPR_COMPILED` holding a reference-counted,
immutable `CompiledFunction` (a bare pointer in the `Expr` union, so no node-size
growth). The VM is a register machine with three banks — scalar `R`, array handle
`V`, strip-mining tile `T` — and a reusable per-program frame (no per-call
`malloc`). It shares the machine-precision kernel registry (`ndkernels.c`) with
NDSolve/Plot/NIntegrate, so every special function with a machine kernel lowers
through the same `KERN_*` opcodes.

**Complexity / limits.** Applying the object to numeric arguments runs the
bytecode with no `Expr` allocation; a symbolic argument, or a body outside the
compilable subset, transparently falls back to the interpreter, so the value is
always what the uncompiled body would give. Integer arithmetic is exact — an
`int64` overflow abandons the compiled call and the interpreter promotes to a
bignum — unless `RuntimeOptions -> "Speed"` opts into wrapping. `WorkingPrecision`
and `"BigIntegers"` opt into MPFR / GMP arithmetic; both are off by default, so
the machine path is unchanged when neither is given.
