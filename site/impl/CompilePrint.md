---
source: src/compile/compiled_function.c
---
**Algorithm.** `CompilePrint[cf]` is `Protected` but — deliberately — **not**
`HoldAll`, unlike `Compile` and `CompileDiagnostics`: the argument has to evaluate
down to the `EXPR_COMPILED` atom, so both `CompilePrint[Compile[…]]` and
`f = Compile[…]; CompilePrint[f]` work. `builtin_compile_print` requires one
argument that is an `EXPR_COMPILED`; anything else returns `NULL` and is left
unevaluated. For a compiled object it calls `compiled_function_disassemble`, prints
the returned text to stdout, frees it, and returns `Null`.

**Data structures.** The disassembler walks the `CompiledFunction`'s bytecode and
emits a header (signature, per-argument register and type, result register, the
three register-bank sizes, instruction/CSE/parallel-loop counts) followed by one
line per instruction with raw operands on the left and a readable rendering on the
right. Registers are named by bank (`R` scalar, `V` array, `T` tile); machine
kernels are resolved back to their symbol names and callee programs and parallel
loops are numbered, so no addresses appear and two versions of a body diff cleanly.

**Complexity / limits.** Linear in the program size. Where `CompileDiagnostics`
reports *how much* code there is, `CompilePrint` shows *which* — the only way to
confirm the optimiser folded constants (`ADD_RK`/`MUL_RK` rather than `CONST`),
that an array chain fused, or that a map fanned out. An object whose body did not
compile has no bytecode, so it instead prints the bail reason and offending
subexpression — the same answer `CompileDiagnostics` gives about a body, asked
about an object.
