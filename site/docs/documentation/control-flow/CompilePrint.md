# CompilePrint

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`CompilePrint[cf] prints the bytecode of the CompiledFunction cf: its argument and result registers with their types, the scalar/array/tile register banks, and one line per instruction giving both the raw operands and a readable rendering. For an object whose body did not compile it reports the bail reason instead. Returns Null.`**

## Examples (2)

Every input below was run against the current Mathilda build and its output recorded.

### Applications (2)

Prints the bytecode to the terminal; the call itself returns Null

```mathematica
In[1]:= CompilePrint[Compile[{{x, _Real}}, x^2 + 2.5 x + 1]]
```

Anything that is not a CompiledFunction is left unevaluated

```mathematica
In[2]:= CompilePrint[42]
Out[2]= CompilePrint[42]
```

## Implementation notes

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

**Attributes:** `Protected`.

## References

**See also:** [CompileDiagnostics](../../control-flow/CompileDiagnostics/), [HoldAll](../../expression-information/HoldAll/), [Compile](../../control-flow/Compile/), [Attributes](../../expression-information/Attributes/)

- Source: [`src/compile/compiled_function.c`](https://github.com/stblake/mathilda/blob/main/src/compile/compiled_function.c)
- Specification: [`docs/spec/builtins/control-flow.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/control-flow.md)
- Tests: [`tests/test_compile_assoc.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compile_assoc.c)
- Tests: [`tests/test_compiledfunction.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compiledfunction.c)

## Notes & additional examples

### Notes

`CompilePrint[cf]` prints the bytecode of a `CompiledFunction`: the signature,
each argument's register and declared type, the result register, the scalar /
array / tile register banks, and one line per instruction giving both the raw
operands and a readable rendering. It returns `Null`.

Unlike `Compile` and `CompileDiagnostics` it is deliberately **not** `HoldAll`, so
its argument evaluates down to the compiled object — both `CompilePrint[Compile[
…]]` and `f = Compile[…]; CompilePrint[f]` work. Where `CompileDiagnostics` says
how *much* code there is, this says *which* — the way to confirm the optimiser
folded constants or an array chain fused. An object whose body did not compile has
no bytecode, so it prints the bail reason instead.

A printed example (from `CompilePrint[Compile[{{x, _Real}}, x^2 + 2.5 x + 1]]`):

```text
Signature   CompiledFunction[{x : Real}, x^2 + 2.5 x + 1]
Arguments   1
              R0   : Real         x
Result      R1 : Real
Registers   3 scalar, 0 array, 0 tile   (frame 3 slots)
Program     5 instructions, 0 CSE, all-Real fast path

    0  POWI_R     R1, R0, 2                       R1 = R0^2
    1  MUL_RK     R2, R0, 2.5                     R2 = R0 * 2.5
    2  ADD_R      R1, R1, R2                      R1 = R1 + R2
    3  ADD_RK     R1, R1, 1                       R1 = R1 + 1
    4  RET        R1                              return R1
```
