### Worked examples

```mathematica
In[1]:= CompilePrint[Compile[{{x, _Real}}, x^2 + 2.5 x + 1]]  (* prints the bytecode to the terminal; the call itself returns Null *)
```

```mathematica
In[1]:= CompilePrint[42]  (* anything that is not a CompiledFunction is left unevaluated *)
```

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
