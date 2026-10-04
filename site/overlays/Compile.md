### Worked examples

```mathematica
In[1]:= f = Compile[{{x, _Real}}, x^2 + 1]; f[3.0]  (* a machine-number function applied to a Real *)
```

```mathematica
In[1]:= f[a]  (* a symbolic argument transparently falls back to the interpreter *)
```

```mathematica
In[1]:= Compile[{{n, _Integer}}, n^3][3000000]  (* integer arithmetic stays exact past the int64 range *)
```

```mathematica
In[1]:= Compile[{{z, _Complex}}, z^2][1.0 + 2.0 I]  (* complex machine arithmetic *)
```

```mathematica
In[1]:= Compile[{{x, _Real}}, If[x > 0, 1., -1.], RuntimeAttributes -> Listable][{1., -2., 3.}]  (* the object threads over a list *)
```

### Notes

`Compile[{arg, …}, expr]` builds a `CompiledFunction` that evaluates `expr` over
raw machine numbers, bypassing the symbolic evaluator. Argument types are `_Real`
(the default for a bare symbol), `_Integer` or `_Complex`, with an optional rank
for array parameters (`{v, _Real, 1}`). It is `HoldAll`, so the body compiles in
its raw form and the argument symbols stay local.

A symbolic argument, or a body outside the compilable subset, transparently falls
back to the interpreter (the second example) — so a `CompiledFunction` always
returns what the original expression would. Integer arithmetic is exact: an
`int64` overflow re-runs the body through the interpreter and promotes to a
bignum. Use `CompileDiagnostics` to find out whether, and why, a body does not
lower; `WorkingPrecision -> n` and `"BigIntegers" -> True` opt into MPFR/GMP
arithmetic.
