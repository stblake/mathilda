### Worked examples

```mathematica
In[1]:= CompileDiagnostics[{{x, _Real}}, Sin[x] + x^2]  (* a body that lowers fully *)
```

```mathematica
In[1]:= CompileDiagnostics[{{x, _Real}}, Sin[x] + y]  (* y is not a declared argument and holds no machine value *)
```

```mathematica
In[1]:= CompileDiagnostics[{{x, _Real}}, Sin[x] + Integrate[x, x]]  (* names the innermost head with no machine lowering *)
```

### Notes

`CompileDiagnostics[argspec, expr]` reports whether `expr` compiles for the given
`Compile[]` argument specification, and if not, **which subexpression stopped
it**. `argspec` is exactly what `Compile` takes, and it is `HoldAll`.

This exists because a bail is otherwise invisible: the caller quietly interprets,
the answer is still correct, and the only symptom is being much slower — and the
cost is not proportional, because the compilable subset is a *cliff*, so one
unsupported head costs the entire body. On success the result also carries the
`"ResultType"` and the instruction count with and without the optimiser; on
failure it gives the `"Reason"` and the innermost `"Subexpression"`. It accepts
the same `WorkingPrecision -> n` / `"BigIntegers" -> True` options, so it can also
report whether the arbitrary-precision subset lowers.
