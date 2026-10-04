# CompileDiagnostics

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`CompileDiagnostics[argspec, expr] reports whether expr compiles for the given Compile[] argument specification, and if not, the innermost subexpression that could not be lowered. Accepts the same WorkingPrecision -> n / "BigIntegers" -> True options as Compile[], so it also reports whether the arbitrary-precision subset lowers (ResultType MPFRReal/MPFRComplex/BigInteger). For a compiled body it also gives the result type and the instruction count with and without the optimiser.`**

## Examples (3)

Every input below was run against the current Mathilda build and its output recorded.

### Applications (3)

A body that lowers fully

```mathematica
In[1]:= CompileDiagnostics[{{x, _Real}}, Sin[x] + x^2]
Out[1]= {"Compiled" -> True, "ResultType" -> "Real", "Instructions" -> 4, "CommonSubexpressions" -> 0, "InstructionsUnoptimized" -> 4}
```

Y is not a declared argument and holds no machine value

```mathematica
In[2]:= CompileDiagnostics[{{x, _Real}}, Sin[x] + y]
Out[2]= {"Compiled" -> False, "Reason" -> "symbol is not a declared argument and holds no machine value", "Subexpression" -> "y"}
```

Names the innermost head with no machine lowering

```mathematica
In[3]:= CompileDiagnostics[{{x, _Real}}, Sin[x] + Integrate[x, x]]
Out[3]= {"Compiled" -> False, "Reason" -> "no machine lowering for this head at these argument types", "Subexpression" -> "Integrate[x, x]"}
```

## Implementation notes

**Algorithm.** `CompileDiagnostics[argspec, expr]` is `HoldAll, Protected` and
takes exactly what `Compile` takes. `builtin_compile_diagnostics` strips the same
trailing `WorkingPrecision`/`"BigIntegers"` options (so it can report whether the
*arbitrary-precision* subset lowers), parses the argspec, and calls the same
`compile_expr_prec` the compiler uses. On success it returns a rule list carrying
`"Compiled" -> True`, the `"ResultType"` (from `ct_name`: `Real`, `Integer`,
`Complex`, `Boolean`, or `MPFRReal`/`MPFRComplex`/`BigInteger` under the managed
options), the `"Instructions"` count, the `"CommonSubexpressions"` the optimiser
hoisted, and `"InstructionsUnoptimized"` — obtained by recompiling with
`COMPILE_NO_OPT`, because the optimiser rewrites in place and the pre-pass count
is otherwise gone. On failure it returns `"Compiled" -> False`, a `"Reason"`
(`compiled_bail_reason`) and the `"Subexpression"` (`compiled_bail_expr`) — the
*innermost* node that could not be lowered, not the construct that contains it.

**Data structures.** A small fixed `Expr*` array of `Rule[key, value]` pairs built
by the `diag_rule` helper and wrapped in a `List`; no compiled program survives
the call (both the optimised and unoptimised programs are freed after their counts
are read).

**Complexity / limits.** It exists because a bail is otherwise invisible: the
caller quietly interprets and the answer stays correct but an order of magnitude
slower, and the compilable subset is a *cliff* — one unsupported head costs the
whole body. A malformed argspec leaves the call unevaluated, exactly as `Compile`
does. `MATHILDA_COMPILE_DIAG=1` prints the same report to stderr whenever an
auto-compiled numeric builtin falls back.

**Attributes:** `HoldAll`, `Protected`.

## References

**See also:** [Compile](../../control-flow/Compile/), [HoldAll](../../expression-information/HoldAll/), [Plot](../../graphics/Plot/), [NIntegrate](../../numerical-calculus/NIntegrate/), [NSum](../../numerical-calculus/NSum/), [ContourPlot](../../graphics/ContourPlot/)

- Source: [`src/compile/compiled_function.c`](https://github.com/stblake/mathilda/blob/main/src/compile/compiled_function.c)
- Specification: [`docs/spec/builtins/control-flow.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/control-flow.md)
- Tests: [`tests/test_compile_arbprec.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compile_arbprec.c)
- Tests: [`tests/test_compile_assoc.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compile_assoc.c)
- Tests: [`tests/test_compile_linalg.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compile_linalg.c)
- Tests: [`tests/test_compile_transforms.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compile_transforms.c)

## Notes & additional examples

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
