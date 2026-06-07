# DifferenceDelta

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
DifferenceDelta[f, i] gives the forward difference (f /. i -> i+1) - f, the discrete analogue of D. It is the left inverse of indefinite Sum.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= DifferenceDelta[i^2, i]
Out[1]= 1 + 2 i

In[2]:= DifferenceDelta[Sum[k k!, k], k]
Out[2]= -Factorial[k] + Factorial[1 + k]
```

## Implementation notes

**Algorithm.** `DifferenceDelta[f, i]` computes the forward difference
`(f /. i -> i+1) - f`, the discrete analogue of `D` and the left inverse of
indefinite `Sum`. `builtin_differencedelta` (src/sum/sum_gosper.c) requires the
second argument to be a symbol, substitutes `i -> i+1` into `f` (`shift_var`,
implemented via `ReplaceAll`), subtracts the original `f`, and returns the
`Expand`-ed result. Returns NULL (unevaluated) unless the variable is a symbol.

**Data structures.** Plain `Expr*` tree manipulation built on the existing
`ReplaceAll`, subtraction, and `Expand` builtins (`shift_var`, `sum_sub`,
`sum_eval`). No closed-form machinery — it is a thin structural operator that
lives alongside Gosper's summation because the two are inverse operations.

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/sum/sum_gosper.c`](https://github.com/stblake/mathilda/blob/main/src/sum/sum_gosper.c)
- Specification: [`docs/spec/builtins/calculus.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/calculus.md)
