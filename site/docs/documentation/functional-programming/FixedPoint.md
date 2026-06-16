# FixedPoint

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
FixedPoint[f, expr]
    starts with expr and applies f repeatedly until the result no longer
    changes, returning the final value.
FixedPoint[f, expr, n]
    stops after at most n applications of f, returning the last value
    obtained even if a fixed point has not been reached.
FixedPoint[f, expr, SameTest -> s]
FixedPoint[f, expr, n, SameTest -> s]
    uses the binary predicate s instead of SameQ to test successive pairs.

FixedPoint[f, expr] gives the last element of FixedPointList[f, expr].
Throw can be used inside f to exit early.
```

## Examples

_No verified examples yet for this function._

## Implementation notes

**Algorithm.** `builtin_fixedpoint` calls `fixedpoint_impl(res, as_list=false)`. It seeds an `ExprBuf` history with a copy of the start value and drives the generic iterator `iter_run` with the `fixedpoint_step` callback. Each step applies `f` to the most recent value (`apply_unary`) and compares the new value to the previous one: by default with structural equality (`expr_eq`), or with a user `SameTest -> g` option (applying `g` and testing for `True`). When they agree the iteration halts (`ITER_STEP_HALT_ADD`); otherwise it continues. A `MaxIterations` count (integer or `Infinity`) is parsed by `parse_fp_opts`; unbounded runs are capped at `ITER_SAFETY_CAP` (exceeding it returns NULL). Throw/Abort/Quit/Return markers produced by `f` are propagated out early. `ebuf_finalize(..., as_list=false)` returns just the final value; the sibling `FixedPointList` returns the whole history.

**Data structures.** `ExprBuf` (a growable `Expr*` vector) is the iteration history; `FixedPointCtx` carries `f`, the optional same-test, and the throw-propagation flag.

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/funcprog.c`](https://github.com/stblake/mathilda/blob/main/src/funcprog.c)
- Specification: [`docs/spec/builtins/functional-programming.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/functional-programming.md)

## Notes & additional examples

### Worked examples

```mathematica
In[1]:= FixedPoint[Floor[#/2] &, 100]
Out[1]= 0
```

```mathematica
In[1]:= FixedPoint[Function[x, (x + 2/x)/2], 1.0]
Out[1]= 1.41421
```

```mathematica
In[1]:= FixedPoint[1 + 1/# &, 1.0]
Out[1]= 1.61803
```

### Notes

`FixedPoint[f, expr]` applies `f` repeatedly until the result stops changing. The second example is Newton's iteration converging to `Sqrt[2]`; the third iterates `x -> 1 + 1/x`, whose attracting fixed point is the golden ratio `(1 + Sqrt[5])/2`. Convergence relies on exact equality of successive results, so use machine reals (or a `SameTest`) for numeric iterations.
