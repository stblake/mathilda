# RootMeanSquare

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
RootMeanSquare[list] gives the root mean square of values in list.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= RootMeanSquare[{a, b, c, d}]
Out[1]= 1/2 Sqrt[a^2 + b^2 + c^2 + d^2]

In[2]:= RootMeanSquare[{{1, 2}, {5, 10}, {5, 2}, {4, 8}}]
Out[2]= {1/2 Sqrt[67], Sqrt[43]}

In[3]:= RootMeanSquare[{1, 2, 3, 4}]
Out[3]= Sqrt[15/2]

In[4]:= RootMeanSquare[{Pi, E, 2}]
Out[4]= Sqrt[1/3 (4 + E^2 + Pi^2)]

In[5]:= RootMeanSquare[{1., 2., 3., 4.}]
Out[5]= 2.73861
```

## Implementation notes

**Algorithm.** `builtin_rootmeansquare` computes `Sqrt[Mean[x^2]]`. It reduces matrices column-wise via `apply_columnwise` and requires a `List`. For data containing a real, it sums squares in `double` and returns `expr_new_real(sqrt(sum_sq/n))`. For exact/symbolic data it builds `Plus[x_i^2 ...]`, then carefully distributes the square root to keep results exact: if the summed result is non-numeric and `n` is a perfect square it factors out `1/Sqrt[n]`; if the mean-square is a rational with a perfect-square denominator it pulls that root out before applying `Power[..., 1/2]`. Otherwise it returns `Power[meanSq, 1/2]` for the evaluator to simplify. `ATTR_PROTECTED`.

- `Protected`.
- Gives the square root of the second sample moment.
- For a list `{x1, x2, ...}`, it computes `Sqrt[1/n Total[{x1^2, x2^2, ...}]]`.
- Handles both numerical and symbolic data.
- Works column-wise on matrices.

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/stats.c`](https://github.com/stblake/mathilda/blob/main/src/stats.c)
- Specification: [`docs/spec/builtins/statistics.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/statistics.md)
