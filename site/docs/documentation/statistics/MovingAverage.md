# MovingAverage

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
MovingAverage[list, r]
    gives the moving average of list, computed by averaging runs of r elements.
MovingAverage[list, {w_1, w_2, ..., w_r}]
    gives the weighted moving average of list with weights w_i (effective weights w_i / Sum[w_i]).
MovingAverage returns a list of length Length[list] - r + 1, and stays unevaluated when r < 1 or r > Length[list].
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= MovingAverage[{1, 5, 7, 3, 6, 2}, 3]
Out[1]= {13/3, 5, 16/3, 11/3}

In[2]:= MovingAverage[{1.2, 5.2, 3.4, 4.5, 2.3, 4.5}, 3]
Out[2]= {3.26667, 4.36667, 3.4, 3.76667}

In[3]:= MovingAverage[{a, b, c, d, e}, 2]
Out[3]= {1/2 (a + b), 1/2 (b + c), 1/2 (c + d), 1/2 (d + e)}

In[4]:= MovingAverage[{a, b, c, d, e}, {1, 2}]
Out[4]= {1/3 a + 2/3 b, 1/3 b + 2/3 c, 1/3 c + 2/3 d, 1/3 d + 2/3 e}

In[5]:= MovingAverage[{2^100, 2^101, 2^102, 2^103}, 2]
Out[5]= {1901475900342344102245054808064, 3802951800684688204490109616128, 7605903601369376408980219232256}

In[6]:= MovingAverage[{1, 2, 3, 4, 5}, 6]
Out[6]= MovingAverage[{1, 2, 3, 4, 5}, 6]
```

## Implementation notes

**Algorithm.** `builtin_moving_average` takes `(list, r)` where `r` is a positive integer window (`EXPR_INTEGER`/`EXPR_BIGINT`) or a `List` of weights. Output length is `n - r + 1`, and the call stays unevaluated unless `1 <= r <= n`. The unweighted form slides a window of `r` elements, builds a sublist, and delegates to `Mean` per window — so it inherits Mean's exact-rational / real / symbolic handling. The weighted form computes `wsum = Plus[w_k]`, the coefficients `w_k / wsum`, and for each window emits `Plus[Times[coef_k, x_{i+k}], ...]`, letting the evaluator simplify. All intermediate trees are built and reduced with `eval_and_free`. `ATTR_PROTECTED`.

- `Protected`.
- Output length is `Length[list] - r + 1`.
- Stays unevaluated when `r < 1`, when `r > Length[list]`, when the second argument is non-integer / non-list, or when the first argument is not a `List`.
- Exact rational arithmetic for integer / rational data; bignums (arbitrary-precision integers) handled natively. Real-valued data or weights yield approximate output. Symbolic data and weights are supported.
- The unweighted form delegates to `Mean` for each window, so it inherits `Mean`'s exact / numeric / symbolic dispatch.

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/stats.c`](https://github.com/stblake/mathilda/blob/main/src/stats.c)
- Specification: [`docs/spec/builtins/statistics.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/statistics.md)

## Notes & additional examples

### Worked examples

```mathematica
In[1]:= MovingAverage[{1, 2, 3, 4, 5}, 2]
Out[1]= {3/2, 5/2, 7/2, 9/2}
```

```mathematica
In[1]:= MovingAverage[{1, 2, 3, 4, 5, 6}, 3]
Out[1]= {2, 3, 4, 5}
```

```mathematica
In[1]:= MovingAverage[Table[k^2, {k, 1, 6}], 3]
Out[1]= {14/3, 29/3, 50/3, 77/3}
```

```mathematica
In[1]:= MovingAverage[{a, b, c, d}, {1, 2, 1}]
Out[1]= {1/4 a + 1/2 b + 1/4 c, 1/4 b + 1/2 c + 1/4 d}
```

### Notes

`MovingAverage[list, r]` slides a window of `r` consecutive elements across the
list, returning their averages; the output has length `Length[list] - r + 1`.
Averages are exact rationals, so smoothing the first six squares with a width-3
window gives `{14/3, 29/3, 50/3, 77/3}` rather than decimals. The list-of-weights
form `MovingAverage[list, {w1, ..., wr}]` performs a weighted moving average with
effective weights `wi / Sum[wj]`; with symbolic data and weights `{1, 2, 1}` it
produces the exact binomial smoothing kernel `a/4 + b/2 + c/4`, the discrete
analogue of a triangular filter. The call is left unevaluated when `r < 1` or
`r > Length[list]`.
