# ExponentialMovingAverage

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
ExponentialMovingAverage[list, alpha]
    gives the exponential moving average of list with smoothing constant alpha.
Defined by the recurrence y_1 = x_1, y_{i+1} = y_i + alpha (x_{i+1} - y_i).
The output has the same length as list. The smoothing constant alpha is typically a number between 0 and 1, but may be any expression; ExponentialMovingAverage handles both numerical (machine and arbitrary precision) and symbolic data.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= ExponentialMovingAverage[Range[10], 1/3]
Out[1]= {1, 4/3, 17/9, 70/27, 275/81, 1036/243, 3773/729, 13378/2187, 46439/6561, 158488/19683}

In[2]:= ExponentialMovingAverage[N[{1, 5, 7, 3, 6, 2}], 1/2]
Out[2]= {1.0, 3.0, 5.0, 4.0, 5.0, 3.5}

In[3]:= ExponentialMovingAverage[{a, b, c, d}, 0]
Out[3]= {a, a, a, a}

In[4]:= ExponentialMovingAverage[{a, b, c, d}, 1]
Out[4]= {a, b, c, d}

In[5]:= ExponentialMovingAverage[{a, b}, x]
Out[5]= {a, a + (-a + b) x}

In[6]:= ExponentialMovingAverage[{2^100, 2^200}, 1/2]
Out[6]= {1267650600228229401496703205376, 803469022129495137770981046171215126561215611592144769253376}
```

## Implementation notes

**Algorithm.** `builtin_exponential_moving_average` takes `(list, alpha)` and applies the recurrence `y[1] = x[1]`, `y[i+1] = y[i] + alpha*(x[i+1] - y[i])`; the output has the same length as the input. It chooses between two paths. The **fast path** (taken when at least one of the list elements or `alpha` is `EXPR_REAL` and all of them are real-valued numerics — no complex, no symbolic, bignums excluded) runs the recurrence in C using `double`s, allocating only the output, returning `EXPR_REAL` elements. The **symbolic / exact path** builds the recurrence out of `Plus`/`Times` nodes per step (via `eval_and_free`), letting the evaluator do exact-rational, bignum, and symbolic arithmetic for arbitrary (including symbolic) `alpha`. Empty or non-`List` first argument leaves the call unevaluated. `ATTR_PROTECTED`.

- `Protected`.
- Output has the same length as `list`.
- Two evaluation strategies: a fast O(n) double-precision path activates when at least one element of `list` or `alpha` itself is a machine-precision real and every other entry is a real-valued numeric (Integer, Real, Rational); otherwise a symbolic recurrence path is taken so exact rationals, bignums (arbitrary-precision integers), and symbolic alpha all work.
- The smoothing constant `alpha` is typically a number between 0 and 1 but may be any expression; with `alpha = 0` the output is constant at `x_1`, and with `alpha = 1` the output equals the input.
- Stays unevaluated when the first argument is not a `List`, when the list is empty, or when the call has the wrong arity.

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/stats.c`](https://github.com/stblake/mathilda/blob/main/src/stats.c)
- Specification: [`docs/spec/builtins/statistics.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/statistics.md)

## Notes & additional examples

### Worked examples

```mathematica
In[1]:= ExponentialMovingAverage[{1, 2, 3, 4, 5}, 1/2]
Out[1]= {1, 3/2, 9/4, 25/8, 65/16}
```

```mathematica
In[1]:= ExponentialMovingAverage[{1, 2, 3, 4, 5}, 0.2]
Out[1]= {1.0, 1.2, 1.56, 2.048, 2.6384}
```

```mathematica
In[1]:= N[ExponentialMovingAverage[{10, 20, 30, 40}, 1/3], 20]
Out[1]= {10.0, 13.3333333333333333334, 18.888888888888888889, 25.9259259259259259259}
```

```mathematica
In[1]:= ExponentialMovingAverage[{a, b, c}, alpha]
Out[1]= {a, a + alpha (-a + b), a + alpha (-a + b) + alpha (-a - alpha (-a + b) + c)}
```

### Notes

`ExponentialMovingAverage[list, alpha]` applies the recurrence
`y_1 = x_1`, `y_{i+1} = y_i + alpha (x_{i+1} - y_i)`. Exact rationals stay
exact, inexact data evaluates at the requested precision, and a symbolic
smoothing constant `alpha` produces the unrolled closed form term by term.
