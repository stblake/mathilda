# StandardDeviation

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`StandardDeviation[data] gives the standard deviation estimate of the elements in data.`**

## Examples (8)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (3)

```mathematica
In[1]:= Median[<|"a" -> 1, "b" -> 3, "c" -> 5|>]
Out[1]= 3

In[2]:= Variance[<|"a" -> 2, "b" -> 4, "c" -> 6|>]
Out[2]= 4

In[3]:= StandardDeviation[<|"a" -> 2, "b" -> 4, "c" -> 6|>]
Out[3]= 2
```

### Applications (5)

```mathematica
In[1]:= StandardDeviation[{1, 2, 3, 4, 5}]
Out[1]= Sqrt[5/2]
```

Results are kept exact; here the (unbiased, n − 1) estimate of a small sample:

```mathematica
In[1]:= StandardDeviation[{2, 4, 4, 4, 5, 5, 7, 9}]
Out[1]= 4 Sqrt[2/7]
```

Numericalize the same estimate to 40 digits:

```mathematica
In[1]:= N[StandardDeviation[{2, 4, 4, 4, 5, 5, 7, 9}], 40]
Out[1]= 2.1380899352993950774764278470380281724321
```

It relates to `Variance` as its square root:

```mathematica
In[1]:= Variance[{1, 2, 3, 4, 5}]
Out[1]= 5/2
```

A constant sample has zero spread:

```mathematica
In[1]:= StandardDeviation[{1, 1, 1, 1}]
Out[1]= 0
```

## Implementation notes

`builtin_standard_deviation` is essentially `Sqrt[Variance[data]]`. It reduces matrices column-wise via `apply_columnwise`. For an all-real numeric vector (`n > 1`) it evaluates `Variance[data]` and, if that returns an `EXPR_REAL`, returns `expr_new_real(sqrt(...))` directly. Otherwise it evaluates `Variance[data]` and raises it to the `1/2` power via a `Power[var, Rational[1,2]]` node, letting the evaluator produce an exact or symbolic radical. `ATTR_PROTECTED`. Inherits Variance's sample (`n-1`) convention.

**Attributes:** `Protected`.

## See also

[Median](../../data-structures/Median/), [Variance](../../data-structures/Variance/), [Mean](../../data-structures/Mean/)

## References

- Source: [`src/stats.c`](https://github.com/stblake/mathilda/blob/main/src/stats.c)
- Specification: [`docs/spec/builtins/data-structures.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/data-structures.md)
- Tests: [`tests/test_association.c`](https://github.com/stblake/mathilda/blob/main/tests/test_association.c)
- Tests: [`tests/test_compiledfunction.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compiledfunction.c)
- Tests: [`tests/test_ndarray_reduce.c`](https://github.com/stblake/mathilda/blob/main/tests/test_ndarray_reduce.c)
- Tests: [`tests/test_stats.c`](https://github.com/stblake/mathilda/blob/main/tests/test_stats.c)

## Notes & additional examples

### Notes

`StandardDeviation[data]` returns the sample (unbiased, divide-by-`n - 1`)
standard deviation, i.e. `Sqrt[Variance[data]]`. Exact inputs give exact
radical output, which `N[..., d]` evaluates to arbitrary precision. A list of
length 1 or a constant list yields `0`.
