# RootMeanSquare

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`RootMeanSquare[list] gives the root mean square of values in list.`**

## Examples (10)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (5)

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

### Applications (5)

```mathematica
In[6]:= RootMeanSquare[{3, 4}]
Out[6]= 5/Sqrt[2]

In[7]:= RootMeanSquare[{a, b}]
Out[7]= Sqrt[1/2 (a^2 + b^2)]

In[8]:= RootMeanSquare[Range[10]]
Out[8]= Sqrt[77/2]

In[9]:= N[RootMeanSquare[{1, 2, 3, 4, 5}], 30]
Out[9]= 3.316624790355399849114932736672

In[10]:= N[RootMeanSquare[Table[Sin[n], {n, 1, 1000}]], 20]
Out[10]= 0.707242937053949660224
```

## Implementation notes

**Algorithm.** `builtin_rootmeansquare` computes `Sqrt[Mean[x^2]]`. It reduces matrices column-wise via `apply_columnwise` and requires a `List`. For data containing a real, it sums squares in `double` and returns `expr_new_real(sqrt(sum_sq/n))`. For exact/symbolic data it builds `Plus[x_i^2 ...]`, then carefully distributes the square root to keep results exact: if the summed result is non-numeric and `n` is a perfect square it factors out `1/Sqrt[n]`; if the mean-square is a rational with a perfect-square denominator it pulls that root out before applying `Power[..., 1/2]`. Otherwise it returns `Power[meanSq, 1/2]` for the evaluator to simplify. `ATTR_PROTECTED`.

- `Protected`.
- Gives the square root of the second sample moment.
- For a list `{x1, x2, ...}`, it computes `Sqrt[1/n Total[{x1^2, x2^2, ...}]]`.
- Handles both numerical and symbolic data.
- Works column-wise on matrices.

**Attributes:** `Protected`.

## References

- Source: [`src/stats.c`](https://github.com/stblake/mathilda/blob/main/src/stats.c)
- Specification: [`docs/spec/builtins/statistics.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/statistics.md)
- Tests: [`tests/test_compiledfunction.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compiledfunction.c)
- Tests: [`tests/test_ndarray_reduce.c`](https://github.com/stblake/mathilda/blob/main/tests/test_ndarray_reduce.c)
- Tests: [`tests/test_stats.c`](https://github.com/stblake/mathilda/blob/main/tests/test_stats.c)

## Notes & additional examples

### Notes

`RootMeanSquare[list]` returns `Sqrt[Mean[list^2]]` and stays **exact** on exact
input — a list of two symbols yields the closed form `Sqrt[(a^2 + b^2)/2]`. The
final example is a numerical demonstration of an analytic fact: the RMS of
`Sin` sampled over many integer arguments approaches `1/Sqrt[2] ≈ 0.70711`,
the continuous RMS of a sinusoid.
