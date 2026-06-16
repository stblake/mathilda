# Quartiles

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
Quartiles[data]
    gives the {q_1/4, q_2/4, q_3/4} quantile estimates of the elements in data.
Quartiles[data,{{a,b},{c,d}}]
    uses the quantile definition specified by parameters a, b, c, d.
Quartiles[dist]
    gives the {q_1/4, q_2/4, q_3/4} quantiles of the distribution dist.
```

## Examples

_No verified examples yet for this function._

## Implementation notes

**Algorithm.** `builtin_quartiles` returns `{Q1, Q2, Q3}` using a parameterised quantile estimator. For a matrix (list of equal-length lists) it `Transpose`s and recurses column-wise. For a flat list it requires all entries to be real (`is_real_numeric`, else emits `Quartiles::rectn`), sorts via the evaluator's `Sort`, and for each `q` in `{1/4, 1/2, 3/4}` computes the interpolated order statistic with the four quantile parameters `{{a, b}, {c, d}}` — defaulting to `{{1/2, 0}, {0, 1}}` (Mathematica's default). The index is `h = a + (n + b)·q`; it clamps to the endpoints when `h ≤ 1` or `h ≥ n`, takes `j = Floor[h]` and the fractional `g = h − j`, and linearly interpolates `sorted[j] + g·(c + d·g type adjustment)` between adjacent sorted elements. All arithmetic is done by building and `eval_and_free`-ing `Plus`/`Times`/`Floor` expression nodes, so exact (rational/symbolic-numeric) inputs stay exact. Non-list input returns the call unevaluated. Attribute: `ATTR_PROTECTED`.

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/stats.c`](https://github.com/stblake/mathilda/blob/main/src/stats.c)
- Specification: [`docs/spec/builtins/statistics.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/statistics.md)

## Notes & additional examples

### Worked examples

```mathematica
In[1]:= Quartiles[{1, 2, 3, 4, 5, 6, 7, 8}]
Out[1]= {5/2, 9/2, 13/2}

In[2]:= Quartiles[{6, 7, 15, 36, 39, 40, 41, 42, 43, 47, 49}]
Out[2]= {81/4, 40, 171/4}
```

The estimates stay perfectly exact on rational data, so the interquartile range
(`q3 - q1`) of an arithmetic progression comes out in closed form:

```mathematica
In[1]:= q = Quartiles[Range[1, 1000]]
Out[1]= {501/2, 1001/2, 1501/2}

In[2]:= q[[3]] - q[[1]]
Out[2]= 500
```

An alternative quantile convention can be selected with a parameter list:

```mathematica
In[1]:= Quartiles[{1, 2, 3, 4, 5, 6, 7, 8}, {{1/2, 0}, {0, 1}}]
Out[1]= {5/2, 9/2, 13/2}
```

### Notes

`Quartiles[data]` gives the three quartile estimates `{q1, q2, q3}` (lower,
median, upper) of the data. The default uses Mathematica's standard quantile
definition; an optional parameter list `{{a, b}, {c, d}}` selects an
alternative quantile convention.
