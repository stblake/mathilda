# Moment

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
Moment[data, r]
    gives the r-th raw (power) moment of data, (1/n) Sum[x_i^r].
Moment[data, {r_1, ..., r_m}]
    gives the multivariate raw moment of data. For a matrix or array the moment is taken columnwise over the first axis.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= Moment[{1, 2, 3, 4}, 2]
Out[1]= 15/2

In[2]:= Moment[{1., 2., 3., 4.}, 2]
Out[2]= 7.5

In[3]:= Moment[{Pi, E, 2}, 1]
Out[3]= 1/3 (2 + E + Pi)

In[4]:= Moment[{{1, 2}, {3, 4}, {5, 6}}, 3]
Out[4]= {51, 96}

In[5]:= Simplify[Moment[{{a, b}, {c, d}}, {1, 2}]]
Out[5]= 1/2 (a b^2 + c d^2)
```

## Implementation notes

- `NHoldAll`, `Protected`.
- The raw moment is `CentralMoment` without the mean subtraction; `Moment[data, 1]` is `Mean[data]`, and `Moment[data, 0]` is `1`.
- For a matrix or array the moment is taken columnwise over the first axis (equivalent to `ArrayReduce[Moment[#, r]&, x, 1]`); because there is no mean to subtract, `Mean[data^r]` threads correctly at every rank.
- Exact input yields exact output; approximate input yields approximate output; symbolic data is handled symbolically.
- Fast path on `NDArray`/packed real buffers (`ndred_moment`); an integer buffer degrades to the exact `Rational` `List` result, like `Variance`.
- Lowerable inside `Compile[]` for a real vector and integer order (participates in auto-compilation).

**Attributes:** `NHoldAll`, `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/statistics.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/statistics.md)
