# CentralMoment

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
CentralMoment[data, r]
    gives the r-th central moment (moment about the mean) of data, (1/n) Sum[(x_i - Mean[data])^r].
CentralMoment[data, {r_1, ..., r_m}]
    gives the multivariate central moment of data. For a matrix or array the moment is taken columnwise over the first axis.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= CentralMoment[{1, 2, 3, 4}, 4]
Out[1]= 41/16

In[2]:= CentralMoment[{1., 2., 3., 4.}, 2]
Out[2]= 1.25

In[3]:= CentralMoment[{{1, 2}, {3, 4}, {5, 6}}, 2]
Out[3]= {8/3, 8/3}

In[4]:= Simplify[CentralMoment[{{a, b}, {c, d}}, {2, 2}]]
Out[4]= 1/16 (a - c)^2 (b - d)^2
```

## Implementation notes

- `Protected`.
- A central moment is `Variance` without the $n/(n-1)$ bias correction: it divides by `n` (not `n-1`), raises to the power `r` (not a square), and needs only `n >= 1`.
- For a matrix or array the moment is taken columnwise over the first axis (equivalent to `ArrayReduce[CentralMoment[#, r]&, x, 1]`).
- Exact input yields exact output; approximate input yields approximate output; symbolic data is handled symbolically.
- Fast path on `NDArray`/packed real buffers (`ndred_central_moment`); an integer buffer degrades to the exact `Rational` `List` result, like `Variance`.
- Lowerable inside `Compile[]` for a real vector and integer order (participates in auto-compilation).

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/statistics.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/statistics.md)
