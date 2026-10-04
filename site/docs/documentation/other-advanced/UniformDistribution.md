# UniformDistribution

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`UniformDistribution[{lo, hi}] represents a continuous uniform distribution; UniformDistribution[] is uniform on {0, 1}.`**

## Examples (3)

Every input below was run against the current Mathilda build and its output recorded.

### Applications (3)

An inert distribution object

```mathematica
In[1]:= Head[UniformDistribution[{0, 1}]]
Out[1]= UniformDistribution
```

Density 1/(hi - lo) inside the support

```mathematica
In[2]:= PDF[UniformDistribution[{0, 1}], 1/2]
Out[2]= 1.0
```

Zero strictly outside [lo, hi]

```mathematica
In[3]:= PDF[UniformDistribution[{2, 6}], 10]
Out[3]= 0.0
```

## Implementation notes

**Definition.** `UniformDistribution[{lo, hi}]` is a **distribution object** representing the
continuous uniform distribution on `[lo, hi]`; `UniformDistribution[]` is uniform on
`{0, 1}`. It is an inert, `Protected` head, not a function: `UniformDistribution[{0, 1}]`
evaluates to itself, and its meaning is supplied by the consumers in `src/ml/dist.c` (`PDF`,
`RandomVariate`), which parse it into an internal `MlDist` record.

**Representation.** The object stays an `EXPR_FUNCTION` with head `UniformDistribution`, so
`Head[UniformDistribution[{0, 1}]]` is `UniformDistribution`. When a consumer reads it, the
parser (`ml_parse_distribution` in `dist.c`) sets the internal kind to `ML_D_UNIFORM` with
endpoints `a = lo`, `b = hi`; the no-argument form defaults to `a = 0, b = 1`. A spec whose
upper endpoint does not exceed the lower (`b > a` fails) is rejected, so the object never
yields a degenerate density.

**Usage & limits.** `PDF[UniformDistribution[{lo, hi}], x]` returns `1/(hi - lo)` for
`lo <= x <= hi` (endpoints included, Wolfram's convention for a continuous uniform) and `0`
strictly outside; `RandomVariate[UniformDistribution[{lo, hi}]]` draws a sample. These
numeric consumers evaluate at machine precision, so `PDF[...]` comes back as a machine real.
Symbolic moment functions such as `Mean`, `Variance` and `CDF` are **not** yet wired to this
object and are left unevaluated. It sits beside `NormalDistribution` in the same module.

**Attributes:** `Protected`.

## References

- Source: [`src/ml/dist.c`](https://github.com/stblake/mathilda/blob/main/src/ml/dist.c)
- Specification index: [`Mathilda_spec.md`](https://github.com/stblake/mathilda/blob/main/Mathilda_spec.md)
- Tests: [`tests/test_ml_dist.c`](https://github.com/stblake/mathilda/blob/main/tests/test_ml_dist.c)

## Notes & additional examples

### Notes

`UniformDistribution[{lo, hi}]` represents the continuous uniform distribution on
`[lo, hi]`, and `UniformDistribution[]` the standard uniform on `{0, 1}`. It is an inert,
`Protected` object — it evaluates to itself, `Head` is `UniformDistribution` — whose meaning
is supplied by the consumers `PDF` and `RandomVariate`.

`PDF[UniformDistribution[{lo, hi}], x]` is `1/(hi - lo)` for `lo <= x <= hi` (endpoints
included) and `0` outside; these numeric consumers work at machine precision, so the density
comes back as a machine real. Symbolic moment functions such as `Mean`, `Variance` and `CDF`
are not yet wired to this object and are left unevaluated.
