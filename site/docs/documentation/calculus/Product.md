# Product

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
Product[f, {i, imax}]
    gives the product of f for i from 1 to imax.

Product[f, {i, imin, imax}], Product[f, {i, imin, imax, di}] and Product[f, {i, {i1, i2, ...}}] use the standard iterator forms; multiple iterators give nested products (an inner bound may depend on an outer index). Product[f, i] gives the indefinite product (anti-quotient). The index is localised (HoldAll). Finite ranges are multiplied out directly; symbolic, indefinite and convergent infinite products are evaluated in exact closed form (n!, Pochhammer, Gamma ratios, base^k, QPochhammer, BarnesG) via a Method polyalgorithm.

Options: Method (Automatic | "Telescoping" | "Rational" | "Geometric" | "QProduct"), VerifyConvergence (default True; a divergent infinite product gives Product::div), GenerateConditions, Assumptions. N[Product[...]] routes to NProduct.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= Product[k, {k, 1, n}]
Out[1]= Factorial[n]

In[2]:= Product[k + a, {k, 1, n}]
Out[2]= Pochhammer[1 + a, n]

In[3]:= Product[2^k, {k, 1, n}]
Out[3]= 2^(1/2 n (1 + n))

In[4]:= Product[1 + 1/k^2, {k, 1, Infinity}]
Out[4]= Sinh[Pi]/Pi

In[5]:= Product[1 - a q^k, {k, 0, n - 1}]
Out[5]= QPochhammer[a, q, n]

In[6]:= Product[k^k, {k, 1, n}]
Out[6]= Hyperfactorial[n]
```

## Implementation notes

**Attributes:** `HoldAll`, `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/calculus.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/calculus.md)
