# HypergeometricPFQ

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
HypergeometricPFQ[{a1, ...}, {b1, ...}, z]
    is the generalized hypergeometric function pFq(a;b;z), the series
Sum (prod_i Pochhammer[a_i, k] / prod_j Pochhammer[b_j, k]) z^k / k!.
Common upper/lower parameters cancel; a non-positive integer upper
parameter terminates the series to a polynomial. Evaluates to machine,
arbitrary-precision (MPFR), and complex numbers by direct summation in
the convergent regime (p<=q for all z; p==q+1 for |z|<1).
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= HypergeometricPFQ[{a1, a2, a3}, {b1, b2, b3}, 0]
Out[1]= 1

In[2]:= HypergeometricPFQ[{}, {}, z]
Out[2]= E^z

In[3]:= HypergeometricPFQ[{a, b, c}, {a, d, e}, z]
Out[3]= HypergeometricPFQ[{b, c}, {d, e}, z]

In[4]:= HypergeometricPFQ[{1, 1}, {3, 3, 3}, 2.]
Out[4]= 1.07893

In[5]:= HypergeometricPFQ[{1, 2, 3, 4}, {5, 6, 7}, {0.1, 0.3, 0.5}]
Out[5]= {1.01164, 1.03627, 1.06296}

In[6]:= N[HypergeometricPFQ[{1, 1, 1}, {3/2, 3/2, 3/2}, 10], 50]
Out[6]= 530.19188827362590438855961685444087792733053398358

In[7]:= D[HypergeometricPFQ[{a1, a2}, {b1, b2, b3}, x], x]
Out[7]= (a1 a2 HypergeometricPFQ[{1 + a1, 1 + a2}, {1 + b1, 1 + b2, 1 + b3}, x])/(b1 b2 b3)
```

## Implementation notes

- Attributes `NumericFunction`, `Protected`.
- `HypergeometricPFQ[a, b, 0]` is `1`; threads over a `List` third argument.
- Parameters common to the upper and lower lists cancel; a non-positive

**Attributes:** `NumericFunction`, `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/special-functions.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/special-functions.md)

## Notes & additional examples

### Worked examples

```mathematica
In[1]:= HypergeometricPFQ[{1, 1}, {2}, z]
Out[1]= -Log[1 - z]/z
```

Common upper and lower parameters cancel, and any of the specialised forms
(`0F1`, `1F1`, `2F1`) is just a particular shape of `HypergeometricPFQ`. A
higher `3F2` evaluates numerically by direct summation:

```mathematica
In[1]:= N[HypergeometricPFQ[{1, 2, 3}, {4, 5}, 1/2], 40]
Out[1]= 1.1898747542564229318256831180919799547257
```

A non-positive integer upper parameter terminates the series to a polynomial:

```mathematica
In[1]:= HypergeometricPFQ[{-3, 1}, {1}, z]
Out[1]= 1 - 3 z + 3 z^2 - z^3
```

### Notes

`HypergeometricPFQ[{a1, ...}, {b1, ...}, z]` is the generalized hypergeometric function `pFq`, the series `Sum[(Product[Pochhammer[a_i, k]] / Product[Pochhammer[b_j, k]]) z^k / k!]`. It converges for all `z` when `p <= q`, and for `|z| < 1` when `p == q + 1`; a non-positive integer upper parameter truncates it to a polynomial. The specialised heads `Hypergeometric0F1`, `Hypergeometric1F1`, and `Hypergeometric2F1` are convenience wrappers.
