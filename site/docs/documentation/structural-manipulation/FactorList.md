# FactorList

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
FactorList[poly] gives a list of the irreducible factors of poly together
with their exponents, as {factor, exponent} pairs.  A thin wrapper over
Factor: options (GaussianIntegers -> True, Extension -> {a1, ...}) are
forwarded verbatim.  The first element is always the overall numerical
factor {c, 1} (it is {1, 1} when there is none); denominator factors of a
rational function appear with negative exponents.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= FactorList[x^2 - 1]
Out[1]= {{1, 1}, {-1 + x, 1}, {1 + x, 1}}

In[2]:= FactorList[2 x^3 + 2 x^2 - 2 x - 2]
Out[2]= {{2, 1}, {-1 + x, 1}, {1 + x, 2}}

In[3]:= FactorList[(x^3 + 2 x^2)/(x^2 - 4 y^2) - (x + 2)/(x^2 - 4 y^2)]
Out[3]= {{1, 1}, {-1 + x, 1}, {1 + x, 1}, {2 + x, 1}, {x - 2 y, -1}, {x + 2 y, -1}}

In[4]:= FactorList[x^2 + 1, Extension -> I]
Out[4]= {{1, 1}, {-I + x, 1}, {I + x, 1}}
```

## Implementation notes

- `Listable`, `Protected`.
- A thin wrapper over `Factor`: it factors via `Factor[poly, opts...]` (options are forwarded verbatim) and splits the product into `{factor, exponent}` pairs.
- The first element is always the overall numerical factor `{c, 1}` — it is `{1, 1}` when there is no numerical factor.
- Denominator factors of a rational function appear with negative exponents.
- `Times @@ Power @@@ FactorList[poly]` reconstructs `poly` (up to `Factor`'s normal form).
- Because it delegates to `Factor`, the factorisation depth and normal form (factor ordering, `Extension`/`GaussianIntegers` handling) are exactly `Factor`'s.

**Attributes:** `Listable`, `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/structural-manipulation.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/structural-manipulation.md)
