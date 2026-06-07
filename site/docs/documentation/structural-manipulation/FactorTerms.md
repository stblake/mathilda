# FactorTerms

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
FactorTerms[poly]
    pulls out any overall numerical factor in poly.
FactorTerms[poly, x]
    pulls out any overall factor in poly that does not depend on x.
FactorTerms[poly, {x1, x2, ...}]
    pulls out any overall factor in poly that does not depend on any of the xi, then progressively factors with respect to smaller subsets {x1, ..., x_{k-1}}.
FactorTerms[poly, x] extracts the content of poly with respect to x.
FactorTerms automatically threads over lists, equations, inequalities and logic functions.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= FactorTerms[3 + 6x + 3x^2]
Out[1]= 3 (1 + 2 x + x^2)

In[2]:= FactorTerms[3 + 3a + 6 a x + 6 x + 12 a x^2 + 12 x^2, x]
Out[2]= 3 (1 + a) (1 + 2 x + 4 x^2)

In[3]:= FactorTerms[12 a^4 + 9 x^2 + 66 b^2]
Out[3]= 3 (4 a^4 + 22 b^2 + 3 x^2)

In[4]:= FactorTerms[7 x + (14 y + 21)/z]
Out[4]= (7 (3 + 2 y + x z))/z

In[5]:= FactorTerms[{5 x^2 - 15, 7 x^4 - 77, 8 x^8 - 24}]
Out[5]= {5 (-3 + x^2), 7 (-11 + x^4), 8 (-3 + x^8)}

In[6]:= FactorTerms[1 < 77 x^3 - 21 x + 35 < 2]
Out[6]= 1 < 7 (5 - 3 x + 11 x^3) < 2

In[7]:= f = 2 x^2 y z + 2 x^2 y + 4 x^2 z + 4 x^2 + 4 y^2 z^2 + 4 z y^2
Out[7]= 4 x^2 + 2 x^2 y + 4 x^2 z + 2 x^2 y z + 4 y^2 z + 4 y^2 z^2

In[8]:= FactorTerms[f, {x, y}]
Out[8]= 2 (1 + z) (2 x^2 + x^2 y + 2 y^2 z)
```

## Implementation notes

- `Protected`.
- Auto-threads over `List`, `Equal`, `Unequal`, `Less`, `LessEqual`,

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/structural-manipulation.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/structural-manipulation.md)
