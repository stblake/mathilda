# QPochhammer

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
QPochhammer[a, q, n]
    gives the q-Pochhammer symbol prod_{k=0}^{n-1} (1 - a q^k).
QPochhammer[a, q] gives the infinite q-Pochhammer (a;q)_Inf for |q|<1. The finite form is exact/symbolic for a non-negative integer n; the infinite form evaluates for machine-real a, q. Listable, NumericFunction.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= QPochhammer[a, q, 3]
Out[1]= (1 - a) (1 - a q) (1 - a q^2)
```

## Implementation notes

**Attributes:** `Listable`, `NumericFunction`, `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/special-functions.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/special-functions.md)
