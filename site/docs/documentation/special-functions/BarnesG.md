# BarnesG

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
BarnesG[z]
    gives the Barnes G-function.
G(z+1) = Gamma[z] G(z) with G(1)=G(2)=1; for a positive integer n, G(n+1) = prod_{k=1}^{n-1} k! (exact via GMP), and G(m)=0 for non-positive integer m. Non-integer orders stay symbolic. Listable, NumericFunction.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= BarnesG[5]
Out[1]= 12
```

## Implementation notes

**Attributes:** `Listable`, `NumericFunction`, `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/special-functions.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/special-functions.md)
