# StieltjesGamma

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
StieltjesGamma[n]
    gives the n-th Stieltjes constant gamma_n, the Laurent coefficients of
Zeta about s = 1. StieltjesGamma[0] is EulerGamma; higher constants are
inert (they stay symbolic) and appear in Series expansions of Zeta. Listable.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= StieltjesGamma[0]
Out[1]= EulerGamma
```

## Implementation notes

**Attributes:** `Listable`, `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/special-functions.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/special-functions.md)
