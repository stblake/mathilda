# FLINT`FactorSquareFree

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

FLINT\`FactorSquareFree\[p\] gives the squarefree factorisation of the polynomial p over the rationals, computed directly via FLINT (fmpq\_mpoly\_factor\_squarefree). Returns unevaluated if out of scope.

## Examples (1)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (1)

```mathematica
In[1]:= FLINT`FactorSquareFree[(x - 1)^2 (x + 1)]
Out[1]= (1 + x) (-1 + x)^2
```

## Implementation notes

**Attributes:** `Protected`.

## References

- Source: [`src/poly/flint_bridge.c`](https://github.com/stblake/mathilda/blob/main/src/poly/flint_bridge.c)
- Specification: [`docs/spec/builtins/flint.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/flint.md)
