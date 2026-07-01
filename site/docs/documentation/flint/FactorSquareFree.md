# FLINT`FactorSquareFree

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
FLINT`FactorSquareFree[p] gives the squarefree factorisation of the polynomial p over the rationals, computed directly via FLINT (fmpq_mpoly_factor_squarefree). Returns unevaluated if out of scope.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= FLINT`FactorSquareFree[(x - 1)^2 (x + 1)]
Out[1]= (1 + x) (-1 + x)^2
```

## Implementation notes

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/poly/flint_bridge.c`](https://github.com/stblake/mathilda/blob/main/src/poly/flint_bridge.c)
- Specification: [`docs/spec/builtins/flint.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/flint.md)
