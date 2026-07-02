# FLINT`Factor

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
FLINT`Factor[p] gives the irreducible factorisation of the polynomial p over the rationals, computed directly via FLINT (fmpq_mpoly_factor), as Times[const, factor^exp, ...]. Multivariate. Returns unevaluated if p is not a polynomial over Q.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= FLINT`Factor[x^4 - 1]
Out[1]= (-1 + x) (1 + x) (1 + x^2)
```

## Implementation notes

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/poly/flint_bridge.c`](https://github.com/stblake/mathilda/blob/main/src/poly/flint_bridge.c)
- Specification: [`docs/spec/builtins/flint.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/flint.md)
