# FLINT`Resultant

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

FLINT\`Resultant\[a, b, x\] gives the resultant of the polynomials a and b eliminating the variable x, over the rationals, computed directly via FLINT (fmpq\_mpoly\_resultant). Other variables are treated as coefficients. Returns unevaluated if out of scope.

## Examples (1)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (1)

```mathematica
In[1]:= FLINT`Resultant[x^2 - 1, x - 2, x]
Out[1]= 3
```

## Implementation notes

**Attributes:** `Protected`.

## References

- Source: [`src/poly/flint_bridge.c`](https://github.com/stblake/mathilda/blob/main/src/poly/flint_bridge.c)
- Specification: [`docs/spec/builtins/flint.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/flint.md)
