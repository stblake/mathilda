# FLINT`PolynomialGCD

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

FLINT\`PolynomialGCD\[a, b\] gives the monic greatest common divisor of the polynomials a and b over the rationals, computed directly via FLINT (fmpq\_mpoly\_gcd). Multivariate. Returns unevaluated if an argument is not a polynomial over Q.

## Examples (2)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (2)

```mathematica
In[1]:= FLINT`PolynomialGCD[x^2 - 1, x^2 - x]
Out[1]= -1 + x

In[2]:= FLINT`PolynomialGCD[x^2 - y^2, x^2 + 2 x y + y^2]
Out[2]= x + y
```

## Implementation notes

**Attributes:** `Protected`.

## References

- Source: [`src/poly/flint_bridge.c`](https://github.com/stblake/mathilda/blob/main/src/poly/flint_bridge.c)
- Specification: [`docs/spec/builtins/flint.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/flint.md)
