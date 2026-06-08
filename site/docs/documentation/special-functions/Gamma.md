# Gamma

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
Gamma[z]
    is the Euler gamma function Gamma(z).
Gamma[a, z]
    is the upper incomplete gamma function Gamma(a, z).
Gamma[a, z0, z1]
    is the generalized incomplete gamma Gamma(a, z0) - Gamma(a, z1).
Integer and half-integer arguments reduce to exact values ((z-1)!,
and rational multiples of Sqrt[Pi]); non-positive integers give
ComplexInfinity. Machine and arbitrary-precision (MPFR) real inputs
evaluate numerically, as do machine-precision complex inputs. Listable.
```

## Examples

_No verified examples yet for this function._

## Implementation notes

- Exact reductions for `Gamma[z]`:
  - Positive integers: `Gamma[n] = (n-1)!` (exact, with GMP BigInt for large `n`).
  - Non-positive integers are poles: `Gamma[0]`, `Gamma[-n]` → `ComplexInfinity`.
  - Half-integers reduce to rational multiples of `Sqrt[Pi]`, e.g.

**Attributes:** `Listable`, `NumericFunction`, `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/special-functions.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/special-functions.md)
