# HurwitzZeta

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
HurwitzZeta[s, a]
    is the Hurwitz zeta function zeta(s, a) = Sum_{k>=0} (k + a)^-s.
Identical to Zeta[s, a] for Re(a) > 0, but built on the principal-branch
power (k + a)^-s, so it differs from Zeta for non-positive real a and has
poles at a = 0, -1, -2, ... . HurwitzZeta[s, 1] is Zeta[s], HurwitzZeta[s, 1/2]
is (2^s - 1) Zeta[s], and a positive integer a reduces to Zeta[s] minus a
finite power sum. A non-positive integer a gives ComplexInfinity for positive
integer s and the Bernoulli-polynomial value for non-positive integer s.
Real, complex, machine and arbitrary-precision (MPFR) arguments evaluate
numerically via an Euler-Maclaurin kernel. Listable.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= HurwitzZeta[s, 1/2]
Out[1]= (-1 + 2^s) Zeta[s]

In[2]:= HurwitzZeta[3, -3.5]
Out[2]= 0.0307784
```

## Implementation notes

**Attributes:** `Listable`, `NumericFunction`, `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/special-functions.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/special-functions.md)
