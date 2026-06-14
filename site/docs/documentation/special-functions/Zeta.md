# Zeta

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
Zeta[s]
    is the Riemann zeta function zeta(s) = Sum_{k>=1} k^-s.
Zeta[s, a]
    is the Hurwitz zeta function zeta(s, a) = Sum_{k>=0} (k+a)^-s.
Even positive integers give rational multiples of Pi^(2n), negative
integers give rationals, Zeta[0] is -1/2, and Zeta[1] is ComplexInfinity;
odd positive integers stay symbolic. Hurwitz zeta at a positive integer a
reduces to Zeta[s] minus a finite power sum. Real, complex, machine and
arbitrary-precision (MPFR) numeric arguments evaluate numerically via
mpfr_zeta (real Riemann) or an Euler-Maclaurin kernel. Listable.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= Zeta[2]
Out[1]= 1/6 Pi^2

In[2]:= Series[Zeta[x], {x, 1, 2}] // Normal
Out[2]= EulerGamma + 1/(-1 + x) - StieltjesGamma[1] (-1 + x) + 1/2 StieltjesGamma[2] (-1 + x)^2
```

## Implementation notes

**Attributes:** `Listable`, `NumericFunction`, `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/special-functions.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/special-functions.md)
