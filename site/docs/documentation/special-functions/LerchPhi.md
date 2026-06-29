# LerchPhi

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
LerchPhi[z, s, a]
    is the Lerch transcendent Phi(z, s, a) = Sum_{k>=0} z^k/(k + a)^s.
It generalizes Zeta, HurwitzZeta and PolyLog: LerchPhi[1, s, a] is
Zeta[s, a] and z LerchPhi[z, s, 1] is PolyLog[s, z]. Exact reductions
cover z = 0 (a^-s), s = 0 (1/(1-z)), z = +-1, positive integer a (a
PolyLog form) and negative integer s (a rational function of z). The
options DoublyInfinite -> True (sum k from -Infinity to Infinity) and
IncludeSingularTerm -> True (keep the k + a = 0 term) are supported.
Inexact arguments with |z| < 1 evaluate numerically at machine or
arbitrary (MPFR) precision; |z| > 1 stays symbolic. Listable.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= LerchPhi[z, s, 1]
Out[1]= PolyLog[s, z]/z

In[2]:= LerchPhi[0.5, 3, 2.5]
Out[2]= 0.0794983
```

## Implementation notes

**Attributes:** `Listable`, `NumericFunction`, `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/special-functions.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/special-functions.md)
