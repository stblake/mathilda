# LerchPhi

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`LerchPhi[z, s, a]`**

is the Lerch transcendent Phi(z, s, a) = Sum\_{k\>=0} z^k/(k + a)^s.

**`Zeta[s, a] and z LerchPhi[z, s, 1] is PolyLog[s, z]. Exact reductions`**

<details>
<summary>Notes</summary>

It generalizes Zeta, HurwitzZeta and PolyLog: LerchPhi\[1, s, a\] is cover z = 0 (a^-s), s = 0 (1/(1-z)), z = +-1, positive integer a (a PolyLog form) and negative integer s (a rational function of z). The options DoublyInfinite -\> True (sum k from -Infinity to Infinity) and IncludeSingularTerm -\> True (keep the k + a = 0 term) are supported. Inexact arguments with |z| \< 1 evaluate numerically at machine or arbitrary (MPFR) precision; |z| \> 1 stays symbolic. Listable.

</details>

## Examples (2)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (2)

```mathematica
In[1]:= LerchPhi[z, s, 1]
Out[1]= PolyLog[s, z]/z

In[2]:= LerchPhi[0.5, 3, 2.5]
Out[2]= 0.0794983
```

## Algorithm

Mathilda -- the Lerch transcendent LerchPhi.

```text
  LerchPhi[z, s, a]   Phi(z, s, a) = Sum_{k>=0} z^k / (k + a)^s
                      (|z| < 1; analytic continuation elsewhere, branch cut
                       z in [1, Infinity))
```

For Re(a) < 0 the principal value uses the symmetric power ((k+a)^2)^(-s/2),

```text
and any term with k + a = 0 is excluded.  LerchPhi is the common
```

generalization of Zeta, HurwitzZeta and PolyLog:

```text
  Phi(1, s, a) = Zeta[s, a],      z Phi(z, s, 1) = PolyLog[s, z].
```

The evaluator routes each kind of argument to the cheapest exact or fastest numeric path (mirroring src/special_functions/{zeta,hurwitzzeta,polylog}.c):

```text
  exact reductions (any z, s, a):
      z = 0                  ->  a^-s                        (the k = 0 term)
      s = 0                  ->  1/(1 - z)                   (geometric sum)
      z = 1                  ->  Zeta[s, a]
      z = -1                 ->  2^-s (Zeta[s,a/2] - Zeta[s,(a+1)/2])
      a positive integer m   ->  z^-m (PolyLog[s,z] - Sum_{j<m} z^j j^-s)
      s negative integer -n  ->  (z d/dz + a)^n [1/(1-z)]    (rational in z)
  options:
      IncludeSingularTerm->True at a non-positive integer a  ->  ComplexInfinity
      DoublyInfinite->True   ->  Phi(z,s,a) + z^-1 Phi(1/z, s, 1-a)
  numeric (>= 1 inexact operand):
      |z| < 1                ->  complex-MPFR power series
      |z| > 1                ->  diverges, stays symbolic (no continuation)
  everything else            ->  stays symbolic (return NULL)
```

Attributes: Listable, NumericFunction, Protected.

## Implementation notes

**Attributes:** `Listable`, `NumericFunction`, `Protected`.

## See also

[Zeta](../../special-functions/Zeta/), [HurwitzZeta](../../special-functions/HurwitzZeta/), [PolyLog](../../special-functions/PolyLog/)

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/special-functions.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/special-functions.md)
- Tests: [`tests/test_compile.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compile.c)
- Tests: [`tests/test_lerchphi.c`](https://github.com/stblake/mathilda/blob/main/tests/test_lerchphi.c)
- Tests: [`tests/test_numeric_stress.c`](https://github.com/stblake/mathilda/blob/main/tests/test_numeric_stress.c)
