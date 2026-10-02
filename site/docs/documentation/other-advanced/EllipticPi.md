# EllipticPi

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`EllipticPi[n, m] is the complete elliptic integral of the third kind, Integrate[1/((1 - n Sin[t]^2) Sqrt[1 - m Sin[t]^2]), {t, 0, Pi/2}], and EllipticPi[n, phi, m] the incomplete one, with upper limit phi. The parameter argument is m = k^2, not the modulus k. EllipticPi[0, m] is EllipticK[m] and EllipticPi[0, phi, m] is EllipticF[phi, m]. For n > 1 the path crosses the pole at Sin[t]^2 == 1/n and the value is the Cauchy principal value.`**

## Examples

_No verified examples yet for this function._

## Implementation notes

**Attributes:** `Listable`, `NumericFunction`, `Protected`.

## References

- Source: [`src/special_functions/elliptic.c`](https://github.com/stblake/mathilda/blob/main/src/special_functions/elliptic.c)
- Specification index: [`Mathilda_spec.md`](https://github.com/stblake/mathilda/blob/main/Mathilda_spec.md)
- Tests: [`tests/test_elliptic.c`](https://github.com/stblake/mathilda/blob/main/tests/test_elliptic.c)
