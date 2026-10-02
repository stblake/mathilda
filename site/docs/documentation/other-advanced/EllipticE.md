# EllipticE

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`EllipticE[m] is the complete elliptic integral of the second kind, Integrate[Sqrt[1 - m Sin[t]^2], {t, 0, Pi/2}], and EllipticE[phi, m] the incomplete one, with upper limit phi. The parameter argument is m = k^2, not the modulus k. EllipticE[0] is Pi/2, EllipticE[1] is 1, and EllipticE[phi, 1] is Sin[phi].`**

## Examples

_No verified examples yet for this function._

## Implementation notes

**Attributes:** `Listable`, `NumericFunction`, `Protected`.

## References

- Source: [`src/special_functions/elliptic.c`](https://github.com/stblake/mathilda/blob/main/src/special_functions/elliptic.c)
- Specification index: [`Mathilda_spec.md`](https://github.com/stblake/mathilda/blob/main/Mathilda_spec.md)
- Tests: [`tests/test_elliptic.c`](https://github.com/stblake/mathilda/blob/main/tests/test_elliptic.c)
