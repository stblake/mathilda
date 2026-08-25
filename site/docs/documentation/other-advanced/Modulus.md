# Modulus

!!! note "Status: Experimental"
    present and registered, but lightly documented and not yet covered by dedicated tests.

## Description

Modulus is an option for Solve.  Solve\[poly == 0, x, Modulus -\> p\] solves a single-variable polynomial equation over the finite ring Z/pZ by residue enumeration, returning {{x -\> r}, ...} with r ascending in \[0, p).  Supported for 2 \<= p \<= 100000; systems, multivariable specs, non-polynomial equations, or an out-of-range modulus leave Solve unevaluated.

## Examples

_No verified examples yet for this function._

## Implementation notes

**Attributes:** none registered.

## References

- Source: [`src/solve/solve.c`](https://github.com/stblake/mathilda/blob/main/src/solve/solve.c)
- Specification index: [`Mathilda_spec.md`](https://github.com/stblake/mathilda/blob/main/Mathilda_spec.md)
