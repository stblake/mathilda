# CoefficientArrays

!!! note "Status: Experimental"
    present and registered, but lightly documented and not yet covered by dedicated tests.

## Description

**`CoefficientArrays[polys, vars] gives the coefficient arrays of the polynomial system polys in the variables vars, grouped by total degree: element d + 1 holds every degree-d coefficient, so {b, m} = CoefficientArrays[eqs, vars] is the constant vector and coefficient matrix of a linear system. A degree-d coefficient sits at the index tuple naming its variables in non-decreasing order, every other permutation being 0; a list of polynomials carries the equation index as the leading axis. CoefficientArrays[poly, vars] treats a non-List first argument as a single polynomial. Option Modulus -> p reduces the coefficients. The arrays are dense Lists rather than SparseArrays.`**

## Examples

_No verified examples yet for this function._

## Implementation notes

**Attributes:** `Protected`.

## References

- Source: [`src/poly/monomials.c`](https://github.com/stblake/mathilda/blob/main/src/poly/monomials.c)
- Specification index: [`Mathilda_spec.md`](https://github.com/stblake/mathilda/blob/main/Mathilda_spec.md)
