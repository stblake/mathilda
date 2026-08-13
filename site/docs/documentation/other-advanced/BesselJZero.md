# BesselJZero

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`BesselJZero[n, k] gives the k-th positive zero of BesselJ[n, x]. Stays symbolic for symbolic arguments.`**

## Examples

_No verified examples yet for this function._

## Algorithm

besseljzero.c -- BesselJZero[n, k], the k-th positive zero of BesselJ[n, x].

Currently a symbolic placeholder: it stays unevaluated for all arguments so that the Hadamard-product recogniser (Product`BesselZero) can match the canonical infinite product

```text
  Product[1 - x^2/BesselJZero[n,k]^2, {k,1,Inf}] = Gamma[n+1] (2/x)^n BesselJ[n,x].
```

A future numeric path can evaluate BesselJZero[n, k] for exact numeric n and positive-integer k via McMahon asymptotics + Newton refinement on BesselJ.

Memory contract: takes ownership of res but must not free it; returns NULL (leave unevaluated) or an owned closed form.

## Implementation notes

**Attributes:** `Listable`, `NumericFunction`, `Protected`, `ReadProtected`.

## References

- Source: [`src/special_functions/besseljzero.c`](https://github.com/stblake/mathilda/blob/main/src/special_functions/besseljzero.c)
- Specification index: [`Mathilda_spec.md`](https://github.com/stblake/mathilda/blob/main/Mathilda_spec.md)
- Tests: [`tests/test_compile_coverage.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compile_coverage.c)
- Tests: [`tests/test_sum_product_families.c`](https://github.com/stblake/mathilda/blob/main/tests/test_sum_product_families.c)
