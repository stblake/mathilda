# BesselJZero

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`BesselJZero[n, k] gives the k-th positive zero of BesselJ[n, x]. Stays symbolic for symbolic arguments.`**

## Examples (3)

Every input below was run against the current Mathilda build and its output recorded.

### Applications (3)

No numeric path yet -- stays symbolic

```mathematica
In[1]:= BesselJZero[0, 1]
Out[1]= BesselJZero[0, 1]
```

A bare, unevaluated head

```mathematica
In[2]:= Head[BesselJZero[n, k]]
Out[2]= BesselJZero
```

The Hadamard product it exists to name

```mathematica
In[3]:= Product[1 - x^2/BesselJZero[n, k]^2, {k, 1, Infinity}]
Out[3]= Gamma[1 + n] BesselJ[n, x] (2/x)^n
```

## Algorithm

besseljzero.c -- BesselJZero[n, k], the k-th positive zero of BesselJ[n, x].

Currently a symbolic placeholder: it stays unevaluated for all arguments so that the Hadamard-product recogniser (Product`BesselZero) can match the canonical infinite product

```text
  Product[1 - x^2/BesselJZero[n,k]^2, {k,1,Inf}] = Gamma[n+1] (2/x)^n BesselJ[n,x].
```

A future numeric path can evaluate BesselJZero[n, k] for exact numeric n and positive-integer k via McMahon asymptotics + Newton refinement on BesselJ.

Memory contract: takes ownership of res but must not free it; returns NULL (leave unevaluated) or an owned closed form.

## Implementation notes

**Algorithm.** `builtin_besseljzero` is a **symbolic placeholder**: for every argument
it returns `NULL`, so `BesselJZero[n, k]` stays unevaluated. It carries no numeric
path yet. Its purpose is to be the canonical name for the `k`-th positive zero of
`BesselJ[n, x]` so that the infinite-product recogniser in the `Product` subsystem
(`Product\`BesselZero`) can match the Hadamard product

```
Product[1 - x^2/BesselJZero[n, k]^2, {k, 1, Infinity}]
    = Gamma[n + 1] (2/x)^n BesselJ[n, x].
```

A future numeric path could evaluate `BesselJZero[n, k]` for exact numeric `n` and
positive-integer `k` via McMahon asymptotics plus Newton refinement on `BesselJ`, as
the source header notes.

**Data structures.** None of its own — it is a bare `EXPR_FUNCTION` head that the
evaluator leaves intact. The product identity above is applied by the `Product` family,
not by `BesselJZero`.

**Complexity / limits.** `O(1)` (immediate decline). Attributes: `Listable`,
`NumericFunction`, `Protected`. No numeric evaluation is implemented, so even
`BesselJZero[0, 1]` remains symbolic; the symbol is useful today only as the product
index that produces the closed Bessel form.

**Attributes:** `Listable`, `NumericFunction`, `Protected`.

## References

- Source: [`src/special_functions/besseljzero.c`](https://github.com/stblake/mathilda/blob/main/src/special_functions/besseljzero.c)
- Specification index: [`Mathilda_spec.md`](https://github.com/stblake/mathilda/blob/main/Mathilda_spec.md)
- Tests: [`tests/test_compile_coverage.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compile_coverage.c)
- Tests: [`tests/test_sum_product_families.c`](https://github.com/stblake/mathilda/blob/main/tests/test_sum_product_families.c)

## Notes & additional examples

### Notes

`BesselJZero[n, k]` names the `k`-th positive zero of `BesselJ[n, x]`. It is currently
a **symbolic placeholder**: it carries no numeric evaluator, so every call — even
`BesselJZero[0, 1]` — stays unevaluated.

Its reason for existing today is the infinite-product recogniser in the `Product`
subsystem, which matches the canonical Hadamard product

```
Product[1 - x^2/BesselJZero[n, k]^2, {k, 1, Infinity}] = Gamma[1 + n] (2/x)^n BesselJ[n, x]
```

and so reproduces the closed Bessel form directly. `BesselJZero` is `Listable`,
`NumericFunction` and `Protected`.
