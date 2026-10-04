---
source: src/special_functions/besseljzero.c
---
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
