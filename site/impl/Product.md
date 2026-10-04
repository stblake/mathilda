---
source: src/product/product.c
---
**Algorithm.** `builtin_product` is the multiplicative analogue of `Sum`, and is
`HoldAll` so the iterator is not evaluated against an outer binding. It strips
trailing options (`Method`, `VerifyConvergence`), rewrites multiple iterators
`Product[f, s1, ..., sk]` into nested single-spec products, then dispatches on
the one spec in `product_one_spec`. A finite numeric range or an explicit list
is multiplied out directly: `expand_list`/`expand_range` bind the index and fold
the evaluated terms with `Times` (an empty product is `1`). Before the slow
evaluate-per-factor loop, a unit-step integer range whose body is inexact at the
first index takes a compiled machine multiply–accumulate (`product_try_compiled`,
via auto-compilation); an exact body (e.g. `n!`) stays on the interpreter so it
keeps exact bignum arithmetic.

**Data structures.** `Expr` trees, the `IterSpec` lattice parser
(`iter_spec_parse_lattice`), and the iterator shadow/restore pair
(`iter_spec_shadow`/`iter_spec_restore`) that localises the bound index. For
symbolic bounds, `Infinity`, or the indefinite form `Product[f, i]`, it runs a
`Method` cascade over context-qualified sub-builtins
(`Product\`Telescoping`, `Product\`Rational`, `Product\`Geometric`,
`Product\`QProduct`, then `Product\`Special`, `Product\`Cantor`,
`Product\`Viete`, `Product\`EulerPrime`, `Product\`RationalInfinite`,
`Product\`BesselZero`, `Product\`Infinite`, `Product\`LogSum`), each returning a
closed form or coming back unevaluated to signal "fall through"
(`result_is_unresolved`). The cascade is ordered cheapest/most-specific first so
the nicest closed form wins.

**Complexity / limits.** A closed-form stage is independent of the span width; a
finite enumeration is linear in the number of factors (guarded by
`PRODUCT_MAX_FINITE_TERMS`, 10^8). The closed-form stages assume a unit step, so
a non-unit step (`{i, 1, n, 2}`) with symbolic/over-wide bounds is left held
rather than given a wrong step-1 form. A finite range whose body hides an
index-dependent predicate (`EvenQ[k]`, `PrimeQ[k]`, ...) is forced to enumerate,
since symbolic evaluation would collapse the predicate and telescope the wrong
factor. When every stage falls through, `Product[...]` is returned unevaluated.
