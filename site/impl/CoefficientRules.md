---
references:
  - "D. Cox, J. Little and D. O'Shea, *Ideals, Varieties, and Algorithms*, 4th ed. (Springer, 2015), ch. 2 — monomial orderings."
source: src/poly/monomials.c
---
**Algorithm.** `builtin_coefficientrules` is the `{exponent-vector -> coefficient}`
rendering of the shared monomial core (see [MonomialList](MonomialList.md)): it
calls the same `parse_and_build` to resolve the variables, reduce by an optional
`Modulus -> m`, `Expand` and split into terms, decompose each term into an
integer exponent vector and coefficient, merge like monomials, and sort by the
requested monomial order (same named orders and explicit weight matrices as
`MonomialList`, default `"Lexicographic"`). It then emits one
`Rule[{e1, ..., ek}, coeff]` per surviving monomial, wrapped in a `List`. Because
the core expands first, it works whether or not `poly` is already expanded.
`FromCoefficientRules` is the exact inverse.

**Data structures.** The same `Mono { int* exps; Expr* coeff; int64_t* key; }`
array as `MonomialList`; the only difference is the output shape — here the
exponent vector becomes a `List` of `Integer`s and is paired with the (moved,
not copied) coefficient under a `Rule`. When FLINT is available and there is no
modulus the monomials are read directly off the packed `fmpq_mpoly`/field
polynomial, bypassing the generic per-term walk.

**Complexity / limits.** Dominated by the expansion/FLINT read-off and the
`O(n log n)` monomial sort. A symbolic structural head — exact exponent vectors
and coefficient trees, so no packed/NDArray or `Compile[]` path. `Protected`.
