---
references:
  - "D. Cox, J. Little and D. O'Shea, *Ideals, Varieties, and Algorithms*, 4th ed. (Springer, 2015), ch. 2 — monomial orderings."
source: src/poly/monomials.c
---
**Algorithm.** `MonomialList`, `CoefficientRules` and `FromCoefficientRules`
share the same core (`parse_and_build` → `build_monomials`). For
`MonomialList`: resolve the variable list (explicit, `All`, or default
`Variables[poly]`); reduce coefficients modulo an optional `Modulus -> m`;
`Expand` and split into additive terms; decompose each term into an integer
exponent vector plus a coefficient w.r.t. the variables; merge like monomials
and drop zero coefficients; and sort by a monomial order. Every named order
(`"Lexicographic"` default, `"DegreeLexicographic"`,
`"DegreeReverseLexicographic"`, and their `"Negative"` variants) is a special
case of descending lexicographic order of the weighted exponent vectors `w.v`,
so a single weight matrix built by `gb_build_order_matrix` — or an explicit user
weight matrix — drives one comparator (`cmp_key_desc`). `builtin_monomiallist`
then renders each `(expvec, coeff)` as `coeff * x1^e1 * x2^e2 * ...` via
`internal_times`/`internal_power` (`Times` drops the leading `1`);
`builtin_coefficientrules` renders the same monomial as `expvec -> coeff`.

**Data structures.** Each monomial is a `Mono { int* exps; Expr* coeff;
int64_t* key; }` — an owned integer exponent vector, an owned coefficient
expression, and a weighted sort key filled just before the `qsort`. When FLINT
is available and there is no modulus, the monomials are read straight off the
packed `fmpq_mpoly` (`build_monomials_polyQ`) or field polynomial
(`build_monomials_field`), skipping the generic `Expand` + per-term walk; the
modular path reduces integer coefficients into `[0, m)` first.

**Complexity / limits.** Dominated by the expansion/FLINT read-off and the
`O(n log n)` monomial sort over `n` terms. A symbolic structural head: the
exponent vectors and coefficients are exact `Expr` trees, so there is no
packed/NDArray or `Compile[]` path. `Protected`.
