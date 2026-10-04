# MonomialList

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`MonomialList[poly] gives the list of monomials of poly, using Variables[poly].`**

**`MonomialList[poly, {x1, x2, ...}] uses the given variables; MonomialList[poly, vars, order]`**

<details>
<summary>Notes</summary>

sorts by order. order is "Lexicographic" (default), "DegreeLexicographic", "DegreeReverseLexicographic", "NegativeLexicographic", "NegativeDegreeLexicographic", "NegativeDegreeReverseLexicographic", or an explicit weight matrix. Modulus -\> m reduces coefficients modulo m. vars may be All (equivalent to Variables\[poly\]).

</details>

## Examples (6)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (2)

```mathematica
In[1]:= MonomialList[(x + y)^3]
Out[1]= {x^3, 3 x^2 y, 3 x y^2, y^3}

In[2]:= MonomialList[x^2 y^2 + x^3, {x, y}, "DegreeLexicographic"]
Out[2]= {x^2 y^2, x^3}
```

### Options (1)

```mathematica
In[3]:= MonomialList[(x + 1)^5, x, Modulus -> 2]
Out[3]= {x^5, x^4, x, 1}
```

### Applications (3)

The individual monomials, with coefficients

```mathematica
In[4]:= MonomialList[x^2 + 2 x y + y^2, {x, y}]
Out[4]= {x^2, 2 x y, y^2}
```

Expands first, so an unexpanded form works

```mathematica
In[5]:= MonomialList[(x + y)^2, {x, y}]
Out[5]= {x^2, 2 x y, y^2}
```

Highest total degree first

```mathematica
In[6]:= MonomialList[1 + x y + x^3, {x, y}, "DegreeLexicographic"]
Out[6]= {x^3, x y, 1}
```

## Implementation notes

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

- `Protected`.
- `MonomialList[poly]` is equivalent to `MonomialList[poly, Variables[poly]]`.
- Works whether or not `poly` is given in expanded form (it is expanded internally).
- `order` is `"Lexicographic"` (default), `"DegreeLexicographic"`,
  `"DegreeReverseLexicographic"`, `"NegativeLexicographic"`,
  `"NegativeDegreeLexicographic"`, `"NegativeDegreeReverseLexicographic"`, or an
  explicit weight matrix `w` (monomials ranked by the lexicographic order of `w.v`).
- `"NegativeLexicographic"` is `Sort` of the exponent vectors; `"Lexicographic"` is
  its reverse.
- `Modulus -> m` reduces coefficients modulo `m` (dropping monomials whose
  coefficient vanishes).
- `vars` may be `All` (equivalent to `Variables[poly]`).
- `Plus @@ MonomialList[poly, vars]` reconstructs the expanded polynomial.

**Attributes:** `Protected`.

## References

**See also:** [Sort](../../data-structures/Sort/)

- D. Cox, J. Little and D. O'Shea, *Ideals, Varieties, and Algorithms*, 4th ed. (Springer, 2015), ch. 2 — monomial orderings.
- Source: [`src/poly/monomials.c`](https://github.com/stblake/mathilda/blob/main/src/poly/monomials.c)
- Specification: [`docs/spec/builtins/structural-manipulation.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/structural-manipulation.md)
- Tests: [`tests/test_coefficient_rules.c`](https://github.com/stblake/mathilda/blob/main/tests/test_coefficient_rules.c)

## Notes & additional examples

### Notes

`MonomialList[poly, {x1, ..., xk}]` gives the list of monomials of `poly`, each
with its coefficient, so their sum is `poly`. The variables default to
`Variables[poly]` (or `All`), and an optional order argument sorts the monomials:
the default `"Lexicographic"`, the degree orders
(`"DegreeLexicographic"`/`"DegreeReverseLexicographic"`), their `"Negative"`
variants, or an explicit weight matrix; `Modulus -> m` reduces coefficients. The
polynomial is expanded first, so an unexpanded form such as `(x + y)^2` is
accepted. `CoefficientRules` is the same decomposition rendered as `expvec ->
coeff` rules.
