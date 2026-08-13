# FactorTerms

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`FactorTerms[poly]`**

pulls out any overall numerical factor in poly.

**`FactorTerms[poly, x]`**

pulls out any overall factor in poly that does not depend on x.

**`FactorTerms[poly, {x1, x2, ...}]`**

pulls out any overall factor in poly that does not depend on any of the xi, then progressively factors with respect to smaller subsets {x1, ..., x\_{k-1}}.

**`FactorTerms[poly, x] extracts the content of poly with respect to x.`**

<details>
<summary>Notes</summary>

FactorTerms automatically threads over lists, equations, inequalities and logic functions.

</details>

## Examples (11)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (8)

```mathematica
In[1]:= FactorTerms[3 + 6x + 3x^2]
Out[1]= 3 (1 + 2 x + x^2)

In[2]:= FactorTerms[3 + 3a + 6 a x + 6 x + 12 a x^2 + 12 x^2, x]
Out[2]= 3 (1 + a) (1 + 2 x + 4 x^2)

In[3]:= FactorTerms[12 a^4 + 9 x^2 + 66 b^2]
Out[3]= 3 (4 a^4 + 22 b^2 + 3 x^2)

In[4]:= FactorTerms[7 x + (14 y + 21)/z]
Out[4]= (7 (3 + 2 y + x z))/z

In[5]:= FactorTerms[{5 x^2 - 15, 7 x^4 - 77, 8 x^8 - 24}]
Out[5]= {5 (-3 + x^2), 7 (-11 + x^4), 8 (-3 + x^8)}

In[6]:= FactorTerms[1 < 77 x^3 - 21 x + 35 < 2]
Out[6]= 1 < 7 (5 - 3 x + 11 x^3) < 2

In[7]:= f = 2 x^2 y z + 2 x^2 y + 4 x^2 z + 4 x^2 + 4 y^2 z^2 + 4 z y^2 + 8 z^2 y + 2 z y - 6 y - 12 z - 12; FactorTerms[f, x]
Out[7]= 2 (2 + y + 2 z + y z) (-3 + x^2 + 2 y z)

In[8]:= FactorTerms[f, {x, y}]
Out[8]= 2 (2 + y) (1 + z) (-3 + x^2 + 2 y z)
```

### Applications (3)

```mathematica
In[9]:= FactorTerms[6 x^2 + 4 x]
Out[9]= 2 (2 x + 3 x^2)

In[10]:= FactorTerms[2 x^2 + 4 x + 2]
Out[10]= 2 (1 + 2 x + x^2)

In[11]:= FactorTerms[3 x^2 y + 6 x y^2, x]
Out[11]= 3 y (x^2 + 2 x y)
```

## Implementation notes

**Algorithm.** `builtin_factorterms` (in `src/poly/facpoly_factorterms.inc`, compiled into `facpoly.c`) pulls out the content of a polynomial without factoring the polynomial part. It threads over `List`/equation/inequality/logic heads (`ft_is_threading_head`), then calls the shared engine `ft_compute_list` and multiplies the resulting factor list back into a single `Times`.

`ft_compute_list` first `Together`-normalises and splits into numerator/denominator. It collects and sorts the numerator's variables, then: (1) extracts the **numerical content** via `ft_content_wrt_set` (content with respect to *all* variables over an empty ground ring, i.e. the integer GCD of coefficients) and divides it out with `ft_divide_out` (exact polynomial division, falling back to symbolic `Times[poly, content^{-1}]`); (2) for a 2-arg call `FactorTerms[poly, {x_1,…,x_k}]`, peels the content with respect to progressively smaller variable subsets, where `ft_content_wrt_set` recursively computes the multivariate `poly_gcd_internal` of the coefficients of each monomial in the chosen variables over the shrinking ground ring; (3) appends the final residue (re-multiplied by `1/den` to round-trip rational inputs).

**Data structures.** `Expr*` polynomials throughout; `poly_gcd_internal` (multivariate polynomial GCD) is the workhorse for symbolic content; variable lists are `Expr**` sorted with `compare_expr_ptrs`.

- `Protected`.
- Auto-threads over `List`, `Equal`, `Unequal`, `Less`, `LessEqual`,
  `Greater`, `GreaterEqual`, `And`, `Or`, `Not`, `Xor`.
- Together-normalises rational inputs before extracting content from the
  numerator, then divides the denominator back through, so rational
  functions round-trip.
- Numerical content over `Z` is computed via the integer GCD of monomial
  coefficients. Gaussian-integer content (e.g. `5 I` from `5 I x^2 + ...`)
  is not extracted as a Gaussian unit; the integer GCD `5` is returned
  instead, with the leading factor of `I` left inside the residue. The
  resulting factorization is mathematically equivalent.

**Attributes:** `Protected`.

## References

**See also:** [List](../../other-advanced/List/), [Equal](../../comparisons/Equal/), [Unequal](../../comparisons/Unequal/), [Less](../../comparisons/Less/), [LessEqual](../../comparisons/LessEqual/), [Greater](../../comparisons/Greater/), [GreaterEqual](../../comparisons/GreaterEqual/), [I](../../mathematical-constants/I/)

- Source: [`src/poly/facpoly_factorterms.inc`](https://github.com/stblake/mathilda/blob/main/src/poly/facpoly_factorterms.inc)
- Specification: [`docs/spec/builtins/structural-manipulation.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/structural-manipulation.md)
- Tests: [`tests/test_factor_terms.c`](https://github.com/stblake/mathilda/blob/main/tests/test_factor_terms.c)

## Notes & additional examples

### Notes

`FactorTerms[poly]` pulls out the overall numerical content (the integer GCD of the coefficients) without touching the polynomial structure — unlike `Factor`, it never splits `1 + 2 x + x^2` into `(1 + x)^2`. With a second argument `FactorTerms[poly, x]` extracts the content with respect to `x`, i.e. the factor that does not depend on `x`: here `3 y` is the part free of `x`, leaving `x^2 + 2 x y`.
