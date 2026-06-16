# FactorSquareFree

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
FactorSquareFree[poly]
    writes poly as a product of pairwise-coprime square-free factors,
    collecting repeated factors into powers.
Computed via the Yun / Musser square-free decomposition using
polynomial GCDs of poly with its derivative; cheaper than full Factor
and sufficient when only multiplicities are needed.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= FactorSquareFree[x^5 - x^3 - x^2 + 1]
Out[1]= (-1 + x)^2 (1 + 2 x + 2 x^2 + x^3)

In[2]:= FactorSquareFree[x^4 - 9x^3 + 29x^2 - 39x + 18]
Out[2]= (-3 + x)^2 (2 - 3 x + x^2)

In[3]:= FactorSquareFree[x^5 - x^3 y^2 - x^2 y^3 + y^5]
Out[3]= (-x + y)^2 (x^3 + 2 x^2 y + 2 x y^2 + y^3)

In[4]:= FactorSquareFree[{(x^2 - 1)(x - 1), (x^4 - 1)(x^2 - 1)}]
Out[4]= {(1 + x) (-1 + x)^2, (1 + x^2) (-1 + x^2)^2}
```

## Implementation notes

**Algorithm.** `builtin_factorsquarefree` (in `src/poly/facpoly_squarefree.inc`, compiled into `facpoly.c`) computes the square-free decomposition `p = ∏ p_i^i` with each `p_i` square-free and pairwise coprime, using **Yun's algorithm**. After a `Rationalize`/`Numericalize` round-trip for inexact inputs, the dispatcher `factor_square_free_dispatcher` expands and routes to `factor_square_free_poly`.

`factor_square_free_poly` works recursively on the variable list. It first extracts the polynomial content (`poly_content`) and recurses on it; on the primitive part `pp` it runs Yun's recurrence: set `A = pp`, `B = gcd(A, A')`, `C = A/B`, `D = A'/B − C'`, then iterate `P_i = gcd(C, D)` (the degree-`i` square-free factor), `C ← C/P_i`, `D ← D/P_i − (C/P_i)'`, accumulating `P_i^i`. Derivatives are computed by the local `poly_deriv` and GCDs by `poly_gcd_internal`; exact divisions by `exact_poly_div`. An F4 "cheap squarefree" pre-check (`sqfree_cheap_check`, a univariate-substitution probe) short-circuits the expensive `gcd(pp, pp')` when `pp` is provably square-free. A final `exact_poly_div` recovers any leading/missing scalar factor so the product round-trips to the input.

**Data structures.** `Expr*` polynomials manipulated through the generic evaluator (`eval_and_free`) and the poly helpers (`get_degree_poly`, `poly_gcd_internal`, `exact_poly_div`); factors collected in a growable `Expr**` buffer with their multiplicities encoded as `Power[P_i, i]`.

**Complexity / limits.** Dominated by the GCD chain; the cheap pre-check avoids one full multivariate GCD on the common square-free case. Square-free decomposition is a prerequisite stage for full `Factor`.

- `Listable`, `Protected`.
- Automatically threads over lists, as well as equations, inequalities and logic functions.
- Works on both univariate and multivariate polynomials.
- Multivariate inputs use a cheap squarefree pre-check (F4 Stage 1): after content extraction in the main variable, `sqfree_cheap_check` substitutes integer values from `{1, -1, 2, -2, 3, -3, 4}` for the other variables and tests `gcd(image, image')` over `Z[x]`.  If any image is squarefree at an alpha that preserves the leading-x degree, the pre-check proves `pp` is squarefree in x and the expensive multivariate `gcd(pp, pp')` is skipped.  Soundness comes from content extraction guaranteeing any repeated factor of `pp` involves the main variable nontrivially.  Measured 6.6× speedup on 4-variable squarefree inputs (6.27 s → 0.95 s); non-squarefree inputs fall through to the original Yun loop with negligible overhead.
- The cheap pre-check's univariate `gcd(image, image')` runs through `zupoly_gcd` (subresultant PRS, GMP `mpz_t` coefficients).  The previous implementation used `poly_gcd_internal` (Knuth-style primitive PRS at the Expr level) which suffers exponential coefficient growth on the intermediate pseudo-remainders; on a degree-31 univariate image (e.g. `Factor[Expand[x^2 (z^13 - x^12)(z^4 + 3 x^9 - y^13)(17 - 5 y - z^14)]]`) it ran for >120 s.  Routing the same gcd through subresultant PRS keeps coefficient sizes polynomially bounded and runs in sub-millisecond time on the same input, bringing the full Factor call to under 1 s.

**Attributes:** `Listable`, `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- D. Y. Y. Yun, "On square-free decomposition algorithms", SYMSAC 1976.
- K. O. Geddes, S. R. Czapor, G. Labahn, *Algorithms for Computer Algebra* (Kluwer, 1992).
- Source: [`src/poly/facpoly_squarefree.inc`](https://github.com/stblake/mathilda/blob/main/src/poly/facpoly_squarefree.inc)
- Specification: [`docs/spec/builtins/structural-manipulation.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/structural-manipulation.md)

## Notes & additional examples

### Worked examples

```mathematica
In[1]:= FactorSquareFree[x^5 - x^4 - x + 1]
Out[1]= (-1 + x)^2 (1 + x + x^2 + x^3)
```

```mathematica
In[1]:= FactorSquareFree[(x^2+1)^3 (x-1)^2]
Out[1]= (-1 + x)^2 (1 + x^2)^3
```

```mathematica
In[1]:= FactorSquareFree[x^8 + 4 x^6 + 6 x^4 + 4 x^2 + 1]
Out[1]= (1 + x^2)^4
```

### Notes

`FactorSquareFree` groups a polynomial into pairwise-coprime square-free factors carrying their multiplicities, using the Yun / Musser decomposition (GCDs of the polynomial with its derivative). It does not split factors that are square-free but reducible: `1 + x + x^2 + x^3` is left intact in the first example even though it factors as `(1+x)(1+x^2)`. The last example recognises `(1 + x^2)^4` as a perfect fourth power without ever calling full `Factor`, which is the point — it is cheaper than factorisation when only multiplicities are needed.
