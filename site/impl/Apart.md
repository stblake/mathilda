---
references:
  - "K. O. Geddes, S. R. Czapor and G. Labahn, *Algorithms for Computer Algebra* (Kluwer, 1992), ch. on partial-fraction decomposition."
source: src/parfrac.c
---
**Algorithm.** `builtin_apart` strips an optional `Extension -> α` (with `Automatic` running `extension_autodetect` for a single algebraic generator) and dispatches to `apart_impl`. Inexact coefficients route through `internal_rationalize_then_numericalize`. The core does an undetermined-coefficients partial-fraction decomposition over `Q` (or `Q(α)`):

1. Thread over `List`/relational/logical heads when the argument is one.
2. For a radical generator (fractional rational exponents), substitute `u -> g^m`, recurse, and back-substitute (`poly_find_radical_gen` / `poly_subst_radical_*`).
3. Combine to a single fraction with `Together`, pick the partial-fraction variable (explicit 2nd arg, else the lexicographically-last collected variable), and split into numerator `N` and denominator `D` (bailing back to the `Together`'d form if either is not polynomial in the variable, checked by `PolynomialQ`).
4. Polynomial-divide `N` by `D` (`PolynomialQuotient`/`PolynomialRemainder`) to get the polynomial part `Q` and proper remainder `R`.
5. `Factor` the denominator (over the extension if given), separate constant factors into a scalar `C`, and for each irreducible base `b_i` of multiplicity `k_i` set up unknown numerators of degree `deg(b_i)-1` over each power `b_i^j`.
6. Build the linear system by expanding each basis-times-`var^r` term, reading coefficients with `get_coeff`, and solving via `RowReduce`; assemble the sum `Q + sum A_ij / b_i^j`, factoring each solved coefficient.

**Data structures.** Everything is `Expr` trees driven through `eval_and_free`/`expr_expand`. The linear system is a dense `Expr***` augmented matrix (`S` rows by `S+1` columns, `S = deg(D)`) of coefficient expressions handed to `RowReduce`; factored denominator bases and their exponents are kept in parallel `Expr**`/`int64*` arrays.

**Complexity / limits.** Dominated by `Factor[D]` and the `O(S^3)` `RowReduce` over symbolic entries. Partial fractions are undefined when `N`/`D` are not polynomial in the chosen variable (handled by returning the combined fraction), and the matrix path assumes integer/rational coefficient structure.
