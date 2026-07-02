# Apart

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
Apart[expr] rewrites a rational expression as a sum of terms with minimal denominators.
Apart[expr, var] treats all variables other than var as constants.
Option Extension -> alpha factors the denominator over Q(alpha) before
decomposition, splitting (x^2 - 2) into (x - Sqrt[2])(x + Sqrt[2]) under
Extension -> Sqrt[2] and producing the corresponding linear-factor
partial fractions.  Default Extension -> None decomposes over Q.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= Apart[1/((1+x)(5+x))]
Out[1]= -1/4/(5 + x) + 1/4/(1 + x)

In[2]:= Apart[(x^5-2)/((1+x+x^2)(2+x)(1-x))]
Out[2]= 2 - x - 34/9/(2 + x) + 1/9/(-1 + x) + (-1 - 1/3 x)/(1 + x + x^2)

In[3]:= Apart[(x+y)/((x+1)(y+1)(x-y)), x]
Out[3]= -(-1 + y)/((1 + x) (1 + y)^2) + 2 y/((1 + y)^2 (x - y))

In[4]:= Apart[1/(y^(2/3) - 1/y^(1/3))]
Out[4]= 1/3/(-1 + y^(1/3)) + (1/3 - 1/3 y^(1/3))/(1 + y^(1/3) + y^(2/3))

In[5]:= Apart[1/(x^2 - 2), x, Extension -> Sqrt[2]]
Out[5]= -1/2 1/(Sqrt[2] (Sqrt[2] + x)) + 1/2 1/(Sqrt[2] (-Sqrt[2] + x))
```

## Implementation notes

**Algorithm.** `builtin_apart` strips an optional `Extension -> α` (with `Automatic` running `extension_autodetect` for a single algebraic generator) and dispatches to `apart_impl`. Inexact coefficients route through `internal_rationalize_then_numericalize`. The core does an undetermined-coefficients partial-fraction decomposition over `Q` (or `Q(α)`):

1. Thread over `List`/relational/logical heads when the argument is one.
2. For a radical generator (fractional rational exponents), substitute `u -> g^m`, recurse, and back-substitute (`poly_find_radical_gen` / `poly_subst_radical_*`).
3. Combine to a single fraction with `Together`, pick the partial-fraction variable (explicit 2nd arg, else the lexicographically-last collected variable), and split into numerator `N` and denominator `D` (bailing back to the `Together`'d form if either is not polynomial in the variable, checked by `PolynomialQ`).
4. Polynomial-divide `N` by `D` (`PolynomialQuotient`/`PolynomialRemainder`) to get the polynomial part `Q` and proper remainder `R`.
5. `Factor` the denominator (over the extension if given), separate constant factors into a scalar `C`, and for each irreducible base `b_i` of multiplicity `k_i` set up unknown numerators of degree `deg(b_i)-1` over each power `b_i^j`.
6. Build the linear system by expanding each basis-times-`var^r` term, reading coefficients with `get_coeff`, and solving via `RowReduce`; assemble the sum `Q + sum A_ij / b_i^j`, factoring each solved coefficient.

**Data structures.** Everything is `Expr` trees driven through `eval_and_free`/`expr_expand`. The linear system is a dense `Expr***` augmented matrix (`S` rows by `S+1` columns, `S = deg(D)`) of coefficient expressions handed to `RowReduce`; factored denominator bases and their exponents are kept in parallel `Expr**`/`int64*` arrays.

**Complexity / limits.** Dominated by `Factor[D]` and the `O(S^3)` `RowReduce` over symbolic entries. Partial fractions are undefined when `N`/`D` are not polynomial in the chosen variable (handled by returning the combined fraction), and the matrix path assumes integer/rational coefficient structure.

- `Protected`, `Listable`.
- Writes `expr` as a polynomial in `var` together with a sum of ratios of polynomials with minimal denominators.
- If `var` is not specified, intelligently selects the main polynomial variable natively.
- Implements exact undetermined coefficients algebraically leveraging row-reduced identity expansions over algebraic inputs avoiding recursive fractional losses natively.
- When `Together[expr]` produces a numerator or denominator that is not polynomial in the chosen variable (e.g. fractional-power inputs whose Together'd form is `y^(1/3)/(y - 1)`), the matrix-of-coefficients algorithm cannot apply; Apart returns the `Together` form unchanged rather than synthesising a spurious zero.
- **Option `Extension -> alpha`** (Phase 0 of the Integrate plan) factors the denominator over `Q(alpha)` before partial-fraction decomposition runs, splitting reducible-over-extension factors (e.g. `x^2 - 2` into `(x - Sqrt[2])(x + Sqrt[2])` under `Extension -> Sqrt[2]`) and producing the corresponding linear-factor partial fractions. The pre-`Together` step is also extension-aware so any algebraic-number cancellations in numerator/denominator fire before splitting.

**Attributes:** `Listable`, `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Geddes, Czapor & Labahn, "Algorithms for Computer Algebra" (1992), on partial fraction decomposition.
- von zur Gathen & Gerhard, "Modern Computer Algebra", on partial fractions and the extended Euclidean algorithm.
- K. O. Geddes, S. R. Czapor and G. Labahn, *Algorithms for Computer Algebra* (Kluwer, 1992), ch. on partial-fraction decomposition.
- Source: [`src/parfrac.c`](https://github.com/stblake/mathilda/blob/main/src/parfrac.c)
- Specification: [`docs/spec/builtins/algebra.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/algebra.md)

## Notes & additional examples

### Worked examples

```mathematica
In[1]:= Apart[1/(x (x+1))]
Out[1]= 1/x - 1/(1 + x)
```

```mathematica
In[1]:= Apart[(x+2)/(x^2 - 1)]
Out[1]= -1/2/(1 + x) + 3/2/(-1 + x)
```

```mathematica
In[1]:= Apart[1/(x^2 (x+1))]
Out[1]= 1/x^2 - 1/x + 1/(1 + x)
```

```mathematica
In[1]:= Apart[(x^3 + 1)/(x^2 - 1)]
Out[1]= x + 1/(-1 + x)
```

```mathematica
In[1]:= Apart[(2 x + 3)/((x+1)^2 (x^2+1))]
Out[1]= 1/2/(1 + x)^2 + 3/2/(1 + x) + (1 - 3/2 x)/(1 + x^2)
```

```mathematica
In[1]:= Apart[1/(x^2 - 2), Extension -> Sqrt[2]]
Out[1]= -1/2 1/(Sqrt[2] (Sqrt[2] + x)) + 1/2 1/(Sqrt[2] (-Sqrt[2] + x))
```

### Notes

Apart computes the partial-fraction decomposition with respect to the (sole)
variable, the inverse of `Together`. It factors the denominator and splits the
quotient into a sum of terms whose denominators are the factors, including
repeated-factor terms: `1/(x^2 (x+1))` decomposes into `1/x^2 - 1/x + 1/(1+x)`,
recovering both the `1/x^2` and the lower-order `1/x` contributions. Rational
residues appear when the factors are linear over the rationals, as in the
`-1/2` and `3/2` coefficients for `(x+2)/(x^2-1)`. When the numerator degree
meets or exceeds the denominator's, Apart first divides out a polynomial part:
`(x^3 + 1)/(x^2 - 1)` becomes `x + 1/(-1 + x)`. Irreducible quadratic factors
over `Q` are kept intact (the `(1 - 3/2 x)/(1 + x^2)` term), while
`Extension -> Sqrt[2]` splits `x^2 - 2` into the conjugate linear factors
`x +- Sqrt[2]` and produces the corresponding partial fractions over
`Q(Sqrt[2])`.
