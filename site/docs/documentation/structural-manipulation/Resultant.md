# Resultant

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
Resultant[p, q, var]
    gives the resultant of p and q as polynomials in var: the unique
    integer / polynomial scalar that vanishes iff p and q share a
    root in var.  Computed via a Sylvester-matrix determinant or, in
    the exact path, a subresultant pseudo-remainder sequence.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= Resultant[x^2 - 2x + 7, x^3 - x + 5, x]
Out[1]= 265

In[2]:= Resultant[x^3 - 5x^2 - 7x + 14, x^3 - 8x^2 + 9x + 58, x]
Out[2]= 0
```

## Implementation notes

- `Protected`, `Listable`.
- Computes the resultant of polynomials `poly1` and `poly2` with respect to the variable `var`.
- The resultant is independent of common roots and vanishes exactly when the polynomials have roots in common.
- Default algorithm is Bronstein's subresultant PRS (Symbolic Integration I, p.24): a linear chain of pseudo-remainders with scalar exact divisions in the coefficient ring, avoiding the (n+m)x(n+m) Sylvester matrix construction and its O(n^3) Bareiss reduction.  For Z/Q coefficients this is materially faster than the matrix path, and on inputs with symbolic coefficients it sidesteps the O(n!) Laplace expansion that the matrix path falls back to when Bareiss exact-division certification fails.
- Inputs containing algebraic-number coefficients (e.g. `Sqrt[N]`, cube roots — any `Power[X, Rational[a,b]]` with `b > 1`) are routed to the Sylvester+Det path instead, because the subresultant chain bloats geometrically when `Power[base, k/m]` forms can't be combined with their `Times[base^q, Sqrt[base]]` equivalents by `Plus` alone.
- A size-budget guard inside the subresultant path falls back to Sylvester+Det for any pathological input where chain elements exceed ~30x the input leaf-count.
- Automatically preserves multiplicativity (e.g., $Res(A \cdot B, Q) = Res(A, Q) Res(B, Q)$ and $Res(A^k, Q) = Res(A, Q)^k$).

**Attributes:** `Listable`, `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- von zur Gathen & Gerhard, "Modern Computer Algebra" (3rd ed.), Ch. 6 (resultants and the Sylvester matrix).
- Geddes, Czapor & Labahn, "Algorithms for Computer Algebra" (1992), Ch. 7 (subresultant PRS).
- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/structural-manipulation.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/structural-manipulation.md)

## Notes & additional examples

### Worked examples

```mathematica
In[1]:= Resultant[x^2 - 1, x^2 - 4, x]
Out[1]= 9
```

```mathematica
In[1]:= Resultant[x^2 - 2, x^2 - 3, x]
Out[1]= 1
```

```mathematica
In[1]:= Resultant[x^2 + a, x + b, x]
Out[1]= a + b^2
```

```mathematica
In[1]:= Resultant[x^2 - y, x^2 + y, x]
Out[1]= 4 y^2
```

### Notes

`Resultant[p, q, x]` returns a scalar in the remaining variables that
vanishes exactly when `p` and `q` share a common root in `x`. For
`x^2 - 1` and `x^2 - 4` (roots `±1` and `±2`) the value is the nonzero `9`,
confirming no shared root, whereas eliminating `x` from `x^2 - y` and
`x^2 + y` gives `4 y^2`, which is zero only at the shared-root locus `y = 0`.
The computation uses either a Sylvester-matrix determinant or, on the exact
path, a subresultant pseudo-remainder sequence; coefficients may themselves
be polynomials in the other variables, as in `a + b^2`.
