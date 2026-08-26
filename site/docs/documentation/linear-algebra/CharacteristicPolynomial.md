# CharacteristicPolynomial

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`CharacteristicPolynomial[m, x]`**

gives the characteristic polynomial Det\[m - x I\] of the square matrix m, as a polynomial in x.

**`CharacteristicPolynomial[{m, a}, x]`**

gives the generalized characteristic polynomial Det\[m - x a\].

<details>
<summary>Notes</summary>

m may have integer, rational, machine- or arbitrary-precision real, complex, or symbolic entries.  The ordinary case is computed by the Faddeev-LeVerrier method (O(n^4)); the generalized case by Laplace expansion of m - x a.

</details>

## Examples (10)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (3)

```mathematica
In[1]:= CharacteristicPolynomial[{{1, 2}, {3, 4}}, x]
Out[1]= -2 - 5 x + x^2

In[2]:= CharacteristicPolynomial[{{a, b}, {c, d}}, x]
Out[2]= -b c + a d - a x - d x + x^2

In[3]:= CharacteristicPolynomial[{{{1, 1, 1}, {1, 0, 1}, {0, 0, 1}}, {{0, 1, 1}, {0, 1, 1}, {1, 0, 0}}}, x]
Out[3]= -1 - x + x^2
```

### Worked examples (1)

```mathematica
In[4]:= CharacteristicPolynomial[IdentityMatrix[3], x]
Out[4]= 1 - 3 x + 3 x^2 - x^3
```

### Applications (6)

```mathematica
In[5]:= CharacteristicPolynomial[{{1, 2}, {3, 4}}, x]
Out[5]= -2 - 5 x + x^2

In[6]:= CharacteristicPolynomial[{{a, b}, {c, d}}, x]
Out[6]= -b c + a d - a x - d x + x^2

In[7]:= CharacteristicPolynomial[IdentityMatrix[3], x]
Out[7]= 1 - 3 x + 3 x^2 - x^3

In[8]:= CharacteristicPolynomial[{{1/3, 1/2, 3/5}, {1/2, 4/5, 1}, {3/5, 1, 9/7}}, x]
Out[8]= 1/10500 - 239/2100 x + 254/105 x^2 - x^3

In[9]:= CharacteristicPolynomial[{{{1, 2}, {5, 4}}, {{4, 3}, {6, 4}}}, x]
Out[9]= -6 + 7 x - 2 x^2

In[10]:= CharacteristicPolynomial[{{{1, 1, 1}, {1, 0, 1}, {0, 0, 1}}, {{0, 1, 1}, {0, 1, 1}, {1, 0, 0}}}, x]
Out[10]= -1 - x + x^2
```

## Algorithm

CharacteristicPolynomial -- the characteristic polynomial of a square matrix.

```text
  CharacteristicPolynomial[m, x]      == Det[m - x I]
  CharacteristicPolynomial[{m, a}, x] == Det[m - x a]   (generalized)
```

This is exactly the polynomial Eigenvalues solves for, so it reuses the eigen module's char-poly machinery (declared in eigen_internal.h):

```text
  - Ordinary case: Faddeev-Leverrier-Souriau (O(n^4)), which returns
    det(lambda I - m).  Wolfram's CharacteristicPolynomial is
    Det[m - x I] = (-1)^n det(x I - m), so the odd-n result is negated.
    This keeps the 100x100 machine-matrix case sub-second: the naive
    Expand[Det[m - x I]] would face an O(n!) Laplace expansion of a
    symbolic-in-x matrix.
  - Generalized case: build (m - lambda a) and take its Laplace determinant,
    which is det(m - lambda a) with the correct sign directly.  Generalized
    inputs are small in practice (<= 3x3), matching the existing generalized
    Eigenvalues path.
```

The polynomial is built in a private internal lambda symbol and the user's variable is substituted in at the end, so the second argument may be a symbol, a number, or any expression (the char-poly coefficient identity makes

```text
det(lambda I - m)|_{lambda -> x} valid even when m itself contains x).  The
```

result is Expand-ed (Det / Laplace do not multiply out) to reach the collected polynomial form.

This head returns a symbolic Plus expression, not a machine buffer, so it is

```text
exempt from the packed/NDArray/Compile surfaces.  It consumes packed input
```

safely: the NDArray guard delists a visible NDArray, faddeev pack_unpacks its own input, and the generalized path materialises m / a before indexing them.

## Implementation notes

**Algorithm.** `builtin_characteristicpolynomial` (`src/linalg/charpoly.c`) is a thin builtin over the eigen module's characteristic-polynomial machinery — the same polynomial `Eigenvalues` solves for. `eigen_extract_matrix_pair` classifies the argument as an ordinary square matrix `m` or a generalised pencil `{m, a}` and validates squareness.

*Ordinary case* `CharacteristicPolynomial[m, x]`: the polynomial is formed by the **Faddeev–LeVerrier–Souriau** recurrence (`eigen_char_poly_faddeev`), which builds the coefficients in `O(n^4)` matrix multiplications — far cheaper than the naïve `Det[m - x I]`, which would face an `O(n!)` Laplace expansion of a symbolic-in-`x` matrix (so a 100×100 machine matrix is sub-second). That routine returns `det(λI − m)`, so for odd `n` the result is negated to match `Det[m − x I] = (−1)^n det(x I − m)` — hence odd-degree characteristic polynomials are monic-negative.

*Generalised case* `CharacteristicPolynomial[{m, a}, x]`: `eigen_build_lambda_matrix` forms the entry-wise matrix `m − λa` and `eigen_compute_det` takes its determinant by Laplace expansion, which already carries the correct sign. A null space shared by `m` and `a` drops the leading term(s), so an infinite generalised eigenvalue shows as a degree deficit.

The polynomial is built in a private internal lambda symbol; the user's second argument (a symbol, a number, or any expression) is then substituted for it via `ReplaceAll` and the result is expanded (`Det`/Laplace do not multiply out). Entries may be integer, rational, machine- or arbitrary-precision real, complex, or symbolic.

**Complexity / limits.** Faddeev–LeVerrier is `O(n^4)` for the ordinary case; the generalised case is Laplace expansion (`O(n!)`, in practice used only for small pencils, matching the generalised `Eigenvalues` path). The result is a symbolic `Plus`, so the head is exempt from the packed-array/`Compile` surfaces, but it accepts a packed or visible `NDArray` matrix as input (the buffer is materialised before the symbolic recurrence runs).

- `Protected`.
- Entries may be integer, rational, machine- or arbitrary-precision real,
  complex, or symbolic. The result is an expanded polynomial in `x`.
- The ordinary case reuses the eigen module's Faddeev–LeVerrier–Souriau fast
  path (`O(n^4)` matrix multiplications), so the characteristic polynomial of a
  large numeric matrix is computed in polynomial time — the naïve
  `Det[m - x I]` would face an `O(n!)` Laplace expansion of a symbolic-in-`x`
  matrix. The leading `(-1)^n` sign of `Det[m - x I]` is applied so odd-`n`
  polynomials are monic-negative (`CharacteristicPolynomial[IdentityMatrix[3],
  x]` → `1 - 3 x + 3 x^2 - x^3`).
- The generalised case is `Det[m - x a]` via Laplace expansion. A shared null
  space of `m, a` drops the leading term(s): an infinite generalised eigenvalue
  shows as a degree deficit (`CharacteristicPolynomial[{a, b}, x]` of degree
  `< n`).
- The second argument may be a symbol, a number, or any expression (the value
  is substituted for the polynomial variable).
- Called with other than two arguments, emits `CharacteristicPolynomial::argrx`
  and stays unevaluated; a non-square matrix stays unevaluated.

**Attributes:** `Protected`.

## References

- G. H. Golub and C. F. Van Loan, *Matrix Computations*, 4th ed., Johns Hopkins University Press, 2013 — the characteristic polynomial and the eigenvalue problem.
- D. K. Faddeev and V. N. Faddeeva, *Computational Methods of Linear Algebra*, W. H. Freeman, 1963 — the Faddeev–LeVerrier–Souriau recurrence.
- Gene H. Golub, Charles F. Van Loan, *Matrix Computations*, 4th ed. (Johns Hopkins University Press, 2013).
- D. K. Faddeev and V. N. Faddeeva, *Computational Methods of Linear Algebra* (W. H. Freeman, 1963) — the Faddeev–LeVerrier–Souriau recurrence.
- Source: [`src/linalg/charpoly.c`](https://github.com/stblake/mathilda/blob/main/src/linalg/charpoly.c)
- Specification: [`docs/spec/builtins/linear-algebra.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/linear-algebra.md)
- Tests: [`tests/test_characteristicpolynomial.c`](https://github.com/stblake/mathilda/blob/main/tests/test_characteristicpolynomial.c)

## Notes & additional examples

### Notes

`CharacteristicPolynomial[m, x]` is `Det[m - x I]` and `CharacteristicPolynomial[{m, a}, x]` is `Det[m - x a]`, returned as an expanded polynomial in `x`. It is the polynomial whose roots are the (generalised) eigenvalues of `m`, so it pairs naturally with `Eigenvalues`.

The ordinary case is computed by the Faddeev–LeVerrier recurrence in `O(n^4)`, so the characteristic polynomial of a large numeric matrix is found in polynomial time rather than through an `O(n!)` symbolic determinant. For an odd-order matrix the polynomial is monic-negative (leading term `-x^n`), matching `Det[m - x I]`. In the generalised case a shared null space of `m` and `a` lowers the degree — the missing leading term corresponds to an infinite generalised eigenvalue, so `Eigenvalues[{m, a}]` returns `Infinity` for each degree drop. The second argument may also be a number or an expression, in which case the polynomial is evaluated at that value.
