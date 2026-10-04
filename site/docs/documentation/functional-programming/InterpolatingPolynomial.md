# InterpolatingPolynomial

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`InterpolatingPolynomial[{f1, f2, ...}, x]`**

gives the single polynomial in x reproducing the values fi at x = 1, 2, ..., in nested (Horner) form. With n values the degree is n-1.

**`InterpolatingPolynomial[{{x1, f1}, {x2, f2}, ...}, x]`**

interpolates the values fi at the abscissae xi (arbitrary real, complex, or -- in 1-D -- symbolic).

**`InterpolatingPolynomial[{{{x1, y1, ...}, f1}, ...}, {x, y, ...}]`**

gives the multidimensional interpolating polynomial of lowest total degree.

**`InterpolatingPolynomial[{{xi, fi, dfi, ...}, ...}, x]`**

reproduces derivatives as well as values (the n-th derivative in m-D is a tensor shaped like D\[f, {{x, ...}, n}\]).

<details>
<summary>Notes</summary>

A value or derivative given as Automatic is filled in from the other conditions. The option Modulus -\> n finds the polynomial over the integers modulo n. Exact data give an exact polynomial.

</details>

## Examples (10)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (4)

```mathematica
In[1]:= InterpolatingPolynomial[{1, 4, 9, 16}, x]
Out[1]= 1 + (-1 + x) (1 + x)
```

Value 8, slope 0 at x=4

```mathematica
In[2]:= InterpolatingPolynomial[{4, 7, 2, {8, 0}, 9}, x]
Out[2]= 4 + (-1 + x) (3 + (-2 + x) (-4 + (-3 + x) (19/6 + (-4 + x) (-107/36 + 109/72 (-4 + x)))))
```

```mathematica
In[3]:= Expand[InterpolatingPolynomial[ {{{0, 0}, 1}, {{1, 0}, 7}, {{0, 1}, 10}, {{2, 1}, 40}, {{3, 3}, 151}, {{1, 2}, 47}}, {x, y}]]
Out[3]= 1 + 2 x + 4 x^2 + 3 y + 5 x y + 6 y^2

In[4]:= Expand[InterpolatingPolynomial[ {{-1, Automatic, 0}, {0, 1, 1}, {1, Automatic, 0}}, x]]
Out[4]= 1 + x - 1/3 x^3
```

### Options (1)

```mathematica
In[5]:= InterpolatingPolynomial[{1, 4, 9, 16}, x, Modulus -> 7]
Out[5]= x^2
```

### Applications (5)

Values at x = 1, 2, 3, in nested Newton form

```mathematica
In[6]:= InterpolatingPolynomial[{1, 4, 9}, x]
Out[6]= 1 + (-1 + x) (1 + x)
```

The cubes expand to x^3

```mathematica
In[7]:= Expand[InterpolatingPolynomial[{1, 8, 27, 64}, x]]
Out[7]= x^3
```

Explicit {abscissa, value} pairs

```mathematica
In[8]:= InterpolatingPolynomial[{{0, 0}, {1, 1}, {2, 4}}, x]
Out[8]= x^2
```

Symbolic data gives a symbolic polynomial

```mathematica
In[9]:= InterpolatingPolynomial[{a, b, c}, x]
Out[9]= a + (-1 + x) (-a + b + 1/2 (a - 2 b + c) (-2 + x))
```

The same fit over Z/7Z

```mathematica
In[10]:= InterpolatingPolynomial[{1, 4, 9, 16, 25}, x, Modulus -> 7]
Out[10]= x^2
```

## Algorithm

interp.c

InterpolatingFunction --- piecewise-polynomial interpolation of tabulated data on a regular (tensor-product) grid, plus the Interpolation[] builder. Modelled on Mathematica's InterpolatingFunction object.

```text
  InterpolatingFunction[domain, table]
  InterpolatingFunction[domain, table, ders]
  InterpolatingFunction[domain, table, ders, orders]
  InterpolatingFunction[domain, table, ders, orders, method]

    domain = {{x1min, x1max}, ...}   -- one interval per dimension; the
             number of intervals m is the dimensionality.
    table  = {{coord, val}, ...}                     -- value-only data, or
             {{coord, val, grad, hess, ...}, ...}     -- derivative-supplied.
             coord is a scalar (1-D value-only) or an {x1,...,xm} list.
             grad = D[f,{vars,1}] (length-m vector), hess = D[f,{vars,2}]
             (m x m matrix), etc.
    ders   = {d1, ..., dm}   -- (optional) derivative-of-interpolant orders.
    orders = {o1, ..., om}   -- (optional) interpolation order per dimension.
    method = "Spline" | "Hermite"   -- (optional) interpolation method.
```

Methods (all evaluate the ders-th mixed derivative so D[ifun[..],..] composes):

```text
  default  : sliding-window Newton divided-difference (order min(3,n-1) or the
             requested InterpolationOrder), per dimension, tensor product.
  "Spline" : natural cubic spline (C2; second derivative 0 at the ends),
             tensor product over the full grid.
  "Hermite": tensor-product piecewise cubic Hermite with node slopes estimated
             by 3-point finite differences.
  supplied : derivative-annotated data is interpolated by tensor-product
             Hermite of per-dimension order k = max(K,1) where K is the highest
             supplied derivative order.  Mixed partials that are not supplied
             are filled by central finite differences across the grid.
```

Precision: machine (double) by default; if the data/argument carry MPFR arbitrary precision the MPFR kernels (interp_mpfr.c) are used instead and an EXPR_MPFR is returned.

Builtin ownership: interp_apply / the Interpolation builtin return a fresh Expr* (or NULL to stay unevaluated); inputs are borrowed.

## Implementation notes

**Algorithm.** `builtin_interpolatingpolynomial` builds the *single* polynomial
reproducing a set of values (unlike `InterpolatingFunction`, which is a piecewise
spline). Two engines share one condition model. **Engine A** — 1-D with every
value present — runs Newton confluent divided differences and emits the nested
Newton–Horner form `f1 + (x - x1)(f[1,2] + (x - x2)(...))`; the arithmetic is the
CAS heads (`ip_add`/`ip_mul`/`ip_sub`/`ip_div`), so exact or symbolic data give an
exact or symbolic polynomial, with a `double` fast path for inexact data.
**Engine B** — multivariate, or 1-D with an `Automatic` condition — solves for the
minimal-total-degree polynomial over a graded monomial basis by exact `Expr`
Gauss–Jordan, looping the degree `d` upward until `C(d + m, m) >= N` conditions
can be met and using *consistency* (not full column rank) as the success gate, so
the `noipf`/`poised` messages fire at the Mathematica-correct degree.

Abscissae default to `1, 2, ...`; the `{{xi, fi}, ...}` form gives explicit
points; a multidimensional point list with a variable list `{x, y, ...}` selects
Engine B; and trailing derivative entries `{{xi, fi, dfi, ...}, ...}` add Hermite
(derivative) conditions. `Modulus -> n` runs the same graded solve in **Z/nZ**
(`ip_solve_modular`), reducing each exact value and coordinate mod `n` and
inverting denominators there. A leading `NDArray` value tensor is materialised to
a nested list first.

**Data structures.** `IPData` is the shared condition model — owned coordinate
`Expr`s and `double`s, per-condition multi-index `alpha[]`, and value `Expr`s —
carrying flags for `any_automatic` (routes 1-D to Engine B) and
`coords_numeric`. Engine B uses a dense `Expr`/`mpz_t` augmented matrix over the
graded monomial basis.

**Complexity / limits.** Engine A is `O(n²)` divided differences; Engine B is
`O(N³)` Gauss–Jordan at the resolved degree. Points that are not *poised* (no
interpolant of the attempted total degree) raise `poised`/`noipf` and the call
declines. `Modulus` requires exact integer or rational data with invertible
denominators.

**Attributes:** `Protected`.

## References

**See also:** [Interpolation](../../functional-programming/Interpolation/), [InterpolatingFunction](../../functional-programming/InterpolatingFunction/), [N](../../arithmetic/N/), [NDArray](../../linear-algebra/NDArray/)

- J. Stoer and R. Bulirsch, *Introduction to Numerical Analysis*, 3rd ed. (Springer, 2002), §2.1 — Newton's divided differences and the Horner form.
- K. O. Geddes, S. R. Czapor and G. Labahn, *Algorithms for Computer Algebra* (Kluwer, 1992) — exact/symbolic interpolation over a ring.
- Source: [`src/interp.c`](https://github.com/stblake/mathilda/blob/main/src/interp.c)
- Specification: [`docs/spec/builtins/functional-programming.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/functional-programming.md)
- Tests: [`tests/test_interp_poly.c`](https://github.com/stblake/mathilda/blob/main/tests/test_interp_poly.c)

## Notes & additional examples

### Notes

`InterpolatingPolynomial[{f1, f2, ...}, x]` gives the one polynomial reproducing
the values `fi` at `x = 1, 2, ...`, with `n` values yielding degree `n - 1`. The
default output is the nested **Newton–Horner** form — compact and numerically
stable — so `{1, 4, 9}` prints as `1 + (-1 + x) (1 + x)` rather than `x^2`; wrap in
`Expand` to see the monomial form.

Abscissae can be given explicitly as `{{xi, fi}, ...}`, and the values may be
exact, inexact, or symbolic — exact data give an exact polynomial. The
`{{{x1, y1, ...}, f1}, ...}` form with a variable list builds a multidimensional
interpolant of lowest total degree, and `Modulus -> n` performs the fit over the
integers mod `n`.
