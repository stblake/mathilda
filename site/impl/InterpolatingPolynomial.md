---
references:
  - "J. Stoer and R. Bulirsch, *Introduction to Numerical Analysis*, 3rd ed. (Springer, 2002), §2.1 — Newton's divided differences and the Horner form."
  - "K. O. Geddes, S. R. Czapor and G. Labahn, *Algorithms for Computer Algebra* (Kluwer, 1992) — exact/symbolic interpolation over a ring."
source: src/interp.c
---
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
