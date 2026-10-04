# NSolve

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`NSolve[expr, vars]`**

gives numerical approximations to the solutions of the equation or system expr for the variables vars, as a list of replacement-rule lists. NSolve\[expr, vars, Reals\] restricts to real solutions; the default domain is the complexes. vars may be a single variable or a list; NSolve\[{e1, e2, ...}, vars\] is the conjunction e1 && e2 && .... A working precision may be given as a trailing positional argument or via WorkingPrecision. Results: {} no solutions, {{x-\>s,...},...} the solutions (univariate roots are repeated by multiplicity), {{}} the universal solution. A univariate polynomial equation is solved with NRoots; square zero-dimensional polynomial systems use a Groebner-basis multiplication-matrix eigenvalue method (Method -\> "Symbolic" uses lexicographic elimination); other equations fall back to Solve or FindRoot seeding. Integer, real, and complex coefficients are handled at machine and arbitrary precision.

<details>
<summary>Notes</summary>

Options: MaxRoots, Method (Automatic | "EndomorphismMatrix" | "Homotopy" | "Symbolic"), WorkingPrecision, AccuracyGoal (default MachinePrecision, forwarded to NRoots), PrecisionGoal, VerifySolutions, RandomSeeding.

</details>

## Examples (11)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (5)

```mathematica
In[1]:= NSolve[x^5 - 2 x + 3 == 0, x, Reals]
Out[1]= {{x -> -1.42361}}

In[2]:= NSolve[{x^2 + y^2 == 1, x^3 - y^3 == 2}, {x, y}]
Out[2]= {{x -> -1.09791 + 0.839887*I, y -> 1.09791 + 0.839887*I}, {x -> -1.09791 - 0.839887*I, y -> 1.09791 - 0.839887*I}, {x -> 1.22333 - 0.0729987*I, y -> 0.125423 + 0.712005*I}, {x -> 1.22333 + 0.0729987*I, y -> 0.125423 - 0.712005*I}, {x -> -0.125423 + 0.712005*I, y -> -1.22333 - 0.0729987*I}, {x -> -0.125423 - 0.712005*I, y -> -1.22333 + 0.0729987*I}}

In[3]:= NSolve[{x^2 + y^3 == 1, 2 x + 3 y == 4}, {x, y}, Reals]
Out[3]= {{x -> 7.93641, y -> -3.95761}}

In[4]:= NSolve[x + 2 y + 3 z == 4 && 3 x + 4 y + 5 z == 6 && 6 x + 7 y + 8 z == 0, {x, y, z}]
Out[4]= {}

In[5]:= NSolve[E^x - x == 7, x, Reals]
Out[5]= {{x -> -6.99909}, {x -> 2.22154}}
```

### Options (1)

```mathematica
In[6]:= NSolve[{x^2 + y^2 == 1, x^3 - y^3 == 2}, {x, y}, WorkingPrecision -> 25]
Out[6]= {{x -> -1.0979116727228235764163996 + 0.83988692161565920362280281*I, y -> 1.0979116727228235764163996 + 0.83988692161565920362280281*I}, {x -> -1.0979116727228235764163996 - 0.83988692161565920362280281*I, y -> 1.0979116727228235764163996 - 0.83988692161565920362280281*I}, {x -> 1.2233348984131033766895813 - 0.072998738390442569855466144*I, y -> 0.12542322569027980027318178 + 0.71200452485314764855498901*I}, {x -> 1.2233348984131033766895813 + 0.072998738390442569855466144*I, y -> 0.12542322569027980027318178 - 0.71200452485314764855498901*I}, {x -> -0.12542322569027980027318178 + 0.71200452485314764855498901*I, y -> -1.2233348984131033766895813 - 0.072998738390442569855466144*I}, {x -> -0.12542322569027980027318178 - 0.71200452485314764855498901*I, y -> -1.2233348984131033766895813 + 0.072998738390442569855466144*I}}
```

### Worked examples (1)

```mathematica
In[7]:= NSolve[Sqrt[x] + 3 x^(1/3) == 5, x]
Out[7]= {{x -> 1.80863}}
```

### Applications (4)

A univariate polynomial, solved through NRoots

```mathematica
In[8]:= NSolve[x^2 - 2 == 0, x]
Out[8]= {{x -> -1.41421}, {x -> 1.41421}}
```

The Reals domain keeps only the real root

```mathematica
In[9]:= NSolve[x^3 - 1 == 0, x, Reals]
Out[9]= {{x -> 1.0}}
```

A linear system

```mathematica
In[10]:= NSolve[{x + y == 3, x - y == 1}, {x, y}]
Out[10]= {{x -> 2.0, y -> 1.0}}
```

A zero-dimensional nonlinear system

```mathematica
In[11]:= NSolve[{x^2 + y^2 == 1, y == x}, {x, y}]
Out[11]= {{x -> 0.707107, y -> 0.707107}, {x -> -0.707107, y -> -0.707107}}
```

## Algorithm

nsolve.c — NSolve[expr, vars, dom, prec, opts]

```text
Numerical equation solver.  NSolve returns approximate solutions of an
```

equation or system of equations as a list of replacement-rule lists:

```text
    {}                          no solutions
    {{x -> r1}, {x -> r2}, ...} one rule list per solution
    {{}}                        universal solution (every point satisfies)
```

Strategy (two specialists, matching the Wolfram Language's "Symbolic" idea):

```text
  1. Univariate polynomial equations  ->  NRoots.
     When the input reduces to a single polynomial equation lhs == rhs in a
     single variable, NSolve calls NRoots (the state-of-the-art Aberth /
     companion-matrix / Jenkins–Traub engine) and repackages its disjunction
     var==r1 || var==r2 || ...  as the rule-list form.  This covers integer,
     real, and complex coefficients, multiple roots (repeated by
     multiplicity), machine and arbitrary working precision, and the Reals
     domain (by discarding the complex roots).

  2. Everything else  ->  Solve, then numericalise.
     Linear systems, radical and inverse-function equations, etc. are solved
     symbolically by Solve and the exact result is rounded to the requested
     working precision.  This is the "Symbolic" method.  Inputs Solve cannot
     handle (e.g. genuine nonlinear polynomial systems) leave NSolve
     unevaluated.

Options:  MaxRoots, Method, WorkingPrecision, VerifySolutions, RandomSeeding,
          PrecisionGoal, MaxIterations.  (Method and the verification/seeding
          options are accepted for compatibility; the polynomial engine is
          always the NRoots default.)

Positional grammar:  NSolve[expr [, vars [, dom [, prec]]], opts...].
  dom  in {Reals, Complexes, Integers}; default Complexes.
  prec a number giving the working precision in decimal digits.
```

Memory contract (builtin): takes ownership of `res`; returns a fresh Expr* on success (the evaluator frees `res`) or NULL to leave NSolve unevaluated.

## Implementation notes

**Algorithm.** `builtin_nsolve` reads `NSolve[expr [, vars [, dom [, prec]]]]`,
peels options from the tail, defaults `vars` to the collected non-constant
symbols and `dom` to Complexes, and dispatches:

1. **Univariate polynomial** → `NRoots` directly (NRoots never frees its
   argument), forwarding `PrecisionGoal`/`AccuracyGoal` so its polishing and
   accuracy contract govern the roots; the disjunction is repackaged into
   `{{x -> r1}, …}` with the `Reals` filter and `MaxRoots` cap applied. A huge
   literal exponent is guarded (`NSolve::deg`) before any machinery allocates.
2. **Square zero-dimensional polynomial system** (`nsolve_system.c`) → the
   **eigenvalue / multiplication-matrix (Möller–Stetter) method**: a greVlex
   Gröbner basis (`gb_buchberger`) gives the quotient ring `A = Q[x]/I`; its
   standard-monomial basis is enumerated; rational multiplication matrices
   `M_{x_i}` are built by normal-form reduction; a generic linear form `M_l = Σ
   c_i M_{x_i}` (deterministic seeded coefficients) is formed, and the
   eigenvalues/eigenvectors of `M_l` at MPFR precision
   (`eigen_all_eigenvectors_real_mpfr`) give each coordinate as `x_i(p) =
   (M_{x_i} v)[j]/v[j]`. Every candidate is verified against the original
   residuals. `Method -> "Symbolic"` instead does lexicographic **elimination**
   (solve the univariate generator with `NRoots`, back-substitute, recurse,
   verify); `"Homotopy"` currently routes to the same eigenvalue engine.
3. **Fallback** → symbolic `Solve` then numericalisation, dropping provably
   extraneous roots; a univariate non-polynomial last resort seeds `FindRoot`
   from a real grid (plus `±2i` unless `Reals`), verified and deduplicated.

Results are a list of rule-lists: `{}` no solutions, `{{}}` the universal
solution. Both the univariate path and both system solvers call `builtin_nroots`
directly, so the NRoots engines (Aberth / companion / Jenkins–Traub) are
NSolve's numeric backbone.

**Data structures.** The Gröbner engine works over `Q` (`GBPoly`); the
multiplication matrices are `mpq_t` rationals, the linear-form matrix and its
per-variable companions `mpfr_t`, and the eigen buffers and recovered
coordinates MPFR/`ncpx`. The standard-monomial basis is a flat `int[d·nvar]`.
`want_machine` holds when no precision digit count is given; otherwise the system
runs at `target_bits + max(32, target_bits/2)` bits and emits MPFR values.

**Complexity / limits.** The eigenproblem is `O(d³)` in the quotient-ring
dimension `d`; hard caps `NSYS_MAX_DIM = 256`, `NSYS_MAX_BOX = 200000`, and a
per-generator total-degree gate `NSYS_MAX_TDEG = 60` make a too-large or
positive-dimensional system fall back / stay unevaluated. The univariate degree
guard is `NSOLVE_MAX_POLY_DEGREE = 10000`. Options: `MaxRoots`, `Method`
(`Automatic` | `"EndomorphismMatrix"` | `"Homotopy"` | `"Symbolic"`),
`WorkingPrecision` (also a trailing positional digit count), `AccuracyGoal`
(default `MachinePrecision`, forwarded to NRoots), `PrecisionGoal`,
`VerifySolutions` (default on), `RandomSeeding` (seed for the generic linear
form, default 1234). `NSolve[expr, vars, Reals]` filters to real values; the
`Integers` domain is left to `Solve`. Diagnostics route through `mth_message`.

**Attributes:** `Protected`.

## References

**See also:** [NRoots](../../numerical-calculus/NRoots/), [Solve](../../solutions-of-equations/Solve/), [VerifySolutions](../../solutions-of-equations/VerifySolutions/), [ConditionalExpression](../../control-flow/ConditionalExpression/), [AccuracyGoal](../../other-advanced/AccuracyGoal/), [PrecisionGoal](../../other-advanced/PrecisionGoal/), [FindRoot](../../calculus/FindRoot/), [Exists](../../solutions-of-equations/Exists/)

- H. M. Möller and H. J. Stetter, *Multivariate polynomial equations with multiple zeros solved by matrix eigenproblems*, Numer. Math. **70** (1995) 311–329.
- D. A. Cox, J. Little and D. O'Shea, *Using Algebraic Geometry*, 2nd ed. (Springer, 2005), ch. 2 — the eigenvalue method for zero-dimensional ideals.
- Source: [`src/numerical_roots/nsolve.c`](https://github.com/stblake/mathilda/blob/main/src/numerical_roots/nsolve.c)
- Specification: [`docs/spec/builtins/numerical-calculus.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/numerical-calculus.md)
- Tests: [`tests/test_nsolve.c`](https://github.com/stblake/mathilda/blob/main/tests/test_nsolve.c)
- Tests: [`tests/test_nsolve_stress.c`](https://github.com/stblake/mathilda/blob/main/tests/test_nsolve_stress.c)

## Notes & additional examples

### Notes

`NSolve[expr, vars]` returns numerical solutions as a list of replacement-rule
lists; `{}` means no solutions and `{{}}` the universal solution. `vars` may be a
single variable or a list, and `NSolve[{e1, e2, ...}, vars]` is the conjunction
`e1 && e2 && ...`. `NSolve[expr, vars, Reals]` restricts to real solutions; the
default domain is the complexes.

A univariate polynomial equation is handed to `NRoots` (roots repeated by
multiplicity), so `NRoots`'s `PrecisionGoal`/`AccuracyGoal` contract governs the
answer. A square, zero-dimensional polynomial **system** uses a Groebner-basis
multiplication-matrix (Moeller-Stetter) eigenvalue method: a Groebner basis gives
the quotient ring, and the eigenvalues/eigenvectors of the multiplication maps
yield each coordinate, every candidate verified against the original residuals.
`Method -> "Symbolic"` instead uses lexicographic elimination (solve the
univariate generator, back-substitute, recurse). Other equations fall back to
symbolic `Solve` then numericalisation, with a univariate `FindRoot` grid as a
last resort.

A working precision may be given as a trailing positional argument or via
`WorkingPrecision`; integer, real, and complex coefficients are handled at
machine and arbitrary precision. `MaxRoots`, `VerifySolutions`, and
`RandomSeeding` (the seed for the generic linear form) are also accepted.
