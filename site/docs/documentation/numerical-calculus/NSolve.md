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

## Examples (6)

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

**Attributes:** `Protected`.

## See also

[NRoots](../../numerical-calculus/NRoots/), [Solve](../../solutions-of-equations/Solve/), [VerifySolutions](../../solutions-of-equations/VerifySolutions/), [ConditionalExpression](../../control-flow/ConditionalExpression/), [AccuracyGoal](../../other-advanced/AccuracyGoal/), [PrecisionGoal](../../other-advanced/PrecisionGoal/), [FindRoot](../../calculus/FindRoot/)

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/numerical-calculus.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/numerical-calculus.md)
- Tests: [`tests/test_nsolve.c`](https://github.com/stblake/mathilda/blob/main/tests/test_nsolve.c)
- Tests: [`tests/test_nsolve_stress.c`](https://github.com/stblake/mathilda/blob/main/tests/test_nsolve_stress.c)
