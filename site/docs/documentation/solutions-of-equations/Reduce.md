# Reduce

!!! warning "Status: Partial"
    implemented with documented limitations or caveats; some argument forms fall through to symbolic/unevaluated output.

## Description

**`Reduce[expr, vars]`**

Reduces the statement expr -- a logical combination (&&, ||, !, Implies, Xor) of equations (==, !=) and inequalities (\<, \<=, \>, \>=) -- to a complete, quantifier-free description of its solution set for the variables vars.  The default domain is Complexes, or Reals when expr contains an ordering inequality (ordering is undefined over the complexes), so e.g. Reduce\[-5 \< 3x+7 \<= 22, x\] is solved over the Reals.

**`Reduce[expr]`**

Infers the variables from expr (all free symbols, excluding named constants such as Pi), so Reduce\[x^2 \< 1\] is Reduce\[x^2 \< 1, {x}\].

**`Reduce[expr, dom]`**

A domain symbol (Reals, Complexes, Integers, Rationals) in the second slot infers the variables and reduces over dom, so Reduce\[x^2 \< 1, Reals\] is Reduce\[x^2 \< 1, {x}, Reals\].

**`Reduce[expr, vars, dom]`**

Reduces over the domain dom: Complexes, Reals, Integers, or Rationals.

<details>
<summary>Notes</summary>

A list {e1, e2, ...} in the expr slot is taken as the conjunction e1 && e2 && ... , so Reduce\[{x + y == 3, x - y == 1}, {x, y}\] is the same as Reduce\[x + y == 3 && x - y == 1, {x, y}\]. Where Solve returns the generic solution of a set of equations as a list of replacement rules and silently drops the degenerate cases, Reduce returns an And/Or tree of equations and inequalities that describes the WHOLE solution set -- every parametric and boundary case kept -- and it solves inequalities over the reals.  A statement that decides returns True or False; an out-of-reach input is left unevaluated, never guessed. Handled so far: - Complexes: univariate polynomial equations carrying the full leading-coefficient case tree -- Reduce\[a x == b, x\] -\> (a != 0 && x == b/a) || (a == 0 && b == 0) -- parametric linear systems by symbolic Gaussian elimination, and zero-dimensional nonlinear systems solved exactly (finitely many points). - Reals, one variable: any Boolean combination of polynomial and rational-function equations and inequalities, solved as a union of intervals and points on an exact real-algebraic sign diagram (denominator roots are breakpoints and are excluded as poles); plus statements built from Abs, real radicals, Log, bounded inverse-trig, Floor/Ceiling/Round/Mod/IntegerPart, Min/Max, and the piecewise heads (Piecewise, Sign, UnitStep, Ramp, Clip, HeavisideTheta, Boole, ...). - Reals, several variables: linear systems by Fourier-Motzkin elimination and nonlinear systems (conics and beyond) by Cylindrical Algebraic Decomposition (McCallum projection), including multivariate Abs/Min/Max/piecewise selectors and square-root radical rationalization.  A zero-dimensional system (equations pinning finitely many points, e.g. three circles with sign constraints) is solved exactly and its branches filtered by the inequalities via the algebraic-number oracle. - Integers / Rationals: the Solve Diophantine engine, reformatted as an Or of Ands with Element\[C\[k\], dom\] for a free parameter. Options (Options\[Reduce\]):

```text
Backsubstitution -> False     the linear-system output is fully
                            solved (grafted); accepted/echoed.
Cubics -> False               emit radicals (True) or Root[]
Quartics -> False             (False) for irreducible cubic /
                            quartic equations.
GeneratedParameters -> C      head of the free parameters C[k]
                            for parametric Integers/Rationals
                            solutions.
Method -> Automatic           reserved (Automatic is the only
                            method).
Modulus -> 0                  a nonzero p solves the equations
                            over Z/pZ by residue enumeration,
                            overriding the domain.
```

WorkingPrecision -\> Infinity  numeric-fallback tolerance for transcendental sign decisions; Infinity keeps the exact-first path. Reduce is sound over complete: an undecidable sign, an unsupported construct, or a positive-dimensional nonlinear system over Complexes (no CAD there yet) leaves the input unevaluated rather than risk a wrong formula.

</details>

## Examples (19)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (7)

```mathematica
In[1]:= Reduce[a x == b, x]
Out[1]= a != 0 && x == b/a || a == 0 && b == 0

In[2]:= Reduce[x^2 < 1]
Out[2]= -1 < x < 1

In[3]:= Reduce[1/x < 1, x, Reals]
Out[3]= x < 0 || x > 1

In[4]:= Reduce[Sin[x] == 1/2 && 0 < x < 2 Pi, x]
Out[4]= x == 1/6 Pi || x == 5/6 Pi

In[5]:= Reduce[x y z > 0, {x, y, z}, Reals]
Out[5]= x < 0 && (y < 0 && z > 0 || y > 0 && z < 0) || x > 0 && (y < 0 && z < 0 || y > 0 && z > 0)

In[6]:= Reduce[x + y + z == 1 && x > 0 && y > 0 && z > 0 && x y + x z + y z < 9 x y z]
Out[6]= False

In[7]:= Reduce[2 x + 3 y == 1, {x, y}, Integers]
Out[7]= Element[C[1], Integers] && x == -1 + 3 C[1] && y == 1 - 2 C[1]
```

### Options (1)

```mathematica
In[8]:= Reduce[x^2 == 2, x, Modulus -> 7]
Out[8]= x == 3 || x == 4
```

### Worked examples (1)

```mathematica
In[9]:= Reduce[x^2 < 1 && y^2 < 1 && z^2 < 1, {x, y, z}, Reals]
Out[9]= -1 < x < 1 && -1 < y < 1 && -1 < z < 1
```

### Applications (10)

The degenerate branch a == 0 is kept, unlike Solve

```mathematica
In[10]:= Reduce[a x == b, x]
Out[10]= a != 0 && x == b/a || a == 0 && b == 0
```

Variables inferred; the ordering inequality selects the Reals

```mathematica
In[11]:= Reduce[x^2 < 1]
Out[11]= -1 < x < 1
```

The pole at 0 is an excluded breakpoint of the sign diagram

```mathematica
In[12]:= Reduce[1/x < 1, x, Reals]
Out[12]= x < 0 || x > 1
```

Periodic families filtered against a bounding region

```mathematica
In[13]:= Reduce[Sin[x] == 1/2 && 0 < x < 2 Pi, x]
Out[13]= x == 1/6 Pi || x == 5/6 Pi
```

Two-variable CAD with algebraic breakpoints

```mathematica
In[14]:= Reduce[x^2 + y^2 < 1 && y > x^2, {x, y}, Reals]
Out[14]= -Sqrt[1/2 (-1 + Sqrt[5])] < x < Sqrt[1/2 (-1 + Sqrt[5])] && x^2 < y < 1/2 Sqrt[4 - 4 x^2]
```

CAD in three variables, split by the sign of each coordinate

```mathematica
In[15]:= Reduce[x y z > 0, {x, y, z}, Reals]
Out[15]= x < 0 && (y < 0 && z > 0 || y > 0 && z < 0) || x > 0 && (y < 0 && z < 0 || y > 0 && z > 0)
```

Exact sum-of-squares certificate proves the region empty

```mathematica
In[16]:= Reduce[x + y + z == 1 && x > 0 && y > 0 && z > 0 && x y + x z + y z < 9 x y z]
Out[16]= False
```

An embedded quantifier is eliminated

```mathematica
In[17]:= Reduce[Exists[y, x^2 + y^2 < 1], {x}, Reals]
Out[17]= -1 < x < 1
```

Linear Diophantine family with an integer parameter

```mathematica
In[18]:= Reduce[2 x + 3 y == 1, {x, y}, Integers]
Out[18]= Element[C[1], Integers] && x == -1 + 3 C[1] && y == 1 - 2 C[1]
```

Honest decline: a positive-dimensional system over the complexes stays unevaluated

```mathematica
In[19]:= Reduce[x^2 + y^2 == 1, {x, y}]
Out[19]= Reduce[x^2 + y^2 == 1, {x, y}]
```

## Algorithm

reduce.c

```text
`Reduce` -- front-end and dispatch skeleton (docs/notes/REDUCE_PLAN.md, Phase 0).
```

Mirrors builtin_solve's argument handling: positional `expr [, vars [, dom]]`, a Solve::ivar-style bad-variable diagnostic, and a True/False short-circuit. The input logical combination is normalised into the internal DNF layer (reduce_form.h); Phase 0 returns True/False when the statement decides and

```text
leaves everything else unevaluated (NULL).  The per-domain solving engines
```

(equational, sign-diagram, Fourier-Motzkin, CAD, integer, quantifier) are wired on top of this skeleton in later phases.

## Implementation notes

**Algorithm.** `builtin_reduce` wraps `reduce_impl` in an internal-function
message mute (`mth_msg_ifun_suppress`), so the speculative probes the engine
runs internally do not leak diagnostics. `reduce_impl` peels trailing option
`Rule`s, validates the variable spec, and rewrites a `List` in the `expr` slot
into an `And`. In the one-argument form, and when the second slot holds a domain
symbol (`reduce_is_domain_symbol`: `Reals`, `Complexes`, `Integers`,
`Rationals`), `reduce_infer_var_list` takes the variables to be every free
symbol, minus named constants, domain symbols and reserved literals, sorted
alphabetically. So `Reduce[x^2 < 1]` and `Reduce[x^2 < 1, Reals]` reach the
same engine as `Reduce[x^2 < 1, {x}, Reals]`. Two heads short-circuit the pipeline: a top-level `Exists`/`ForAll`
routes to `reduce_qe_dispatch` (quantifier elimination), and a `Modulus -> p`
option routes to `reduce_modular` (residue enumeration over `Z/pZ`). Otherwise the
statement is normalised to a DNF `RForm` by `reduce_form_from_expr` and simplified
(`rform_simplify`); the domain is chosen (default `Complexes`, but `Reals` when an
ordering `<`/`<=` is present with no explicit domain, or when a univariate
real-function / multivariate piecewise/radical statement is force-routed and
pre-rewritten by `reduce_realfn_preprocess` / `reduce_piecewise_preprocess`), and
the per-domain dispatch runs.

The dispatch is: **Integers/Rationals** → `reduce_integers`, which reuses the
`Solve` Diophantine engine and reformats its answer as an `Or` of `And`s carrying
`Element[C[k], dom]` for each free parameter. **Complexes** (equations only) → a
single univariate polynomial equation goes to `reduce_eq_univariate` (which carries
the full leading-coefficient case tree) or, when the residual is transcendental, to
`reduce_eq_transcendental` via `Solve`; a linear system to `reduce_eq_system`; and
anything else (nonlinear, or equations plus `!=`) to the zero-dimensional engine
`reduce_zerodim`. **Reals, one variable** → the exact real-algebraic sign-diagram
engines `reduce_univar` / `reduce_univar_general`, with trig/hyperbolic/exponential
equations pre-empted through `Solve` (the sign diagram is unsound on them) and
periodic-equation-plus-bounds shapes handled by `reduce_periodic_region` /
`reduce_trig_ineq_region`. **Reals, several variables** → `reduce_fm`
(Fourier-Motzkin for linear systems), then `reduce_sos`, then `reduce_cad`
(Cylindrical Algebraic Decomposition, McCallum projection), and finally
`reduce_zerodim` for a zero-dimensional system with irrational fibres. A
rational atom `p/q REL 0` with a variable denominator is cleared first to its
sign-exact polynomial form (`p q < 0`, plus `q != 0` guards), in
`reduce_piecewise_preprocess`.

`reduce_sos` (`reduce_sos.c`, v0.322) is a Positivstellensatz emptiness
certificate. For a region `{g_i >= 0, h_j == 0, q < 0}` it looks for a Putinar
representation `q = σ₀ + Σ σ_i g_i` modulo the equality ideal, with each
`σ` a sum of squares. The Gram matrices come from a numeric SDP solved by
alternating projection, through LAPACK. The interior-zero case uses facial
reduction about a rational zero found by minimisation. The matrices are then
rounded to rationals, and the identity and positive semidefiniteness are
re-verified in exact `mpq` arithmetic by a rational LDLᵀ. Only a verified
certificate returns `False`, and any other shape declines quickly. It runs
ahead of the CAD because the high-degree regions it targets are exactly the
ones on which the CAD grinds.

The n-variable CAD wraps `cad_build` in `cad_build_augmented`, which
quantifier elimination shares. When a fibre polynomial nullifies, that
factor's coefficients are added to the projection level below and the stack is
rebuilt, which repairs McCallum well-orientedness. The rounds are bounded, and
a nullification that persists declines. Every CAD loop head polls `CAD_POLL()`,
one clock read per 64 iterations and only inside a `TimeConstrained` scope.
When the deadline has passed the engine returns `NULL` through its ordinary
decline path, freeing every local, so `TimeConstrained[Reduce[...]]` yields
promptly without leaking.

**Data structures.** The working form is an `RForm`: a DNF of `RConj`
conjunctions of `RAtom`s, each a polynomial (an `Expr`) paired with a relation code
(`R_EQ`, `R_LT`, `R_LE`, `R_ELEM`). The CAD (`reduce_cad.c`) keeps per-level
`PolySet` projection stacks and a `CADRegion` cell tree, and delegates every
numeric step — `Discriminant`, `Resultant`, `FactorList`, `Solve[..,Reals]`
real-root isolation, and the exact `qqbar` sign oracle — to the evaluator and
`reduce_real_util`.

**Complexity / limits.** CAD cost is doubly exponential in the number of
variables (McCallum projection then lifting). The 2-variable case is fully
realised, and the n-variable path is an iterated projection (`cad_project_out` for
each level) with recursive lifting and partial-CAD pruning. The SOS stage needs
LAPACK and declines cleanly without it. The governing invariant is **soundness over completeness**: an
undecidable sign or ordering (`qqbar` returning `-2`), an unwired case, or an
irrational section all return `NULL`, leaving the input unevaluated rather than
emitting a wrong or partial formula. Consequently a returned formula describes the
entire solution set exactly, and a literal `False` is a genuine proof of emptiness.

**Attributes:** `Protected`.

## References

**See also:** [Unequal](../../comparisons/Unequal/), [Pi](../../mathematical-constants/Pi/), [E](../../mathematical-constants/E/), [Degree](../../mathematical-constants/Degree/), [Implies](../../control-flow/Implies/), [Xor](../../control-flow/Xor/), [Inequality](../../comparisons/Inequality/), [Solve](../../solutions-of-equations/Solve/)

- Collins & Hong, "Partial Cylindrical Algebraic Decomposition for Quantifier Elimination", J. Symbolic Computation 12 (1991).
- McCallum, "An Improved Projection Operation for Cylindrical Algebraic Decomposition", in Quantifier Elimination and Cylindrical Algebraic Decomposition (1998).
- Basu, Pollack & Roy, "Algorithms in Real Algebraic Geometry" (2nd ed., 2006), Ch. 5 & 11 (sign determination, CAD).
- Schrijver, "Theory of Linear and Integer Programming" (1986), §12.2 (Fourier-Motzkin elimination).
- G. E. Collins and H. Hong, *Partial Cylindrical Algebraic Decomposition for Quantifier Elimination*, J. Symbolic Computation **12** (1991) 299-328.
- S. McCallum, *An Improved Projection Operation for Cylindrical Algebraic Decomposition*, in *Quantifier Elimination and Cylindrical Algebraic Decomposition* (Springer, 1998) 242-268.
- S. Basu, R. Pollack and M.-F. Roy, *Algorithms in Real Algebraic Geometry*, 2nd ed. (Springer, 2006), Ch. 5 & 11 (sign determination, CAD).
- A. Schrijver, *Theory of Linear and Integer Programming* (Wiley, 1986), §12.2 (Fourier-Motzkin elimination).
- Source: [`src/solve/reduce.c`](https://github.com/stblake/mathilda/blob/main/src/solve/reduce.c)
- Specification: [`docs/spec/builtins/solutions-of-equations.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/solutions-of-equations.md)
- Tests: [`tests/test_minimize.c`](https://github.com/stblake/mathilda/blob/main/tests/test_minimize.c)
- Tests: [`tests/test_reduce.c`](https://github.com/stblake/mathilda/blob/main/tests/test_reduce.c)
- Tests: [`tests/test_reduce_corpus.c`](https://github.com/stblake/mathilda/blob/main/tests/test_reduce_corpus.c)

## Notes & additional examples

### Notes

- **Reduce and Solve answer different questions.** `Solve[eqns, vars]`
  returns generic solutions as replacement rules and drops the special
  cases. `Reduce[expr, vars]` returns a quantifier-free `And`/`Or` formula
  covering the whole solution set, every parametric and boundary branch
  included. Use `Reduce` when you need the exact conditions, especially
  with symbolic parameters or with inequalities.
- **Domain and variables.** The default domain is `Complexes`. A statement
  with an ordering inequality, or with a real-only function such as `Abs`,
  `Floor` or `Max`, defaults to `Reals`. `Reduce[expr]` and
  `Reduce[expr, dom]` infer the variables, taking every free symbol except
  named constants, in alphabetical order. A third argument (`Complexes`,
  `Reals`, `Integers` or `Rationals`) chooses the domain explicitly.
- **Engines over the reals.** One variable goes to an exact
  real-algebraic sign diagram, which handles poles, radicals, `Abs`,
  `Floor`, piecewise heads and bounded trig regions. Several variables go
  to Fourier–Motzkin for linear systems. Higher-degree regions go to an
  exact SOS/Positivstellensatz emptiness certificate, then an n-variable
  McCallum CAD, then a zero-dimensional solver. Rational atoms in several
  variables are cleared to polynomial form first.
- **Options.** `Cubics`, `Quartics`, `Modulus -> p`,
  `GeneratedParameters`, `Backsubstitution`, `WorkingPrecision` and
  `Method` are accepted. `TimeConstrained[Reduce[...], t]` interrupts the
  CAD cleanly, without leaking.
- **Soundness over completeness.** The engine is designed to return either
  the exact solution set or nothing. An undecidable sign, an unsupported
  construct, or an engine not yet wired (positive-dimensional nonlinear
  systems over the complexes) leaves the call unevaluated.
- `Protected`. Symbolic, with no NDArray or `Compile[]` path.
- See also: `Solve`, `Resolve`, `FindInstance`, `CylindricalDecomposition`,
  `Refine`.
