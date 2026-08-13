# Solve

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Solve[expr, vars]`**

Attempts to solve the equation or system expr for the variables vars.

**`Solve[expr, vars, dom]`**

Solves over the domain dom.  Default Complexes; Reals filters down to real roots via per-degree discriminant and sign tests; Integers further restricts the output to provably concrete integer solutions (Integer / BigInt only -- Rationals, Sqrt\[\], and held Root\[\] objects are dropped).

**`ConditionalExpression[..., Element[C[k], Integers]].  Emits`**

<details>
<summary>Notes</summary>

Options: Cubics              -\> False     (radical form for cubics) Quartics            -\> False     (radical form for quartics) InverseFunctions    -\> Automatic (use inverse-function peel) GeneratedParameters -\> C         (head for parameters C\[k\]) VerifySolutions     -\> Automatic (reserved) Solves single polynomial equalities, radical equations, linear systems, zero-dimensional nonlinear polynomial systems (via a lexicographic Groebner basis and triangular back-substitution; positive-dimensional systems emit Solve::nsdim and stay unevaluated), and -- via the inverse-function specialist -- single-variable equations whose outermost dependence is an elementary invertible head (Log, Exp, Sin/Cos/Tan/Cot/Sec/Csc, their hyperbolic counterparts, the inverse trig/hyperbolic forms, and Power\[g, n\] for integer n \>= 2).  Multi-branch heads introduce an integer parameter C\[k\] wrapped in Solve::ifun the first time inverse functions are used.

</details>

## Examples (18)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (8)

```mathematica
In[1]:= Solve[2 x + 3 == 0, x]
Out[1]= {{x -> -3/2}}

In[2]:= Solve[x^2 - 5 x + 6 == 0, x]
Out[2]= {{x -> 2}, {x -> 3}}

In[3]:= Solve[x^2 + 1 == 0, x]
Out[3]= {{x -> -I}, {x -> I}}

In[4]:= Solve[x^2 + 1 == 0, x, Reals]
Out[4]= {}

In[5]:= Solve[(x-1)^2 == 0, x]
Out[5]= {{x -> 1}, {x -> 1}}

In[6]:= Solve[x^4 - 5 x^2 + 4 == 0, x]
Out[6]= {{x -> -2}, {x -> -1}, {x -> 1}, {x -> 2}}

In[7]:= Solve[x^3 + x + 1 == 0, x]
Out[7]= {{x -> Root[1 + #1 + #1^3 &, 1]}, {x -> Root[1 + #1 + #1^3 &, 2]}, {x -> Root[1 + #1 + #1^3 &, 3]}}

In[8]:= Solve[Sin[x] == 0, x]
Out[8]= {{x -> ConditionalExpression[Pi + 2 C[1] Pi, Element[C[1], Integers]]}, {x -> ConditionalExpression[2 C[1] Pi, Element[C[1], Integers]]}}
```

### Worked examples (1)

```mathematica
In[9]:= Solve[x^3 - 6 x^2 + 11 x - 6 == 0, x, Integers]
Out[9]= {{x -> 1}, {x -> 2}, {x -> 3}}
```

### Applications (9)

```mathematica
In[10]:= Solve[x^2 - 5 x + 6 == 0, x]
Out[10]= {{x -> 2}, {x -> 3}}

In[11]:= Solve[x^2 + 1 == 0, x]
Out[11]= {{x -> -I}, {x -> I}}

In[12]:= Solve[x^2 - 2 == 0, x]
Out[12]= {{x -> -Sqrt[2]}, {x -> Sqrt[2]}}

In[13]:= Solve[{x + y == 3, x - y == 1}, {x, y}]
Out[13]= {{x -> 2, y -> 1}}

In[14]:= Solve[x y == 1 && x + y == 3, {x, y}]
Out[14]= {{x -> 1/2 (3 - Sqrt[5]), y -> 1/2 (3 + Sqrt[5])}, {x -> 1/2 (3 + Sqrt[5]), y -> 1/2 (3 - Sqrt[5])}}

In[15]:= Solve[x y == 6 && x + y == 5, {x, y}, Integers]
Out[15]= {{x -> 3, y -> 2}, {x -> 2, y -> 3}}

In[16]:= Solve[a x^2 + b x + c == 0, x]
Out[16]= {{x -> (1/2 (-b + Sqrt[b^2 - 4 a c]))/a}, {x -> (1/2 (-b - Sqrt[b^2 - 4 a c]))/a}}

In[17]:= Solve[x^4 - 1 == 0, x]
Out[17]= {{x -> -1}, {x -> 1}, {x -> -I}, {x -> I}}

In[18]:= Solve[Sin[x] == 0, x]
Out[18]= {{x -> ConditionalExpression[Pi + 2 C[1] Pi, Element[C[1], Integers]]}, {x -> ConditionalExpression[2 C[1] Pi, Element[C[1], Integers]]}}
```

## Options & behaviour

### Options

- `Cubics -> False`: Emit cubic roots as held `Root[]` objects (default).
  `Cubics -> True` switches to closed-form Cardano radicals.
- `Quartics -> False`: Emit quartic roots as held `Root[]` objects (default).
  `Quartics -> True` switches to closed-form Ferrari radicals (Complexes only).
- `InverseFunctions -> Automatic`: Enables the inverse-function specialist
  (default).  Set to `False` to disable the specialist; equations that can
  only be solved through inversion then return unevaluated.
- `GeneratedParameters -> C`: Head used by the inverse-function specialist
  when minting fresh integer-parameter symbols `C[1], C[2], ...`.  Only the
  bare-symbol form is honoured; the `Function` form is reserved.
- `VerifySolutions -> Automatic`: Reserved.

## Algorithm

solve.c

The `Solve` router: classifies the input equation system, parses

```text
options, and dispatches to a specialist solver.  The only specialist
```

wired up in this initial cut is Solve`SolvePolynomialEquality (src/solvepoly.c) for a single polynomial equality in one variable.

```text
`Solve` does not hold its arguments -- the evaluator delivers
`expr` and `vars` already evaluated, matching Mathematica's
attribute set ({Protected}).  When `vars` has been substituted to
```

a non-symbol (typically because the user previously assigned

```text
`x = 5` and then called `Solve[..., x]`), the router emits
`Solve::ivar` and returns unevaluated.
```

## Performance

Against other systems, from the benchmark suite (same input, results cross-checked for agreement):

| case | Mathilda | Wolfram | Python |
|---|---:|---:|---:|
| NSolve degree 40 | 9.27 s | 0.779 s | 341 s |
| Discriminant of deg 20 | 2.51 s | 0.068 s | 0.182 s |
| Solve factorable sextic | 0.469 s | 0.747 s | 1.52 s |
| Solve sextic, root sum | 0.439 s | 0.755 s | 1.49 s |
| Expand (1+x)^400 | 0.434 s | 0.107 s | 0.003 s |
| Cancel deg-60 over deg-58 | 0.337 s | 0.569 s | 7.37 s |

## Implementation notes

**Algorithm.** `builtin_solve` is a classifier/router, not a solver: it parses options, validates the variable spec, normalises the input, then dispatches to one of five specialists in `src/solvepoly.c`, `src/solverad.c`, `src/solvelinsys.c`, `src/solvetrig.c`, `src/solveinv.c`.

The router first peels trailing `Rule`/`RuleDelayed` options (`Cubics`, `Quartics`, `InverseFunctions`, `GeneratedParameters`, `VerifySolutions`, `Assumptions`, `Method`, `Modulus`) off the end of the positional args via `is_known_option_name`; an unrecognised trailing option name emits `Solve::optx`. Positional args are `expr [, vars [, dom]]`. `is_valid_solve_vars` rejects numeric-literal variables with `Solve::ivar`. Compound variables (`Dt[y]`, `f[a,b]`, `x^2`) are rewritten to fresh internal symbols `Solve$var$N` throughout `expr` by `collect_and_subst_compound_vars`, then restored in the output by `unsubst_compound_vars`, so the specialists only ever see bare symbols. Inexact-coefficient inputs are detected (`common_scan_inexact`), force-rationalised to the minimum bit precision found (`common_rationalize_input`), solved exactly, and numericalised back at the tail (`common_numericalize_result`) — exact-in/exact-out, inexact-in/inexact-out, mirroring Cancel/Together/Integrate. `True`/`False` short-circuit to `{{}}` (tautology) / `{}` (contradiction). `Abs[u]==0` is rewritten to `u==0` (`try_abs_zero_rewrite`).

The dispatch cascade: an `And`/`List` of equations, or a single `Equal` over a ≥2-symbol variable list, routes to `solvelinsys_solve_linear_system` (linear-system specialist, canonicalises each equation to `lhs - rhs`, returns NULL when non-affine). Otherwise the single-variable path tries `solvepoly_solve_polynomial_equality` first (the polynomial specialist, `src/poly/solvepoly.c`, also exposed as `Solve\`SolvePolynomialEquality`); on NULL it falls back in order to `solveinv_solve_inverse_equality` (peels one elementary invertible head — Log/Exp/trig/hyperbolic/inverse-trig and integer Power — introducing `C[k]` integer parameters wrapped in `ConditionalExpression`), then `solvetrig_solve_trig_equality` (multi-trig canonicalisation), then `solverad_solve_radicals_equality` (radical equations). A specialist returning NULL leaves the call unevaluated.

The `dom` third argument selects the solution domain: default `Complexes`; `Reals` filters via per-degree discriminant/sign tests inside the polynomial specialist; `Integers` further drops non-concrete-integer solutions. `Cubics`/`Quartics -> False` (the defaults) return cubic/quartic roots as held `Root[]` objects rather than radical formulas.

**Data structures.** Everything is `Expr*`. Options accumulate into a stack `SolveOpts` bundling `SolvePolyOpts` and `SolveInvOpts`. Compound-var substitutions are tracked in a fixed `SolveVarSub subs[32]` array (`SOLVE_MAX_VAR_SUBS`). Output is the standard `List` of `List` of `Rule` rewrite-rule form `{{x -> ...}, ...}`.

**Complexity / limits.** Dominated by the chosen specialist (polynomial root-finding, linear-system elimination, radical isolation). The router itself is linear in expression size plus the substitution passes. The compound-variable cap is 32 distinct variables per call.

- `Protected`.  Uses the standard attribute set -- arguments are
  evaluated by the evaluator before reaching the router.  When the
  second argument has been substituted to a numeric atom (typically
  because an OwnValue like `x = 5` was previously set, or the user
  literally passed `Solve[..., 5]`), the router emits `Solve::ivar`
  and returns unevaluated.
- **Generalised (compound) variables.**  `vars` may contain any
  non-numeric expression, not only symbols: `Solve[lhs == rhs, Dt[y]]`,
  `Solve[f[a] + b == c, f[a]]`, `Solve[a x^2 + b == 0, x^2]`, and
  multi-var forms like `Solve[{...}, {Dt[x], Dt[y]}]` are all
  accepted.  The router substitutes each non-symbol entry with a
  fresh internal symbol (`Solve$var$N`), runs the standard dispatch,
  and reverses the substitution on the result so the user sees Rule
  LHSes like `Dt[y] -> ...` directly.  The substitution is purely
  structural (literal `expr_eq`); polynomial identifications like
  `x^4 == (x^2)^2` are not yet recognised, so `Solve[x^4 - 1 == 0,
  x^2]` returns the substituted form rather than `{{x^2 -> 1},
  {x^2 -> -1}}`.
- Acts as a router that classifies its input and dispatches to a specialist:
  - Single equality, single variable -> `Solve`SolvePolynomialEquality` (below).
  - Single equality, single variable, polynomial specialist declines because
    the outermost dependence on `var` is an elementary invertible head ->
    inverse-function specialist (`src/solveinv.c`): peels `Log`, `Exp`,
    `Sin`/`Cos`/`Tan`/`Cot`/`Sec`/`Csc`, the hyperbolic counterparts, the
    inverse trig/hyperbolic forms, and `Power[g, n]` for integer `n >= 2`.
    Multi-branch heads introduce a fresh integer parameter `C[k]` and wrap
    each solution in `ConditionalExpression[..., Element[C[k], Integers]]`.
    Emits `Solve::ifun` on first use per call.
  - Single equality, single variable, both specialists above decline (because
    the equation carries `Sqrt[...]` / `x^(p/q)` / nested radicals) ->
    `Solve`SolveRadicalsEquality` (also below).
  - Multi-variable list, or `And`/`List` of equations -> `Solve`SolveLinearSystem`
    (also below).  The linear-system specialist accepts the same input shapes
    that the router uses to decide dispatch; it canonicalises each equation
    `lhs_i == rhs_i` to `lhs_i - rhs_i` and refuses (returns `NULL`) when the
    system is not affine in the variables.
  - When the linear-system specialist declines a non-affine system ->
    `Solve`SolveNonlinearSystem` (also below).  This handles genuinely
    nonlinear polynomial systems whose solution set is zero-dimensional
    (finitely many solutions) via a lexicographic Gröbner basis and
    triangular back-substitution.  Positive-dimensional systems (infinitely
    many solutions) emit `Solve::nsdim` and leave `Solve` unevaluated;
    non-polynomial systems also stay unevaluated.
- Inequalities and multi-equation transcendental systems are reserved for
  future work and currently leave `Solve[...]` unevaluated.  When the
  inverse-function specialist's outermost peel succeeds but the inner
  equation is unsolvable and the peel was over `var` itself, Solve returns
  `{{var -> InverseFunction[head][rhs]}}` under `Solve::ifun`.
- **Approximate-number input**: if the equation contains any inexact numeric
  leaf (`Real` / MPFR), it is force-rationalised via the shared preprocessor
  in `src/common.c` before dispatch (so `1.5` becomes `3/2`, `N[Pi]` becomes
  a bit-exact rational, etc.), then the exact bindings produced by the
  specialist are numericalised on the way out -- same `inexact-in /
  inexact-out` contract `Integrate` and the exact-symbolic builtins
  (`Apart`, `Cancel`, `Together`, `Factor`, ...) follow.  The `vars`
  argument is never rationalised.  The preprocessor also tracks the
  *minimum* precision (in bits) across every inexact leaf and uses it
  both as the rationalisation tolerance and as the output precision, so a
  pure 30-digit-MPFR input flows back out at 30 digits, while a mixed
  Real + MPFR input drops to machine precision (the lower of the two)
  -- matching standard inexact-arithmetic semantics.
- Returns the solution set as a `List` of `List` of `Rule` pairs:
  - `{}` -- no solutions.
  - `{{}}` -- tautology (full-dimensional solution set).
  - `{{x -> v1}, {x -> v2}, ...}` -- one inner list per solution.  Multiplicity
    is preserved (repeated roots appear once per unit of multiplicity).
- **Rational-equality canonicalisation**: both sides are run through `Together`
  to combine into single fractions `N1/D1 == N2/D2`, then cross-multiplied to
  `N1*D2 - N2*D1 == 0` and `Collect`-ed in the solving variable before
  dispatch.  This routes equations like `a/x + b == 0` or `1/(x-1) == 2`
  through the polynomial specialist.  Any candidate root that provably zeroes
  one of the cleared denominators is dropped as extraneous (e.g.
  `Solve[x/(x-1) == 2/(x-1), x]` returns `{{x -> 2}}`, not `{{x -> 1}, {x -> 2}}`).
  Symbolic / undetermined denominator values are kept (parametric inputs like
  `Solve[a/x + b == 0, x]` return `{{x -> -a/b}}`).
- **Hidden-zero coefficient stripping**: after `Collect[Expand[...], var]` the
  per-degree coefficients are tested in turn with `PossibleZeroQ` (top down).
  Coefficients that test as zero but are not structurally zero -- e.g.
  `Sqrt[5 + 2 Sqrt[6]] - Sqrt[3] - Sqrt[2]`, recognised through the Stage-2
  numeric ladder -- are folded out and the polynomial is rebuilt at its true
  degree before the fast-path classifier sees it.  Without this pass the
  quadratic formula would divide by such a hidden-zero leading coefficient
  (`Solve[Sqrt[5 + 2 Sqrt[6]] x^2 - Sqrt[3] x^2 - Sqrt[2] x^2 - x - 1 == 0, x]`
  reduces to the linear `-x - 1 == 0` and returns `{{x -> -1}}`).  A
  hidden-zero constant is treated as a tautology
  (`Solve[Sqrt[5 + 2 Sqrt[6]] - Sqrt[3] - Sqrt[2] == 0, x]` returns `{{}}`).
- Per-degree handling for irreducible factors:
  - Degree 1 / 2: closed-form rules.
  - Quadratic in `Reals`: discriminant-aware.  Δ < 0 → no real roots;
    Δ = 0 → the double root is emitted *twice* (multiplicity preserved
    in step with the `Complexes` path); Δ > 0 → two distinct real
    roots.
  - Binomial `a*x^n + b == 0`: all n complex roots, or the real
    radical(s) in `Reals`.  Odd-`n` real branch: `(−b/a)^(1/n)` when
    `−b/a > 0`, `0` when `−b/a == 0`, and `−((b/a)^(1/n))` when
    `−b/a < 0` -- the last case is the *real* `n`-th root, not the
    principal complex one that `Power[base, 1/n]` produces by default.
    Even-`n`: ±r with `−b/a > 0`, `0` with `−b/a == 0`, `{}` with
    `−b/a < 0`.  Complex roots (no Reals constraint) are emitted as
    `r * (-1)^(2k/n)` for the principal radical `r = (-b/a)^(1/n)` and
    `k = 0..n-1`, then folded by `Power`'s rational-exponent canonicaliser
    so output matches Mathematica's standard form (e.g.
    `Solve[x^5 + 1 == 0, x]` returns
    `{{x -> (-1)^(1/5)}, {x -> (-1)^(3/5)}, {x -> -1},
       {x -> -(-1)^(2/5)}, {x -> -(-1)^(4/5)}}`).
  - n-quadratic `a*x^(2n) + b*x^n + c == 0`: substitution `u = x^n` followed by
    two binomial sub-solves; 2n radical roots regardless of `Cubics` / `Quartics`.
  - Degree 3: held `Root[Function[t, p[t]], k]` objects unless `Cubics -> True`.
  - Degree 4: held `Root[]` objects unless `Quartics -> True`, which emits the
    four roots in closed-form radicals via Ferrari's resolvent-cubic method
    (Complexes only; a `Reals` request still yields `Root[]`).
  - Degree ≥ 5: held `Root[]` objects per irreducible factor.
- `Integers` domain is implemented as a post-pass over the `Reals` output:
  every candidate value is type-checked against `EXPR_INTEGER` /
  `EXPR_BIGINT` and dropped otherwise.  `Rational[p, q]`, irrational
  radicals (`Sqrt[2]`, `Power[2, 1/3]`, ...), held `Root[]` objects, and
  symbolic / parametric residues are *not* trusted to be integer-valued
  and are silently removed.  This means polynomials with one or more
  rational integer roots are returned correctly (`Solve[x^3 - 6 x^2 + 11
  x - 6 == 0, x, Integers]` -> `{{x -> 1}, {x -> 2}, {x -> 3}}` via
  factoring), but polynomials that only have irrational or symbolic
  integer roots return `{}`.  Higher-degree irreducibles default to
  `Root[]` form (`Cubics -> False`, `Quartics -> False`) and therefore
  yield `{}` under `Integers` unless the user opts into radical output.

**Attributes:** `Protected`.

## References

**See also:** [Log](../../elementary-functions/Log/), [Exp](../../elementary-functions/Exp/), [Sin](../../elementary-functions/Sin/), [Cos](../../elementary-functions/Cos/), [Tan](../../elementary-functions/Tan/), [Cot](../../elementary-functions/Cot/), [Sec](../../elementary-functions/Sec/), [Csc](../../elementary-functions/Csc/)

- von zur Gathen & Gerhard, "Modern Computer Algebra" (3rd ed.), Ch. 14 (polynomial roots and resolution).
- Geddes, Czapor & Labahn, "Algorithms for Computer Algebra" (1992), Ch. 9 (solving systems).
- Source: [`src/solve.c`](https://github.com/stblake/mathilda/blob/main/src/solve.c)
- Specification: [`docs/spec/builtins/solutions-of-equations.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/solutions-of-equations.md)
- Tests: [`tests/test_integrate_line.c`](https://github.com/stblake/mathilda/blob/main/tests/test_integrate_line.c)
- Tests: [`tests/test_root_numeric.c`](https://github.com/stblake/mathilda/blob/main/tests/test_root_numeric.c)
- Tests: [`tests/test_solve.c`](https://github.com/stblake/mathilda/blob/main/tests/test_solve.c)
- Tests: [`tests/test_solvenlsys.c`](https://github.com/stblake/mathilda/blob/main/tests/test_solvenlsys.c)

## Notes & additional examples

### Notes

`Solve` returns a list of solution rule-lists, one per solution; each inner
list assigns every requested variable. Complex roots are produced by default,
so `x^2 + 1 == 0` yields the conjugate pair `±I`, and irrational roots come
back in exact radical form (`±Sqrt[2]`). Linear systems are solved directly
and return a single rule-list. Cubic roots are reported using `(-1)^(1/3)`
style radicals; pass `Cubics -> False` / `Quartics -> False` to suppress
explicit radical forms when desired.

Nonlinear polynomial systems with finitely many solutions (a zero-dimensional
ideal) are solved via a lexicographic Gröbner basis and triangular
back-substitution, honouring the `Reals` / `Integers` domain. Systems with
infinitely many solutions (positive-dimensional ideals, e.g.
`Solve[x^2 - y^2 == 0, {x, y}]`) emit `Solve::nsdim` and are left unevaluated.
