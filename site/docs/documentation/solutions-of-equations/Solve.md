# Solve

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Solve[expr, vars]`**

Attempts to solve the equation or system expr for the variables vars.

**`Solve[expr, vars, dom]`**

Solves over the domain dom.  Default Complexes; Reals filters down to real roots via per-degree discriminant and sign tests; Integers finds all integer solutions.  For a single univariate polynomial it filters the roots to concrete integers; a polynomial equation or system with constraints invokes a Diophantine engine -- linear via Hermite Normal Form, Pell and generalised/negative Pell by continued fractions, binary and ternary quadratic forms, Mordell curves, Thue equations (Tzanakis-de Weger), sum-of-three-cubes (Booker; the mod-9 impossibility globally), exponential Diophantine (Catalan, Ramanujan-Nagell), and additive power-sum searches (meet-in-the-middle).  An empty result {} is always a proof of no solution; an out-of-reach input is left unevaluated, never guessed.

**`Log[x]^2-3Log[x]+2), and -- via the inverse-function specialist --`**

**`ConditionalExpression[..., Element[C[k], Integers]].  Emits`**

<details>
<summary>Notes</summary>

Options: Cubics              -\> False     (radical form for cubics) Quartics            -\> False     (radical form for quartics) InverseFunctions    -\> Automatic (use inverse-function peel) GeneratedParameters -\> C         (head for parameters C\[k\]) VerifySolutions     -\> Automatic (True: drop non-verifying) Modulus             -\> 0         (solve over Z/pZ when p\>0) Solves single polynomial equalities, radical equations, linear systems, zero-dimensional nonlinear polynomial systems (via a lexicographic Groebner basis and triangular back-substitution; positive-dimensional systems emit Solve::nsdim and stay unevaluated), a single non-affine equation in several variables (solved for the earliest variable it is polynomial in, e.g. x y == 1 -\> {{x -\> 1/y}}), equations that are a polynomial in one transcendental kernel g(x) (u = g(x); e.g. E^(2x)-3E^x+2 and single-variable equations whose outermost dependence is an elementary invertible head (Log, Exp, Sin/Cos/Tan/Cot/Sec/Csc, their hyperbolic counterparts, the inverse trig/hyperbolic forms, and Power\[g, n\] for integer n \>= 2).  Multi-branch heads introduce an integer parameter C\[k\] wrapped in Solve::ifun the first time inverse functions are used.

</details>

## Examples (40)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (21)

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

In[9]:= Solve[a/x + b == 0, x]
Out[9]= {{x -> -a/b}}

In[10]:= Solve[1/(x-1) == 2, x]
Out[10]= {{x -> 3/2}}

In[11]:= Solve[x/(x-1) == 2/(x-1), x]
Out[11]= {{x -> 2}}

In[12]:= Solve[x^2 - 5 x + 6 == 0, x, Integers]
Out[12]= {{x -> 2}, {x -> 3}}

In[13]:= Solve[x^2 - 2 == 0, x, Integers]
Out[13]= {}

In[14]:= Solve[1.5 x + 3 == 0, x]
Out[14]= {{x -> -2.0}}

In[15]:= Solve[{1.5 x + y == 4.5, x - y == 0.5}, {x, y}]
Out[15]= {{x -> 2.0, y -> 1.5}}

In[16]:= Solve[N[Pi, 50] x == 1, x]
Out[16]= {{x -> 0.318309886183790679458993446443125399383289683168242}}

In[17]:= Solve[x^3 + 1 == 0, x, Reals]
Out[17]= {{x -> -1}}

In[18]:= Solve[x^3 - 6 x^2 + 11 x - 6 == 0, x, Integers]
Out[18]= {{x -> 1}, {x -> 2}, {x -> 3}}

In[19]:= Solve[3 x + 2 y == 11 && x + y == 12, {x, y}]
Out[19]= {{x -> -13, y -> 25}}

In[20]:= Solve[a x + c == 1 && b x - d y == 2, {x, y}]
Out[20]= {{x -> (1 - c)/a, y -> (-2 a + b - b c)/(a d)}}

In[21]:= Solve[3 x + 2 y == 11 && x + y == 12 && 3 x + y == 32, {x, y}]
Out[21]= {}
```

### Worked examples (10)

```mathematica
In[22]:= Solve[x y == 1, {x, y}]
Out[22]= {{x -> 1/y}}

In[23]:= Solve[x^2 + y^2 == 1, {x, y}]
Out[23]= {{x -> -1/2 Sqrt[-4 (-1 + y^2)]}, {x -> 1/2 Sqrt[-4 (-1 + y^2)]}}

In[24]:= Solve[E^(2x)-3E^x+2==0, x]
Out[24]= Solve[2 - 3.0 2.71828^x + E^(2 x) == 0, x]

In[25]:= Solve[Log[x]^2-3Log[x]+2==0, x]
Out[25]= {{x -> E}, {x -> E^2}}

In[26]:= Solve[3 x == 1, x, Modulus -> 7]
Out[26]= {{x -> 5}}

In[27]:= Solve[{x^2 + y^2 == 1, x == y}, {x, y}, Modulus -> 7]
Out[27]= {{x -> 2, y -> 2}, {x -> 5, y -> 5}}

In[28]:= Solve[x^3 - 6 x^2 + 11 x - 6 == 0, x, Integers]
Out[28]= {{x -> 1}, {x -> 2}, {x -> 3}}

In[29]:= Solve[x^2 + 2 y^3 == 3681 && x > 0 && y > 0, {x, y}, Integers]
Out[29]= {{x -> 15, y -> 12}, {x -> 41, y -> 10}, {x -> 57, y -> 6}}

In[30]:= Solve[x^2 - 61 y^2 == 1 && x > 0 && y > 0 && x < 10^10, {x, y}, Integers]
Out[30]= {{x -> 1766319049, y -> 226153980}}

In[31]:= Solve[x^3 + y^3 == z^3 && z > x > y > 0 && x,y,z < 10000, {x,y,z}, Integers]
Out[31]= Solve[x^3 + y^3 == z^3 && z > x > y > 0 && x, y, z < 10000, {x, y, z}, Integers]
```

### Applications (9)

```mathematica
In[32]:= Solve[x^2 - 5 x + 6 == 0, x]
Out[32]= {{x -> 2}, {x -> 3}}

In[33]:= Solve[x^2 + 1 == 0, x]
Out[33]= {{x -> -I}, {x -> I}}

In[34]:= Solve[x^2 - 2 == 0, x]
Out[34]= {{x -> -Sqrt[2]}, {x -> Sqrt[2]}}

In[35]:= Solve[{x + y == 3, x - y == 1}, {x, y}]
Out[35]= {{x -> 2, y -> 1}}

In[36]:= Solve[x y == 1 && x + y == 3, {x, y}]
Out[36]= {{x -> 1/2 (3 - Sqrt[5]), y -> 1/2 (3 + Sqrt[5])}, {x -> 1/2 (3 + Sqrt[5]), y -> 1/2 (3 - Sqrt[5])}}

In[37]:= Solve[x y == 6 && x + y == 5, {x, y}, Integers]
Out[37]= {{x -> 3, y -> 2}, {x -> 2, y -> 3}}

In[38]:= Solve[a x^2 + b x + c == 0, x]
Out[38]= {{x -> (1/2 (-b + Sqrt[b^2 - 4 a c]))/a}, {x -> (1/2 (-b - Sqrt[b^2 - 4 a c]))/a}}

In[39]:= Solve[x^4 - 1 == 0, x]
Out[39]= {{x -> -1}, {x -> 1}, {x -> -I}, {x -> I}}

In[40]:= Solve[Sin[x] == 0, x]
Out[40]= {{x -> ConditionalExpression[Pi + 2 C[1] Pi, Element[C[1], Integers]]}, {x -> ConditionalExpression[2 C[1] Pi, Element[C[1], Integers]]}}
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
- `VerifySolutions -> Automatic`: With `VerifySolutions -> True`, every
  returned solution is back-substituted into the equation(s) and dropped when
  `PossibleZeroQ` proves the residual non-zero; solutions that verify or are
  undecidable (`Root[]`, free parameters, `ConditionalExpression`) are kept.
  The default `Automatic` keeps per-specialist verification (e.g. radicals).
- `Modulus -> 0`: With `Modulus -> p` (`2 <= p <= 100000`), solve a
  single-variable polynomial over `Z/pZ` (see **Modular solving** above).

**Domains** (third positional argument): `Complexes` (default), `Reals`
(discriminant / sign filtering), `Integers` (keep provably-concrete integer
roots), and `Rationals` (keep provably-concrete `Integer`/`Rational` roots --
,
).  `Algebraics`, `Booleans`, and
`Primes` are not yet wired and leave `Solve` unevaluated.

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
  - A **single non-affine equation in several variables** (linear-system
    specialist declined, but the input is one `Equal`, not a conjunction) is
    solved for the earliest-listed variable it is polynomial in, treating the
    rest as symbolic parameters: `Solve[x y == 1, {x, y}]` -> `{{x -> 1/y}}`,
    `Solve[x^2 + y^2 == 1, {x, y}]` -> `{{x -> -Sqrt[1-y^2]}, {x -> Sqrt[1-y^2]}}`.
    This yields explicit rules only (never inequalities or case splits, which
    belong to `Reduce`).
  - When the linear-system specialist declines a genuine multi-equation
    system -> `Solve`SolveNonlinearSystem` (also below).  This handles
    nonlinear polynomial systems whose solution set is zero-dimensional
    (finitely many solutions) via a lexicographic Gröbner basis and
    triangular back-substitution.  Positive-dimensional systems (infinitely
    many solutions) emit `Solve::nsdim` and leave `Solve` unevaluated;
    non-polynomial systems also stay unevaluated.
  - **Polynomial in a single transcendental kernel** `g(x)` (single equation,
    single variable, the peel/trig/radical passes all declined): if
    substituting `u = g(x)` makes the equation a polynomial in `u` free of
    `x`, it is solved in `u` and each root `u0` unwound through `g(x) == u0`.
    Two kernel shapes: exponential `E^(c x)` (`Solve[E^(2x)-3E^x+2==0, x]` ->
    `x = 0, Log[2]` with periodic families in `Complexes`) and generic
    invertible heads `H[x]^k` (`Solve[Log[x]^2-3Log[x]+2==0, x]` ->
    `{{x -> E}, {x -> E^2}}`).  Implemented in `src/solvetrig.c`
    (`solvetrig_solve_poly_in_kernel`), reusing the polynomial and
    inverse-function specialists.
- **Modular solving.** `Solve[poly == 0, x, Modulus -> p]` solves a
  single-variable polynomial equation over the finite ring `Z/pZ` by residue
  enumeration (`src/solvemod.c`), returning `{{x -> r}, ...}` with `r`
  ascending in `[0, p)`: `Solve[x^2 == 2, x, Modulus -> 7]` -> `{{x -> 3},
  {x -> 4}}`, `Solve[3 x == 1, x, Modulus -> 7]` -> `{{x -> 5}}`.  Supported
  for `2 <= p <= 100000` (prime or composite; rational coefficients handled
  via modular inverse).  Non-polynomial equations and out-of-range moduli leave
  `Solve` unevaluated -- the option is never silently ignored.
- **Modular systems.** `Solve[{system}, {vars}, Modulus -> p]` with **prime** `p`
  solves a polynomial system over the finite field `GF(p)`: a finite-field
  Gröbner basis (`src/poly/gbmod.c`) is computed and its lex triangular form is
  walked with per-variable residue enumeration.  `Solve[{x^2 + y^2 == 1, x == y},
  {x, y}, Modulus -> 7]` -> `{{x -> 2, y -> 2}, {x -> 5, y -> 5}}`; an
  inconsistent system (unit ideal) -> `{}`; an under-determined system
  enumerates the free variables over `GF(p)` (`Solve[{x + y == 1}, {x, y},
  Modulus -> 3]` -> the three points).  A **composite** modulus is not a field,
  so systems with composite `p` are refused (unevaluated); a coefficient whose
  denominator is divisible by `p` (no image in `GF(p)`) is likewise refused.
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
  - **`Reals`/`Integers`/`Rationals` reality filter for `Root[]`.** Held
    `Root[]` objects are emitted with the *full* index range for an irreducible
    factor; a post-dispatch filter at `Solve`'s funnel then drops every solution
    whose bound value is a *provably non-real* number (numericalised to a
    `Complex[re, im]` with a concrete `|im| > 1e-9`). So
    `Solve[x^5 - x - 1 == 0, x, Reals]` returns the single real `Root[.., 1]`
    (not all five), and `Solve[x^6 - x - 1 == 0, x, Reals]` returns two. The
    filter is conservative — real `Root[]` objects, concrete reals, and
    symbolic/parametric values that do not numericalise (e.g. `Sqrt[a]`) are
    kept — and covers polynomial *systems* the same way (complex `Root`-tuples
    are dropped over `Reals`). The default `Complexes` domain is untouched.
- For a **single univariate polynomial equation** (no constraints), the
  `Integers` domain is implemented as a post-pass over the `Reals` output:
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
- **Diophantine solving (`Integers` with constraints).** When the input is a
  polynomial equation (or system) conjoined with inequality / ordering /
  disequation constraints, a dedicated pre-pass (`src/solve/`) finds *all*
  integer solutions:
  `Solve[x^2 + 2 y^3 == 3681 && x > 0 && y > 0, {x, y}, Integers]` ->
  `{{x -> 15, y -> 12}, {x -> 41, y -> 10}, {x -> 57, y -> 6}}`.
  The method is bound propagation to a finite box (explicit bounds, ordering
  chains like `0 < x <= y <= z`, **absolute-value ordering chains** like
  `Abs[x] < Abs[y] < Abs[z] < B`, and an interval-positivity rule that turns a
  sign-definite term of `Σ term == constant` into a per-variable bound, both
  above and — for odd powers, deducing the sign — below). An abs-value ordering
  is the natural way to ask for the *ordered representatives* of a symmetric
  solution set (one per permutation orbit) instead of every permutation: the
  magnitude chain propagates a box onto every variable (`|x| < |y| ≤ B ⟹
  |x| ≤ B-1`) and filters the result to the ordered subset, so
  `Solve[x^3 + y^3 + z^3 == 63 && Abs[x] < Abs[y] < Abs[z] < 10000, {x,y,z},
  Integers]` returns the 6 ordered solutions rather than the 36 permutations of
  the box form. A partial chain (`Abs[x] < Abs[y] < 10000 && Abs[z] < 10000`)
  bounds and orders only the named variables. A variable that
  appears only with **even exponents** is sign-symmetric, so even without a
  sign constraint it is bounded on both sides to `[-B, B]` — this makes the
  unconstrained sum of even powers finite, so
  `Solve[x^2 + y^2 == 25, {x, y}, Integers]` returns all 12 signed pairs
  (and `x^2 + y^2 == 0` the origin) rather than the empty set. Then recursive
  elimination that enumerates all but one variable and solves the last
  *exactly* (integer k-th root, quadratic discriminant, or rational-root),
  every candidate re-verified against the original conjunction. A **single
  separable additive equation** (`Σ g_i(x_i) == c`, e.g. sums of powers, the
  taxicab equation) is instead solved by **meet-in-the-middle** in
  ~`N^ceil(n/2)` work. Only necessary conditions tighten a bound, so an
  exhausted finite search returns `{}` as a proof of no solutions; an input
  that cannot be bounded to a finite box (an unbounded Pell orbit, a
  constraint-free Thue equation) is left **unevaluated** rather than answered
  wrongly. `Solve`SolveIntegers[eqns, vars]` is the independently-testable
  entry point.
  - **`Solve::svars` diagnostic.** If the system carries a symbol in an
    inequality/ordering constraint (`… && d > 0 && d < 100000`) that is not among
    the solve variables — almost always a mistyped variable list, e.g. the `d`
    dropped from `{x, y, z, y}` — `Solve` emits `Solve::svars` ("Equations may not
    give solutions for all \"solve\" variables") rather than silently declining.
    Operator heads and named constants are skipped, and a bare parameter in an
    *equation* (`a x == b`) is not flagged, so the warning is low-noise.
- **Divisor-factoring and reciprocal special forms.** Two shapes that
  positivity cannot bound are still finite once the right identity is applied:
  - A single **bilinear** equation `a*u*v + b*u + c*v + d == 0` (reached by
    eliminating unit-coefficient linear equations — the router tries each pair
    of variables to keep) factors as `(a*u + c)(a*v + b) = b*c - a*d`, so the
    integer solutions come from the **divisors** of that constant with no
    enumeration of `u, v`.  This solves the Pythagorean-with-perimeter case
    `x^2 + y^2 == z^2 && x + y + z == 3000 && 0 < x < y && z > 0` ->
    `{{500, 1200, 1300}, {600, 1125, 1275}, {750, 1000, 1250}}` (with `z > 0`;
    the constraint-free system also admits negative-`z` solutions, which are
    returned when not excluded).
  - A sum of unit fractions `sum 1/x_i == R` with an ordering chain
    `x_1 <= ... <= x_k` bounds the smallest variable to
    `[ceil(1/R), floor(k/R)]` and recurses, the last variable determined
    exactly.  This solves the Egyptian-fraction case
    `4/2027 == 1/x + 1/y + 1/z && 0 < x <= y <= z`.  Because the equation is
    fully symmetric in its variables, **no ordering need be supplied**: with only
    positivity, the ascending representatives are found and every distinct
    permutation of each is emitted (re-verified), so
    `4/5 == 1/x + 1/y + 1/z && x > 0 && y > 0 && z > 0` returns the full
    unordered set (all 12 permutations of `{2,4,20}` and `{2,5,10}`).
  - A separable **odd-power sum** (e.g. `x^3 + y^3 + z^3 == 42`) over a box too
    large for the leaf search is solved by the **divisor method**: because
    `e` is odd, `s = x + y` divides `m = x^e + y^e`, and the power sum in terms
    of `s` and `p = x y` is a degree-`e/2` polynomial, so for each divisor `s`
    of `m` the integer roots `p` give `(x, y)`.  Fixing the remaining variables
    turns the `O(N^2)` inner search into `O(N * factoring)` --
    `x^3 + y^3 + z^3 == 42 && Abs[...] < 10^5` is settled in ~7 s (the search
    space is 8x10^15), and `x^3 + y^3 == 1729 && 0 < x <= y` gives Ramanujan's
    `{{1, 12}, {9, 10}}`.  It applies to any odd exponent (`x^5 + y^5 == 1267`
    -> `{{3, 4}}`); higher powers are admitted only over boxes small enough
    that `m` stays in the fast-factoring range.
  - A **Pell** equation `x^2 - D y^2 == +/-1` (D a positive non-square) is
    solved from the continued fraction of `sqrt(D)`: the fundamental unit
    generates the whole orbit, enumerated up to any explicit bound.
    `Solve[x^2 - 61 y^2 == 1 && x > 0 && y > 0 && x < 10^10, {x, y}, Integers]`
    -> `{{x -> 1766319049, y -> 226153980}}`; the negative Pell
    `x^2 - 3 y^2 == -1` correctly returns `{}` (unsolvable).
  - **Multi-leaf staged elimination.** A variable that appears in exactly one
    equation and is univariate-solvable there is *peeled* -- resolved by an
    exact root per free-variable assignment rather than enumerated -- so a
    system with several "determined" variables reduces to a search over only
    the coupled ones.  The **Euler brick**
    `x^2+y^2==a^2 && x^2+z^2==b^2 && y^2+z^2==c^2 && 0<x<y<z<500 && a,b,c>0`
    peels `a,b,c` (three square roots) and walks only `x<y<z`, returning all
    three bricks in ~1 s, the smallest `(44,117,240;125,244,267)`.
  - **Ordered box + int64 fast leaf.** The search-space guard divides the raw
    box by the factorial of the longest ordering chain, and a degree-≤2 leaf
    over small coefficients is solved in machine integers (GMP fallback on
    overflow), so a genuinely ordered four-variable box is enumerated rather
    than declined: `2(x^2+y^2+z^2+w^2)==(x+y+z+w)^2 && 0<x<=y<=z<=w<1000` and
    the Markov-Hurwitz `x1^2+x2^2+x3^2+x4^2==x1 x2 x3 x4 && 0<x1<=...<=x4<=1000`.
  - **Non-polynomial power-leaf.** When one side is a pure power `m^e` of a
    leaf that appears nowhere else and every other variable is bounded, the
    others are enumerated, the remaining side is evaluated through the
    interpreter (so `Factorial`, `Binomial`, ... work), and `m` is solved by an
    exact integer `e`-th root.  **Brocard's problem** `n! + 1 == m^2 && 0<n<100`
    returns the Brown numbers `n = 4, 5, 7`.
  - **Binary-quadratic conic.** `Y^2 == A X^2 + B X + C` with a perfect-square
    leading coefficient `A` completes to a difference of squares
    `(2 p Y)^2 - (2 A X + B)^2 = 4 A C - B^2` and factors that constant over its
    divisors -- exhaustive, so an empty result is a proof.  Euler's
    `n^2 + n + 41 == y^2` -> `{n -> 40, y -> 41}`; `x^2 - y^2 == 15` ->
    `{{4,1},{8,7}}`.  (A non-square `A` is a genuine Pell conic, left to the
    continued-fraction path.)
  - **Definite binary quadratic (ellipse).** A single 2-variable degree-2
    equation with a **negative** discriminant `delta = B^2 - 4AC < 0` is a compact
    ellipse -- finite, but a rotated one (`B != 0`) escapes the interval bounder.
    Solved as a quadratic in `x` for each `y` in the finite interval where the
    `x`-discriminant `delta y^2 + (2BD - 4AE) y + (D^2 - 4AF)` (a downward
    parabola) is `>= 0`, exhaustively -- so `{}` is a proof:
    `x^2 + x y + y^2 == 7` -> 12 points; `x^2 + x y + y^2 == 2` -> `{}`.
    Negative-definite forms are normalised; linear terms and constraints are
    handled.
  - **Factorable binary quadratic (Runge's simplest case).** A single 2-variable
    equation `A x^2 + B x y + C y^2 + D x + E y + F == 0` whose quadratic part has
    a **cross term** and a perfect-square discriminant `δ = B^2 - 4AC > 0` factors
    into two rational linear forms, so it is a hyperbola with finitely many
    integer points.  Completing the square (via `U = 2Ax + By + D`) reduces it to
    a difference of squares `(2 k U)^2 - V^2 = W` (`k = √δ`,
    `W = -(P^2 + 4 k^2 Q)`, `P = 4AE - 2BD`, `Q = 4AF - D^2`) and factors `W` over
    its divisors — exhaustive, so an empty result is a proof:
    `x^2 + x y - 2 y^2 == 4` -> six points, and `(x - y)(x + 2 y) == 15` -> `{}`
    (a mod-3 obstruction, not a decline).  Handles non-unit square coefficients
    (`2 x^2 + 3 x y - 2 y^2 == 7` -> `{(-3,1),(3,-1)}`) that the conic form above
    cannot.  A non-square `δ` (Pell-type) or `δ ≤ 0` (parabola/ellipse) is
    declined here.
  - **Prouhet-Tarry-Escott -> {}.** Two `k`-element groups with equal power sums
    for degrees `1..k` are the same multiset (Newton's identities); with a strict
    ordering inside each group they are forced equal, so a disequation such as
    `a != d` proves the system empty -- e.g. `a+b+c==d+e+f && ...(deg 2)... &&
    ...(deg 3)... && 0<a<b<c && 0<d<e<f && a!=d` -> `{}` though every variable is
    unbounded.
  - **Unbounded Mordell.** `y^2 == x^3 + k` factors as
    `x^3 = (y - sqrt k)(y + sqrt k)` in `Z[sqrt k]`; the cube factors give the
    COMPLETE integer-point set whenever the descent is sound -- `k < 0`, `|k|`
    squarefree, `k = 2,3 (mod 4)` (units `{+/-1}`, and a mod-8 argument forces
    the two factors coprime), and `3` not dividing the class number of
    `Q(sqrt k)` (so an ideal cube is a principal cube).  So `y^2 == x^3 - 2` ->
    `(3, +/-5)`, `y^2 == x^3 - 13` -> `(17, +/-70)`, and `y^2 == x^3 - 5` -> `{}`
    (proved).  A bounded `x` uses the ordinary leaf search; the half-integer ring
    (`k = 1 mod 4`), `3 | h`, and the real-quadratic case (`k > 0`, infinite units
    -- e.g. `y^2 == x^3 + 3`) are left unevaluated.
  - **Sum of three cubes (`x^3 + y^3 + z^3 == k`, Booker method).** A box-bounded
    sum of three cubes with a fixed nonzero integer `k` is solved by the
    cube-root-mod-`d` method of A. R. Booker ("Cracking the problem with 33"),
    rather than linearly enumerating one variable and factoring the ~`B^3` value
    `k - c^3` (the classical divisor path, whose reach is capped at the ~3×10⁵
    outer budget). Since `k - c^3 = (a+b)(a^2-ab+b^2)`, the divisor `d = |a+b|`
    satisfies `c^3 ≡ k (mod d)`; enumerating the SMALL `d` up to `α·B`
    (`α = ∛2 - 1`) and cube-rooting `k mod d` pins `c` to arithmetic progressions,
    with `{a,b} = (sgn(m)·d ± √((4|m|/d - d^2)/3))/2`. This reaches coordinates
    up to ~10⁶ where the classical path declines, e.g.
    `Solve[x^3 + y^3 + z^3 == 2 && -200000 <= x <= 200000 && … , {x,y,z}, Integers]`
    returns all 195 solutions (including `(162001, -161999, -5400)`) in ~0.4 s.
    Because its `O(α·B · roots)` work also beats the leaf search's `O(B^2)` and
    the classical divisor path's `O(B · factoring)` for *small* boxes, it now
    engages for any non-trivial box (`|coord| > 100`; validated to return the
    identical solution set to the leaf/divisor pipeline over a wide `(k, box)`
    sweep, including the `(a,-a,∛k)` family and signed equations). E.g.
    `x^3 + y^3 + z^3 == 63 && Abs[...] < 10000` drops from ~0.46 s to ~0.02 s,
    and the `|coord| < 5000` box from ~3.9 s (leaf) to ~0.01 s; only a trivially
    small box (already sub-millisecond) is left to the leaf/mitm paths.
    The result is complete (small-coordinate solutions and the `(a,-a,∛k)` family
    covered by a divisor sub-search, the rest by Booker's `d < α·|c|` bound across
    the three roles). Restricted to `|k| < ~10⁹`; a box that would hold a huge
    parametric family (> 200 000 tuples) is declined rather than materialised. The
    underlying "all cube roots of `k` mod `d`" primitive is
    `Solve`CubeRootsMod[k, d]`. Divisors are factored through a smallest-prime-
    factor sieve (built once per solve, O(log d) per `d`) and the coordinate
    arithmetic is 128-bit, extending the reach to coordinates ~10⁷ (e.g. the
    point at 5 821 795 for a radius-6×10⁶ box in ~16 s); a box beyond the divisor
    budget or the candidate backstop declines (unevaluated) rather than running
    unbounded. The detector accepts coefficient `±1` on each cube: a `−v^3` is
    normalised away by the substitution `u = −v` (mirroring that variable's box),
    so `± x^3 ± y^3 ± z^3 == k` over a pure box reduces to the same solver — e.g.
    `Solve[x^3 + y^3 − z^3 == 227 && −200000 <= x,y,z <= 200000, {x,y,z}, Integers]`
    returns `(24579, 51748, 53534)` (the classic `227 = 24579^3 + 51748^3 −
    53534^3`). The sign substitution is used only for a pure box (no orderings /
    disequations), so it stays exact.
    A **global mod-9 obstruction** (`si_solve_three_cubes_mod9`) short-circuits
    the *unbounded* case the Booker box search must otherwise decline: every cube
    is `≡ {−1, 0, 1} (mod 9)`, so `± x^3 ± y^3 ± z^3` can never be `≡ 4` or `5
    (mod 9)`; hence `x^3 + y^3 + z^3 == k` with `k ≡ ±4 (mod 9)` has no integer
    solution at all, and `Solve[x^3 + y^3 + z^3 == 4, {x, y, z}, Integers] ->
    {}` is returned as a proof with **no bound required** (checked before the
    Booker engine, independent of the cube signs).
  - **Sum of like powers = a like power (ordering-aware 128-bit MITM).** A
    single **separable** additive equation `Σ cᵢ vᵢ^k == c₀ y^k` over an ordered
    box is solved by a meet-in-the-middle that improves on the plain int64
    `mitm_solve` in two ways: partial sums are `__int128` (so k-th powers past
    2⁶³ do not force a decline), and the variables are split along their ordering
    chain into a contiguous hash prefix and iterate suffix, enumerating only
    **ordered** tuples (combinations, not the Cartesian product). This is what
    the **Lander–Parkin** quintic needs:
    `Solve[x^5 + y^5 + z^5 + w^5 == r^5 && 0 < x < y < z < w < r < 1000, {x,y,z,w,r}, Integers]`
    returns the complete set — `27^5 + 84^5 + 110^5 + 133^5 = 144^5` (the 1966
    counterexample to Euler's conjecture) and its 2×–6× multiples — in ~6.5 s,
    where the plain MITM rejects the box because its *unordered* iterate product
    (`~10^9`) exceeds the node budget while the *ordered* count (`C(1000,3) ≈
    1.6×10^8`) fits. A modular residue sieve (built from the collected hash sums
    mod a structured `M`, e.g. mod 11/25/31/41 for fifth powers) skips the
    binary search for iterate tuples whose complement residue is unreachable.
    Engages only when it adds capability the plain MITM lacks (values overflow
    int64, or an ordering chain makes an otherwise-too-big iterate side fit); a
    small box is left to the existing path unchanged. Exhaustive over the box, so
    it returns the complete set or a proven `{}`.
  - **Sum of three biquadrates = a biquadrate (Frye's search).** `x^4 + y^4 +
    z^4 == w^4` over a box is searched by Frye's method (R. E. Frye, "Finding
    95800⁴ + 217519⁴ + 414560⁴ = 422481⁴ on the Connection Machine", 1988) — the
    minimal counterexample to Euler's conjecture. A box up to `10^6` cannot be
    exhaustively verified interactively, so this is a **witness search**: it finds
    and returns the minimal solution, ascending in `w`, and declines (never a
    spurious `{}`) if the node budget is spent with nothing found. The number
    theory: for a primitive solution the fourth powers mod 5 force exactly one
    summand `C != 0 (mod 5)`, the other two `A, B` to be multiples of 5, and
    `w != 0 (mod 5)`; then `625 | (w^4 − C^4)` and `N = (w^4 − C^4)/625 = a^4 +
    b^4` is decomposed by scanning `a` over `[~0.841 N^{1/4}, N^{1/4}]`. Extra
    moduli coprime to 5 (whose sum-of-two-fourth-powers residue set is a proper
    subset) and Frye's prime-factor constraint (any odd prime `P != 1 (mod 8)`
    dividing `N` must appear to an exponent `≡ 0 (mod 4)`) prune the `(w,C)`
    pairs before the decompose — all sound necessary conditions, so no real
    primitive solution is dropped. `__int128` throughout, plus a mod-2¹⁶ 4th-power
    fast reject. The full
    `Solve[x^4 + y^4 + z^4 == w^4 && 0 < x < y < z < w < 1000000, {x,y,z,w}, Integers]`
    finds `{x → 95800, y → 217519, z → 414560, w → 422481}` in ~11 min on a
    single modern core (Frye needed a 16384-processor Connection Machine for
    ~33 hours in 1988). Engages only for a box too large for the exhaustive MITM
    (`w > 20000`); a tunable node cap (`MATHILDA_FRYE_MAXNODES`) bounds the
    no-solution case. The target `w` need **not** be two-sided bounded by the
    user: `w^4 = x^4+y^4+z^4` gives a sound magnitude bound `|w| <= (Σ hi^4)^{1/4}`
    from the summand box, so a one-sided box (`w < 500000`) or an unconstrained
    `w` also engage — the search runs on `|w|` and emits **both** signs, each
    filtered by the user's constraints (so `0 < w` keeps only `+w`). Engagement
    is not tractability: a cold scan of a wide window is intrinsically the 1988
    computation, so the fast path is a witness window near the known `w`.
  - **Modular-sieved leaf search (large non-separable boxes).** When the
    ordering-pruned leaf box still exceeds `SI_MAX_NODES` (so the ordinary leaf
    search declines), a single polynomial equation is searched **exhaustively**
    by pruning the innermost enumerated variable to the residues (mod a small
    `M`) for which the leaf equation can vanish — a sound necessary condition, so
    the search stays complete: it returns the full set / a proven `{}` when it
    finishes in the raised budget, and declines (never a partial list) otherwise.
    This closes the gap for big cross-term boxes, e.g.
    `Solve[x^2 + x y + y^2 == z^2 && 0 < x < y < z < 15000, {x,y,z}, Integers]`
    (16386 solutions, ~6 s), validated identical to the ordinary engine on the
    sub-box it can already exhaust.
  - **Fermat's Last Theorem.** `x^n + y^n == z^n` (n >= 3, equal coefficient
    magnitudes, no constant) with `x, y, z` all strictly positive has no
    solutions (Wiles 1995), so this returns `{}` *immediately*, before any
    search and with no need for a finite box:
    `Solve[x^3 + y^3 == z^3 && z > x > y > 0 && x,y,z < 10000, {x,y,z}, Integers]`
    -> `{}` in ~0.4 ms (was ~2.6 s), and the unbounded
    `Solve[x^3 + y^3 == z^3 && x>0 && y>0 && z>0, {x,y,z}, Integers]` -> `{}`
    instantly rather than declining. The check fires only on the exact
    `a^n + b^n == c^n` shape with every lower bound `>= 1`; `n = 2` still returns
    Pythagorean triples, `x^3 + y^3 == z^3 + 1` still solves normally, and a box
    admitting `0` / negatives still enumerates the `(0, a, a)` solutions.
  - **Unbounded Pell -> parametric family.** `x^2 - D y^2 == 1 && x>0 && y>0`
    with no bound returns the fundamental-unit family as a
    `ConditionalExpression` on `C[1] >= 1`:
    `x -> ((x1+y1 Sqrt[D])^C[1] + (x1-y1 Sqrt[D])^C[1]) / 2` and the matching
    `y`, using the fundamental solution `(x1, y1)` from the continued fraction.
  - **Unbounded generalised Pell -> a family per class.** `x^2 - D y^2 == N`
    with `x>0 && y>0`, no bound, and **any `N != +1`** (including negative Pell
    `N = -1`) returns one `ConditionalExpression` family per solution class:
    `x, y -> ((a+b√D)(t+u√D)^C[1] ± (a-b√D)(t-u√D)^C[1]) / (2 or 2√D)`, `C[1] >= 0`,
    where `(t, u)` is the fundamental unit and `(a, b)` the class's minimal
    positive representative. The class representatives come from the **Nagell
    bound** `y ≤ u√(|N|/(2(t±1)))` (a finite search), advanced into the positive
    orthant and reduced by `ε⁻¹` to the minimal member so one class yields one
    family. Exhaustive, so an empty result is a proof:
    `x^2 - 2 y^2 == 7` -> two families with fundamentals `(3,1)`, `(5,3)`;
    `x^2 - 2 y^2 == 5` -> `{}` (a mod-8 obstruction); `x^2 - 3 y^2 == -1` -> `{}`
    (`√3` has even CF period). Without the positivity constraints the family is
    declined (unevaluated).
  - **Homogeneous linear system -> parametric ray.** `n-1` homogeneous linear
    equations in `n` positive unknowns have a one-dimensional integer kernel (the
    generalised cross product, via a fraction-free Bareiss determinant); if the
    primitive kernel vector is entirely positive the solutions are
    `{v_i C[1] : C[1] >= 1}`, otherwise the positive orthant meets the kernel only
    at the origin and there is no positive solution.
  - **General linear system -> HNF integer family.** An unconstrained system of
    `m >= 2` linear equations `A x == b` in `n` unknowns is solved completely
    over `Z` via the Hermite normal form (`HermiteDecomposition`, see the
    linear-algebra reference). With `P A^T == R` (row HNF), the substitution
    `x = P^T y` triangularises the system; forward substitution over the pivots
    with an exact-division test yields a particular solution (a divisibility
    failure is a **proof of no integer solution**, e.g.
    `2 x + 2 y == 3 && x - y == 0 -> {}`), and the free columns of `P^T` (the
    integer kernel lattice) become the parameters `C[k]`:
    `Solve[{x + 2 y + 3 z == 10, x - y + z == 2}, {x, y, z}, Integers]` ->
    `{{x -> 18 + 5 C[1], y -> 8 + 2 C[1], z -> -8 - 3 C[1]}}`. This replaced a
    silent wrong `{}` (the Complexes-oriented linear-system dispatch expressed
    the pivots as a *rational* family in the free variable and then discarded it
    as non-integer). A determined system reads off its unique solution or proves
    `{}`; a *bounded* system (any inequality present) still uses the finite leaf
    search, not this path.
- **Exponential Diophantine (variable exponents).** Equations such as
  `x^a - y^b == 1`, where the exponent is a solve variable, are handled before
  the polynomial stage (which cannot represent `x^a`).  A fully bounded box
  (`2^a - 3^b == -23 && 0 < a < 10 && 0 < b < 10` -> `{{a -> 2, b -> 3}}`) is
  enumerated exactly; the **Catalan** shape `x^a - y^b == +/-1` with bases and
  exponents `>= 2` is settled by **Mihailescu's theorem** -- the unique solution
  is `3^2 - 2^3 = 1`, so
  `Solve[x^a - y^b == 1 && 1 < x < 100 && 1 < y < 100 && a > 1 && b > 1,
  {x, y, a, b}, Integers]` -> `{{x -> 3, y -> 2, a -> 2, b -> 3}}` (and `{}` when
  the box excludes it).  The **fixed-base** Pillai form `P^m - Q^n == +/-1`
  (constant bases `P, Q`, variable exponents) is likewise settled by Mihailescu
  plus the exponent-1 cases, so `3^m - 2^n == 1` -> `{(1,1),(2,3)}` and
  `2^n - 3^m == 1` -> `{(1,2)}` even though `m, n` are unbounded.
- **Elliptic / hyperelliptic curves over a box.** `y^m == f(x)` with a bounded
  `x` (Mordell `y^2 = x^3 + k`, hyperelliptic `y^2 = quartic`) is solved by the
  ordinary bounded search -- enumerate `x`, test that `f(x)` is a perfect
  `m`-th power -- so `y^2 == x^3 - 10000 && 0 < x < 10^5 && y > 0` finds
  `{{25, 75}}` and `y^2 == x^3 - 2` gives Fermat's `{{3, 5}}`.  The *unbounded*
  Mordell curve is solved for every imaginary `k = 2,3 (mod 4)` with `|k|`
  squarefree and `3` not dividing the class number (see the `Z[sqrt k]`
  factorisation above); the half-integer ring, `3 | h`, the real-quadratic case
  `k > 0`, and higher-genus hyperelliptic curves need Mordell-Weil / Baker
  methods and are left unevaluated.
- **Thue equations** `F(x, y) == m` (`F` irreducible homogeneous of degree
  `>= 3`, `m` constant) have finitely many solutions, returned as a plain list
  by the **Tzanakis-de Weger** engine (`src/solvethue.c`): build `K = Q(theta)`
  for a root of `F(t, 1)`, reduce to a unit equation, bound the unit exponents
  by **Baker's linear forms in logarithms** (Waldschmidt) + **de Weger LLL**
  reduction, then enumerate and verify each `(x, y)` exactly. Scope today is a
  **monic** form (`|a0| = 1`) over a real field; out-of-scope inputs
  (`|a0| != 1`, precision out of reach, and the cases below) DECLINE rather than
  guess. **Non-monogenic** fields (`Z[theta] != O_K`) are handled by computing
  `O_K` with **Round 2 (Pohst-Zassenhaus)** and searching for units over the
  `O_K` lattice — so `x^3 - 17 y^3 == 1 -> {(1,0),(18,7)}`,
  `x^3 - 20 y^3 == 1 -> {(1,0),(-19,-7)}`, and the non-monogenic quartic
  `x^4 - 12 y^4 == 1` now solve. **Reducible** forms (`F(x,1)` factors into
  `>= 2` coprime irreducibles) are not Thue equations but still finite: they are
  factored and solved by enumerating the divisor assignments across the factors
  (any `m`) -- so `x^3 - x^2 y - 3 x y^2 - y^3 == 1 -> {(-1,2),(0,-1),(1,0)}` and
  `x^4 - y^4 == 15` returns its four points; a pure power of one factor
  (`(x-y)^3 == 1`, infinitely many) DECLINEs. Fundamental units come from an exact
  coefficient-box search
  certified by p-saturation; for **large-regulator complex cubics** whose unit
  exceeds any box (`Q(cbrt 15)` has a coordinate `30`, `Q(cbrt 41)` a 24-digit
  unit, regulator `56.3`), the box falls back to **Voronoi's algorithm** — a
  walk along the chain of relative minima of `O_K`, polynomial in the regulator
  — which proposes the unit for the *same* p-saturation certifier. So
  `x^3 - 15 y^3 == 1`, `x^3 - 41 y^3 == +-1`, `x^3 - 42 y^3 == 1`,
  `x^3 - 97 y^3 == 1` now return their complete sets (each `{(+-1, 0)}`) instead
  of declining. For **rank-2** fields the box likewise misses an
  intrinsically-large fundamental unit: the monogenic quintic `Q(5^{1/5})`
  (`x^5-5y^5`) is reached by a wider search box, and the signature-(2,1) quartic
  `Q(10^{1/4})` (`x^4-10y^4 == +-1`) by a **rank-2 Voronoi minima walk** in the
  relative-unit direction (`nfvoronoi2.c`) paired with the box's subfield unit —
  again certified by p-saturation (its regulator matches PARI's `bnfinit`).
  - **General `m` (`|m| != 1`).** For a monic form `N(x - theta*y) = F(x,y) =
    m`, so `beta = x - theta*y` is a norm-`m` integer `mu * unit`, with `mu`
    over a finite set of **bounded-norm representatives**.  For a rank-1 complex
    cubic these are enumerated (canonical-orbit box from the fundamental domain
    of the unit + the norm constraint), the Baker/de-Weger bound is made
    μ-aware (the linear form's constant gains the `log(mu^(k)/mu^(j))` term),
    and `beta = mu * prod eps^b` is enumerated per μ.  So `x^3 - 2 y^3 == 2 ->
    {(0,-1)}`, `== 3 -> {(1,-1),(-5,-4)}`, `== 10 -> {(2,-1),(4,3)}`, and
    `== 4/5/9/73/100 -> {}` (proven).  Rank-2 totally-real fields with
    `|m| != 1` still DECLINE.
  - **Totally complex fields (`r1 = 0`, any `m`).** When the field has no real
    embedding the Baker/unit machinery has no real type-index, but no factor can
    be small either: every root `theta_i` is non-real, so for real `x, y`,
    `|x - theta_i y| >= |Im(theta_i)| * |y|`, and
    `|m| = prod_i |x - theta_i y| >= |y|^n * prod_i |Im theta_i|` bounds `|y|`
    **elementarily and rigorously** -- no units, torsion, or Baker bound. Each
    `y` is closed by exact univariate root-finding. So the cyclotomic quartic
    `x^4 + x^3 y + x^2 y^2 + x y^3 + y^4 == 1` (over `Q(zeta_5)`) returns its 6
    points, `x^4 + y^4 == {1, 2, 17, 82}` their `{4, 4, 8, 8}`, `== 3 -> {}`, and
    higher cyclotomics (`Phi_7`, `Phi_10`) likewise.
- **Homogeneous ternary quadratic / Legendre** (`si_solve_ternary_quadratic`,
  `src/solve/solveint_ternary.c`). A single homogeneous diagonal degree-2
  equation in exactly 3 variables over an unbounded domain,
  `a x^2 + b y^2 + c z^2 == 0`, is decided by **Legendre's theorem** rather than
  searched: after clearing the gcd and normalising, solvability is the condition
  that `-1` is a quadratic residue modulo every odd prime factor of the
  opposite-sign coefficient (tested per prime with `mpz_legendre`). An
  unsolvable form returns the **trivial solution** `{{x->0,y->0,z->0}}` — a
  proof of no nontrivial solution, *not* `{}` (matching Mathematica) — so
  `Solve[x^2 + y^2 == 3 z^2, {x,y,z}, Integers] -> {{0,0,0}}`. A solvable form
  returns the **complete** integer family: `Solve[x^2 + y^2 == z^2, …]` is the
  Pythagorean parametrisation and `Solve[x^2 + y^2 == 2 z^2, …]` the full family
  including the tangent point `(1,1,1)`. The family is the sign/swap orbit of a
  chord/tangent parametrisation from a sum-of-two-squares witness `P0`, unioned
  with the tangent line `C·P0` (whose odd multiples the chord map alone misses);
  it is validated complete + sound against a brute-force oracle for every
  single-representation solvable case. Scope: the symmetric pattern
  `x^2 + y^2 == k z^2` with `k` squarefree and at most one prime factor `== 1
  (mod 4)`; this integer-exact solver keeps its stronger completeness guarantee
  for those cases. Everything it declines — cross terms, non-symmetric `a != b`
  diagonals, and multi-representation `k` (e.g. `65 = 5*13`) — is picked up by
  the **general** ternary solver below.
- **General homogeneous ternary quadratic** (`si_solve_ternary_general`,
  `src/solve/solveint_ternary_general.c`). Any single homogeneous degree-2
  equation in exactly 3 variables with cross terms and/or general coefficients,
  `a x^2 + b y^2 + c z^2 + d x y + e y z + f z x == 0`, over an unbounded domain.
  A non-degenerate ternary quadratic with one rational point is a **genus-0
  curve**, hence rational — the complete solution is a two-parameter family from
  the **chord construction**. The integer symmetric matrix `M` (2× convention,
  `G(v)=v^T M v`) is congruently **diagonalised over Q** (mpq) only to decide
  Legendre solvability and find one witness: the diagonal form is reduced to a
  squarefree, pairwise-coprime, integer `a X^2+b Y^2+c Z^2` and decided by the
  general three-coefficient **Legendre** conditions (`-bc` a QR mod `a`, etc.).
  Anisotropic (definite, or Legendre fails) → the trivial-only
  `{{x->0,y->0,z->0}}` (a proof). Otherwise a witness is found within **Holzer's
  bound** and mapped back to a primitive integer point `P0` (only that *one*
  point is transformed back — the family is built directly in the original
  coordinates, so it never inherits the diagonalisation's denominators). The
  family, `G(V)·P0 − 2(P0^T M V)V` over `V = C[1] e_i + C[2] e_j` with scale
  `C[3]`, plus the tangent line `C[1]·P0`, is the **complete projective**
  solution — every branch is a solution and every integer solution is a rational
  multiple of a family point, exactly the representation Mathematica / sympy
  return. So `Solve[4 x^2 - 5 y^2 + z^2 == 0, {x,y,z}, Integers]`,
  `Solve[x^2 + 3 x y + 2 y^2 - z^2 == 0, …]` and the multi-representation
  `Solve[x^2 + y^2 == 65 z^2, …]` now solve. Degenerate (rank-deficient) forms,
  a witness beyond the Holzer box, and any linear/constant term decline.
- **Extended ("general") Pythagorean** (`si_solve_general_pythagorean`,
  `src/solve/solveint_pythag.c`). A homogeneous sum of `k >= 3` squares equal to
  a square, `x_1^2 + ... + x_k^2 == y^2` (so `n >= 4` variables; the `k <= 2`
  cases are the binary / ternary solvers). The cone has the rational point
  `(1,0,...,0,1)`, giving the standard stereographic parametric family
  (parameters `C[1..k]`): `x_i -> 2 C[i] C[k]` (`i < k`),
  `x_k -> Σ_{i<k} C[i]^2 - C[k]^2`, `y -> Σ_{i<k} C[i]^2 + C[k]^2`. One family is
  emitted, matching Mathematica / sympy. So `Solve[x^2 + y^2 + z^2 == w^2,
  {x,y,z,w}, Integers]` and `Solve[x^2 + y^2 + z^2 + w^2 == v^2, …]` solve.
  Weighted summand coefficients decline (a documented follow-up).
- **General binary quadratic — parabolic and hyperbolic families**
  (`si_solve_bqf_parametric`, `src/solve/solveint_bqf_parametric.c`). The two
  unbounded conic types the diagonal Pell paths miss because of a cross term
  `B x y`, for `A x^2 + B x y + C y^2 + D x + E y + F == 0` with
  `delta = B^2 - 4 A C`:
  - **Parabolic** (`delta == 0`): the quadratic part is a perfect square, so the
    curve is a parabola with integer points given by a finite union of
    one-parameter families `{x -> quad(t), y -> quad(t)}`, found by the classical
    congruence method (a residue `u mod |_c|` per family, `_c = sb·sqc·D −
    sqa·E`). So `Solve[x^2 - 4 x y + 4 y^2 - 3 x == 0, {x,y}, Integers]` returns
    two families, and the degenerate parabola `Solve[y == x^2, {x,y}, Integers]`
    the single family `x = t, y = t^2`. Exhaustive over residues, so an empty
    result is a proof.
  - **Hyperbolic** (`delta > 0`, not a perfect square, `D == E == 0`,
    `x>0 && y>0`): a Pell conic whose automorphism group is generated by the
    fundamental unit `(t,u)` of `X^2 - delta Y^2 = 1` acting as the integer,
    determinant-1 matrix `M = [[t - B u, -2 C u],[2 A u, t + B u]]`. Reducing
    `A x^2 + B x y + C y^2 + F == 0` to `X^2 - delta Y^2 == N` (`X = 2A x + B y`,
    `Y = y`, `N = -4 A F`), the finitely many orbit-**base** solutions are found
    by scanning `Y` over the **Nagell bound** `Y <= u*sqrt(|N|/(2(t±1)))` — which
    is INDEPENDENT of the (possibly enormous) fundamental-unit size — mapping each
    `(X,Y)` representative back to `(x,y)` and advancing it into the positive
    quadrant. Each base yields one closed-form `ConditionalExpression` family in
    the fundamental unit (`C[1] >= 0`), exactly as the diagonal Pell path does. So
    `Solve[x^2 - 3 x y + y^2 == 1 && x>0 && y>0, {x,y}, Integers]` returns the
    complete set (six families reproducing `(1,3),(3,1),(3,8),(8,3),…`), and the
    large-discriminant `x^2 - 9 x y + 5 y^2 == 1` (`delta = 61`, fundamental unit
    ≈ 3.5×10⁹) and `x^2 - 1001 x y + 500 y^2 == 1` (`delta = 10⁶`) — which a direct
    `(x,y)` search cannot bound — are solved in sub-millisecond. Perfect-square
    `delta > 0` (rational asymptotes, finite) and `delta < 0` (ellipse, finite)
    are handled by the bounded conic / elliptic paths; general linear `D, E` terms
    in the hyperbolic case decline.
- **Ramanujan–Nagell-type exponential** `x^2 + D == 2^n`
  (`si_solve_ramanujan_nagell`, `src/solve/solveint_rn.c`), handled before the
  MPoly stage (variable exponent). For the class-number-1 imaginary-quadratic
  case — `D` squarefree, `D == 7 (mod 8)` (so 2 splits in the half-integer ring
  `O_K = Z[(1+sqrt(-D))/2]`), and `h(Q(sqrt(-D))) == 1`; for base 2 exactly
  `D = 7` — the factorisation `(x+sqrt(-D))/2 = +- alpha^(n-2)` forces
  `U_{n-2} = +-1` for the Lucas sequence `U_m = U_{m-1} - Q U_{m-2}`,
  `Q = (1+D)/4`. The **Bilu–Hanrot–Voutier primitive-divisor theorem** bounds
  `n <= 32`, so an exact perfect-square scan of that finite window returns the
  complete set: `Solve[2^n - 7 == x^2 && n > 0 && x > 0, {n, x}, Integers] ->
  {(3,1),(4,3),(5,5),(7,11),(15,181)}`. A Lucas cross-check guards the answer;
  an empty result within the gate is a proof, and out-of-scope inputs (base `!=
  2`, `D` not `== 7 mod 8`, `h != 1`) DECLINE (the general linear-forms-in-logs
  route B is a documented future extension, never a guess).
- **Linear Diophantine.** A single **linear** equation is solved through its
  solution lattice (gcd staircase, particular solution + `(n-1)`-vector
  homogeneous basis):
  - Unconstrained -> the full **parametric family**
    `{{x_i -> x0_i + sum_j basis[j][i] C[j+1]}}` with `C[k]` integer
    parameters (`Solve[x + y == 10, {x, y}, Integers]` ->
    `{{x -> C[1], y -> 10 - C[1]}}`); an unsolvable equation
    (`gcd(a)` does not divide `b`) gives `{}`.
  - Over a finite box, an unsolvable equation is reported as `{}` from the gcd
    test (so a large no-solution box is instant).  A solvable box is enumerated
    through the **LLL-reduced solution lattice** (`LatticeReduce`): the search
    is over the coefficient box obtained by projecting the value box through
    the lattice pseudoinverse, so it is small exactly when the coefficients are
    large (few solutions) — e.g. `1000003 x + 999983 y == 7 && Abs[x] < 10^9 &&
    Abs[y] < 10^9` returns its 2000-point arithmetic progression.  A box whose
    dense lattice would yield an intractable number of points is left
    unevaluated rather than truncated.

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
- Tests: [`tests/test_solve_corpus.c`](https://github.com/stblake/mathilda/blob/main/tests/test_solve_corpus.c)

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
