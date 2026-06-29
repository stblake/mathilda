# Solutions of Equations

The equation solver and its supporting machinery: `Solve` and `SolveAlways`, the algebraic-number representation `Root` and its radical conversion `ToRadicals`, and the cubic/quartic closed-form helpers (`Cubics`, `Quartics`). Related options documented elsewhere include `GeneratedParameters`, `InverseFunctions`, `VerifySolutions`, and `Eliminate`.

## Solve

Attempts to solve an equation or system of equations for one or more variables.
- `Solve[expr, vars]`: Solve `expr` for `vars` over the complex numbers (default).
- `Solve[expr, vars, dom]`: Solve over the domain `dom`. Supported: `Complexes` (default), `Reals`, `Integers`.

**Features**:
- `Protected`.  Matches Mathematica's attribute set -- arguments are
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
    system is not affine in the variables, in which case the router leaves
    `Solve` unevaluated.
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

**Options**:
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
Out[6]= {{x -> 1}, {x -> -1}, {x -> 2}, {x -> -2}}

In[7]:= Solve[x^3 + x + 1 == 0, x]
Out[7]= {{x -> Root[Function[1 + #1 + #1^3], 1]}, ...}

In[8]:= Solve[Sin[x] == 0, x]
Out[8]= Solve[Sin[x] == 0, x]

In[9]:= Solve[a/x + b == 0, x]
Out[9]= {{x -> -a/b}}

In[10]:= Solve[1/(x-1) == 2, x]
Out[10]= {{x -> 3/2}}

In[11]:= Solve[x/(x-1) == 2/(x-1), x]
Out[11]= {{x -> 2}}             (* x = 1 dropped as extraneous *)

In[12]:= Solve[x^2 - 5 x + 6 == 0, x, Integers]
Out[12]= {{x -> 2}, {x -> 3}}

In[13]:= Solve[x^2 - 2 == 0, x, Integers]
Out[13]= {}                     (* Sqrt[2] is not an Integer *)

In[14]:= Solve[1.5 x + 3 == 0, x]
Out[14]= {{x -> -2.0}}          (* approximate-in / approximate-out *)

In[15]:= Solve[{1.5 x + y == 4.5, x - y == 0.5}, {x, y}]
Out[15]= {{x -> 2.0, y -> 1.5}}

In[16]:= Solve[N[Pi, 50] x == 1, x]
Out[16]= {{x -> 0.31830988618379067...}}   (* 50-digit MPFR result *)

In[14]:= Solve[x^3 + 1 == 0, x, Reals]
Out[14]= {{x -> -1}}            (* real cube root, not (-1)^(1/3) *)

In[15]:= Solve[x^3 - 6 x^2 + 11 x - 6 == 0, x, Integers]
Out[15]= {{x -> 1}, {x -> 2}, {x -> 3}}

In[16]:= Solve[3 x + 2 y == 11 && x + y == 12, {x, y}]
Out[16]= {{x -> -13, y -> 25}}

In[17]:= Solve[a x + c == 1 && b x - d y == 2, {x, y}]
Out[17]= {{x -> (1 - c)/a, y -> (-2 a + b - b c)/(a d)}}

In[18]:= Solve[3 x + 2 y == 11, {x, y}]
Solve::svars: Equations may not give solutions for all "solve" variables.
Out[18]= {{y -> 11/2 - 3 x/2}}     (* x is free; only y has a rule *)

In[19]:= Solve[3 x + 2 y == 11 && x + y == 12 && 3 x + y == 32, {x, y}]
Out[19]= {}                        (* over-determined, inconsistent *)
```

## SolveAlways

Finds values of parameters (the symbols appearing in `eqns` but **not** in
`vars`) that make every equation in `eqns` hold for every value of `vars`.
The reduction is: each `lhs == rhs` is rewritten as the polynomial
`lhs - rhs`; `CoefficientList[lhs - rhs, vars]` exposes every coefficient
as a polynomial in the remaining symbols; every such coefficient must
vanish; the resulting system is fed to `Solve` with the parameters as
unknowns.
- `SolveAlways[eqns, vars]`

**Scope**:
- `eqns` may be `Equal[lhs, rhs]`, a `List[Equal[...], ...]`, or an
  `And[Equal[...], ...]`.  `True`/`False` sentinels arising from the
  evaluator's pre-pass on `==` are folded into the trivial answers.
- `vars` is a symbol or a `List` of symbols.
- The empty-parameter case (every symbol in `eqns` already appears in
  `vars`) returns `{}` — Mathematica's convention that there are no
  parameter values to report regardless of whether the polynomial is
  identically zero.
- `Unequal` (`!=`), `Or`-combinations of equations, radicals
  (`Sqrt[a x] == ...`), and `Series`-with-`O[x]^n` stripping are **not**
  handled in this version; those inputs will produce a `Solve`-level
  result on whatever coefficient system `CoefficientList` produces, which
  may not be the SolveAlways-correct answer.

**Diagnostics**:
- `SolveAlways::argt` — wrong number of arguments.
- `SolveAlways::eqf` — `eqns` contained a non-`Equal` element.
- `SolveAlways::ivar` — `vars` was not a symbol or non-empty list of
  symbols.

**Examples**:

```
In[1]:= SolveAlways[a x + b == 0, x]
Out[1]= {{b -> 0, a -> 0}}

In[2]:= SolveAlways[(a + b) x + (a - b) y == 0, {x, y}]
Out[2]= {{a -> 0, b -> 0}}

In[3]:= SolveAlways[{a x + b == 0, c x + d == 0}, x]
Out[3]= {{b -> 0, a -> 0, d -> 0, c -> 0}}

In[4]:= SolveAlways[(a - b) x == 0, x]
Out[4]= {{b -> a}}
```

## Root


Held symbolic representation of an indexed root of a univariate polynomial.

- `Root[Function[t, p[t]], k]` — the `k`-th root of `p` (1-indexed).

**Canonical index `k`** (matches Mathematica):

1. **Real roots first**, ordered ascending by value.
2. **Complex roots** afterwards, ordered by `Re` ascending; ties broken by
   `|Im|` ascending; within a conjugate pair the negative-`Im` member comes
   first.

This is the convention used by both `Solve`'s emission and `N[Root[..]]`'s
numericalization, so `Root[f, 1]` always refers to the same root regardless
of how it was produced.

**Numericalization** — `N[Root[f, k]]` and `N[Root[f, k], prec]`:

The pipeline is companion-matrix QR → Sturm certificate → canonical sort →
Newton refinement → basin verification. Both real and complex roots are
returned as MPFR values (`Complex[MPFR, MPFR]` for complex). Failure modes:

- `Root::nonint` — polynomial has non-integer coefficients (deferred case).
- `Root::indx`   — `k` is outside `1..deg(p)`.
- `Root::conv`   — QR or Newton did not converge after one precision retry.

Examples:
```
In[1]:= N[Root[Function[#^3 - 2 # - 5], 1], 30]
Out[1]= 2.094551481542326591482386540580

In[2]:= N[Root[Function[#^3 + # + 1], 1], 20]    (* real root first *)
Out[2]= -0.68232780382801932737

In[3]:= N[Root[Function[#^3 + # + 1], 2], 20]    (* conj pair: -Im first *)
Out[3]= 0.34116390191400966368 - 1.1615414252683233453 I

In[4]:= N[Root[Function[#^3 + # + 1], 3], 20]
Out[4]= 0.34116390191400966368 + 1.1615414252683233453 I
```

## ToRadicals


Convert held `Root[Function[poly], k]` objects in an expression into
closed-form radical expressions.

- `ToRadicals[expr]`

**Features**:

- `Protected`.
- Closed-form radicals are always returned when the polynomial has degree
  at most four — linear (trivial), quadratic (`Sqrt`), cubic (Cardano), and
  quartic (Ferrari via the depressed quartic + resolvent cubic).
- Binomial Root objects `Root[Function[a #^n + b], k]` are reduced to
  radicals for any degree `n`, using the principal `n`-th root multiplied
  by `(-1)^(2 (k-1) / n)`.
- Other Root objects of degree ≥ 5 are returned unchanged — the system
  makes no attempt at decomposition or solvable-Galois detection (cf.
  Mathematica's note "ToRadicals cannot find them").
- The k-th radical root is selected to agree with `N[Root[poly, k]]`'s
  canonical ordering (real-first ascending, complex by `Re` / `|Im|` /
  negative-`Im` first) — each formula's natural emission order is
  numerically matched against `root_numericalize` at machine precision.
  When the polynomial carries parametric coefficients (no numericalisation
  possible), the natural per-formula index `k - 1` is used and the result
  is allowed to disagree with `expr` for some parameter values, matching
  Mathematica's `nongen` behaviour.
- Walks its argument recursively, so `Root[..]` nodes inside `List`,
  `Equal`, `Less`, `Greater`, `And`, `Or`, `Not`, `Implies`, ... thread
  automatically — every `Root` anywhere in the tree is processed
  independently and the surrounding structure is preserved.
- Idempotent: `ToRadicals[ToRadicals[expr]] === ToRadicals[expr]`, since a
  successful conversion produces an expression free of `Root[..]` nodes.

```mathematica
In[1]:= ToRadicals[Root[Function[#^2 + 3 # + 5], 1]]
Out[1]= 1/2 (-3 - I Sqrt[11])

In[2]:= ToRadicals[Root[Function[#^2 + 3 # + 5], 2]]
Out[2]= 1/2 (-3 + I Sqrt[11])

In[3]:= ToRadicals[Root[Function[#^5 - 2], 3]]
Out[3]= (-1)^(4/5) 2^(1/5)

In[4]:= With[{r = ToRadicals[Root[Function[#^4 + 3 #^3 - 5 #^2 - 7 # + 9], 1]]},
              Chop[N[r^4 + 3 r^3 - 5 r^2 - 7 r + 9, 30]]]
Out[4]= 0

In[5]:= ToRadicals[Root[Function[#^5 - # - 1], 1]]      (* non-binomial deg 5 *)
Out[5]= Root[#1^5 - #1 - 1 &, 1]

In[6]:= ToRadicals[Root[Function[#^2 - 2], 2] < 3]      (* threading *)
Out[6]= True
```

## Solve`SolveLinearSystem

The linear-system specialist invoked by `Solve` for multi-variable inputs
(`And` / `List` of equations, or a single equation paired with a multi-symbol
variable list).  Reachable directly via its context-qualified name when the
caller has already classified its input.
- `Solve`SolveLinearSystem[eqns, vars]`
- `Solve`SolveLinearSystem[eqns, vars, dom]`

**Features**:
- `Protected`.
- `eqns` may be a single `Equal[lhs, rhs]`, `And[Equal[...], ...]`, or
  `List[Equal[...], ...]`.  `vars` must be a `List` of distinct symbols.
- Each equation is canonicalised to `lhs - rhs` and `Expand`-ed, then
  asserted affine in `vars`: coefficient of each `var` must be free of the
  variables, and the residual after subtracting `sum_j coeff_j * vars[j]`
  must also be free of the variables.  Non-affine systems return `NULL`
  (caller leaves `Solve` unevaluated).
- The m x (n+1) augmented matrix is built with variable columns in
  **reversed order** (M[i][0] is the coefficient of `vars[n-1]`).  This is
  what produces Mathematica's `Solve::svars` convention for under-determined
  systems: left-to-right Gauss--Jordan then naturally pivots on the
  right-most variable first, leaving left-most variables free.
- Gauss--Jordan elimination with symbolic-pivot selection: among non-zero
  candidates in the current column, prefer a concretely-non-zero entry
  (`Integer`, `Rational`, `Real`) over a symbolic one.  A column whose
  entries all simplify to zero (via `Cancel[Together[.]]`) becomes a free
  variable.  After reduction, any zero row whose augmented column is
  non-zero is detected as an inconsistency.
- Output shape:
  - Unique solution: `{{v1 -> e1, v2 -> e2, ...}}` (rules in input order).
  - Inconsistent system: `{}`.
  - Under-determined system: `{{pivot_vars -> exprs in free vars}}`; emits
    `Solve::svars`; free variables produce no rule.
  - Empty equation list (`Solve[True, vars]`): `{{}}` (tautology).
- Domain filtering (post-pass):
  - `Integers`: every emitted rule's RHS must be a concrete `EXPR_INTEGER`;
    otherwise the entire solution is dropped (`{}`).
  - `Reals`: any RHS that syntactically contains a `Complex[_, _]` head is
    treated as non-real and the whole solution is dropped.
  - `Complexes` (default): no filter.

```mathematica
In[1]:= Solve`SolveLinearSystem[{x + y == 3, x - y == 1}, {x, y}]
Out[1]= {{x -> 2, y -> 1}}

In[2]:= Solve`SolveLinearSystem[{x + y == 0}, {x, y}]
Solve::svars: Equations may not give solutions for all "solve" variables.
Out[2]= {{y -> -x}}                (* x free *)

In[3]:= Solve`SolveLinearSystem[
            {x + y + z == 6, 2 x - y + z == 3, x + 2 y - z == 2},
            {x, y, z}]
Out[3]= {{x -> 1, y -> 2, z -> 3}}
```

## Solve`SolvePolynomialEquality

The polynomial-equality specialist invoked by `Solve`.  Reachable directly via
its context-qualified name when the caller has already classified its input.
- `Solve`SolvePolynomialEquality[lhs == rhs, var]`
- `Solve`SolvePolynomialEquality[lhs == rhs, var, dom]`

**Features**:
- `Protected`.
- Same algorithm and output shape as `Solve` for single polynomial equalities
  in one variable.  Does not parse options; the caller supplies them through
  the C-level entry point.

## Solve`SolveInverseFunctions

The inverse-function specialist invoked by `Solve` when the outermost
dependence on `var` is an elementary invertible head.  Reachable directly
via its context-qualified name when the caller has already classified its
input.
- `Solve`SolveInverseFunctions[lhs == rhs, var]`
- `Solve`SolveInverseFunctions[lhs == rhs, var, dom]`

**Features**:
- `Protected`.
- Recognised heads: `Log`, `Exp`, `Sin`, `Cos`, `Tan`, `Cot`, `Sec`, `Csc`,
  the hyperbolic counterparts (`Sinh`, `Cosh`, `Tanh`, `Coth`, `Sech`,
  `Csch`), their inverses (`ArcSin`, `ArcCos`, ..., `ArcCsch`), and
  `Power[g, n]` for integer `n >= 2`.  Also recognises `Power[E, g(x)]`
  as the canonical form of `Exp[g(x)]`.
- Additive-shift isolation pre-pass: equations of the form
  `c * head[g(x)] + free_of_var == 0` are reduced to
  `head[g(x)] == new_rhs` before head dispatch.
- Multi-branch heads introduce a fresh integer parameter `C[k]` (head
  controlled by the parent `Solve`'s `GeneratedParameters` option) and
  wrap each branch in `ConditionalExpression[..., Element[C[k], Integers]]`.
- Inverse heads (`ArcSin`, `ArcCos`, `ArcTan`) use a vertical-strip
  predicate on `Re[a]`/`Im[a]` matching Mathematica's principal-branch
  domain.
- Inner equations are solved by hand-off to
  `Solve`SolvePolynomialEquality` -> `Solve`SolveInverseFunctions`
  (depth-capped at 8) -> `Solve`SolveRadicalsEquality`.  The recursion
  bypasses `Solve` itself, so the parent's inexact-rationalisation pre-
  pass runs only once.
- Emits `Solve::ifun` to stderr on first multi-branch peel per call and
  on the `InverseFunction[head][rhs]` fallback.
- Does not parse options; the caller supplies them through the C-level
  entry point `solveinv_solve_inverse_equality`.  When called via the
  qualified builtin, defaults `InverseFunctions -> True` and
  `GeneratedParameters -> C` are used.

```mathematica
In[1]:= Solve`SolveInverseFunctions[Sin[x] == a, x]
Out[1]= {{x -> ConditionalExpression[Pi - ArcSin[a] + 2 Pi C[1],
                                     Element[C[1], Integers]]},
         {x -> ConditionalExpression[ArcSin[a] + 2 Pi C[1],
                                     Element[C[1], Integers]]}}

In[2]:= Solve`SolveInverseFunctions[Log[x^2 + 1] + 1 == 0, x]
Out[2]= {{x -> -I Sqrt[1 - 1/E]}, {x -> I Sqrt[1 - 1/E]}}
```

## Solve`SolveRadicalsEquality

The radicals-equation specialist invoked by `Solve` when the polynomial
specialist declines (because the equation contains `Sqrt[...]`, fractional
powers `x^(p/q)`, or nested radicals).  Reachable directly via its
context-qualified name when the caller has already classified its input.
- `Solve`SolveRadicalsEquality[lhs == rhs, var]`
- `Solve`SolveRadicalsEquality[lhs == rhs, var, dom]`

**Features**:
- `Protected`.
- Algorithm:
  1. Canonicalise the equation: compute `e = Numerator[Together[lhs - rhs]]`.
     This combines rational-radical inputs (e.g.
     `(x + Sqrt[x])/Sqrt[x] + Sqrt[x]/(x + Sqrt[x]) == 4`) into a single
     polynomial-style residual in `var` and the radicals it contains.
  2. Iteratively locate radical atoms `Power[base, p/q]` (q > 1) anywhere in
     the working system (main equation + accumulated side equations).  For
     each distinct base `g_i`, introduce a fresh generator `u_i` so that
     `u_i = g_i^(1/L_i)`, where `L_i` is the LCM of denominators of *all*
     exponents of `g_i` (so `Sqrt[x]` and `x^(1/4)` share a single
     generator `u = x^(1/4)` with `L = 4`).  Replace every `g_i^(p/q)`
     by `u_i^(p*L_i/q)` in every equation, and append the side equation
     `u_i^L_i - g_i == 0`.  Nested radicals are picked up automatically
     -- a fresh atom inside a previously substituted base becomes its
     own generator on the next iteration.
  3. Eliminate `u_1, u_2, ...` from the main equation by chained
     `Resultant_{u_i}(main, side_eq_i, u_i)` in introduction order, so
     each side equation contributes exactly one fresh generator and the
     end-result is a polynomial in `var` alone.
  4. Hand the eliminated polynomial to `Solve`SolvePolynomialEquality`.
  5. Verify every candidate by back-substitution into the *original*
     equation.  The residual is first evaluated numerically with `N[]`
     and rejected when its magnitude exceeds `1e-9`; only when the
     numerical pass cannot decide -- free parameters, removable-
     singularity `Indeterminate`, etc. -- does the verifier fall back
     to a symbolic `Simplify` pass to catch structural zeros.  This
     ordering matters for candidates with algebraic coefficients
     (e.g. `Sqrt[2]` in the elimination): `Simplify` on the back-
     substituted residual can run for seconds per candidate, while
     `N[]` evaluates the same residual in microseconds.  Candidates
     whose residual still depends on free parameters (and so cannot
     be decided either way) are kept and trigger `Solve::nongen`,
     matching Mathematica's convention.
- Output shape matches `Solve`SolvePolynomialEquality`: a `List` of
  singleton-rule `List`s, plus the empty `List[]` when no candidate
  survives verification.  The `dom` argument flows through to the
  polynomial specialist (so `Reals` filters the candidate polynomial
  via the same per-degree discriminant tests as the polynomial path).
- The substitution introduces fresh generator symbols whose names follow
  the template `$radu<n>$`.  They are local to the call -- they never
  appear in the returned solution list (the resultant elimination
  removes every generator).
- The verifier accepts `Root[poly, k]` candidates without further
  checks: the polynomial elimination is exact, and `Root[]` objects
  describe the unique algebraic root of an irreducible factor that
  is not amenable to back-substitution.  Reflects Mathematica's
  policy of keeping `Root[]`-form solutions when they cannot be
  further simplified.
- The substitution-then-elimination strategy is "complete up to
  verification": every actual solution survives if it is closed-form,
  while extraneous roots from cross-multiplication (Together) or
  L-th-root branching are filtered out at the verifier.

```mathematica
In[1]:= Solve[Sqrt[x] + 3 x == 5, x]
Out[1]= {{x -> 1/18 (31 - Sqrt[61])}}

In[2]:= Solve[Sqrt[x] + 3 == 5, x]
Out[2]= {{x -> 4}}

In[3]:= Solve[x - 8 Sqrt[x] + 15 == 0, x]
Out[3]= {{x -> 9}, {x -> 25}}

In[4]:= Solve[Sqrt[x] + 3 x^(1/4) == 5, x]
Out[4]= {{x -> 1/2 (311 - 57 Sqrt[29])}}

In[5]:= Solve[(x + Sqrt[x])/Sqrt[x] + Sqrt[x]/(x + Sqrt[x]) == 4, x]
Out[5]= {{x -> 2 (2 + Sqrt[3])}}

In[6]:= Solve[Sqrt[x + 1] + Sqrt[x - 1] == 3, x]
Out[6]= {{x -> 85/36}}

In[7]:= Solve[Sqrt[x + 5] + Sqrt[x] == -1, x]
Out[7]= {}

In[8]:= Solve[x + Sqrt[x - 1] == 1, x]
Out[8]= {{x -> 1}}

In[9]:= Solve[Sqrt[a x + c] + 3 x == 5, x]
Solve::nongen: There may be values of the parameters for which some or
              all solutions are not valid.
Out[9]= {{x -> 1/18 (30 + a - Sqrt[60 a + a^2 + 36 c])},
         {x -> 1/18 (30 + a + Sqrt[60 a + a^2 + 36 c])}}
```

## Cubics

Option for `Solve` that controls whether cubic equations are solved via
explicit radical formulas.
- `Cubics -> False` (default): emit held `Root[]` objects.
- `Cubics -> True`: emit closed-form Cardano radicals.

## Quartics

Option for `Solve` that controls whether quartic equations are solved via
explicit radical formulas.
- `Quartics -> False` (default): emit held `Root[]` objects.
- `Quartics -> True`: emit closed-form radicals via Ferrari's resolvent-cubic
  method (Complexes only; a `Reals` request still yields `Root[]`, since a
  faithful real/non-real split of quartic radicals is not attempted).
