# Solutions of Equations

The equation solver and its supporting machinery: `Solve` and `SolveAlways`, the algebraic-number representation `Root` and its radical conversion `ToRadicals`, and the cubic/quartic closed-form helpers (`Cubics`, `Quartics`). Related options documented elsewhere include `GeneratedParameters`, `InverseFunctions`, `VerifySolutions`, and `Eliminate`.

## Solve

Attempts to solve an equation or system of equations for one or more variables.
- `Solve[expr, vars]`: Solve `expr` for `vars` over the complex numbers (default).
- `Solve[expr, vars, dom]`: Solve over the domain `dom`. Supported: `Complexes` (default), `Reals`, `Integers`.

**Features**:
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
- **Diophantine solving (`Integers` with constraints).** When the input is a
  polynomial equation (or system) conjoined with inequality / ordering /
  disequation constraints, a dedicated pre-pass (`src/solveint.c`) finds *all*
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
    `4/2027 == 1/x + 1/y + 1/z && 0 < x <= y <= z`.
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
  of declining.
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
`Solve[x^2 == 4, x, Rationals]` -> `{{x -> -2}, {x -> 2}}`,
`Solve[x^2 == 2, x, Rationals]` -> `{}`).  `Algebraics`, `Booleans`, and
`Primes` are not yet wired and leave `Solve` unevaluated.

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

## Solve`SolveNonlinearSystem

The nonlinear polynomial-system specialist invoked by `Solve` when the
linear-system specialist declines a non-affine system.  Reachable directly via
its context-qualified name.
- `Solve`SolveNonlinearSystem[eqns, vars]`
- `Solve`SolveNonlinearSystem[eqns, vars, dom]`

**Features**:
- `Protected`.
- `eqns` may be a single `Equal[lhs, rhs]`, `And[Equal[...], ...]`, or
  `List[Equal[...], ...]`.  `vars` must be a `List` of distinct symbols.
- Each equation `lhs == rhs` is put over a common denominator (`Together`); the
  **`Numerator`** is taken as the polynomial equation and a non-constant
  **`Denominator`** is recorded so denominators clear for rational systems
  (`1/x + 1/y == 1`).  The numerator is then **`Expand`-ed**, so products and
  powers of sums (`(x - 1)^2 + y^2 == 1`, `(x + y)(x - y) == 0`) distribute to a
  sum of monomials before Gröbner conversion — the GBPoly single-term parser
  cannot ingest `Power[Plus, k]` / `Times[Plus, ...]` directly.  Every numerator
  must then be a polynomial over Q in `vars` (a transcendental head, a radical /
  non-integer power, or a foreign symbol makes the specialist decline -> `NULL`).
  The numeric analogue `NSolve` expands identically.  For a pure polynomial the
  denominator is `1` and this reduces to `Expand[lhs - rhs]`.
- **Spurious-root pruning:** a completed tuple that drives any recorded
  denominator provably to zero is dropped (a den-zero point is a spurious root
  introduced by clearing).  `{(x - y)/(x + y - 2) == 0, x y == 1}` returns only
  `{{x -> -1, y -> -1}}` — the `{1, 1}` candidate zeroes `x + y - 2`;
  `{1/(x - 1) == 1/(y - 1), x + y == 2}` returns `{}` (its only candidate is a
  pole).  `Together` cancels removable factors, so `(x^2 - 1)/(x - 1)` is treated
  as `x + 1`.
- **Parametric systems:** a free symbol that is neither a solve variable nor a
  known constant (`Pi`, `E`, …) is treated as a **parameter**.  When any are
  present, the lex Gröbner basis is computed over the field `Q(parameters)` by
  reusing `GroebnerBasis`'s `CoefficientDomain -> RationalFunctions` engine, and
  the same triangular back-substitution runs with the symbolic univariate solver.
  `{x + y == a, x y == b}` → `{{x -> (a - Sqrt[a^2 - 4 b])/2, …}, …}`;
  `{x^2 == a, x + y == 0}` → `x -> ±Sqrt[a]`.  A basis generator free of all solve
  variables is a coefficient-field constant: a nonzero number → `<1>` → (generic)
  inconsistency → `{}`; a parameter-dependent one is a consistency/nondegeneracy
  condition and is **declined** (unevaluated) rather than guessed — case-splitting
  is reserved for `Reduce`.  No `ConditionalExpression` nondegeneracy guards are
  emitted (matching `Solve`, not `Reduce`).  More than three distinct parameters
  declines.
- A lexicographic Gröbner basis is computed via the Gröbner walk
  (`gb_groebner_walk`).  For a zero-dimensional ideal this basis is
  triangular: the univariate generator in the last variable is solved with
  `Solve`SolvePolynomialEquality` (forwarding the domain and the
  `Cubics` / `Quartics` flags), each root is back-substituted, and the search
  recurses variable by variable.  Every completed tuple is verified against
  the original equations (via `zero_test_decide`) before being accepted;
  duplicate tuples are removed.
- Output shape:
  - Finite solutions: `{{v1 -> a1, ...}, {v1 -> b1, ...}, ...}`.
  - Provably inconsistent system (unit ideal, or no domain points): `{}`.
  - Empty equation list (`Solve[True, vars]`): `{{}}` (tautology).
- Correctness policy: an empty `{}` is returned **only** when the ideal is
  provably inconsistent or a fully-solved zero-dimensional ideal has no
  solutions in the requested domain.  A branch on which no triangular
  generator can be found, or on which the univariate solver cannot reduce the
  generator, marks the search incomplete and the specialist returns `NULL`
  (Solve stays unevaluated) rather than emit a false `{}`.
- **Positive-dimensional ideals** (not every variable owns a pure-power
  leading monomial -> infinitely many solutions) are detected and left
  unevaluated with the advisory `Solve::nsdim`.  Solving in terms of free
  variables is reserved for future work.
- Domain filtering: `dom = Complexes` (default), `Reals`, or `Integers`.
  The domain is forwarded to the univariate step at every level, and a final
  per-tuple pass drops any solution with a non-real `Complex[_, _]` RHS
  (`Reals`) or a non-integer RHS (`Integers`).
- Safety: a per-generator total-degree gate (60) declines systems that could
  make the Gröbner computation blow up.

```mathematica
In[1]:= Solve[x y == 1 && x + y == 3, {x, y}]
Out[1]= {{x -> (1/2) (3 - Sqrt[5]), y -> (1/2) (3 + Sqrt[5])},
         {x -> (1/2) (3 + Sqrt[5]), y -> (1/2) (3 - Sqrt[5])}}

In[2]:= Solve[x^2 + y^2 == 1 && x == y, {x, y}]
Out[2]= {{x -> -1/Sqrt[2], y -> -1/Sqrt[2]}, {x -> 1/Sqrt[2], y -> 1/Sqrt[2]}}

In[3]:= Solve[x y == 6 && x + y == 5, {x, y}, Integers]
Out[3]= {{x -> 3, y -> 2}, {x -> 2, y -> 3}}

In[4]:= Solve[x^2 - y^2 == 0, {x, y}]
Solve::nsdim: The solution set is not zero-dimensional (infinitely many
solutions); Solve returned unevaluated.
Out[4]= Solve[x^2 - y^2 == 0, {x, y}]
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
     matching Mathematica's convention.  **`Root[]` candidates are
     verified the same way, not exempted:** the elimination hands every
     branch of the substituted `u = base^(1/L)` to the polynomial solver,
     but only the branch matching the principal `Sqrt` / `x^(p/q)`
     satisfies the original equation, so `N[Sqrt[Root[..]] + ...]` is
     `~0` for the valid root and `O(1)` for the spurious complex ones.
     This is what makes `Solve[Sqrt[x] + 3 x^(1/3) == 5, x]` return the
     single valid root rather than all three roots of the resultant cubic.
     (When a `Root[]` cannot be numericalised — e.g. a `USE_MPFR=0` build —
     the numeric pass abstains and the candidate is kept.)
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
