# Solve over the Integers — status and roadmap

This is the **development-driving document** for Diophantine solving in
`Solve[eqns && constraints, vars, Integers]`. It is self-contained: a fresh
session should be able to read this file (plus the source pointers it names) and
continue the roadmap without any prior conversation. Tier 1 is done and shipped;
**Tier 2 and Tier 3 are the open work**, specified below in enough detail to
implement against.

Keep this file current: when a Tier item lands, move it to "Done", add its
file/function pointers, and add a held-out case (see §4).

---

## 1. Orientation — how the solver is put together

The integer solver is an **Integers-domain pre-pass** on `Solve`. `builtin_solve`
(`src/solve.c`) calls `solveint_solve_integer(expr, vars, dom)` (`src/solveint.c`)
*before* the ordinary Complexes/Reals dispatch. It returns one of:

- a **finite solution list** `{{x -> …, …}, …}` — concrete integer tuples;
- a **parametric family** — one or more `ConditionalExpression` tuples in `C[k]`
  (unbounded families);
- the **empty list `{}`** — a *proof* of no solutions (only when genuinely
  proven, e.g. an exhausted finite search or a divisor/parity argument);
- **`NULL`** — decline; `Solve` falls through / is left unevaluated.

The evaluator owns `res`; a builtin returns a new `Expr*` or `NULL` (never frees
`res` itself). All the closed-form methods below are C functions over a shared
`SICtx` (parsed variables, equations as sparse integer `MPoly`, bound arrays,
constraint store) and a `SearchState` (candidate accumulator).

### Dispatch order inside `solveint_solve_integer`

Methods are tried in this sequence; the first that handles the input wins
(returns non-NULL). Order matters — specific/cheap/closed-form before general
search.

| # | Function | Handles |
|--:|---|---|
| 1 | `si_solve_exponential` | variable exponents `x^a - y^b == …` (Catalan/Mihailescu) |
| 2 | `si_solve_bounded_powerleaf` | non-poly leaf `n! + 1 == m^2`, … |
| — | `classify_conjunct` (Stage A) → `derive_bounds` / `derive_even_only_bounds` (Stage B) | split eqns/constraints; propagate a finite box |
| 3 | `si_solve_linear_parametric` | one linear equation → gcd-staircase `C[k]` family |
| 4 | `si_solve_pell_parametric` | unbounded `x^2 - D y^2 == 1`, `x,y>0` |
| 5 | `si_solve_genpell_parametric` | unbounded `x^2 - D y^2 == N` (any `N != +1`), `x,y>0` |
| 6 | `si_solve_linear_system_ray` | homogeneous linear system + positivity → ray |
| 7 | `si_solve_linear_system_hnf` | general linear system (m≥2 eqns) via HNF |
| 8 | `si_solve_power_sum_equal` | Prouhet-Tarry-Escott → `{}` via Newton's identities |
| 9 | `si_solve_mordell` | unbounded imaginary Mordell `y^2 == x^3 + k` (class-number) |
| 10 | `si_try_special_forms` → `si_solve_pell` / `si_solve_conic` / `si_solve_factorable_conic` / `si_solve_reciprocal` / `si_solve_linelim_bilinear` | bounded Pell; conic; factorable BQF; Egyptian fractions; bilinear divisor |
| 11 | `si_solve_multileaf` | staged elimination (Euler brick, …) |
| 12 | *decline if `n_unbounded >= 2`* | hand back to later phases / unevaluated |
| 13 | `mitm_solve` / `search_rec` / `si_solve_linear_bounded` / `si_solve_powersum_divisor` | the bounded engine (meet-in-the-middle, leaf search, LLL box, divisor method) |

### The two correctness invariants (do not break these)

1. **Never emit an unproven `{}`.** `{}` means "provably no solution". If a
   method cannot prove emptiness, it returns `NULL` (decline → unevaluated), not
   `{}`. Backstops in `src/solve.c`: `solution_set_is_parametric` converts any
   post-filter empty set that came from a positive-dimensional parametric answer
   back to unevaluated; `solvelinsys.c` declines underdetermined integer systems
   rather than integer-filtering a rational family to `{}`.
2. **A decline is always safe; a wrong answer is not.** Prefer unevaluated over
   guessing. The held-out gate (§4) enforces this: a wrong `{}`/finite/parametric
   answer fails CI; a decline is merely a reported gap.

### Output conventions

- Parameters are `C[k]`: `expr_new_function(mk_sym("C"), (Expr*[]){ mk_int(k) }, 1)`.
- Unbounded families wrap each rule RHS in `ConditionalExpression[val, cond]`
  where `cond` is e.g. `C[1] >= 0`. See `si_solve_genpell_parametric` /
  `si_genpell_family` for the canonical construction (a `Sqrt[D]`-based orbit).
- Positive-orthant convention: unbounded Pell-type families engage only with
  `x > 0 && y > 0` present (matching `si_solve_pell_parametric`), keeping the
  emitted family single-signed and finite-per-`C`.

### Where things live

- `src/solveint.c` — the whole integer engine (entry + every `si_*` method).
- `src/solve.c` — `builtin_solve`, the Integers pre-pass wiring, the spurious-`{}`
  guards.
- `src/solvelinsys.c` — the Complexes/Reals linear-system specialist (declines
  the underdetermined integer case).
- `src/linalg/hnf.c` / `hnf.h` — `linalg_hnf` (row HNF + unimodular transform) and
  the `HermiteDecomposition` builtin. New linalg `.c` files must be added to the
  test `COMMON_SRC` in `tests/CMakeLists.txt`.
- `src/poly/mpoly.h` — sparse integer polynomials (`mpoly_get_coef`,
  `mpoly_total_deg`, `mpoly_deg_var`, term arrays `eq->exps` / `eq->coefs`).
- Tests: `tests/test_solve_integers.c` (string-compare + property-based via
  `run_test`), plus `test_hermite_decomposition` in `tests/test_latticereduce.c`.
- Docs: `docs/spec/builtins/solutions-of-equations.md` (per-method reference),
  weekly changelog `docs/spec/changelog/<Mon>.md`.

---

## 2. Done — Tier 1 (shipped, v0.067–0.069)

| Item | Method(s) | Notes |
|---|---|---|
| **P0a** correctness | `solution_set_is_parametric` (solve.c), solvelinsys decline | underdetermined linear system no longer returns a silent `{}` |
| **P0b** HNF linear systems | `HermiteDecomposition` builtin (`src/linalg/hnf.c`), `si_solve_linear_system_hnf` | `A x = b` over ℤ: particular soln by forward-substitution + divisibility (fail ⇒ `{}` proof), kernel lattice as `C[k]` |
| **P1** Runge / factorable BQF | `si_solve_factorable_conic` | 2-var degree-2, cross term, perfect-square discriminant → difference-of-squares divisor enumeration (exhaustive ⇒ `{}` proof) |
| **P2** general-N Pell | `si_solve_genpell_parametric` (+ `si_genpell_detect/_bases/_family`) | unbounded `x²-Dy²==N`, any `N≠+1` incl. negative Pell; Nagell fundamentals + ε-orbit, one family per class; unsolvable ⇒ `{}` proof |
| **Systemic** | `benchmarks/87/{heldout.py,validate.py}` | the held-out silent-wrong-answer gate (§4) |

Design memories (in the session memory store): `project_solve_integers_linear_system_hnf`,
`project_solve_genpell_parametric`, `feedback_benchmark_overfit_needs_heldout_corpus`.

---

## 3. The two families of coverage, and where the gaps are

Mathematica's `Solve`/`Reduce` over the integers is an algorithmic dispatcher of
~25 engines. Mathilda already owns the **bounded-search + closed-form** half
strongly. Every open item is in the **unbounded / algebraic-number-theory** half.

Coverage of the classical taxonomy:

- Linear: single eq (gcd staircase) ✅; systems (HNF) ✅; LLL for huge coeffs ✅;
  ILP/branch-and-bound — only bounded enumeration, no cutting planes.
- Quadratic: Pell `±1` (CF) ✅; general-N Pell ✅; factorable BQF ✅; conic
  `Y²=AX²+BX+C` ✅. **Gap (Tier 2 D):** a *unified* general binary quadratic
  covering the hyperbolic (cross+linear terms → reduce to Pell) and elliptic
  (Δ<0, bounded) cases in one classified solver. **Gap (Tier 2 E):** ternary
  quadratic / Legendre solvability (Hilbert symbols) with a proof + witness.
- Thue (irreducible homogeneous degree ≥3 = const): **Gap (Tier 3 F)** —
  Tzanakis–de Weger (Baker bounds + LLL). Needs number-field arithmetic.
- Elliptic / Mordell: imaginary Mordell class-number method ✅; bounded leaf ✅.
  **Gap (Tier 3 G):** general integral points via Mordell-Weil + descent +
  elliptic logs. Largest item.
- Systems / higher: staged elimination ✅; finite-field Gröbner ✅. Possible
  future: Gröbner over ℤ (or modular+CRT) for general zero-dimensional integer
  systems. CAD for guaranteed real fences (Tier 3 H, orthogonal).

---

## 4. The validation workflow (drive every change through this)

**`make check-diophantine-heldout`** runs Mathilda cold on
`benchmarks/87-diophantine-integers/heldout.py` (equations from standard
references, none in the developed-against `cases.py`) and cross-checks each
answer against an independent Python brute-force oracle over the same box.
Verdicts: **OK** (matches), **DECLINE** (unevaluated where acceptable — a gap,
not a failure), **WRONG** (a `{}`/finite/parametric answer the oracle
contradicts → nonzero exit). Report: `HELDOUT_REPORT.md`.

**Whenever you add or change a method:**

1. Add held-out cases in `heldout.py` — include at least one *solvable* case, one
   *unsolvable* case (must return `{}`), and one *out-of-scope* case that must
   decline. For boxes too big for a naive grid (elliptic, Pell windows), give the
   case an equation-aware `oracle` (see the `_oracle_*` helpers).
2. Run the gate. A new WRONG is a real bug; a new DECLINE is expected for
   research-grade inputs.
3. Add C unit tests in `tests/test_solve_integers.c`. For parametric output,
   assert **properties**, not the giant FullForm string: class count
   (`Length[...]`), fundamentals (`{x,y} /. sol[[i]] /. C[1] -> 0`), and that the
   equation holds under substitution (`Simplify[(lhs) /. (sol[[i]] /. C[1] -> k)]`).
4. During development, validate completeness against a throwaway Python brute
   force *before* trusting the C — this is how P2's `ε⁻¹`-reduction bug was found.
   Compare as **lists/counts**, not sets: set-equality hides duplicate families.
5. Update `docs/spec/builtins/solutions-of-equations.md` and the weekly changelog.
6. `make check-c99`; valgrind delta vs the macOS startup baseline must be 0.

---

## 5. Tier 2 — unify and broaden the quadratics

### D. Unified general binary quadratic  `A x² + B x y + C y² + D x + E y + F == 0`

**Goal.** One classified solver, `si_solve_binary_quadratic`, that solves *any*
2-variable degree-2 integer equation completely — finite set, full parametric
family, or a proven `{}` — subsuming `si_solve_conic`, `si_solve_factorable_conic`,
`si_solve_pell*`, and `si_solve_genpell_parametric`. This is sympy's
`diop_quadratic` / the Alpern–Matthews method. It closes the case where a
hyperbola carries a cross term *and* linear terms (not yet reduced to Pell) and
the rotated ellipse the interval bounder can't box.

**Algorithm — classify by the discriminant `Δ = B² − 4AC`:**

- **Δ is a positive perfect square** → the form factors into two rational lines:
  already handled by `si_solve_factorable_conic` (hyperbola with rational
  asymptotes, divisor enumeration). Route here.
- **Δ = 0 (parabolic).** New. The form is a perfect square in the quadratic part
  plus linear terms. Standard reduction (sympy `_diop_parabolic`): with
  `g = gcd(A, B, C)` and the leading part `(√A x + √C y)²` (sign chosen so
  `B = 2√(AC)`), substitute `u = √A x + √C y`; the equation becomes quadratic in
  `u` with `u` linear in the other variable → solve a small number of residue
  classes; each yields a **linear** family in one parameter (a parabola has
  `O(√range)` points but they lie on finitely many arithmetic progressions).
  Emit `C[k]` families; `{}` when the residue test fails.
- **Δ > 0 non-square (hyperbolic).** New (the main content). Reduce the general
  BQF to a **Pell normal form** `X² − Δ Y² = 4 a M` and dispatch to the
  generalised-Pell machinery. Concretely (sympy `_diop_DN` route): multiply by
  `4A`, complete the square `X = 2A x + B y + D`, obtaining
  `X² − Δ y² + (…)y + (…) = 0`; a further completion in `y` (mirror
  `si_solve_factorable_conic`'s `U/V/W` reduction but with **non-square** Δ)
  gives `X'² − Δ Y'² = N'`. Then reuse `si_genpell_bases` / `si_genpell_family`
  to produce one `ConditionalExpression` family per class, mapping `(X', Y')`
  back to `(x, y)` by the (rational, integrality-checked) inverse substitution.
  Engage only with the positive-orthant / bound convention already used for Pell.
- **Δ < 0 (elliptic).** ✅ **DONE** (`si_solve_elliptic_bqf`, in
  `si_try_special_forms`). Treats the equation as a quadratic in `x` per fixed
  `y`, over the finite `y`-interval where the `x`-discriminant
  `Δ y² + (2BD−4AE)y + (D²−4AF)` (a downward parabola, Δ<0) is `≥ 0`; solves the
  integer quadratic in `x` exactly. Exhaustive ⇒ `{}` is a proof. Handles
  rotated forms, negative-definite (normalised), linear terms, constraints.
  Validated vs brute force over `|.|≤200`; `test_elliptic_bqf`; held-out cases
  `ellipse-*`. Closed the hole where a rotated ellipse escaped `derive_bounds`.

**Where to hook.** A single `si_solve_binary_quadratic(&c)` called from
`si_try_special_forms` *before* the individual conic/factorable/pell helpers (or
have it delegate to them for the sub-cases they already cover). Keep the existing
helpers; the new function is the classifier/router that also fills the Δ=0, Δ>0
non-square-with-cross-terms, and Δ<0-rotated gaps.

**Tests.** Parabolic (`x² − 4 x y + 4 y² − 6 x + y == 0`-style); hyperbolic with
cross+linear terms reduced to Pell (assert families satisfy the equation, count
classes, brute-force window); rotated ellipse (finite, `{}` when empty);
regression that the already-covered forms are unchanged. Held-out: one of each
class, plus an unsolvable of each.

**References.** sympy `sympy/solvers/diophantine/diophantine.py` (`_diop_quadratic`,
`_diop_DN`, `_diop_parabolic`, `transformation_to_DN`, `find_DN`); D. Alpern,
"Methods to solve `ax² + bxy + cy² + dx + ey + f = 0`"; Andreescu & Andrica,
*Quadratic Diophantine Equations*.

### E. Ternary quadratic / Legendre  `a x² + b y² + c z² == 0`

**Goal.** For a homogeneous ternary quadratic, **decide solvability** (a proof,
not a bounded search) and, when solvable, return a **primitive nontrivial
witness**; when not, return `{}`. Extends to the general ternary form
`Q(x,y,z) == 0` via diagonalisation. Today Mathilda only bounded-searches these
(finds points in a box, or an unproven empty over an unbounded domain).

**Algorithm.**

1. Normalise: make `a, b, c` squarefree and pairwise coprime (absorb square
   factors into the variables), signs not all equal (else only trivial over ℝ).
2. **Legendre's condition (solvability).** `a x² + b y² + c z² = 0` has a
   nontrivial integer solution iff `−bc` is a QR mod `|a|`, `−ca` is a QR mod
   `|b|`, and `−ab` is a QR mod `|c|` (plus the real sign condition). Implement
   with Jacobi/Legendre symbols (`JacobiSymbol` exists in `src/numbertheory/`).
   This is the Hasse–Minkowski / Hilbert-symbol content specialised to ℚ.
3. **Witness (when solvable).** Gauss's reduction / descent — reduce the
   determinant `|abc|` by successive substitutions until a solution is read off;
   or the Cremona–Rusin lattice method (`conic_param`) which uses `LatticeReduce`
   (already available) to return a small solution. Return a primitive `(x,y,z)`.
4. For a *general* (non-diagonal) ternary quadratic, diagonalise over ℚ first
   (symmetric-matrix congruence), solve, map back.

**Output.** `Solve[a x²+b y²+c z² == 0, {x,y,z}, Integers]` over an unbounded
domain: return the parametric family (a conic has a rational parametrisation once
one point is known — 2-parameter integer family), or `{}` proved by Legendre.
With a positivity/box constraint, intersect with the box.

**Hook.** New `si_solve_ternary_quadratic(&c)`, engaged for a single homogeneous
degree-2 equation in exactly 3 variables with no linear terms (after clearing).
Place among the special forms.

**Tests.** `3x²+5y²==7z²` → `{}` (Legendre fails; already a held-out DECLINE that
should become a proven `{}`); a solvable case e.g. `x²+y²==2z²` → witness +
family; diagonalisable general form. Held-out: solvable + unsolvable + boxed.

**References.** J. E. Cremona & D. Rusin, "Efficient solution of rational
conics", *Math. Comp.* 72 (2003); Cohen, *A Course in Computational Algebraic
Number Theory*, §5.3; sympy `diop_ternary_quadratic`, `descent`, `holzer`.

---

## 6. Tier 3 — deep algebraic-number-theory engines (research-grade)

These are **subsystems**, not functions. Each is weeks of work and depends on
number-field / elliptic-curve machinery Mathilda does not yet have. Until built,
the honest behaviour is **decline** (unevaluated) — and the held-out gate now
guarantees they can never silently degrade into wrong answers. Build the shared
prerequisite first.

### Prerequisite: an algebraic-number-field arithmetic layer

Needed by both F and G: ring of integers `O_K`, integral basis, unit group and
**fundamental units**, regulator, real+complex **embeddings**, ideal
factorisation, class group. Mathilda has FLINT (`src/poly/flint_bridge.c`,
`qqbar`, algebraic-tower reduce) as a foundation but no explicit number-field
unit-group / embedding API. Scope this as its own milestone; FLINT/ANTIC provide
much of it (`nf_t`, `nf_elem`), so a bridge may be feasible rather than a
from-scratch build.

> **Status (2026-08-19): built and validated for the monogenic case.**
> `src/numbertheory/numberfield.{c,h}` provides field setup, `arb`/`acb`
> embeddings, exact `disc`, and **Gate 1** — maximal-order certification via
> Dedekind's criterion (monogenic-first; non-maximal orders decline).
> `src/numbertheory/nfunits.{c,h}` provides **Gate 2** — fundamental units +
> regulator, certified by p-saturation (rank test of the mod-`p` character
> matrix). Regulators validated against LMFDB (`Q(2^1/3)`, the cyclic cubic
> `t^3-3t+1`, `Q(2^1/4)`). Round-2 order enlargement (non-monogenic fields) and
> class-group / ideal factorisation remain for G.

### F. Thue equations  `F(x, y) == m`  (F irreducible homogeneous, deg ≥ 3)

**Method (Tzanakis–de Weger).**
1. `K = ℚ(θ)`, `θ` a root of `F(x, 1)`; factor `F(x,y) = a₀ ∏ (x − θⁱ y)`.
2. Reduce the equation to unit equations in `O_K`; parameterise `x − θ y` by the
   fundamental units.
3. **Baker's linear forms in logarithms** → an explicit (astronomically large)
   upper bound on `max(|x|, |y|)`.
4. **LLL reduction** of the associated lattice shrinks the bound to a
   computationally exhaustible range.
5. Enumerate the reduced range; verify.

**Scope / staging.** Start with **cubic Thue** over a real cubic field (one
fundamental unit, simplest linear form). Then general degree. Emits a **finite**
solution set (Thue equations have finitely many solutions), so output is a plain
list — no parametric families. `{}` is a proof.

**Hook.** `si_solve_thue(&c)` for a single homogeneous irreducible 2-var equation
of total degree ≥3 equal to a constant. Currently such inputs decline (e.g.
`x³ − 2 y³ == 1` unbounded); bounded versions already work via the leaf search.

> **Status (2026-08-19): DONE for the monic `|m|=1` monogenic case.**
> `src/solvethue.{c,h}` implements the full Tzanakis–de Weger method (steps
> 1–5): field setup, reduce to unit equations, **Baker (Waldschmidt) initial
> bound + de-Weger LLL reduction** (`thue_exponent_bound`, `arb`/`acb` at 1600
> bits, reusing `lll_reduce_q`), the Q-dependent case (iii) via
> relation-detection + the L-trick, and exact reconstruction/verification.
> `si_solve_thue` (`src/solveint.c`) dispatches, so **`Solve[…, Integers]` now
> returns the complete set** for `x³−2y³=±1`, `x³−7y³=1`, `x³−3xy²+y³=1`,
> `x⁴−2y⁴=−1`, the Thomas cubic, etc. (`tests/test_thue.c`, incl. the automatic
> path). Out-of-scope inputs DECLINE safely: non-monogenic field (`x³−17y³=1`),
> `|m|≠1`, `|a₀|≠1`, degree/precision beyond reach.
> **M1 (2026-08-19): Voronoi units for large-regulator complex cubics.**
> `src/numbertheory/nfvoronoi.c` walks the chain of relative minima of `O_K` to
> propose a fundamental unit whose coordinates exceed the coefficient box
> (`ℚ(∛41)`: 24-digit unit, reg 56), certified by the same p-saturation gate
> (now mpz end-to-end). `x³−{15,41,42,97}y³=±1` now solve; benchmark 88 CORRECT
> 48→56.
> **M2 (2026-08-19): general `m` (`|m|≠1`) via µ-enumeration** (rank-1 complex
> cubic). `β=x−θy` is a norm-`m` integer `µ·unit`; `thue_norm_reps_cubic11`
> enumerates the bounded-norm reps µ, and `thue_exponent_bound` is now µ-aware
> (the linear form's constant gains `log|µ^(k)/µ^(j)|`; `C4`/`Y2p`/`V0` over-
> estimated). `x³−2y³={2,3,10}` etc. now return the full set, `={4,5,9,73,100}`
> proven `{}`; benchmark 88 CORRECT 56→65, 270-case PARI grid 0 WRONG.
> **M3 (2026-08-20): Round-2 maximal order + O_K-basis unit search** (cubics +
> quartics). `src/numbertheory/nfround2.c` computes `O_K` by Pohst–Zassenhaus
> (p-radical Frobenius kernel, ring of multipliers via HNF), returned as an
> integral basis `(1/D)W`; the unit search walks the O_K lattice `L={Σ c_i W[i]}`
> (`|N(v)|=Dⁿ`). Non-monogenic `Q(∛{10,12,17,19,20})`, the Dedekind cubic, and
> `Q(d^{1/4})` (incl. index-16 `Q(12^{1/4})`) now solve; benchmark 88 CORRECT
> 65→81, 130-case non-monogenic PARI grid 0 WRONG.
> **M5 (2026-08-20): totally complex fields (`r1=0`), any `m`.**
> `thue_solve_totally_complex` (`src/solvethue.c`): every root is non-real, so
> `|x−θᵢy| ≥ |Im θᵢ|·|y|` gives the elementary rigorous bound
> `|y| ≤ (|m|/∏|Im θᵢ|)^{1/n}` — no units/torsion/Baker; each `y` closed by exact
> univariate root-finding. Solves `Q(ζ₅)` cyclotomic quartic (6 pts),
> `x⁴+y⁴={1,2,17,82}`/`=3→{}`, `Φ₇`/`Φ₁₀`; benchmark 88 CORRECT 98→99.
> **Perf (2026-08-20):** the small-|Y| gap-closing brute box now uses exact
> univariate root-finding for a wide x-window (`x³−2y³=100`: 244→21 ms).
> **Validation (2026-08-20):** reproducible randomized grid
> `benchmarks/88-thue-equations/grid.py` (deg 3–6, mixed `m`) vs PARI, with an
> adjudicator that verifies disputed points directly — so it distinguishes a real
> bug from a PARI `thue()` incompleteness (found one on a `Q(ζ₅)` generator).
> Follow-on (the last coverage gap): **M4** — rank-≥2 units for `Q(10^{1/4})` and
> `x⁵−5y⁵`. **Full completion roadmap (algorithms, files, order, verification):**
> [`docs/design/thue_completion_plan.md`](docs/design/thue_completion_plan.md).
> Stress-tested vs PARI/GP `thue()` in `benchmarks/88-thue-equations/`.

**References.** N. Tzanakis & B. de Weger, "On the practical solution of the Thue
equation", *J. Number Theory* 31 (1989); Bilu & Hanrot; PARI/GP `thue`.

### G. Integral points on elliptic curves  `y² == x³ + a x + b`

**Method (largest item).**
1. **Mordell–Weil group** `E(ℚ) = ℤʳ ⊕ E_tors`: torsion via Nagell–Lutz /
   Mazur; **rank `r` and generators** via descent (2-descent: compute the
   2-Selmer group; this is the hard, sometimes-conditional part).
2. Bound integral points by **elliptic logarithms** + David's lower bound
   (linear forms in elliptic logs), then **LLL** to shrink, then enumerate.

**Scope.** Depends on a working rank/descent implementation — the deepest
dependency in this document. The imaginary-Mordell class-number method
(`si_solve_mordell`) already covers a slice for free; a full solver generalises
it to all `k` and to general Weierstrass curves. Output is a finite list; `{}` a
proof.

**Hook.** Extend `si_solve_mordell` / a new `si_solve_elliptic(&c)` for
`y² == cubic(x)` unbounded.

**References.** Stroeker & Tzanakis, "Solving elliptic Diophantine equations by
estimating linear forms in elliptic logarithms", *Acta Arith.* 67 (1994);
Cohen, *Number Theory Vol. I* §8; Cremona, *Algorithms for Modular Elliptic
Curves*; Sage `E.integral_points()`.

### H. Cylindrical Algebraic Decomposition (orthogonal)

CAD at the `Reduce` level yields **guaranteed real fences** around the real
solution set → integer search boxes for systems the interval bounder can't box.
Large, and mostly a *real*-solving feature that Solve/Integers would piggyback
on. Note as future; not on the critical path for the items above.

---

## 7. Recipe — adding a closed-form method

1. Write `static bool si_solve_<name>(SICtx* c, SearchState* st)` (emits into
   `st`, returns handled) **or** `static Expr* si_solve_<name>(SICtx* c)`
   (returns the family/`{}`/`NULL`) depending on whether output is a candidate
   set or a symbolic family. Read coefficients from `c->eq[q]` via `eq->exps` /
   `eq->coefs` (see `si_pell_detect`, `si_linear_detect` for the term-walk
   idiom); use `mpz_*` throughout.
2. Gate precisely: exact shape, and *decline* (`false`/`NULL`) on anything you
   cannot solve completely. Never emit an unproven `{}`.
3. Wire into `solveint_solve_integer` (or `si_try_special_forms`) at the right
   point in the dispatch order (§1).
4. Manage memory manually: every `mpz_init` paired with `mpz_clear`, every owned
   `Expr*` freed or returned. Reuse `mk_int/mk_mpz/mk_fn1/mk_fn2/mk_rule/mk_list`
   and `eval_and_free`.
5. Tests + held-out cases + docs + changelog + `check-c99` + valgrind (§4).

---

## 8. Suggested order of work

1. **Tier 2 D** (unified BQF) — highest value; finishes the quadratic story and
   reuses the P1/P2 machinery. Do the four Δ-classes incrementally (parabolic and
   Δ<0-rotated are self-contained; the Δ>0-hyperbolic reduction leans on P2).
2. **Tier 2 E** (ternary/Legendre) — self-contained, needs only Jacobi symbols +
   `LatticeReduce`, no number-field layer.
3. **Number-field arithmetic layer** — the shared Tier-3 prerequisite (evaluate a
   FLINT/ANTIC `nf_t` bridge first).
4. **Tier 3 F** (Thue, cubic first), then **G** (elliptic), then **H** (CAD) if
   ever.

Each step is independently shippable and independently gated by
`make check-diophantine-heldout`.
