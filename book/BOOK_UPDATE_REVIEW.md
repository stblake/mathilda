# The Mathilda Book — Update Review (git history since the last book pass)

> **STATUS: IMPLEMENTED 2026-09-25.** All backlog items below were written into the
> book and verified (`make examples`/`check-links`/`usage`/`pdf` clean, 240 pages).
> New sections §4.3 Calculus+DSolve, §4.4 Linear Algebra, §4.7 Special Functions,
> §4.9 Computational Geometry, §4.10 Graphs, and Ch. 6/7/9; UPDATEs to §4.2/§4.5/
> §4.6/§4.8. This file is kept as the record of what was done. (Two kernel bugs
> were found during the pass — `Simplify` on a polar Laplacian, and `DownValues`
> display — and are flagged, not fixed, being out of scope for a book update.)

This is a **triage backlog**, not prose. It maps every *book-relevant* capability
that landed in Mathilda **after the book's last substantive pass** onto the
existing writing campaigns in `ROADMAP.md`, so each campaign has a ready list of
what it must now cover. Read `CONTEXT.md` and `ROADMAP.md` first; this file is
scoped to *the delta*.

## Anchor and window

- **Book anchor:** commit `8c4332ba` (2026-08-28), Mathilda **v0.111/0.112** —
  "Statistics sections (§3.11 tour + §4.8 deep dive) + ChineseRemainder /
  CylindricalDecomposition reference pages". This is the last commit that touched
  `book/`.
- **Window reviewed:** `8c4332ba..HEAD` — **163 `src/` commits**, Mathilda
  **v0.113 → v0.189**, 2026-08-30 → 2026-09-25.
- **Sources of record:** the weekly changelogs `docs/spec/changelog/2026-08-24.md`
  (entries *above* its `## Book — Statistics` marker are post-anchor) through
  `2026-09-21.md`; the git log; cross-checked against `docs/spec/builtins/` and
  the book's current `\B{}` / `\usagebox` coverage.
- **Verification done for this file:** all 65 heads named below were confirmed
  registered in the live `./Mathilda` binary (v0.189) — no phantom functions — and
  every "not covered" claim was confirmed by `grep` over `book/chapters/`
  (all cited heads currently return 0 references).

## How to read the classification

| Tag | Meaning | Where it lands |
|-----|---------|----------------|
| **INCLUDE** | New capability, no book coverage. | A *Planned/stub* section — write it as part of that section's campaign. |
| **UPDATE** | Extends a chapter already marked **Verified** — its prose must be broadened. | §4.1, §4.2, §4.5, §4.6, §4.8 (the written math sections). |
| **REVIEW** | A behaviour/semantics change, not a new head. Transcripts auto-update on rebuild; surrounding *prose* must be re-checked. | Any Verified chapter whose claims the change touches. |

Effort: **S** = a few `.m` examples + a subsection · **M** = a full subsection or
several · **L** = a large multi-subsection section · **XL** = a chapter-scale
campaign in its own right.

---

## Priority-ranked summary

| # | Capability | Target campaign | Class | Effort | Landed |
|---|-----------|-----------------|-------|--------|--------|
| 1 | **DSolve** — the whole symbolic ODE/PDE subsystem | §4.3 Calculus *(needs new subsection — see note)* | INCLUDE | **XL** | v0.113–0.189 (M0–M59) |
| 2 | **`Integrate`** ParallelMixedTower + Cherry di/polylog + DiffUnderInt + Fermi–Dirac | §4.3 Calculus (§4.3.5 Risch story) | INCLUDE | **L** | v0.114–0.189 |
| 3 | **Number-field family** (`ToNumberField`, `AlgebraicNumber*`, `NumberFieldIntegralBasis`, `AlgebraicIntegerQ`) | §4.2 Algebra | UPDATE | M | 2026-09-14 → 09-21 (v0.16x–0.184) |
| 4 | **Gröbner / polynomial views** (`PolynomialReduce`, `MonomialList`, `CoefficientRules`, `FromCoefficientRules`, `FactorSquareFreeList`, order fix) | §4.2 Algebra | UPDATE | M | 2026-09-14 |
| 5 | **Reduce transcendental solving + QE/CAD extensions + `Simplify` over `Root`** | §4.2 Algebra (solving) | UPDATE | M | v0.113/0.124/0.125/0.165, 08-31 wk |
| 6 | **`JordanDecomposition` / `SchurDecomposition`** + `RowReduce`/`NullSpace` `ZeroTest` | §4.4 Linear Algebra | INCLUDE | M | v0.118 / v0.180–0.181 |
| 7 | **Statistics `Quantile` family** (`Quantile`, `InterquartileRange`, `MeanDeviation`, `MedianDeviation`) | §4.8 Statistics | UPDATE | **S** | 08-24 wk (post-anchor) |
| 8 | **Computational geometry** (`Area`, `Perimeter`, `RegionCentroid`, `RegionMember`, `ConvexHullRegion`) | *No home — propose* | INCLUDE | M | 08-24 wk (post-anchor) |
| 9 | **Graph additions** (`StarGraph`, `FindVertexColoring`, `EdgeWeight`, `WeightedAdjacencyMatrix`, weighted `FindShortestPath`/`GraphDistance`, `RandomGraph` count form) | *No home — propose* | INCLUDE | M | 08-24 / 08-31 wks |
| 10 | **I/O stream layer** (`Read`, `OpenRead/Write/Append`, `Write`, `WriteString`, `Close`, `Streams`, `StreamPosition`, `SetStreamPosition`, `ReadList`) | Ch. 9 Data I/O | INCLUDE | M | 2026-09-14 |
| 11 | **`LegendreQ`** | §4.7 Special Functions | INCLUDE | S | 2026-08-31 |
| 12 | **List/array ops** (`ArrayReshape`, `ArrayPad`, `ListGradient`, `FirstPosition` levelspec/default) | Ch. 6 Data Structures | INCLUDE | S | 08-24 / 08-31 / 09-14 |
| 13 | **Interpolation** (`ListInterpolation`, `Interpolation` default Options) | §4.5 Numerical Calculus *(or Ch. 6)* | UPDATE | S | 2026-08-31 |
| 14 | **`Inactive` / `Activate`** evaluation-control primitives | Ch. 7 Programming | INCLUDE | S | v0.175 |
| 15 | **Packed-array surface** (`PackedArrayQ`, `ToPackedArray`, `ToNDArray`) | Ch. 6 Data Structures | INCLUDE | S | 2026-08-31 |
| 16 | **`BitLength`** (seeds a Bit\* family; new `src/bitwise/`) | *No home — propose §4.1 or Ch. 6* | INCLUDE | S | 2026-09-14 |
| 17 | **Behaviour/semantics changes** (A1–A10 divergence fixes, `Function` closure, `Apart` 2-arg, `DSolve` not HoldAll) | REVIEW written chapters | REVIEW | S | 08-31 / 09-21 wks |

---

## Detail by campaign

### §4.3 Calculus — *stub (5 lines); the single biggest gap* — INCLUDE

The entire symbolic-ODE subsystem and the marquee integration engines were built
*after* the book. The section is a bare stub, so this is a from-scratch campaign;
DSolve alone is large enough to be its own multi-week effort.

**DSolve (items #1).** Confirmed live: `DSolve[y'[x]==y[x],y[x],x]` →
`{{y[x] -> C[1] E^x}}`. Method families to cover, all new (M0–M59):
- *First order:* separable, linear, exact, homogeneous, Bernoulli, Riccati, Chini,
  Abel, Clairaut, Lagrange/d'Alembert parametric, `FirstOrderSubstitution`,
  `AutonomousReduction`, linearizable `u=φ(y)`, integrating factors / implicit
  first integrals, `SolvableForY`/`SolvableForX` (the `y=G(x,y')` method).
- *Second order & higher, linear:* constant-coefficient (undetermined
  coefficients, variation of parameters, Dirac/step forcing), Cauchy–Euler,
  variable-coefficient via NormalForm + **Kovacic** (Cases 1 & 2) + Frobenius /
  power series, `ReductionOfOrder`, `OperatorFactor`/`DFactor`, `ExactODE`,
  higher-order autonomous reduction.
- *Special-function recognizers:* Bessel, Kummer/Gauss `₂F₁`, confluent
  Whittaker/`₁F₁`, generalized power-potential, Airy, Pöschl–Teller, Legendre/
  Chebyshev (via `DSolve\`SpecialFunctionForm`).
- *Systems:* general/defective/singular/triangular linear systems, coupled linear
  systems, variable-coefficient systems, 2-D autonomous nonlinear systems.
- *Lie point-symmetry engine (M10):* heuristic symmetries, `chi` heuristic,
  `function_sum`, `abaco1_product`/`abaco2_*` quadrature ansätze,
  `SecondOrderSymmetry` for nonlinear 2nd-order, `PolynomialShiftSubstitution`.
- *PDEs (M6 Phase 2):* first-order linear (characteristics), quasilinear Lagrange
  + Clairaut, Charpit (nonlinear), heat-equation Cauchy (heat kernel),
  wave-equation IVP (d'Alembert), separation of variables, 2nd-order
  constant-coefficient linear PDE, `PDEClassify` discriminant classifier.
- *BVP / eigenvalue (M11):* formal BVP, Sturm–Liouville eigenvalues.

> **ROADMAP flag.** The current §4.3 scope outline (4.3.1 Derivatives … 4.3.6
> Definite integration) has **no DSolve entry.** Add a **§4.3.7 Differential
> Equations (DSolve)** — or, given the size, break DSolve out as its own section.
> `CONTEXT.md`'s dual mandate means this can't be a syntax tour: it should teach
> the *methods* (the classification cascade, Kovacic, Lie symmetry) in
> `Theory`/`Under the hood` callouts.

**Integrate engines (item #2), for §4.3.5 (the Risch story):**
- `ParallelMixedTower` method — parallel Risch–Norman over a mixed radical tower;
  native number-field arithmetic; emits `Integrate::nonelem` non-elementarity
  certificates and non-torsion divisor certificates for genus ≥ 1.
- **Cherry** engines — exponential-tower dilogarithm and general-weight
  polylogarithm ladder integrals.
- `Integrate\`DiffUnderInt` — finite-domain Feynman families (power-log,
  secant-radical, tangent-power).
- High-fugacity (degenerate) Fermi–Dirac half-line integrals.

*(These are `Integrate` **methods/engines**, not standalone builtin heads — write
them as capabilities of `Integrate`, not as new functions.)*

### §4.2 Algebra — *Verified; predates all of this* — UPDATE

The written §4.2 covers factoring "over algebraic number fields" in prose but
predates the actual number-field API, and predates the new polynomial views and
transcendental solving.

**Number-field family (item #3).** `ToNumberField[Sqrt[2]]` →
`AlgebraicNumber[Sqrt[2], {0, 1}]` (confirmed). Cover: `AlgebraicNumber` as an
element of ℚ(θ) in the power basis; `ToNumberField`; field arithmetic under
`+ * / ^`; `AlgebraicNumberNorm`, `AlgebraicNumberTrace`,
`AlgebraicNumberDenominator`, `AlgebraicNumberPolynomial`;
`NumberFieldIntegralBasis` + `AlgebraicIntegerQ` (maximal order / integrality).
FLINT-`qqbar`-backed — good `Under the hood` material.

**Gröbner / polynomial views (item #4):** `PolynomialReduce` (multivariate
division with cofactors / ideal membership — confirmed `PolynomialReduce[x^2+1,{x},{x}]`
→ `{{x},1}`); `MonomialList`, `CoefficientRules`, `FromCoefficientRules` (sparse
views over all six monomial orders); `FactorSquareFreeList` + the new `Extension`
option on `FactorSquareFree`; `GroebnerBasis` now honours `DegreeLexicographic`
and the `Negative*` order family (was silently falling back to Lexicographic).

**Solving extensions (item #5):** `Reduce` now solves trig/hyperbolic equations &
inequalities over bounded regions, periodic trig equations, forward circular-trig
`f[A(x)]==c`, generalized two-argument trig/hyperbolic equations, and general
log/exp transcendental equations. QE/CAD gained multi-variable & alternating
quantifier elimination (v0.124), McCallum well-orientedness (v0.125), and
real-algebraic-coefficient fibre isolation (v0.113). `Simplify` over `Root[...]`
objects no longer hangs and collapses algebraically (v0.165). `Apart`'s
two-argument form now decomposes radicals of the given variable.

> **Note.** The base quantifier heads `Exists`/`ForAll`/`Resolve`/`LogicalExpand`
> (v0.095) actually *predate* the anchor but were never written into §4.2 either —
> so the whole quantifier-elimination story is an open UPDATE, not just the
> post-anchor extensions.

### §4.4 Linear Algebra — *stub* — INCLUDE

- `JordanDecomposition` (Jordan canonical form, v0.116) and `SchurDecomposition`
  (Schur + generalized/QZ, v0.118). Schur is already in the ROADMAP §4.4 scope
  ("LU/QR/Cholesky/SVD/Schur"); **add Jordan to that list.**
- `RowReduce` / `NullSpace` gained a user-facing `ZeroTest` option (v0.180–0.181)
  for predicate-driven exact reduction — worth a subsection on exact vs. machine
  reduction.

*(The many packed/NDArray and FLINT fast-path commits behind these are perf — see
the out-of-scope appendix — but the `Performance` callouts in §4.4 should mention
the machine-LAPACK vs exact/MPFR split, which is already in the ROADMAP scope.)*

### §4.7 Special Functions — *stub* — INCLUDE

- `LegendreQ` — Legendre function of the second kind (v0.113 wk). Small; fold into
  the §4.7 campaign alongside the Gamma/Bessel/hypergeometric families already
  scoped there.

### §4.8 Statistics — *Verified* — UPDATE (quick win)

- `Quantile` (confirmed `Quantile[{1,2,3,4,5},1/2]` → `3`), `InterquartileRange`,
  `MeanDeviation`, `MedianDeviation` — landed ~2 days *after* the §4.8 commit, so
  the Verified chapter is missing exactly this family. Natural home: the
  five-number-summary / spread subsections (§4.8 already discusses `Quartiles`,
  `Variance`, `MeanDeviation`-adjacent robustness). Smallest self-contained UPDATE.

### Ch. 6 Data Structures / Ch. 7 Programming — *stubs* — INCLUDE

- **Arrays/lists (item #12):** `ArrayReshape`, `ArrayPad`, `ListGradient`
  (numpy.gradient port), `FirstPosition` (now takes a levelspec + default).
- **Packed arrays (item #15):** `PackedArrayQ`, `ToPackedArray`, `ToNDArray` — the
  user-facing edge of the packed/NDArray substrate that Ch. 6 is scoped to explain.
- **Evaluation control (item #14):** `Inactive` / `Activate` — general
  evaluation-control heads (DSolve uses them for inert first integrals). Fits Ch. 7
  (or an expression-control subsection).

### §4.5 Numerical Calculus — *Verified* — UPDATE (small)

- `ListInterpolation` (values on a regular grid) and default `Options` on
  `Interpolation`. §4.5 already covers `InterpolatingFunction`; add the grid form.
  (Could equally live in Ch. 6 — writer's call.)

### Ch. 9 Data I/O — *stub, ROADMAP caveat now outdated* — INCLUDE

- A real stream layer landed (`src/io/`): `Read`, `OpenRead`/`OpenWrite`/
  `OpenAppend`, `Write`, `WriteString`, `Close`, `Streams`, `StreamPosition`,
  `SetStreamPosition`, and `ReadList` (type-directed file reading). ROADMAP's Ch. 9
  note ("*mostly not implemented yet … write last*") is now partly wrong — there is
  real, documentable I/O. Update that caveat and scope the chapter to what exists.

### No current ROADMAP home — *propose placement*

- **Computational geometry (item #8):** `Area`, `Perimeter`, `RegionCentroid`,
  `RegionMember`, `ConvexHullRegion` — exact GMP-rational coordinates over
  `Polygon`/point sets. No chapter or §4.x section covers geometry. Propose a new
  math section (e.g. **§4.9 Computational Geometry**) or a Ch. 6 sub-topic.
- **Graphs (item #9):** `StarGraph`, `FindVertexColoring` (exact minimal
  colouring), `EdgeWeight` / `WeightedAdjacencyMatrix`, weighted `FindShortestPath`
  / `GraphDistance` (Dijkstra), `RandomGraph[{n,m},k]`. The graph subsystem
  (`src/graph/`, `graphs.md` in the spec) has **no book campaign at all** — propose
  a Graphs section/chapter; these additions are the trigger to open one.
- **Bitwise (item #16):** `BitLength` (confirmed `BitLength[255]` → `8`) + the new
  `src/bitwise/` module that seeds a Bit\* family. Propose a short §4.1 Arithmetic
  bitwise subsection or a Ch. 6 topic.

### REVIEW — behaviour/semantics changes (no new head)

Verified transcripts auto-update on rebuild, but re-check the surrounding *prose*:

- **A1–A10 Mathematica-divergence fixes** (v0.503-era, 09-21 wk): packed-list
  iterators (`Do`/`Table`/`Sum`/`Product`); `Modulus->p` for `PolynomialQuotient`/
  `PolynomialRemainder`/`Mod`/`ExtendedGCD`; ragged `PadRight`/`Transpose`/
  `Lookup`; `Discriminant` deg ≤ 1 = 1; `OptionValue` under `OptionsPattern[other]`.
  → Re-check **§4.1** (Mod/Quotient), **§4.2** (Discriminant, polynomial mod), and
  the Ch. 3 tour.
- **`Function` lexical closure** over an enclosing `Function`'s parameter (v0.16x).
  → Re-check the Ch. 3 functional-programming tour (and Ch. 7 when written).
- **`Apart` two-argument radical form**, **`DSolve` no longer `HoldAll`** — minor;
  note if any written prose asserts the old behaviour.
- *Minor correctness fixes* not worth prose unless a written example relies on them:
  `TrigToExp[Coth]` sign fix, `Factor`/`Simplify` dropped-leading-term fix,
  `PossibleZeroQ` exponential-combining, `Equal`/`Unequal` numeric decision. (The
  `Degree`-in-exact-eval, bare-`N`→machine, and `Eigenvalues` conjugate-order fixes
  predate the anchor and are already in the book era.)

---

## Explicitly out of scope for the book (triaged, not missed)

Internal / performance / build / churn — no reader-facing surface:

- **Perf:** Orderless-canonicalisation symbol-set memo & per-node cache
  (v0.185/0.188), field-aware `Can` (v0.189), FLINT exact-rational `RowReduce` fast
  path (v0.186), native `nf_elem` `RowReduce` over `AlgebraicNumber` matrices
  (v0.179), `RootReduce` `Root→qqbar` memoization (v0.169), FLINT `nmod_poly`
  modular fast path, matrix-decomposition packed/NDArray buffer paths
  (v0.119/0.120), Interpolation packed-buffer table build.
- **Churn (do NOT document as features):** `RischNorman` / `RischNormanBlake`
  methods were added (v0.121–0.123) then **removed** (subsumed by
  `ParallelMixedTower`); `TowerCRE` is dormant infra.
- **Build / packaging / stability:** `make install`/`uninstall` layout (#77),
  `USE_MPFR=0`/`USE_LAPACK=0` clean builds, GCC-16 miscompile & `-Wmaybe-uninitialized`
  fixes, `TimeConstrained` per-call leak fix, `ReadProtected` attribute removed,
  TeXForm `C[k]→c_k` rendering, per-call `Print` flush.
- **Testing/hardening:** the vast DSolve corpus milestones (M15, M21–M55: the
  12000.org §2.1.2/§2.2.\* harnesses) drove the DSolve *features* above but are not
  themselves book content.

---

## ROADMAP status-impact note (for a future writer — not yet applied)

When the campaigns above are opened, `ROADMAP.md` should be updated:

1. **§4.3 Calculus** — add a **Differential Equations (DSolve)** subsection to the
   scope outline (currently absent); consider splitting DSolve into its own
   section given its size. Add ParallelMixedTower/Cherry/DiffUnderInt to the §4.3.5
   Risch story.
2. **§4.4 Linear Algebra** — add **Jordan** to the decompositions list; add the
   `ZeroTest` exact-reduction option.
3. **New section** — add **Computational Geometry** (§4.9?) and open a **Graphs**
   section/chapter; neither exists today.
4. **Ch. 9 Data I/O** — soften the "mostly not implemented, write last" caveat: the
   stream layer now exists.
5. Statuses to revisit once written: §4.2 and §4.8 remain **Verified** only until
   their UPDATEs land.
