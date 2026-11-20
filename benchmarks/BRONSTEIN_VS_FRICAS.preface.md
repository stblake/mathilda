## Summary (2026-11-21, Mathilda v0.6646, FriCAS 1.3.6)

The comparison runs on three corpora. Benchmarks 256 and 257 are the integration shootout
corpora; benchmark 259 holds the C completeness campaign's hunt families. On the Mathilda
side there are three runs:

- **Mathilda**: Bronstein's algorithm as Mathilda implements it. That is RischTranscendental
  (*Symbolic Integration I*) for rational and transcendental integrands, Trager for algebraic
  ones, and the Bronstein package for the mixed algebraic and transcendental cases.
- **PMT**: ParallelMixedTower, the parallel Risch (Risch–Norman) method over a mixed tower.
- **Automatic**: the full cascade.

The first version of this report (2026-11-20, v0.6641) led to an improvement plan, items
C0–C5 in `src/internal/bronstein/BRONSTEIN_PLAN.md`. Releases v0.6642–v0.6646 implement
C0, C1, C2 (in part), C3.2 and C4 (in part). This run measures the result.

**Completeness: how often each system decides.** Deciding means giving a verified closed
form, or proving that no elementary integral exists.

| corpus | Mathilda | PMT | Automatic | FriCAS |
|---|---|---|---|---|
| 256 (1266) | **96%** (was 93%) | 77% | **96%** (was 91%) | 75% |
| 257 (1000, all elementary) | **100%** (was 99%) | 86% | **100%** | 90% |
| 259 (640, transcendental) | **98%** (was 96%) | 66% | **98%** (was 92%) | 95% |

- **Mathilda decides more than FriCAS on every corpus.** It is ahead on algebraic integrands
  (256 class A: 100% against 82%; Automatic 91% -> 99%) and far ahead on a radical over a
  log or exp tower (class M: 100% against 44%).
- **Class B (a radical at the bottom) is the one class where FriCAS still leads:** 93%
  against 72%, which was 54%. Several independent square roots are now integrated by
  Galois components (C3.2, theory §E8). The other cause, Risch over an algebraic base
  (a radical inside a log or exp with further monomials above it), is Phase 5 (C3.1),
  not started.
- **The cascade now keeps what its methods find (C1).** Trager's certificates are no longer
  hidden behind an `Inactive[Integrate]`. A large tower integrand reaches the Risch stage
  before the stages that grind on it. The decision's verified antiderivative is the
  integrator's answer. Automatic's timeouts on 259 went from 34 to 6.
- **On the transcendental corpus** Mathilda now decides 625 and Automatic 627, against
  FriCAS's 610. FriCAS still *closes* more cases with special functions (polylog, Ei, erf)
  where Mathilda certifies non-elementarity.

**Soundness.** No Mathilda method gives a false "no" or a wrong answer on any corpus.
The Automatic wrong answer of the first run (the PolyLog stage, `Log[-I + Tan[2 x]]`) is
fixed (C0, v0.6642).

- **FriCAS: 14 false "no".** 11 are on integrands elementary by construction (`gauss-elem`);
  3 more have an elementary answer from Mathilda that verifies (`gauss-rand`). All are
  exponentials differing by a constant over `Q(i)`. The minimal case is
  `integrate(exp(1+2*%i*x)/(1+exp(%i*x)), x)`, whose integral is
  `-%i*%e*(exp(%i*x) - log(1+exp(%i*x)))`. A draft report for FriCAS is in
  `benchmarks/FRICAS_UPSTREAM_DRAFT.md` (not submitted).
- **FriCAS: 40 answers fail their derivative check** at the real points. Rechecked in FriCAS,
  18 of them are formally correct (`normalize(D(r, x) - f) = 0`, a branch choice); for the
  other 22, `normalize` does not reach 0 on the nested radicals, so they are unresolved.
- **One lower-layer soundness bug was found and fixed while doing C2.** `Cancel` over `Q(i)`
  left `(-E + I E u)/(I + u)` uncancelled, and the base Risch DE turned that into a false
  "no". The hunts caught it before release (v0.6644).

**Speed.** The comparison uses the cases both systems close, with each system's own
integration time (start-up excluded).

- **Mathilda is faster on the typical case.** Its medians are 16–71 ms against FriCAS's
  70–110 ms, and it is faster on 1345 of 1927 cases. FriCAS's timer has a floor of about
  60 ms per call, so the medians favour Mathilda.
- **FriCAS has no slow tail; Mathilda has one.** 51 of Mathilda's common cases take over
  1 s, against none of FriCAS's. Mathilda's 90th percentile is 226–341 ms against
  110–160 ms. Most of the tail is the Bronstein package's `.m` code on radical-on-top
  families; a profile puts the evaluator's scoping substitution (`Module`/`With` renaming)
  ahead of the algebra. Fixing that needs a per-node symbol summary in `Expr`, a core
  design decision, recorded under C4.
- **PMT is the slowest of the three,** with a geometric-mean ratio of 0.9–1.8 against FriCAS.
  It still decides some cases the dedicated methods do not, so it stays as the fallback.

**Fairness notes.** FriCAS's unevaluated `integral(...)` counts as its proof of
non-elementarity, except on the elementary-by-construction cases, where it is counted as
false. Its answers in Ei, li, erf and erfi are verified exactly, because `complexNumeric`
cannot evaluate them. The FriCAS runs on 256 and 257 are the cached runs of
2026-11-12/13; 259 ran on 2026-11-20.

**What follows from this.** The improvement plan and its status are in
[`src/internal/bronstein/BRONSTEIN_PLAN.md`](../src/internal/bronstein/BRONSTEIN_PLAN.md),
under "Plan from the FriCAS comparison". What is left: C2 (the Liouvillian primitive's
bottom coefficient, a few undecided tower cases), C3.1 (Phase 5), the C4 evaluator change,
and C5 (special-function answers where Mathilda certifies; submitting the FriCAS report).
