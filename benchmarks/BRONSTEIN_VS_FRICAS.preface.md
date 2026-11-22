## Summary (2026-11-22, Mathilda v0.6688, FriCAS 1.3.6)

The comparison runs on four corpora:
- 256 and 257, the integration shootout corpora;
- 259, the C completeness campaign's hunt families;
- 260, a radical at the bottom (Campaign D0).

On the Mathilda side there are three runs:
- **Mathilda**: Bronstein's algorithm as Mathilda implements it. That is RischTranscendental
  (*Symbolic Integration I*) for rational and transcendental integrands, Trager for algebraic
  ones, and the Bronstein package for the mixed algebraic and transcendental cases.
- **PMT**: ParallelMixedTower, the parallel Risch (Risch–Norman) method over a mixed tower.
- **Automatic**: the full cascade.

The previous run (2026-11-21, v0.6646) led to Campaign D in
`src/internal/bronstein/BRONSTEIN_PLAN.md`. Releases v0.6647–v0.6688 contain:
- D0 and D1 (support groups);
- D2 (the free square-root split);
- D3 M1 (Phase 5: Risch over the curve, two monomials);
- D4 and D5;
- D6's evaluator change;
- D7 in part.

This run measures the result.

**Completeness: how often each system decides.** Deciding means giving a verified closed
form, or proving that no elementary integral exists.

| corpus | Mathilda | PMT | Automatic | FriCAS |
|---|---|---|---|---|
| 256 (1266) | **100%** (was 96%) | 77% | **99%** (was 96%) | 75% |
| 257 (1000, all elementary) | **100%** | 86% | **100%** (999: one timeout under load, 14 s alone) | 90% |
| 259 (640, transcendental) | **99%** (was 98%) | not run | **100%** (was 98%) | 95% |
| 260 (420, radical at the bottom) | **99%** (was 85%) | 46% | **98%** (was 89%) | 95% |

- **Mathilda now decides more than FriCAS on every corpus and every class.** Class B (a
  radical at the bottom), the one class where FriCAS led, is 169 of 169 on 256 (was 72%,
  FriCAS 93%) and 99% on 260 (FriCAS 95%). Three things closed it:
  - the support groups (D1);
  - the free square-root split (D2);
  - transcendental Risch over `Q(x)(y)` (Phase 5 M1).
- **On the transcendental corpus** Mathilda decides 636 and Automatic 637, against
  FriCAS's 610.
- **Not yet in these numbers:** Phase 5 M2 (any number of monomials over the curve, v0.6689,
  on a branch during this run). The comparison is rerun after it lands.

**Soundness.** No Mathilda method gives a false "no" or a wrong answer on any corpus.
The Automatic wrong answer of the first run (the PolyLog stage, `Log[-I + Tan[2 x]]`) is
fixed (C0, v0.6642).

- **FriCAS: 14 false "no".** 11 are on integrands elementary by construction (`gauss-elem`);
  3 more have an elementary answer from Mathilda that verifies (`gauss-rand`). All are
  exponentials differing by a constant over `Q(i)`. The minimal case is
  `integrate(exp(1+2*%i*x)/(1+exp(%i*x)), x)`, whose integral is
  `-%i*%e*(exp(%i*x) - log(1+exp(%i*x)))`. A draft report for FriCAS is in
  `benchmarks/FRICAS_UPSTREAM_REPORT.md`, submitted as https://github.com/fricas/fricas/issues/249.
- **FriCAS: 42 answers fail their derivative check** at the real points (40 at v0.6646). Rechecked in FriCAS,
  18 of them are formally correct (`normalize(D(r, x) - f) = 0`, a branch choice); for the
  other 22, `normalize` does not reach 0 on the nested radicals, so they are unresolved.
- **One lower-layer soundness bug was found and fixed while doing C2.** `Cancel` over `Q(i)`
  left `(-E + I E u)/(I + u)` uncancelled, and the base Risch DE turned that into a false
  "no". The hunts caught it before release (v0.6644).

**Speed.** The comparison uses the cases both systems close, with each system's own
integration time (start-up excluded). **This run's Mathilda timings are not comparable with the
previous ones.** It shared the machine with a 4-process stress run and a test-suite build, so
its tail (58 cases over 1 s) is inflated. The paragraphs below are from the 2026-11-21 run. The
speed item (D6) is re-measured on an idle machine.

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
2026-11-12/13; 259 ran on 2026-11-20 and 260 on 2026-11-21 (cached).

**What follows from this.** The improvement plan and its status are in
[`src/internal/bronstein/BRONSTEIN_PLAN.md`](../src/internal/bronstein/BRONSTEIN_PLAN.md),
under "Plan from the FriCAS comparison". What is left: C2 (the Liouvillian primitive's
bottom coefficient, a few undecided tower cases), C3.1 (Phase 5), the C4 evaluator change,
and C5 (special-function answers where Mathilda certifies; submitting the FriCAS report).
