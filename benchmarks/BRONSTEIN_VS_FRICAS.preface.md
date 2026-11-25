## Summary without the easy integrals (2026-11-25, FriCAS 1.3.6)

From 2026-11-25 every corpus leaves out its **easy** integrals (`easy.tsv`; see the
note above). 3191 of the 4876 cases of 256-261 are easy: rational functions, integrands
a substitution makes rational, and those `DerivativeDivides` closes. 257 keeps 18 of
1000 cases and 261 keeps 58 of 1550. 262 (1000) and 263 (1850) are integrals generated
outside the easy classes. 4535 cases are left; the counts below are of those only. 263
ran at v0.6733, 262 at v0.6732, 261 at v0.6714, the others at v0.6691.

**Decided** (a verified closed form, or a proof that none is elementary):

| corpus | n | Mathilda | PMT | Automatic | FriCAS |
|---|---|---|---|---|---|
| 256 | 777 | **100%** (775) | 73% (568) | 99% (771) | 62% (481) |
| 257 | 18 | 94% (17) | 67% (12) | **100%** (18) | 78% (14) |
| 259 | 561 | 99% (557) | not run | **100%** (560) | 96% (536) |
| 260 | 271 | **98%** (266) | 35% (96) | 97% (264) | 93% (252) |
| 261 | 58 | 98% (57) | 55% (32) | **100%** (58) | 24% (14) |
| 262 | 1000 | **98%** (976) | 66% (660) | 97% (970) | 64% (640) |
| 263 | 1850 | 98% (1817) | 72% (1339) | **99%** (1826) | 69% (1276) |
| all | 4535 | **98%** (4465) | 68% (2707 of 3974) | **98%** (4467) | 71% (3213) |

- **Leaving out the easy cases widens the gap.** FriCAS decides 71% of what is left,
  against 95% of 261 and 90% of 257 with the easy cases in. Mathilda's share is
  94-100% on every corpus (94% only on 257's 18 cases).
- **262 by class.** FriCAS closes no M case (0 of 216) and 51 of 118 B cases; Mathilda
  closes all of both. On A, Trager closes 312 of 336 and FriCAS 266. On T both are near
  complete (330 and 323).
- **Soundness on 262.** No Mathilda method gives a wrong answer or a false "no". All 7
  FriCAS answers that failed their own check are zero by `PossibleZeroQ` in Mathilda.
- **263: FriCAS gives 11 false "no".** On B4 integrands elementary by construction,
  derivatives of `E^y R` with y the root of a cubic or quartic, FriCAS returns
  `integral(...)`. Minimal case in FriCAS 1.3.6:
  `integrate(D(exp(sqrt(x^3+1))/(x+1), x), x)`, which should be
  `exp(sqrt(x^3+1))/(x+1)` (with `1/x` it integrates). FriCAS closes none of 263's
  400 M cases.
- **Automatic below the dedicated method on 262: 14 cases** (11 B, 3 A). The cascade
  stops at about 20 s, and Bronstein or Trager needs 10-26 s on these sixth-root radicals
  (`cB-201-559`, `Log[x]/(x^6 (3 + 2 x^6)^(1/6))`). It is a speed gap, not a decline.
  263 has the same gap for Trager: 22 A cases it closes in 11-15 s, where Automatic
  gives up at 20-27 s.
- **Speed on 262** (623 cases both close): median 328 ms against FriCAS's 100 ms, and
  154 Mathilda cases over 1 s. The hard cases are where the `.m` Bronstein and Trager
  code is slow. The run shared the machine with other jobs (load average 50-66), which
  inflates both systems' times. On 263 (1259 cases both close, rerun at low load) the
  medians are 222 ms against 170 ms, but Mathilda's 90th percentile is 3.6 s against
  250 ms. The tail is Trager on class A.

The sections below are from earlier runs. Their counts **include the easy integrals**.

## Benchmark 261 (2026-11-24, Mathilda v0.6714, FriCAS 1.3.6; easy integrals included)

261 is 1550 new random integrands, each of whose antiderivative is more than four
times its LeafCount. They are spread over five classes: R 50, T 250, A 250, M 550 and
B 450. 261 ran at v0.6714; the four corpora below ran at v0.6691.

- **Decided:** RischTranscendental, Trager or Bronstein 1547 (100%), Automatic 1549 (100%),
  FriCAS 1472 (95%). FriCAS falls behind on M, a radical over a Log/Exp tower: 481 of 550
  (87%) against Mathilda's 548 and Automatic's 549. Neither system gave a special-function
  answer.
- **One rule for "wrong".** An answer that fails its own system's check, or that its own
  check cannot evaluate, is translated into Mathilda. It is then judged by
  `PossibleZeroQ[D[F, x] - f]`, the test Mathilda's answers get.
  - FriCAS: 97 of its 100 such answers are zero, so FriCAS has no wrong answers on 261.
    The point check had failed on a branch choice of FriCAS's own `complexNumeric`; the
    same answers pass the point check at 40 digits in Mathilda.
  - Axiom: 156 answers are wrong by both checks.
- **Soundness:** no Mathilda method gives a wrong answer or a false "no".

## Summary (2026-11-22, Mathilda v0.6691, FriCAS 1.3.6; easy integrals included)

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

Since the v0.6646 run, releases v0.6647–v0.6691 contain:
- Campaign D: D0–D2, D3 (Phase 5 M1 and M2: transcendental Risch over `Q(x)(y)` with any number
  of monomials), D4, D5, L1 (gap I2), D7 in part;
- D6's evaluator change;
- the v0.6691 cascade fixes found by the intermediate v0.6688 run.

**Completeness: how often each system decides.** Deciding means giving a verified closed
form, or proving that no elementary integral exists.

| corpus | Mathilda | PMT | Automatic | FriCAS |
|---|---|---|---|---|
| 256 (1266) | **100%** (1264; was 96%) | 77% | **99%** (1259; was 96%) | 75% |
| 257 (1000, all elementary) | **100%** | 86% | **100%** | 90% |
| 259 (640, transcendental) | **99%** (636; was 98%) | not run | **100%** (639; was 98%) | 95% |
| 260 (420, radical at the bottom) | **99%** (415; was 85%) | 47% | **98%** (412; was 89%) | 95% |

- **Mathilda decides more than FriCAS on every corpus and every class.** Class B (a radical at
  the bottom), the one class where FriCAS led at v0.6646 (93% against 72%), is now 169 of 169 on
  256 and 99% on 260.
- **The cascade fixes of v0.6691** close cases where Automatic gave a special-function form
  for an elementary integrand, or timed out where the Risch stage answers in under a second (259
  `tan-elem-26`: 60 s to 0.95 s). They also add Erf/Ei/PolyLog answers FriCAS gives (D7: one case
  of the 18 left).
- **Not in these numbers:** v0.6692 (the homogeneous-RDE witness, L2) and v0.6693 (the
  exact-support split: 260 `two-elem-6`, 48 s to 2.2 s), on a branch during this run.

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
integration time (start-up excluded). This run had the machine largely to itself.
- **Medians:** Mathilda's are now 21–60 ms against FriCAS's 70–130 ms. It is faster on 1539 of
  2109 shared cases.
- **Tail:** 45 cases over 1 s across the four corpora (10, 20, 11 and 4), against 51 on three
  corpora at v0.6646.
- **90th percentiles:** 201–332 ms against 110–170 ms.
- **Shape of the gap:** FriCAS still has no slow tail. The p90 gap is the remaining D6 work, so
  the C4 criterion (p90 within 2x) is not yet met on 257 and 259.

The older paragraphs below are from the 2026-11-21 run.

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
