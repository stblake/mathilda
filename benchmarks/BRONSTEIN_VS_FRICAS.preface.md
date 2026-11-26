## Summary: easy integrals left out, every answer judged by PossibleZeroQ (2026-11-26, FriCAS 1.3.6)

**What is counted.**
- **Constant factors are removed from the integrands** (`normalize_constants.py`):
  Integrate[c f, x] is c Integrate[f, x], so the -1/2 of -1/2 (x^2 - 3x - 2)/Sqrt[g],
  or the Sqrt[2] in Sqrt[4 + 2 E^(-2x) + 8x + 4x^2], is no part of the problem. Every
  corpus was rewritten factor by factor, with each reference divided by the same
  constant and checked by PossibleZeroQ. An integrand that became a duplicate of another was
  left out (on 264 it was replaced by a spare). Every system was then run again.
- **The easy integrals are left out** of every corpus (`easy.tsv`; see the note above):
  rational functions, integrands a substitution makes rational, and those
  `DerivativeDivides` closes. 3129 of the 4706 cases of 256-261 are easy; 257 keeps 11
  of 985 cases and 261 keeps 43 of 1512. 15 cases of 256 are undecided by the test
  (it runs out of time) and are kept.
- **262 (958), 263 (1801) and 264 (500)** were generated outside the easy classes; 3
  and 37 of them turned out to be easy once DerivativeDivides no longer depended on a
  constant factor (v0.6760). 4796 cases are left.
- **Every answer of every system is judged by one test**, `PossibleZeroQ[D[F, x] - f]`
  in Mathilda (`pzq_verify.py`), with the other systems' answers translated first. The
  standard is an antiderivative valid for every complex x. Each system's own
  three-point check no longer decides anything.
- **Versions.** Every run is from 2026-11-26 with Mathilda v0.6750, and so are the
  verdicts. The easy labels are from v0.6760.

**Decided** (a verified closed form, or a proof that none is elementary):

| corpus | n | Mathilda | PMT | Automatic | FriCAS |
|---|---|---|---|---|---|
| 256 | 699 | 99% (695) | 72% (506) | **100%** (696) | 60% (421) |
| 257 | 11 | 91% (10) | 82% (9) | **100%** (11) | 64% (7) |
| 259 | 555 | **99%** (552) | not run | **99%** (552) | 95% (528) |
| 260 | 269 | **99%** (265) | 35% (94) | **99%** (265) | 97% (260) |
| 261 | 43 | 98% (42) | 67% (29) | **100%** (43) | 19% (8) |
| 262 | 955 | 98% (932) | 67% (641) | **98%** (937) | 65% (625) |
| 263 | 1764 | 98% (1735) | 73% (1285) | **100%** (1761) | 71% (1255) |
| 264 | 500 | **100%** (500) | 0% (0) | **100%** (500) | 0% (0) |
| all | 4796 | 99% (4731) | 60% (2564 of 4241) | **99%** (4765) | 65% (3104) |

- **No Mathilda method gives a wrong answer or a false "no"** under PossibleZeroQ, on
  any corpus.
- **Mathilda's unverified answers:** 7, all answers PossibleZeroQ does not decide
  (logxy-1-1 and logxy-3-6 on 256, by Automatic and Bronstein; tani-elem-26 and
  gauss-elem-23 on 259). The answers are no longer cut when they are recorded.
- **264: 500 class M integrals, seven systems.** Bronstein and Automatic close all of
  them. Of FriCAS, Axiom, Maxima, sympy, Wolfram (wolframscript, Wolfram Language
  14.0.0) and Reduce (PSL, revision 7327), only Wolfram closes one. FriCAS stops with
  "implementation incomplete (has polynomial part)" (274) or "(constant residues)" (226),
  Axiom with "(constant residues)" (500). Wolfram returns 496 unevaluated and times out
  on 3 (median 4.1 s), Reduce 446 unevaluated and 54 timeouts.
- **The other benchmarks keep their systems:** 256 and 257 ran Mathilda, FriCAS,
  Maxima, sympy and Wolfram; 259 and 260 Mathilda and FriCAS; 261-263 Mathilda,
  FriCAS, Axiom, Maxima and sympy.
- **FriCAS's error records were partly the harness's.** A fault in the adapter's
  output file recorded successful FriCAS runs as errors. Re-run with the fault fixed,
  FriCAS's counts rose on 256, 259, 260, 262 and 263.
- **Answers rejected by PossibleZeroQ:** FriCAS 8 (7 on 259, 1 on 256). They are
  antiderivatives on part of the plane, nonzero at a point such as x = -5/2 or 3i.
  Maxima 28 (6 on 256, 3 on 262, 19 on 263), sympy 43 (15 on 262, 28 on 263), Axiom 6
  (2 on 262, 4 on 263) and Wolfram 3 (on 256). Maxima's are `abs()` and
  `asinh(.../(x+1))` forms valid only on a half-line.
- **FriCAS gives 21 false "no"**: 10 on 259 (gauss-elem) and 11 on 263. On 263's B4
  integrands, elementary by construction as derivatives of `E^y R` with y the root of
  a cubic or quartic, FriCAS returns `integral(...)`. Minimal case in FriCAS 1.3.6:
  `integrate(D(exp(sqrt(x^3+1))/(x+1), x), x)`, which should be
  `exp(sqrt(x^3+1))/(x+1)` (with `1/x` it integrates).
- **Automatic decides at least as many as the dedicated method** on every corpus: 937
  against 932 on 262 and 1761 against 1735 on 263. Since v0.6744 an expired cascade
  budget resumes Trager or Bronstein once the later stages decline.
- **Speed** (the cases both close, each system's own time):
  - 262 (608 cases): median 248 ms against FriCAS's 100 ms; Mathilda is faster on 300.
  - 263 (1239 cases): median 59 ms against 140 ms; Mathilda is faster on 880, and its
    90th percentile is 1.3 s against 230 ms.
  - The tail is the `.m` Trager code on class A.

The sections below are from earlier runs. Their counts **include the easy integrals**, and their answers were verified by each system's own check.

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
