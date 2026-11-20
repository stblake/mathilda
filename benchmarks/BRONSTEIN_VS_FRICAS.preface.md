## Summary (2026-11-20, Mathilda v0.6641, FriCAS 1.3.6)

The comparison runs on three corpora. Benchmarks 256 and 257 are the integration
shootout corpora; benchmark 259 is new and holds the families of the C completeness
campaign. On the Mathilda side there are three runs:

- **Mathilda**: Bronstein's algorithm as Mathilda implements it. That is
  RischTranscendental (*Symbolic Integration I*) for rational and transcendental
  integrands, Trager for algebraic ones, and the Bronstein package for the mixed
  algebraic and transcendental cases.
- **PMT**: ParallelMixedTower, the parallel Risch (Risch–Norman) method over a
  mixed tower, run on every class.
- **Automatic**: the full cascade.

**Completeness: how often each system decides.** Deciding means giving a verified
closed form, or proving that no elementary integral exists.

| corpus | Mathilda | PMT | Automatic | FriCAS |
|---|---|---|---|---|
| 256 (1266) | **93%** | 77% | 91% | 75% |
| 257 (1000, all elementary) | 99% | 86% | **100%** | 90% |
| 259 (640, transcendental) | **96%** | 66% | 92% | 95% |

- **Where Mathilda decides more.** It is ahead on algebraic integrands (256 class
  A: 100% against 82%) and decisively ahead on a radical over a log or exp tower
  (class M: 100% against 44%). FriCAS answers "integrate failed" or "implementation
  incomplete" on most of the M cases.
- **Where FriCAS decides more.** It is ahead when the radical is at the bottom of
  the tower (class B: 93% against 54%), on the perturbed `bgenus`, `bcurve` and
  `bbottom` families. Bronstein's own scope check names two causes:
  - several independent radicals, a non-cyclic Kummer extension (31 cases);
  - transcendental Risch over an algebraic base, i.e. a log or exp of a radical
    with further monomials above it (47 cases).
- **On the transcendental corpus** the two decide about the same share (616 and
  610 of 640). FriCAS *closes* more (429 against 405, and 21 cases that no Mathilda
  method closes): it answers in polylog, Ei and erf where Mathilda certifies
  non-elementarity. FriCAS decides 3 elementary cases (`tan-elem`, `gauss-elem`,
  `f12-elem`) that no Mathilda method does.
- **PMT** decides 19 cases the Bronstein-family methods leave open. FriCAS also
  fails on only 2 of those. 6 of them are certificates in class B, the class where
  Mathilda trails. Everywhere else PMT is a subset of the dedicated methods, and
  slower.
- **Since 2026-11-12** the dedicated methods decide 1140 -> 1182 of 256 and close
  939 -> 991 of 257, on corpora the campaign was never tuned on.

**Soundness.** No Mathilda method gives a false "no" on any corpus.

- **FriCAS: 11 false "no".** All are in `gauss-elem` (259), integrands that are
  elementary by construction. FriCAS returns an unevaluated `integral` when the
  exponentials differ by a constant factor, such as `E^(1 + 2 I x)` beside
  `E^(I x)` or `E^(I (2 + x))`. This is the same class of error that Mathilda's
  v0.6638 fixed in Mathilda.
- **FriCAS: 40 answers fail their derivative check** at the real points (1 in
  256, 39 in 257). Rechecked in FriCAS, 18 of them are formally correct
  (`normalize(D(r, x) - f) = 0`): the radicals are evaluated on a different branch.
  For the other 22, `normalize` does not reach 0; it cannot denest radicals like
  `12^(1/4)`, so those are unresolved rather than proven wrong.
- **Mathilda: 1 wrong answer**, `I x^2 Log[-I + Tan[2 x]]` (259, `tani-rand-9`),
  from Automatic's special-function stage. It predates the campaign. The stage
  rewrites `Log[-I + Tan[2 x]]` as logarithms of exponentials, an identity valid on
  one branch interval only, so the PolyLog answer is wrong past the pole at
  `x = Pi/4`. RischTranscendental certifies the case non-elementary, correctly.

**Speed.** The comparison uses the cases both systems close, with each system's
own integration time (start-up excluded).

- **The Bronstein-family methods are faster on the typical case.** Their medians
  are 19–67 ms against FriCAS's 70–110 ms, and they are faster on 1354 of 1923
  cases. FriCAS's timer has a floor of about 60 ms per call, so the medians favour
  Mathilda.
- **FriCAS has no slow tail; Mathilda has one.** None of FriCAS's common cases
  takes over 1 s, against 50 of Mathilda's (Bronstein on radical-on-top families;
  a few Tan and Tanh towers). So Mathilda's 90th percentile is slower (256–343 ms
  against 110–160 ms), and so are its totals (287 s against 179 s).
- **PMT is slower than FriCAS.** Its geometric-mean ratio is 0.9–2.0 against
  FriCAS, and 169 of its cases take over 1 s.
- **The cascade costs time and decisions.** On 259, Automatic times out (30 s) on
  34 cases that RischTranscendental decides; in the one profiled, the
  DerivativeDivides stage (12, before Risch) used 57.7 s. On 256 class A, Automatic
  returns an `Inactive[Integrate]` rewrite on 35 integrands that Trager certifies
  non-elementary, so the certificate is lost.

**Fairness notes.** FriCAS's unevaluated `integral(...)` counts as its proof of
non-elementarity, except on the elementary-by-construction cases, where it is
counted as false. Its answers in Ei, li, erf and erfi are verified exactly,
because `complexNumeric` cannot evaluate them. The FriCAS runs on 256 and 257 are
the cached runs of 2026-11-12/13 (the corpora are unchanged); 259 ran on
2026-11-20.

**What follows from this.** The improvement plan built on these results is in
[`src/internal/bronstein/BRONSTEIN_PLAN.md`](../src/internal/bronstein/BRONSTEIN_PLAN.md),
under "Plan from the FriCAS comparison (2026-11-20)".
