## Summary (2026-11-20, Mathilda v0.6641, FriCAS 1.3.6)

The comparison runs on three corpora. Benchmarks 256 and 257 are the integration
shootout corpora; benchmark 259 is new and holds the families the C completeness
campaign was built against. Mathilda here means Bronstein's algorithm as Mathilda
implements it:
- RischTranscendental, from *Symbolic Integration I*, for rational and
  transcendental integrands;
- Trager, for algebraic integrands;
- the Bronstein package, for the mixed algebraic and transcendental cases.

**Completeness: how often each system decides.** Deciding means giving a verified
closed form, or proving that no elementary integral exists.

| corpus | Mathilda | FriCAS | Mathilda, full cascade |
|---|---|---|---|
| 256 (1266) | **93%** | 75% | 91% |
| 257 (1000, all elementary) | **99%** | 90% | 100% |
| 259 (640, transcendental) | **96%** | 95% | 92% |

- **Where Mathilda decides more.** It is ahead on algebraic integrands (256 class
  A: 100% against 82%) and decisively ahead on a radical over a log or exp tower
  (class M: 100% against 44%). FriCAS answers "integrate failed" or "implementation
  incomplete" on most of the M cases.
- **Where FriCAS decides more.** It is ahead when the radical is at the bottom of
  the tower (class B: 93% against 54%). The difference is the perturbed `bgenus`,
  `bcurve` and `bbottom` families: FriCAS proves about 60 of them non-elementary,
  and Mathilda's Bronstein package leaves them undecided. This is the clearest
  remaining completeness gap.
- **On the transcendental corpus** the two decide about the same share (617 and
  606 of 640). FriCAS *closes* more cases (428 against 406), because it answers in
  polylog, Ei and erf where Mathilda certifies non-elementarity. It decides 3
  elementary cases that Mathilda leaves undecided or times out on.
- **Since 2026-11-12** the dedicated methods decide 1140 -> 1182 of 256 and close
  939 -> 991 of 257. This is the campaign's gain on corpora it was never tuned on.

**Soundness.** Mathilda has no false "no" on any corpus.
- **FriCAS: 11 false "no".** All of them are in `gauss-elem` (259), integrands
  that are elementary by construction. FriCAS returns an unevaluated `integral`
  when the exponentials differ by a constant factor, such as `E^(1 + 2 I x)`
  beside `E^(I x)` or `E^(I (2 + x))`, and treats them as independent monomials.
  This is the same class of error that Mathilda's v0.6638 fixed in Mathilda.
- **FriCAS: 40 answers fail their derivative check** at the real points (1 in
  256, 39 in 257). Rechecked in FriCAS, 18 of them are formally correct
  (`normalize(D(r, x) - f) = 0`): the radicals are evaluated on a different branch.
  For the other 22, `normalize` does not reach 0; it cannot denest radicals like
  `12^(1/4)`, so those are unresolved rather than proven wrong.
- **Mathilda: 1 wrong answer**, `I x^2 Log[-I + Tan[2 x]]` (259, `tani-rand-9`),
  from Automatic's special-function stage. It predates the campaign (the same at
  v0.6639). The stage rewrites `Log[-I + Tan[2 x]]` as logarithms of
  exponentials, an identity valid on one branch interval only, so the PolyLog
  answer is wrong past the pole at `x = Pi/4`. Even `Integrate[Log[-I + Tan[2 x]], x]`
  is off by `2 Pi I` there. RischTranscendental certifies the case
  non-elementary, correctly.

**Speed.** The comparison uses the cases both systems close, with each system's
own integration time (start-up excluded).
- **Mathilda is faster on the typical case.** Its medians are 18–75 ms against
  FriCAS's 70–100 ms, and it is faster on 1297 of 1921 cases. FriCAS's timer has a
  floor of about 60 ms per call, so the medians favour Mathilda.
- **FriCAS has no slow tail.** None of its common cases takes over 1 s; 53 of
  Mathilda's do (`expcube`, compact radical-on-top integrands, Tan and Tanh
  towers). So Mathilda's 90th percentile is slower (259–383 ms against 110–150 ms),
  and so are its totals (310 s against 175 s).
- **The cascade costs time.** On 259, Automatic times out (30 s) on 27
  elementary cases that RischTranscendental closes in under 0.1 s. In the one
  profiled, the DerivativeDivides stage (12, before Risch) used 57.7 s of a 60 s
  budget.

**Fairness notes.** FriCAS's unevaluated `integral(...)` counts as its proof of
non-elementarity, except on the elementary-by-construction cases, where it is
counted as false. Its answers in Ei, li, erf and erfi are verified exactly,
because `complexNumeric` cannot evaluate them. The FriCAS runs on 256 and 257 are
the cached runs of 2026-11-12/13 (the corpora are unchanged); 259 ran on
2026-11-20.
