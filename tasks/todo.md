# M63 — DSolve corpus §2.2.36 (3501–3600) + the parameter-as-indvar converter repair

Next hundred in the `DSOLVE_PLAN.md` campaign: upstream §2.1.36
(`Ch2.S1.SS36.htm`), Problems 3501–3600, 100 records / 100 scalar / 9 IVPs.

Converting it turned up a bug with reach far beyond the new section: the
converter promotes a lone **parameter** letter to the independent variable, so
`3570` (`y'' - 2a y' + a² y == 0`) became an ODE *in `a`*. Measured over all 36
committed corpora, **17 records are silent wrong equations of that class** — one
(`§2.2.16-1534`) already written down in `README.md` as an unexplained residue.
A mis-transcribed record scores PASS or UNEVAL against an equation the book never
asked, so this is a correctness bug in the gated corpora.

Two commits, each bumped and tagged (user's call).

## Commit 1 — converter fix + the 17-record repair

- [x] `detect_symbols`: restrict the step-3 lone-letter adoption to `y` (+ the
      existing Greek path) — every other lone Latin letter is a parameter
- [x] `detect_symbols`: under a `_missing_x` classification, take the indvar only
      from a function-argument position, else the fresh-letter fallback (needs
      `classif` plumbed through `convert_row` ← `main`)
- [x] Verify in isolation: old-vs-new converter on the SAME fetched page, **19**
      sections. Exactly 13 records move, all known victims; the other 10
      sections (incl. §2.2.1/20/21/24/25/29/30/35 as negative controls) are
      byte-identical
- [x] Patch the **17** record lines into the **9** committed `DE_examples_*.m`
      (record lines only — keep upstream LaTeX drift out of the diff). The 17th
      is §2.2.28-2789, found by the new audit, not by the hand scans
- [x] Re-measure those 9 sections before/after, same machine, back to back:
      **+3, zero lost, 0 FAIL in all 18 runs**; gates + `reports/*` + `STATUS.md` updated
- [x] `tools/check_corpus_indvar.py` + `make check-corpus-indvar`: green on the
      repaired tree, 24 findings on the pre-fix copies
- [x] `README.md`: the stale `Ch2.S2.SSN.htm` fetch URL, the "cannot be
      regenerated" note, and an audit-first instruction
- [x] v0.255, changelog, `DSOLVE_PLAN.md` M63a entry, tagged `v0.255`

## Commit 2 — the §2.2.36 wave

- [x] Generate `DE_examples_2236.m`; validate (100 records, all parse, trap greps
      clean, indvar census x/z/t only, `check-corpus-indvar` green)
- [x] `add_test(dsolve_corpus_2_2_36_tests)` at baseline 1 + `STATUS.md` block + `README.md` row
- [x] Baseline **96/100**, twice per-case identical; **99/100** after the fix, also
      twice identical; bucket report regenerated
- [x] The `3521`/`3527`/`3599` root cause, **re-diagnosed by measurement twice**.
      (i) Not a missing-answer case: `DSolve\`Linearizable` already solves them
      from the original equation in 0.11 s, three cascade slots after
      `Separable` — which eats the whole 8 s budget first on an `Integrate` of
      its own SAMPLED integrand. (ii) And not a normalisation case either: a
      `TrigExpand`-before-sampling rewrite was built, measured and **discarded**
      (faster at 0.08 s, but it claims the records for Separable's implicit twin,
      whose relation carries the sampling artefact `Cot[2]` where the cascade
      returns the explicit `ArcCos` Mathematica gives; and it half-expands
      multiple angles into something worse). What landed is a **6 s wall-clock
      deadline SHARED by both Separable entries** — a per-integral bound cannot
      work: 1 s loses §2.1.2-1134 (whose answer legitimately carries an
      unevaluated `Integrate` needing ~5 s to be DECIDED), 2 s leaves it on the
      boundary, and raising it walks the repaired records toward the wall because
      each path pays its own
- [x] `tests/test_dsolve_m63_stress.c`, five families, negative controls + a
      latency bound (the fix is a latency property, so an answer-only test would
      pass before AND after). The generator found three PRE-EXISTING unbounded
      steps on its own: `DSolve\`Exact`, `DSolve\`LieSymmetry` (~10 s) and
      `DSolve\`Homogeneous` (16 s on a scaled mixed-angle RHS), all identical on
      a pre-fix binary — named in `DSOLVE_PLAN.md`, not fixed here
- [x] Regression: **17 sections, zero regressions, one gain** (§2.2.16-1581); §2.1.2
      632 PASS vs 635 pre-fix, 0 FAIL both, all five movers solving on BOTH
      binaries over interleaved reps (boundary-band noise, not a regression);
      `check-c99` / `check-messages` / `check-corpus-indvar` green; valgrind
      13,440 B vs the pre-fix 13,632 B on identical input
- [x] v0.257 (0.256 taken by a concurrent notebook commit), docs, `STATUS.md`
      incl. M62's missing wave-history bullet, tagged `v0.257`

## Review

**Landed as two commits.** `v0.255` — the converter repair: eighteen corpus records
across nine sections were silently the wrong equation, verified against themselves by
the harness and therefore invisible. Two gates in `detect_symbols`, a new assert-empty
`make check-corpus-indvar` (which found the eighteenth victim on its first run, a
*system* record every hand scan had skipped), +3 PASSes with zero lost. `v0.257` — the
§2.2.36 wave: 96 → 99/100, one root-cause fix, residue 3579 alone.

**What this wave cost in wrong turns, and what each one taught.** The solver fix went
through three designs and measurement killed the first two. The symptom said Separable
could not see the split; it could not, but `Linearizable` had the answer three slots
later all along, so no capability was missing. The obvious repair then was to normalise
before sampling, which worked and was 75x faster — and was discarded because it made
the answer *worse* (an implicit relation carrying a `Cot[2]` sampling artefact in place
of the explicit `ArcCos` Mathematica gives). What landed is a shared 6 s deadline, whose
shape was also forced by a loss rather than by reasoning: a per-integral bound cannot
serve both §2.1.2-1134 (needs ~5 s to get an integral *decided*) and 3521 (needs
Separable to quit). Three process lessons are in `tasks/lessons.md`; the most expensive
was that a batch-ordered A/B showed a phantom +0.3 s regression that interleaving
removed entirely.

**Honest residue and named-not-done.** 3579 (`y=_G(x,y')`, `sympy=False`). Three
PRE-EXISTING unbounded steps the stress generator found on its own — `DSolve\`Exact`,
`DSolve\`LieSymmetry`, and a 16 s `DSolve\`Homogeneous` spin on a scaled mixed-angle
RHS, all identical on a pre-fix binary. The corpora also carry real pre-existing drift:
§2.2.4's gate was aspirational (3 cases at 5.7–11 s against an 8 s bound, identical on
both binaries), and §2.2.12-1133 / §2.2.13-1201/1219 regressed on main at some earlier
point. Tags are local — the branch also carries a concurrent session's v0.256, so
pushing is not mine to do.

---

# M64 — is the Integrate cascade ordered wrong?

Hypothesis (user): expensive stages may run before cheap ones in
`Integrate`'s Automatic cascade. Worth testing, but note two things first:

- The cascade is **already** deliberately ordered, and not by cost: *confidence
  first* — deterministic correct-by-construction methods ahead of
  search-and-verify ones and ahead of Risch-Norman's complex-log forms. That
  rule came from a user correction ("Weierstrass should run before deriv-divides
  as it's a domain-specific algorithm guaranteed to succeed without needing to
  check the result"), so sorting by cost would partly invert it.
- Reordering is correctness-neutral but **performance-significant**: moving
  `try_goursat` last in v0.191 exposed latent grind in the stages that had been
  short-circuiting it, and needed `has_pseudoelliptic_radical` gates on
  `linearity`, `derivdivides` and PMT before it was safe.

And the case that prompted this is NOT an ordering problem: `try_weierstrass`
(`integrate.c:1577`) already runs immediately before `try_derivdivides` (1578).
The 20 s goes to derivdivides because Weierstrass *declines in microseconds* —
its gate `wj_has_kernel_in_denominator` needs a literal `Power[base, negative]`
and the canonicaliser writes `Cos[x−a]/Sin[x]` as `Times[Csc[x], Cos[…]]`.

## Plan

## Result (measured 2026-10-02) — the order is NOT the problem; one stage is

**The hypothesis does not survive the measurement.** Decline corpus: of 150
candidate integrands (Charlwood 50 + the overintegration 100), the Automatic
cascade closes **125**, declines 9 in 0.5–6.3 s, and aborts 16 at the 20 s wall.
Timing every pinnable stage on each of those 25 non-closing integrands, in its own
process under a 15 s bound:

| pos | stage | total (s) | mean | capped | share |
|---|---|---:|---:|---:|---:|
| 02 | BronsteinRational | 0.11 | 0.004 | 0 | 0.0% |
| 03 | LinearRadicals | 0.00 | 0.000 | 0 | 0.0% |
| 04 | QuadraticRadicals | 0.02 | 0.001 | 0 | 0.0% |
| 05 | LinearRatioRadicals | 0.01 | 0.000 | 0 | 0.0% |
| 06 | ChebychevAlgebraic | 0.00 | 0.000 | 0 | 0.0% |
| 10 | Weierstrass | 0.06 | 0.002 | 0 | 0.0% |
| **11** | **DerivativeDivides** | **200.65** | **8.72** | **13/23** | **87.0%** |
| 12 | RischTranscendental | 15.02 | 0.65 | 1 | 6.5% |
| 13 | CRCTable | 0.78 | 0.034 | 0 | 0.3% |
| 14 | ParallelMixedTower | 6.13 | 0.27 | 0 | 2.7% |
| 15 | GoursatAlgebraic | 0.05 | 0.002 | 0 | 0.0% |
| 16 | ParallelMixedSpecial | 7.81 | 0.34 | 0 | 3.4% |

**There is no expensive-before-cheap pair to fix.** The six stages ahead of
position 11 cost **0.19 s in total** across all 25 integrands — they are already
first and they are already free. The stages *after* position 11 are cheap too
(Goursat 0.002 s mean, CRCTable 0.034, PMT 0.27, PMS 0.34). Reordering has nothing
to work with, and sorting by cost would have inverted the confidence-first rule for
no gain.

**`DerivativeDivides` is 87% of the entire decline cost** — 200.65 s of 230.6 s,
mean 8.72 s per integrand it cannot close, and it **fails to finish inside 15 s on
13 of 23**. That is the same defect shape as M63's `DSolve\`Separable`: an
unbounded Eliminate/Solve search in a cascade stage, charged to every input it
declines. `RischTranscendental` is a distant second (one case of 23 caps).

### CORRECTION — the pinned profile measured the wrong thing

The table above is the PINNED surface (`Method -> "DerivativeDivides"`), and
attributing it to the cascade was wrong twice over. Fixed by building
`MATHILDA_INTEGRATE_PROFILE=1` (src/calculus/integrate.c), which times each stage
from inside the cascade at the outermost frame, and attributes a stage killed by
the caller's `TimeConstrained` (a siglongjmp skips the normal accounting, so the
most expensive stage on every aborting integrand was invisible).

In-cascade, over the same 25-integrand decline corpus:

| pos | stage | total (s) | mean | share |
|---|---|---:|---:|---:|
| 01-08 | the eight cheap exact stages | 0.14 | ~0.001 | **0.0%** |
| **09** | **Linearity** | **120.39** | 4.82 | **43.5%** |
| 10 | Weierstrass | 0.04 | 0.002 | 0.0% |
| **11** | **DerivativeDivides** | **130.86** | 6.23 | **47.3%** |
| 12 | RischTranscendental | 19.97 | 1.54 | 7.2% |
| 13-16 | CRCTable / PMT / Goursat / PMS | 5.48 | - | 2.0% |

**`try_linearity` is 43.5% and has NO pinned surface**, so no outside measurement
could see it. With DerivativeDivides that is 90.8% in two stages. Both are
expensive for the same reason -- speculative recursion into the FULL cascade,
per `Plus` term and per candidate kernel respectively (the cascade switch runs at
every depth). **The ordering hypothesis is answered: the order is not the
problem.** The eight stages ahead of the hot pair cost 0.0%.

Two further corrections to the record: the cascade DOES run the Eliminate/Solve
search (`try_derivdivides` sends only a *pseudo-elliptic* integrand down the
direct-only path), and the test suites print `FAIL:` lines in their body while
still printing "All tests passed!" and exiting 0 -- reading the tail line is not
reading the result.

### What was tried, measured, and NOT shipped

1. **The recursion gate** (the user's proposal: recurse only when the reduced
   integrand is polynomial, rational, or algorithmically convertible). Corpus
   effect was good -- closes 125 -> 128, aborts 16 -> 9, five aborts becoming
   closes at 3.6-10.8 s. But it LOSES FOUR: O079/O080 of the overintegration
   corpus and, in the test suite, `Integrate[x ArcSin[x]/Sqrt[1-x^2], x,
   Method -> "DerivativeDivides"]` plus its ArcCos twin, which go from a verified
   answer to unevaluated. The predicate rejects reduced integrands that are
   transcendental in u, and those are this method's core business. **Disabled at
   its two call sites; the predicate is kept (MATHILDA_MAYBE_UNUSED) with this
   note.** The idea is right; the class needs deriving from the actual reduced
   integrands, case by case, not guessed.
2. **A stage-boundary `TimeConstrained` wrapper**, which DID make TimeConstraint a
   hard bound (3 s measured 3.00 s, 1 s measured 1.00 s). Removed: it cost an
   answer -- an integrand closing in 0.066 s became a 3 s decline. Rationale left
   as a source note in integrate.c.

### Open leads

- **Reaching the stage through its registered head is ~680x faster on one
  integrand**: `Integrate`DerivativeDivides[Sqrt[Sin[x]]/(Sin[x]^2+1), x]` closes
  in 0.066 s where the direct C call `integrate_derivdivides_full` takes 45 s and
  aborts. Either a real optimisation in the eval path (memoisation /
  canonicalisation) or the two routes are not doing equivalent work. Understand
  this BEFORE attempting a hard bound again.
- `try_linearity`, the larger of the two hot stages, is untouched.
- The gate's class predicate (above).

### Superseded: the bound decided from the pinned profile

Timing `Method -> "DerivativeDivides"` on the 125 integrands the cascade closes
makes the case worse than the decline table alone showed, and the fix free:

- On those 125 it spends **138.6 s**, of which only **33.8 s** is on its own 45
  closes — **104.8 s is spent on 80 integrands it declines or aborts** and a later
  stage then closes. Adding the 200.65 s from the non-closing set: of roughly
  **339 s total, ~306 s (90%) is spent on integrands it does not close.**
- Its close cost is tight: min 0.003 s, **median 0.40 s, p90 1.72 s**, and only
  seven of 45 closes exceed 1 s.
- A **3 s** bound therefore loses exactly ONE close of 45 — `O069` at 11.31 s —
  and saves 65.1 s on six cases in the closing set alone, plus capping the 13
  cases that currently do not finish in 15 s.
- **And that one close is free to lose:** `O069` is closed by
  `ParallelMixedTower` (position 14) in **0.21 s**, i.e. 54x faster, three stages
  later. So a 3 s bound on `DerivativeDivides` costs **zero closes**.

- [x] Decide the bound by measuring the close cost — 3 s, zero closes lost.
- [ ] Bound the **cascade call site**, not the shared routine — the documented rule
      from v0.191: `Method -> "DerivativeDivides"` is a deliberate request for the
      full search and must keep it. Mirror M63's shape (`sep_remaining`) if a
      deadline is wanted rather than a per-call `TimeConstrained`.
- [ ] A/B the Integrate suites: `integrals_tests`, `integrate_derivdivides_tests`,
      `integrate_dispatch_tests`, the CRC corpus, and the Charlwood 50 — the bound
      can only turn a slow close into a decline, so the suites are the gate.
- [ ] Only then revisit `wj_has_kernel_in_denominator` (separate change, alters the
      *spelling* of existing answers).

## Original plan (superseded by the result above)

- [x] Build a **decline corpus** — a stage is charged for what it DECLINES, not
      what it closes. Candidates: `mixed/charlwood_mathilda.json` (50,
      `integrand_wl`) + `mixed/overintegration_corpus.json` (100, `f_plus_g`).
      Classify each: closes / declines / times out.
- [ ] For each declining integrand, time every cascade-relevant pinned
      `Integrate\`<Method>` in its own process under a bound. Bias to state: the
      pinned surface keeps its full search where the cascade call site may be
      gated, so pinned decline time ≥ cascade decline time.
- [ ] Table: stage × total decline seconds × cascade position. An expensive
      stage ahead of a cheap one is the finding; no such pair means the order is
      already right and the real lever is the gates.
- [ ] Only if the table shows a genuine inversion: instrument in-cascade
      (`MATHILDA_INTEGRATE_PROFILE=1`, mirroring the `MATHILDA_PACK_DIAG`
      convention — lazily-read env var, per-stage accumulators, atexit dump) to
      confirm before moving anything.
- [ ] Separately: the `wj_has_kernel_in_denominator` widening is its own change
      (it alters the *spelling* of existing answers — `Integrate[Csc[x], x]`
      moves to `Log[Tan[x/2]]`), so it needs its own A/B against the Integrate
      suites.
