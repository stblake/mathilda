# NMinimize / NMaximize Benchmarking Implementation Plan

**Ticket**: DEMO-1
**Research**: `thoughts/shared/research/2026-08-21-DEMO-1-nminimize-nmaximize-benchmarking.md`

## Overview

Add two benchmark experiments covering `NMinimize`/`NMaximize` against the strongest
available competing implementations — SciPy 1.17.1 and Mathematica (`wolframscript`,
present on this host) — and fold them into `make bench-gap`.

## Current State Analysis

- `NMinimize` has 7 applied cases in `benchmarks/63-global-optimization/`; the four newer
  methods have `79`–`86`. **`NMaximize` appears on one line in the entire tree.**
- No experiment measures the default `Method -> Automatic` path (`nm_driver.c:394`: → DE).
- No experiment uses the canonical academic corpus (Rastrigin/Ackley/Griewank/…).
- `benchmarks/REPORT.md` and the last `history.jsonl` record are dated **2026-08-06**;
  `NMinimize` landed **2026-08-14**. `ABSENT.md:26` still says `_declared, not yet a row_`.

## Desired End State

Two new experiment folders, discovered automatically by `run_all.py`'s glob
(`run_all.py:100-123`, no manifest edit), producing rows in `REPORT.md`:

- `benchmarks/89-nminimize-nmaximize/` — the **timed race**, fair-comparison envelope
  applied. Feeds speed ratios.
- `benchmarks/90-nminimize-testbed/` — the **full corpus**, including landscapes where
  divergence is expected. Measures robustness, not speed.

This mirrors the existing `79/80`, `81/82`, `83/84`, `85/86` experiment/testbed pairing.

### Key Discoveries
- Discovery is a pure glob on `^(\d{2})-(.+)$`; `89`/`90` are free and land in group
  "D uncovered subsystems" via `group_of()` (`run_all.py:126-158`). No registry edit.
- Wire format is tab-separated: `BENCH\t<label>\t<ms>`, `CHECK\t<label>\t<value>`
  (`run_all.py:165`). `<label>` is the join key and must be **byte-identical** across
  `.m` and `.py`.
- `CHECK-FAIL` **discards timings** (`run_all.py:521-523`) — the value gate outranks
  every timing.
- The **same `.m` file** is fed to Mathilda and `wolframscript -file`
  (`run_all.py:400-407`), so it must be valid Wolfram Language as well.
- Mathilda's PRNG is a seedable SplitMix64, so `"RandomSeed"` makes runs reproducible —
  which is the only reason a `CHECK` gate is viable for a stochastic method at all.

## What We're NOT Doing

- Not modifying `NMinimize`/`NMaximize` source. This is measurement, not optimization work.
- Not modifying `run_all.py` or the harness. If the harness cannot express something,
  that is a finding, not a licence to change the scorer.
- Not regenerating the stale `REPORT.md` from 2026-08-06 as part of this ticket — a full
  `make bench-gap` is a ~20-minute job across all 90 experiments and is a separate concern.
  We run `--only 89,90` to validate, and note the staleness.
- Not benchmarking the eight methods individually; `79`–`86` already do that.

## Implementation Approach

**Resolved design decision (was research Open Question 4).** Split the work: `89` applies
the fair-comparison envelope so every reported ratio compares like with like; `90` runs the
full corpus and uses the `CHECK` channel to carry *solution quality* rather than the
objective value, so a robustness divergence surfaces as a real signal instead of a
discarded row.

**The `90` trick.** Instead of `check[label, objective]`, experiment 90 emits
`check[label, Boole[Abs[fbest - fstar] < tol]]` — "did you find the known global optimum?"
Then `CHECK-FAIL` means precisely *the systems disagree about whether they solved it*,
which is the robustness gap made legible. When all systems solve it, the timing is
legitimate and is reported normally.

**Engine matching in `89`.** Two families, so neither question is fudged:
- `D*` cases pin `Method -> "DifferentialEvolution"` in all three systems — a genuine
  engine-matched race against `scipy.optimize.differential_evolution`.
- `A*` cases use `Method -> Automatic` — what a real caller actually writes, measuring
  Mathilda's default against Mathematica's default and SciPy's closest analogue.

**Rounding.** `CHECK` uses `Round[10^4 f]`, not the `10^6` used by `15-optimization`.
Basins differ by far more than `1e-4`, so `10^4` still catches a wrong basin, while
`10^6` would CHECK-FAIL on last-digit differences between three different local-polish
implementations. This is a deliberate loosening and is documented in the README.

---

## Phase 1: Experiment 89 — the timed race

### Changes Required

**File**: `benchmarks/89-nminimize-nmaximize/nminimize_nmaximize.m` (new)
**File**: `benchmarks/89-nminimize-nmaximize/nminimize_nmaximize.py` (new)
**File**: `benchmarks/89-nminimize-nmaximize/README.md` (new)

Case families:
- `D1`–`D6` engine-matched DE: Rastrigin 2D/5D, Ackley 10D, Griewank 5D, Levy 5D,
  Rosenbrock 5D.
- `A1`–`A3` default path: Branin 2D, Six-hump camel 2D, Hartmann 6D.
- `M1`–`M3` **NMaximize** — the coverage gap. Includes a wrapper-overhead pair: the same
  problem posed as `NMinimize[f]` and `NMaximize[-f]`, whose ratio isolates the cost of
  the negation wrapper (`nm_driver.c:555-605`).
- `C1`–`C3` constrained: inequality, equality, mixed-integer (`Element[x, Integers]`).

Preamble per `benchmarks/README.md:243-247`: `Get["../harness.m"]` on the `.m` side,
`sys.path.insert` off `__file__` on the `.py` side. `require[{"NMinimize","NMaximize"}]`.

### Success Criteria

#### Automated Verification
- [ ] `./Mathilda -file benchmarks/89-nminimize-nmaximize/nminimize_nmaximize.m` exits 0
      and emits a `BENCH` and a `CHECK` line per case
- [ ] `python3 benchmarks/89-nminimize-nmaximize/nminimize_nmaximize.py` likewise
- [ ] Label sets match: `python3 benchmarks/run_all.py --only 89 --check-labels` exits 0
- [ ] `python3 benchmarks/run_all.py --only 89` produces zero `INCOMPLETE` rows
- [ ] Every case classifies `AHEAD` or `SLOWER` — no `CHECK-FAIL` (that is what the
      envelope is *for*; a `CHECK-FAIL` here means a case belongs in 90 instead)

#### Manual Verification
- [ ] Every excluded case is listed in the README's `Fair-comparison envelope` section
      with the basin each system found — the convention set by `81`/`85`
- [ ] Ratios are plausible against the `81`/`83`/`85` precedent (those report Mathilda
      ahead by 1.5×–650×); a wildly different number means a measurement bug, not a win

---

## Phase 2: Experiment 90 — the full-corpus testbed

### Changes Required

**File**: `benchmarks/90-nminimize-testbed/nminimize_testbed.m` (new)
**File**: `benchmarks/90-nminimize-testbed/nminimize_testbed.py` (new)
**File**: `benchmarks/90-nminimize-testbed/README.md` (new)

Full corpus including the landscapes deliberately excluded from 89: Schwefel 5D,
Eggholder 2D, Bukin N.6, Michalewicz 5D, Drop-wave 2D, Cross-in-tray 2D, Styblinski-Tang
5D, Rastrigin 10D. `CHECK` carries `Boole[Abs[fbest - fstar] < tol]` against the published
global optimum for each function.

### Success Criteria

#### Automated Verification
- [ ] Both halves run and emit matching labels: `run_all.py --only 90 --check-labels`
- [ ] Zero `INCOMPLETE` rows
- [ ] Each `CHECK` value is `0` or `1` on every system

#### Manual Verification
- [ ] Any `CHECK-FAIL` row corresponds to a *real* solve/no-solve disagreement, verified
      by reading the reported optima — not to a tolerance artefact
- [ ] README states the known global optimum and source for every function

---

## Phase 3: Integration

- [ ] `python3 benchmarks/run_all.py --only 89,90` completes; writes `REPORT.partial.md`
      only, leaving canonical `REPORT.md`/`history.jsonl` untouched (`run_all.py:1377`)
- [ ] Confirm `group_of(89)`/`group_of(90)` → "D uncovered subsystems"
- [ ] Note in the READMEs that `REPORT.md` is stale as of 2026-08-06 and that a full
      `make bench-gap` is required to surface these rows in the canonical report

## Testing Strategy

The harness *is* the test: `--check-labels` catches join failures, the `CHECK` gate
catches fast-but-wrong, and `INCOMPLETE` catches a head returning unevaluated. Beyond
that, the three-way agreement between Mathilda, SciPy and Mathematica on the objective
value is itself a correctness check on `NMinimize` that no unit test currently performs.

## Performance Considerations

Stochastic global optimizers are timed as min-of-3 after a warm-up (`harness.m:53-90`).
Seeds are pinned so the three timed reps solve the identical problem; without that, the
variance across reps would exceed the differences being measured.

## References

- Research: `thoughts/shared/research/2026-08-21-DEMO-1-nminimize-nmaximize-benchmarking.md`
- Harness contract: `benchmarks/run_all.py:100-123`, `:165`, `:514-589`
- Envelope precedent: `benchmarks/81-dual-annealing/README.md`,
  `benchmarks/85-basin-hopping/README.md`
- Implementation: `src/numerical_calculus/nm_driver.c`, `nm_de.c`, `findmin_nm_common.c`
