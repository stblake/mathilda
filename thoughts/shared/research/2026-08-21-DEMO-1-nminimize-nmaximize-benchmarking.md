---
date: 2026-08-21T00:00:00-04:00
researcher: Michael Sollami
git_commit: b614d1ed79fc84e9b944a21830a8a42b86ab1c99
branch: claude/where-are-you-5f1147
repository: mathilda
topic: "Benchmarking NMinimize/NMaximize against the strongest competing implementations"
tags: [research, codebase, nminimize, nmaximize, benchmarks, optimization, DEMO-1]
status: complete
last_updated: 2026-08-21
last_updated_by: Michael Sollami
---

# Research: Benchmarking NMinimize / NMaximize

**Date**: 2026-08-21
**Researcher**: Michael Sollami
**Git Commit**: b614d1ed79fc84e9b944a21830a8a42b86ab1c99
**Branch**: claude/where-are-you-5f1147
**Repository**: mathilda (github.com/stblake/mathilda)
**Ticket**: DEMO-1

## Research Question

How do `NMinimize`/`NMaximize` work, how does the `make bench-gap` suite ingest new
experiments, and what would a fair benchmark against the strongest competing
implementations look like?

## Summary

Three findings reframe the task.

**1. `NMinimize` is already benchmarked; `NMaximize` is not.** `benchmarks/63-global-optimization/`
carries seven cases (A1–A7), all *applied/industrial* problems (refinery pooling, VWAP
tracking, risk parity). Eight further folders (`79`–`86`) benchmark the four newer methods
against their `scipy.optimize` counterparts. But `NMaximize` appears on exactly **one line**
across the entire `benchmarks/` tree. The maximization wrapper — including the sign handling
that is its only real logic — is effectively untested for performance.

**2. Nothing benchmarks the default path a user actually hits.** Every existing optimization
experiment pins an explicit `Method`. No experiment measures `NMinimize[f, vars]` with
`Method -> Automatic`, which is what `Automatic → DifferentialEvolution` (`nm_driver.c:394`)
resolves to and what any real caller writes. The canonical academic corpus (Rastrigin,
Ackley, Griewank, Schwefel, Levy, Michalewicz, Hartmann, Shekel, Branin, Six-hump camel,
Styblinski-Tang, Eggholder) is absent — the applied cases in `63` are not a substitute,
because they do not probe multimodality in a controlled way.

**3. The aggregate report has never seen any of this.** `benchmarks/REPORT.md` and the last
`history.jsonl` record are both dated **2026-08-06**. `NMinimize` landed **2026-08-14**
(`docs/spec/changelog/2026-08-10.md:1110`), and the four newer methods on **2026-08-17**.
`benchmarks/ABSENT.md:26` still reads `` `NMinimize` | _declared, not yet a row_ ``. So a
full `make bench-gap` has not run in 15 days, and every optimization benchmark added since
exists on disk but has never appeared in a report.

## Detailed Findings

### NMinimize / NMaximize implementation

- Registration: `src/numerical_calculus/findmin.c:12-30`. Both `ATTR_PROTECTED` only —
  deliberately **not** `HoldAll` (`findmin.c:17-25` explains: holding broke `Method`
  sub-option values like `"RandomSeed" -> s`). Note the changelog entry at
  `docs/spec/changelog/2026-08-10.md:1110` still claims `HoldAll, Protected`; that was
  corrected later in the same week. **The changelog is stale relative to the code.**
- Driver: `src/numerical_calculus/nm_driver.c`, shared core `findmin_nm_common.c` (1517 lines).
- Eight methods (`findmin_internal.h:205-206`), dispatched at `nm_driver.c:423-433`:
  `DifferentialEvolution` (default), `NelderMead`, `RandomSearch`, `SimulatedAnnealing`,
  `SHGO`, `DualAnnealing`, `DIRECT`, `BasinHopping`.
- `Automatic` resolves to DE (`nm_driver.c:394`). DE is DE/rand/1/bin, population
  `NP = 10n` clamped `[15,200]`, `F=0.6`, `CR=0.9`, bounce-back on bound violation
  (`nm_de.c:66-139`), followed by multi-start BFGS polish on up to `min(2n,50)` members
  (`nm_de.c:180-251`).
- Constraints use **Deb's feasibility rules** (`nm_better`, `findmin_nm_common.c:287-293`),
  not penalty-weight tuning. Supports `==`, inequalities, chained `Inequality`, `And`, `Or`
  disjunctions, and `Element[x, Integers|Reals]`.
- **`NMaximize` is a thin negation wrapper** (`nm_driver.c:555-605`): negates the objective
  to `Times[-1, f]`, builds a synthetic `NMinimize` call, negates the returned optimum
  (handling both `EXPR_REAL` and `EXPR_MPFR`). This is exactly why it is worth benchmarking
  separately — the wrapper allocates a rewritten expression tree per call, and that cost is
  currently unmeasured.
- Determinism: SplitMix64 PRNG, seedable via `Method -> {..., "RandomSeed" -> n}`. Results
  are reproducible, which is what makes a `CHECK` gate viable at all for a stochastic method.
- Objective auto-compilation: `nm_driver.c:326-360` compiles objective and constraints to
  bytecode for the machine-precision trial loop, falling back to the interpreter on
  non-lowerable bodies. MPFR `WorkingPrecision` bypasses this entirely (`nm_driver.c:458-459`).

### Benchmark harness contract

- `make bench-gap` → `python3 benchmarks/run_all.py` (`makefile:576-577`).
- Discovery is a **pure directory glob**, no manifest: `EXP_RE = ^(\d{2})-(.+)$` over
  `benchmarks/`, taking the first `.m` and first `.py` alphabetically
  (`run_all.py:100-123`). Adding `benchmarks/89-<slug>/` requires **zero registry edits**.
- `group_of()` (`run_all.py:126-158`) buckets 1–10 A, 11–20 B, 21–28 C, 53–62 E, everything
  else → "D uncovered subsystems". `89` lands in D, alongside the whole 63–88 series.
- Wire format is tab-separated tagged lines on stdout, parsed by one regex
  (`run_all.py:165`): `BENCH\t<label>\t<ms>`, `CHECK\t<label>\t<value>`,
  `REQUIRE\t<head>\tpresent|absent`, `SKIP\t<label>\t<head>`.
- **`<label>` is the join key and must be byte-identical across `.m` and `.py`.** Mismatches
  are reported as a symmetric-difference "label mismatch" (`run_all.py:1358-1363`).
- Classification (`classify()`, `run_all.py:514-589`): `ABSENT` (Mathilda emitted `SKIP`) →
  `INCOMPLETE` (no Mathilda `BENCH` line) → `CHECK-FAIL` (systems disagree; **timings are
  discarded**) → `SLOWER`/`AHEAD` on `mathilda_ms / min(others)` with `SLOWER_AT = 1.5`.
- `CHECK` agreement is `CHECK_RTOL = 1e-6` relative to max magnitude (`run_all.py:489-507`).
- Timing: one untimed warm-up then **min of 3** (`$BenchReps`), `AbsoluteTiming` /
  `perf_counter`, never `Timing[]` (`docs/design/performance.md`, `harness.m:53-90`).
- `wolframscript` is autodetected and optional (`run_all.py:70-81`). **It is present on this
  machine** (`/usr/local/bin/wolframscript`), so a genuine three-way race is available.
- The **same `.m` file** is fed to both Mathilda and `wolframscript -file`
  (`run_all.py:400-407`). The script must therefore be valid Wolfram Language too.

### Prior methodology conventions

- `docs/design/performance.md` — wall clock via `AbsoluteTiming`, min of several reps after
  a warm-up, max recorded as a caching tripwire; **median not mean** when aggregating.
- `benchmarks/REPORT.md` methodology block — `SLOWER` (carries a ratio) and `ABSENT` (never
  carries a ratio) are **never pooled**, explicitly so a missing feature cannot be laundered
  into a fake speed multiplier.
- `benchmarks/81-dual-annealing/README.md` and `85-basin-hopping/README.md` each define a
  **"Fair-comparison envelope"**: cases are *excluded from the timed race* when the two
  implementations land in different basins, on the grounds that a 1e-6 check race would then
  be comparing different answers rather than different speeds. **This convention is
  per-experiment and hand-curated — it is not implemented in the harness.**
- `plans/HPC_IMPROVEMENT_PLAN.md` — "N4: Parity is a gate, not a report"; warns against
  changing the exactness contract to win a benchmark.

## Code References

- `src/numerical_calculus/findmin.c:12-30` — registration and attributes
- `src/numerical_calculus/nm_driver.c:394` — `Automatic` → DifferentialEvolution
- `src/numerical_calculus/nm_driver.c:555-605` — `NMaximize` negation wrapper
- `src/numerical_calculus/nm_driver.c:326-360` — objective auto-compilation
- `src/numerical_calculus/findmin_nm_common.c:287-293` — Deb feasibility rule
- `src/numerical_calculus/nm_de.c:66-251` — DE parameters and multi-start polish
- `benchmarks/run_all.py:100-123` — glob discovery, no manifest
- `benchmarks/run_all.py:514-589` — classification
- `benchmarks/harness.m:53-90` / `benchmarks/harness.py:37-129` — emission helpers
- `benchmarks/63-global-optimization/global_optimization.m:18` — `require[{"NMinimize","NMaximize"}]`

## Architecture Insights

The harness's value gate (`CHECK` before timing) is the load-bearing design decision, added
after an incident where `RandomReal[{}, dims]` made 8 cases "run" in 0.1 ms while computing
nothing (`run_all.py:37-40`). For a **stochastic global optimizer** that gate does more work
than it was designed for: it is the only thing preventing "gave up early" from reading as
"fast". But it is a binary agreement test at 1e-6, so it cannot express the thing that
actually matters for this class of algorithm — *cost to reach a given solution quality*.
That is the gap the `81`/`85` READMEs papered over by hand-curating a fair-comparison subset.

## Related Research

- `thoughts/shared/research/2026-08-18-unitbox-builtin.md` — unrelated (UnitBox), but the
  only prior artifact of this workflow in the repo.

## Open Questions

1. **RESOLVED** — Experiment number: `89` is the next free integer (`88-thue-equations` is
   the current max). Lands in group "D uncovered subsystems".
2. **RESOLVED** — Competitor set for the timed race: `scipy.optimize.differential_evolution`
   is the engine-matched competitor for Mathilda's default path, and Mathematica's own
   `NMinimize`/`NMaximize` runs the identical `.m` file. Both are available here.
3. **RESOLVED** — `NMaximize` needs first-class coverage; it currently has one line.
4. **OPEN — blocks planning.** *Should experiment 89 adopt the hand-curated
   "fair-comparison envelope" convention from `81`/`85` (drop cases where Mathilda and SciPy
   converge to different basins, so every reported ratio compares like with like), or should
   it report the full corpus and let `CHECK-FAIL` rows stand as honest evidence of a
   robustness gap?*

   These are not cosmetic variants. The envelope produces a flattering, methodologically
   clean speed table that silently omits the hardest landscapes. The full corpus produces
   `CHECK-FAIL` rows whose timings the harness **discards** (`run_all.py:521-523`), so the
   robustness gap shows up as an absence of data rather than as a measured result — and
   `REPORT.md` will surface those as check disagreements against Mathematica, which reads
   as a correctness defect rather than a known property of stochastic search.

   There is no repo convention that settles this: `63-global-optimization` does neither
   (it uses applied problems with single basins), while `81`/`85` chose the envelope without
   recording a rationale that generalizes. **Requires a human decision.**
