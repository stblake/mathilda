---
ticket: STATS-1
created: 2026-08-27T18:00:00-04:00
researcher: msollami (beta run, auto mode)
source_sha: 15b088da
branch: main
repository: mathilda
topic: "Core probability & statistics gaps: what exists, what is missing, smallest coherent vertical to add"
tags: [research, codebase, statistics, quantile, src-stats]
subsystems: []
type: research
lifecycle: active
status: complete
last_updated: 2026-08-27
last_updated_by: msollami
---

# Research: Core probability & statistics gaps in mathilda

**Date**: 2026-08-27T18:00:00-04:00
**Researcher**: msollami (beta run, auto mode)
**Git Commit**: 15b088da
**Branch**: main
**Repository**: mathilda

## TL;DR
*(Human, ≤80 words)*
Asked which core Wolfram-style statistics builtins mathilda lacks and which gap is
smallest to close coherently. Found a mature `src/stats/` subsystem (16 builtins,
one-per-file) whose `Quartiles` already implements the full Wolfram 4-parameter
quantile engine — but general `Quantile` itself is absent and tracked as a known gap.
Recommendation: add the Quantile/dispersion family (Quantile, InterquartileRange,
MeanDeviation, MedianDeviation). Unresolved: nothing blocking; scope choices flagged
as assumptions below.

## Summary
*(Human, ≤200 words)*
Runtime probes (baseline binary, built clean at 15b088da) and source inventory agree:
Mean, Median, Variance, StandardDeviation, Quartiles, Covariance, Correlation,
Moment/CentralMoment, Skewness, Kurtosis, RootMeanSquare, Moving* exist and are exact
(Mean[{1,2,3,4}] → 5/2). Absent with zero occurrences in src/: **Quantile** (tracked in
benchmarks/ABSENT.md:30), InterquartileRange, MeanDeviation, MedianDeviation,
GeometricMean, HarmonicMean, Mode, TrimmedMean, CDF (PDF exists — the distribution
triad is one-third present).

The decisive discovery: `src/stats/quartiles.c:118-176` already implements Wolfram's
parameterized quantile `{{a,b},{c,d}}` definition (h = a+(n+b)q, clamp, floor,
interpolate), with q hardcoded to {1/4,1/2,3/4} and defaults {{1/2,0},{0,1}}. A general
`Quantile[data, q]` (Wolfram default `{{0,0},{1,0}}`), `Quantile[data, {q...}]`, and
`Quantile[data, q, {{a,b},{c,d}}]` is mostly a parameterization of proven code.
InterquartileRange, MeanDeviation, MedianDeviation then compose from existing builtins.
This is the smallest coherent vertical; it stays entirely inside `src/stats/`
conventions and closes a gap the repo itself declares.

## Open Questions
*(Human, no cap — a question list, not prose)*

### Unresolved
_None._

### Resolved
- [x] Which missing vertical to build? — Quantile + InterquartileRange + MeanDeviation +
  MedianDeviation (assumed — beta test; evidence: Quantile is the repo's own declared gap
  at benchmarks/ABSENT.md:30, its engine already exists in quartiles.c, and the other
  three are one-file compositions over existing Mean/Median/Quartiles machinery; the
  distribution/CDF vertical touches the inert-head ml/dist.c design and special
  functions — larger blast radius).
- [x] Which Quantile forms to support? — `Quantile[list, q]`, `Quantile[list, {q..}]`,
  `Quantile[list, q, {{a,b},{c,d}}]`, exact arithmetic, Wolfram default parameters
  {{0,0},{1,0}} (left-continuous: result = sorted[[Ceiling[n q]]], clamped). Matches
  Wolfram; Quartiles keeps its own {{1/2,0},{0,1}} default. Verified against Wolfram
  semantics: Quantile[{1,2,3,4},1/2] = 2 while Median = 5/2 — the two MUST differ.
- [x] Distribution arguments (Quantile[dist, q])? — Out of scope; the two shipped
  distributions are inert heads consumed by RandomVariate/PDF via strcmp
  (src/ml/dist.c:230,240). Decline (return NULL → unevaluated), like Median[dist] today.
- [x] Packed/NDArray fast path for the new heads? — Correctness path only: new heads NOT
  added to pack.c AWARE[] (packed args get boxed by the gate — correct, slower) and NOT
  to Compile ND_REDS (Quartiles precedent: AWARE but not ND_REDS, quartiles returns a
  list). Intentional omissions recorded via the audit tools' baseline/exempt mechanism
  if `tools/check_packed_aware.py` / `nd_fastpath_sweep.py` flag them (assumed — beta
  test; matches CLAUDE.md's "intentional omissions are declared, never silent").
- [x] Non-numeric elements? — Follow Median/Quartiles exactly: emit
  `Head::rectn: Rectangular array of real numbers is expected at position 1 in ...` via
  printf and return `expr_copy(res)` (median.c:53-57).
- [x] stats_tests never runs under ctest (add_executable at tests/CMakeLists.txt:1654-1655,
  no add_test anywhere) — pre-existing repo bug. Fix in-scope with one
  `add_test(NAME stats_tests COMMAND stats_tests)` line, since this change's tests extend
  the stats family and a suite that never runs would false-green the V phase (assumed —
  beta test; recorded for plan-deviation tracking).

## Requires Approval
*(Human, ≤100 words)*
_None._ — LOW risk tier (qualify.md: 10/100 projected), no architecture impact, no new
dependencies. Scope selections above are flagged "(assumed — beta test)" and are the
kind of call a human owner would normally ratify at gate 1.

---

## Research Question
"MISSION: fill core probability & statistics gaps (Wolfram's ProbabilityAndStatistics
guide as inspiration). Research the codebase first — mathilda likely already has
Mean/Variance-class basics; enumerate what exists, then pick the smallest coherent
missing vertical your research supports (candidates: Median/Quantile family;
Correlation/Covariance; a distribution object with PDF/CDF/RandomVariate for
Normal+Uniform; StandardDeviation edge semantics)."

## Detailed Findings

### What exists (runtime-verified on the baseline binary, then source-confirmed)
- Registered in `stats_init()` (src/stats/stats.c:13-46, called from core.c:807):
  Mean, RootMeanSquare, Median, Quartiles, Variance, Moment (+ATTR_NHOLDALL),
  CentralMoment, Skewness, Kurtosis, StandardDeviation, Covariance, Correlation,
  MovingAverage, MovingMedian, ExponentialMovingAverage — all ATTR_PROTECTED, none
  Listable (vector functions; threading would be wrong).
- Exact semantics: `Mean[{1,2,3,4}]` → `5/2`; `Skewness[{1,2,3,10}]` → `18/25 Sqrt[2]`;
  int64 accumulation uses the `ci_*_i64` family with symbolic fallthrough on overflow
  (src/stats/mean.c:98-108 documents a real wrong-answer bug fixed that way).
- List-layer relatives exist: Total, Accumulate, Commonest, Tally, MinMax, Differences,
  Rescale, Standardize; RandomReal/RandomInteger family; Erf/InverseErf.
- Distributions: `NormalDistribution`/`UniformDistribution` are INERT heads (attributes +
  docstring only, src/ml/dist.c:887,897) consumed by `RandomVariate` and `PDF` via
  strcmp dispatch (dist.c:230,240). `CDF` does not exist at all.

### What is absent (zero occurrences in src/, runtime-unevaluated)
- **Quantile** — declared-but-absent in benchmarks/ABSENT.md:30 and required by
  benchmarks/19-statistics. GeometricMean, HarmonicMean, MeanDeviation,
  MedianDeviation, InterquartileRange, Mode, TrimmedMean, WeightedData,
  CDF/InverseCDF/SurvivalFunction, Probability, Expectation, Binomial/Poisson/
  ExponentialDistribution, BinCounts, HistogramList.

### The quantile engine already exists
- src/stats/quartiles.c:99-102 — defaults a=1/2, b=0, c=0, d=1 when no parameter matrix
  given; :83-96 parses an optional `{{a,b},{c,d}}` second argument (already!).
- :118-176 — for each q: h = a + (n+b)q (exact Expr arithmetic via SYM_Plus/SYM_Times),
  edge clamps h<=1 → first, h>=n → last, j = Floor[h] clamped to [1, n-1],
  result = A[j] + (c + d(h-j))(A[j+1] - A[j]). This IS Wolfram's general definition;
  only the hardcoded q_vals {1/4,1/2,3/4} (:120-123) and the default parameter matrix
  distinguish it from Quantile.
- Sorting: `Sort[data]` evaluated via **pack_eval_plain, never evaluate** (quartiles.c:104-110,
  median.c:59-64) because Sort returns a PACKED list for large machine-number input and
  the subsequent `.data.function.args` walk would break. Comparison order is
  expr_compare (src/sort.c:324), the canonical total order incl. GMP rationals.
- Matrix input: columnwise via Transpose (quartiles.c:33-65 hand-rolls it to carry the
  extra parameter; median.c:34-39 uses stats_apply_columnwise(name, matrix) from
  stats_common.c, which builds Map[name, Transpose[m]] — works for any registered head
  by name string; expr_new_symbol accepts a plain C string).

### Conventions that will bite (verified in-source)
- Builtin contract: `Expr* builtin_x(Expr* res)`; return NULL to decline (never free
  res); fresh Expr* on success (docs/extending.md:10-131, six-step recipe).
- Message + unevaluated return: printf `"Head::rectn: ...\n"` then `return expr_copy(res)`
  (median.c:52-57) — NOT NULL — for wrong-typed elements.
- Docstrings: central `symtab_set_docstring` block in src/info.c (Median :4150,
  Quartiles :4152); terse, no inline examples. Examples go in
  docs/spec/builtins/statistics.md; weekly changelog docs/spec/changelog/<Monday>.md
  (this week: 2026-08-24.md) + a row in Mathilda_spec.md's table if the file is new.
- sym_names: SYM_* are interned `const char*` (sym_names.h:637,655; sym_names.c:1515,1533).
  New heads only need entries if referenced by pointer; expr_new_symbol("Quantile")
  works (stats_common.c stats_apply_columnwise passes string literals).
- Tests: CMake only; gcc-16 enforced (tests/CMakeLists.txt:20-25 FATAL_ERROR on clang).
  Three manual steps per new test: (1) add any new src file to the COMMON_SRC stats
  block :475-490 (explicit list — the link trap), (2) add_executable +
  target_include_directories + add_test 3-liner (pattern at :1024-1026),
  (3) assert_eval_eq(input, expected, is_fullform) from test_utils.h:16.
  **stats_tests is built but never registered with add_test** (:1654-1655, :3216 only)
  — pre-existing bug; README's `for t in *_tests` loop is what actually runs it today.
- Main binary picks up src/stats/*.c by WILDCARD (makefile:362) — no makefile edit
  needed for new stats files; only tests/CMakeLists.txt is manual.
- Fast-path surfaces (CLAUDE.md treats as correctness): pack.c AWARE[] stats block
  :553-566 (existing stats heads listed; deliberately NOT INT64_OK — exact results are
  Rationals); Compile ND_REDS at compile_ndtables.c:232 (Quartiles is AWARE but not in
  ND_REDS — returns a list); audit tools in tools/: check_packed_aware.py,
  nd_fastpath_sweep.py, nd_surface_audit.py, compile_coverage.py, check_array_exactness.py,
  check_c99_portability.py. Intentional omissions belong in the tools' exempt/baseline
  lists with a reason.
- C99: M_PI needs #ifndef fallback; feature-test macros before any include; gate is
  `make check-c99` + Linux CI. (No UTF-8-in-comments prohibition observed — em-dashes
  appear throughout src/, contrary to the mission briefing.)
- Multi-statement `-file` scripts mis-parse (runtime observation: 25 Print lines → 8
  garbled outputs); C-level assert_eval_eq tests are the trustworthy medium.

## Code References
- `src/stats/stats.c:13-46` — stats_init registration + attributes
- `src/core.c:807` — stats_init() call site
- `src/stats/quartiles.c:83-96` — existing {{a,b},{c,d}} parameter parsing
- `src/stats/quartiles.c:99-102` — Quartiles default parameters {1/2,0,0,1}
- `src/stats/quartiles.c:104-110` — Sort via pack_eval_plain (packed-return trap)
- `src/stats/quartiles.c:118-176` — the general quantile interpolation engine
- `src/stats/median.c:15-70` — decline/message/columnwise/sort conventions
- `src/stats/stats_common.h` — stats_is_real_numeric, stats_apply_columnwise, helpers
- `src/stats/mean.c:98-108` — exact int64 accumulation + overflow fallthrough
- `src/info.c:4150-4152` — docstring block for Median/Quartiles
- `tests/CMakeLists.txt:475-490` — COMMON_SRC stats block (link trap)
- `tests/CMakeLists.txt:1024-1026` — canonical 3-line test registration pattern
- `tests/CMakeLists.txt:1654-1655,3216` — stats_tests missing add_test (pre-existing bug)
- `tests/test_stats.c` — 622 lines, existing stats tests, assert_eval_eq style
- `src/pack.c:553-566` — AWARE[] stats block
- `src/compile/compile_ndtables.c:232` — ND_REDS
- `benchmarks/ABSENT.md:30` — Quantile declared absent
- `docs/extending.md:10-131` — the six-step builtin recipe
- `docs/spec/builtins/statistics.md` — examples home (369 lines)

## Architecture Insights
- One-builtin-per-TU under src/stats/, filename = lowercased head; stats.c is
  registration-only. New files are picked up by the main build wildcard automatically
  but must be listed explicitly in tests' COMMON_SRC.
- The subsystem prefers building Expr trees through the evaluator (SYM_Plus/SYM_Times +
  eval_and_free) over C doubles for exactness; only h's numeric value is read as double
  for clamping/floor. Quantile should do the same — exact in, exact out.
- Stats heads are conspicuously NOT Listable: Quantile[list, {q1,q2}] must thread over
  the SECOND argument manually (Wolfram behavior), not via ATTR_LISTABLE.
- The repo audits its own fast-path coverage; silent omission is the failure mode the
  tools exist to catch, so any deliberate skip must land in their baselines.

## Historical Context (from thoughts/)
- thoughts/shared/{plans,research,tickets} exist but are empty of statistics content;
  no subsystem docs (`thoughts/shared/subsystems/` absent) — subsystem lookup: none
  apply. tasks/*.md (repo's own planning notes) contain no quantile work.

## Related Research
- thoughts/shared/tickets/STATS-1/qualify.md — QRISPY Q receipt (10/100 LOW, projected).
