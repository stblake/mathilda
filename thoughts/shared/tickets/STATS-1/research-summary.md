---
ticket: STATS-1
created: 2026-08-27
type: research
lifecycle: active
status: complete
full_research: thoughts/shared/tickets/STATS-1/research.md
---

# STATS-1 research summary

**Question**: which core Wolfram-style statistics builtins does mathilda lack, and what
is the smallest coherent vertical to add?

**Answer**: the Quantile/dispersion family. `src/stats/` already has 16 exact-arithmetic
builtins, and `Quartiles` (src/stats/quartiles.c:118-176) contains Wolfram's full
4-parameter quantile engine with q hardcoded to {1/4,1/2,3/4}. General `Quantile` is
absent and tracked as a gap in benchmarks/ABSENT.md:30. `InterquartileRange`,
`MeanDeviation`, `MedianDeviation` are absent one-file compositions over existing
machinery.

**Key constraints**: builtin ownership contract (NULL declines, never free res);
message-then-`expr_copy(res)` for bad element types; Sort via `pack_eval_plain`;
COMMON_SRC explicit list in tests/CMakeLists.txt (link trap); docstrings central in
src/info.c; examples in docs/spec/builtins/statistics.md + weekly changelog; fast-path
audits require declared exemptions; gcc-16 only, C99, `make check-c99`.

**Found bug (pre-existing)**: `stats_tests` is built but never `add_test`-registered —
it never runs under ctest.

**Open questions**: none unresolved; scope choices flagged "(assumed — beta test)" in
the full document.
