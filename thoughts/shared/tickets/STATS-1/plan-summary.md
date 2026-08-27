---
ticket: STATS-1
created: 2026-08-27
type: plan
lifecycle: active
status: draft
full_plan: thoughts/shared/tickets/STATS-1/plan.md
---

# STATS-1 plan summary — Quantile & dispersion family

**Recommendation**: add `Quantile` (Wolfram 4-parameter form, list-of-qs, matrix
columnwise), `InterquartileRange`, `MeanDeviation`, `MedianDeviation` to `src/stats/`,
by extracting the quantile engine that already lives inside `Quartiles`
(src/stats/quartiles.c:118-176) into `stats_common.c` and parameterizing it.

**Options considered**: CDF completion for the two inert distributions (rejected —
larger blast radius, touches ml/dist.c dispatch design); copying the engine instead of
extracting (rejected — drift risk); C-double arithmetic (rejected — subsystem exactness
discipline).

**Decision criteria**: repo-declared gap (benchmarks/ABSENT.md:30), smallest coherent
vertical, LOW risk tier (qualify.md 10/100 projected), engine reuse over rewrite.

**Decisions**: engine extraction with the integer-h/Ceiling correction (Wolfram
semantics; Quartiles provably unchanged); correctness path only — NDArray/packed input
materialized in-head via pack_unpack, omissions declared to the audit tools (no kernels
yet); fix the pre-existing `stats_tests`-never-runs bug in scope; decline
symbolic/distribution args. Scope calls are "(assumed — beta test)".

**Non-goals**: distribution arguments, CDF family, GeometricMean/HarmonicMean/Mode/
TrimmedMean, packed fast paths, Association support, benchmark rows, any Quartiles
behavior change.

**Open Questions — Unresolved**: _None._ (three resolved with flagged assumptions)

**Requires Approval**: _None._ — LOW tier, no architectural impact.

**Architecture Impact**: none across the board. **Subsystems**: none.

Full detail, AC table (16 rows), phases 1-4: see plan.md (appendix).
