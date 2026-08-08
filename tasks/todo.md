# Pattern-Matcher Stress Test & Refinement

## Plan
- [x] Explore matcher, test infra, edge cases (3 Explore agents)
- [x] Build stress-test harness (fork-per-case corpus runner)
- [ ] Write ~230 asserted + ~25 observational stress cases (documented WL semantics)
- [ ] Run corpus, triage failures (correctness / missing-construct / robustness / perf)
- [ ] Fix: Verbatim, PatternSequence, top-level Longest/Shortest
- [ ] Fix: correctness divergences surfaced by asserted tier
- [ ] Fix: recursion-depth robustness (graceful, no false no-match)
- [ ] Fix: heap-allocate reorder storage (drop MATCH_REORDER_CAP 64); widen eval_guard_true stack
- [ ] Add bench_match.c (doubling-ratio perf gate)
- [ ] Docs: docs/spec/builtins + weekly changelog; re-run pre-existing matcher suites; valgrind

## Review
(to fill in at end)

## Review (2026-08-08)

Built a 228-case asserted conformance corpus + ~20-case observations corpus for
the pattern matcher, ran it, triaged, and fixed everything the asserted tier
surfaced (207/228 -> 228/228):

Correctness fixes (src/match.c):
- `_Rational`/`_Complex`/`_Integer`-bigint head-typed blanks (blank_head_matches)
- nested `Condition` guard eval (eval_guard_true recurses + evals leaves; drops 64 cap)
- `__?test`/`___?test` per-element PatternTest on sequence blanks

New constructs:
- `Verbatim[p]` (literal match), `PatternSequence[...]` named+unnamed,
  top-level `Longest`/`Shortest`. `OrderlessPatternSequence` interned+documented,
  not yet implemented (rare).

Robustness / perf:
- match recursion bounded by $RecursionLimit -> graceful non-match, no SIGSEGV
- reorder storage heap fallback removes the 64-element MATCH_REORDER_CAP
- bench_match.c doubling-ratio gate (all ops ~2.0, O(n))

Verification: asserted 228/228; pre-existing match/patterns/replace/rule_dispatch
suites pass; bench_match + bench_pack pass on a quiet machine; valgrind +48B over
baseline (one-time symbol interning, no per-call leak); make check-c99 clean.

Symbols/attrs/docs: sym_names (+2), attr.c (Protected x3), info.c docstrings,
docs/spec/builtins/pattern-matching.md, changelog 2026-08-03.md.

## Review addendum (2026-08-08, part 2)

Took on both remaining gaps + promoted all cases to gating unit tests:
- Implemented `OrderlessPatternSequence` (ops_unwrap/ops_assign backtracking;
  first-OPS-to-front; named binding; preds + trailing ___). 10 new asserted cases.
- Fixed the Orderless x two-sequence-blank blowup: last pattern element forced to
  k=n_exprs (must consume all) -> Plus[x__,y__] over 200 terms >10s -> ~0.5ms.
- Merged the report-only observations tier into the asserted corpus (deleted
  match_stress_observations.m + its ctest entry). Corpus now 258/258 GATING.

Verified: 258/258; all pre-existing matcher suites pass; bench_match + bench_pack
pass quiet; valgrind byte-identical to baseline over OPS/blowup paths (0 leaks);
docs (pattern-matching.md, changelog) + memory updated.
