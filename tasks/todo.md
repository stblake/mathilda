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

## Review addendum (2026-08-08, part 3) — 12 new user stress cases

Ran 12 fresh adversarial cases; all now match documented WL, all efficient
(<7ms; no exponential blowup). Findings & fixes:
- Cases 8 (Flat + `x_?(Total[{##}]>10&)`) and 11 (top-level `(a|b|c|d)...` over a
  50-elt List) were ALREADY correct at `False` — `##` in the test is the single
  bound value, and a List is one expression, not a bare alternative.
- **Case 5 & the latent root of 11 — top-level `Repeated`/`RepeatedNull`.**
  `MatchQ[a, a..]` was False (should be True). Added an `is_repeated` handler at
  the top of `match_internal`, matching the single subject as a length-1 seq.
  This also fixes nested `(a...)..` and `{(((x:a...)...)..)}` (routed per-repetition
  through the same path). Case 5 -> True.
- **Case 12 — `Unique` was unimplemented** (`Table[Unique["sym"],{n}]` stayed
  head-`Unique`, so `x_Symbol` correctly failed). Implemented `Unique[]` /
  `["x"]` / `[x]` / `[{...}]` in src/modular.c (shared `$ModuleNumber`, fresh via
  `symtab_lookup`, `Temporary`). Case 12 -> True.
- **Perf:** trimmed a dead `subset` malloc+copy for plain unnamed/untyped
  `__`/`___` (only remainder recurses). ~15-20% on the duplicate-search pattern.

Corpus 258 -> 287 (sections 31-36, +29 gating cases). Verified: 287/287; all
pre-existing matcher/replace/dispatch suites pass; bench_match linear (~2.0);
valgrind definitely-lost == macOS libobjc baseline (13,440B/420 blocks), zero
frames from our code. Docs (scoping-constructs.md +Unique, pattern-matching.md
top-level-Repeated note, changelog) + memory updated.
