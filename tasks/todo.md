# Reduce: radical∘Abs composition + univariate domain-gate soundness

Plan: `~/.claude/plans/golden-herding-unicorn.md`

## Design B — univariate domain-gate per-conjunct scoping (fixes `Sqrt[Abs[x]]<1 → x==0`) ✅
- [x] Rework `reduce_univar_general` domain collection to per-conjunct arrays
- [x] Move gate inside `form_truth_general` per-conjunct loop; drop global gate
- [x] Keep breakpoints as a union; keep scanning past undecidable domain (nested-radical fix)
- [x] Free per-conjunct arrays at all exits (success + decline)
- [x] Verified: `Sqrt[Abs[x]]<1 → -1<x<1`, `Log[Abs[x]]<0` correct; nested radicals restored; corpus+units green

## Design A — radical rationalization pass (fixes multivariate silence) ✅
- [x] Add `reduce_stmt_has_radical` detector (+ prototype in `reduce_realfn.h`)
- [x] Implement `rationalize_tree`/`rationalize_relation` (NNF walk, isolation, 6-row table, decline policy)
- [x] Add `rationalize_radical_leaves` step to `reduce_piecewise_preprocess` fixpoint (threaded vars/nv)
- [x] OR the detector into the multivariate dispatch gate (`reduce.c`)
- [x] Verified flagship = correct region; per-relation rows; declines sound; sampled-equivalent; no regressions

## Tests ✅
- [x] Added 8 `"solved"` (mm-sqrt-*) + 2 univariate `"solved"` + 3 `"decline"` rows to `tests/reduce_corpus.m`
- [x] Pinned `Sqrt[Abs[x]]<1` and `Log[Abs[x]]<0` in `tests/test_reduce.c`
- [x] `reduce_corpus_tests`, `reduce_tests`, `solve_corpus_tests`, `solve_radicals_reals_tests` all green
- [x] valgrind: leak profile byte-identical to main; no leak allocated by new code (all pre-existing)

## Docs ✅
- [x] Updated `docs/spec/builtins/solutions-of-equations.md` (multivariate radical bullet + per-conjunct gate)
- [x] Added changelog section to `docs/spec/changelog/2026-08-24.md`
- [x] Added QE/Phase-6b note to `REDUCE_PLAN.md`
- [x] Rebuilt code-review graph
- [x] Recorded memory: radical-rationalization+domain-gate (project) + scan-past-undecidable (feedback)

## Follow-up — And/Or precedence bug (chased down the "pre-existing CAD bug")
- [x] Diagnosed: NOT a CAD bug — the tree was correct (sampling matched). A **parser+printer**
      precedence bug: `And` and `Or` both had precedence 2800, so `And[a,Or[b,c]]` printed as
      `a && b || c` (re-parses to `Or[And[a,b],c]`) and `a || b && c` parsed to `And[Or[a,b],c]`.
- [x] Fixed: lowered `Or` to 2700 (< `And` 2800) in `src/parse.c`, `src/print.c`, `docs/spec/operators.md`.
- [x] Verified roundtrip `ToExpression[ToString[e]] === e`; `x^2-y^2<1` now prints guards correctly.
- [x] Full test suite: 225/228 pass. The 3 failures (moebiusmu, primenu, interp) are PRE-EXISTING
      (fail identically on main; big-number factoring env-dependence + known interp issue) — zero new
      failures. `parse_tests` (incl. `test_unparenthesised_chains_still_chain`) and `boolean_tests` green.
- [x] Docs: `operators.md` split And/Or row; changelog entry added. Memory + graph updated.

## Review
**Outcome:** both original examples solve. `Reduce[Sqrt[Abs[x]]+Abs[y]<1,{x,y},Reals]` now
returns the correct region (was unevaluated); the newly-found univariate wrong answer
`Sqrt[Abs[x]]<1 -> x==0` is fixed to `-1<x<1`, and `Log[Abs[x]]<0 -> False` to the correct set.

**Files changed (engine):** `src/solve/reduce_realfn.{c,h}` (radical rationalization pass +
`reduce_stmt_has_radical` + `reduce_piecewise_preprocess` now takes vars/nv), `src/solve/reduce.c`
(dispatch gate), `src/solve/reduce_realdiag.c` (per-conjunct domain gate). Tests:
`tests/reduce_corpus.m` (+13 rows), `tests/test_reduce.c` (+2 pins).

**Verification:** reduce_corpus / reduce_tests / solve_corpus / solve_radicals_reals /
linearsolve all green; sampled-equivalence 0 mismatches; `make check-c99` clean; valgrind leak
profile byte-identical to main (no leak from new code).

**Key decision:** rationalize-by-square (keeps original dimension, reuses CAD) over aux-variable
purification (blocked on unbuilt Phase 7 QE + Phase 6b algebraic-coeff CAD fibres).

**Self-correction caught:** first per-conjunct gate broke nested radicals by breaking on the first
undecidable domain sign; fixed to keep scanning for a decidable failure (see feedback memory).

**Not done (out of scope / not requested):** git commit; book/ example addition; a pre-existing CAD
emission mis-association on `x^2-y^2<1` (noted, unrelated to radicals/Abs); MEMORY.md compaction.
