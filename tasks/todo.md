# Task: Implement `Refine[expr, assum]` — full Mathematica parity  ✅ COMPLETE

Plan: /Users/user/.claude/plans/let-s-plan-the-implementation-woolly-kitten.md

## Phase 0 — scaffold + wiring
- [x] src/refine.c + src/refine.h (orchestrator: parse, rewrite path, identity)
- [x] core.c refine_init(); sym_names SYM_Refine; info.c docstring; Options; ATTR_PROTECTED
- [x] tests/CMakeLists.txt (COMMON_SRC + refine_tests); build + smoke test
## Phase 1 — predicate decision path (refine.c)
- [x] export element_decide; var collector + Reduce entailment; Element/Equal/inequality dispatch
## Phase 2 — shared-engine gaps (also improve Simplify)
- [x] G1 a^p b^p (rule, not simp_power.c); G2 (x^m)^r + nonneg bucket; G2b (a^b)^c |b|<1; G3 Log[x<0]; trig-shift Cos/Sin/Tan[x+kPi]
- [x] G5/G6 assume_structural_rewrite (Floor/Ceiling/Round/IntegerPart/FractionalPart/Mod)
- [x] G7 prov_re/prov_int compound Element (Gamma/Factorial/LogGamma/Pochhammer, Floor->int, int^nonneg-int); inequality⟹real
## Phase 3 — deep positivity post-pass (refine.c: Sign/Abs/Sqrt via Reduce)
- [x] wire + guard (node cap 24, deadline)
## Phase 4 — docs, tests, housekeeping
- [x] extensive test_refine.c (all pass); valgrind (leak-clean over baseline); docs/spec/builtins/simplification.md + changelog; version 0.196->0.197

## Review
- Full Mathematica parity: every example in the Refine spec reproduces (see
  smoke1/smoke2 + test_refine.c). Architecture held: refine.c is a thin
  orchestrator; all capability gaps landed in the SHARED engine
  (simp_assume.c, simp_assume_rewrite.c) so Simplify benefits too — only the
  Reduce-backed positivity wiring is Refine-local.
- No regressions: simp/simplify/fullsimplify(+corpus)/element/possiblezeroq/
  zero_test/reduce(+corpus)/simp_log/denest/cuberoot/hang suites all PASS.
- Memory: valgrind `definitely lost` identical to the Simplify baseline
  (13,440 B / 420 blocks = macOS libobjc noise); 0 leaks cite refine.c /
  apply_assumption_rules / simp_assume after switching `evaluate(fresh)` ->
  `eval_and_free(fresh)`.
- Bugs caught in-flight: (1) pattern vars written with trailing `_` on a rule
  RHS (`p_`, `base_`, `c_`) leak into the result — RHS must use the bare name;
  (2) a symbol-collector that recurses into function HEADS mis-collects
  Plus/Times/Complex as "variables"; (3) `evaluate()` borrows its arg -> use
  `eval_and_free` for freshly-built exprs. See tasks/lessons.md.
- Left tagging/commit to the user per CLAUDE.md (bump to v0.197 staged in
  src/version.h; tag `v0.197` on commit).

---

# Task: Chapter 10 — "The Internals of Mathilda"  ✅ COMPLETE

Plan: /Users/user/.claude/plans/let-s-plan-the-writing-vectorized-petal.md

## Setup / grounding
- [x] Read book infra: `mathilda.sty`, `08-compilation.tex`, `appendix-architecture.tex`
- [x] Verify exact C for excerpts: `expr.h`, `attr.h`, `match.h`, `symtab.h`
- [x] `./Mathilda` binary present; cross-ref labels confirmed

## Infra
- [x] Add `ccode` listing style to `book/mathilda.sty` (blue, "not build-verified" label)

## Example sessions (`book/examples/10-internals/*.m`) — 18 files, all generated
- [x] expr-fullform, expr-bignum
- [x] eval-trace, eval-fixedpoint, eval-attributes
- [x] nonstandard-hold, nonstandard-scoping
- [x] matcher-basics, matcher-downvalues
- [x] rules-replace, symtab-context
- [x] memory-inuse
- [x] parser-fullform, printer-forms
- [x] repl-history (pipe leaves %/Out/$Line symbolic — used as the layering proof)
- [x] numbers-tower, packed-gate, twotier-rules

## Chapter prose (`book/chapters/10-internals.tex`) — 10 sections
- [x] Opener + pipeline signpost (references appendix fig:pipeline)
- [x] 10.1 Expression trees (+ ExprType/Expr C excerpt)
- [x] 10.2 The evaluator (+ step-order figure, hold C excerpt, attribute table)
- [x] 10.3 Non-standard evaluation and scoping
- [x] 10.4 The pattern matcher (+ MatchEnv/env_rollback excerpt, backtracking figure)
- [x] 10.5 Rules and the symbol table (+ Rule excerpt; corrects the "front" myth)
- [x] 10.6 Memory management (+ expr_copy/expr_free excerpt)
- [x] 10.7 The parser and the printer
- [x] 10.8 The read–eval–print loop and session state
- [x] 10.9 Numbers and the packed-array substrate
- [x] 10.10 Bootstrapping: startup + two-tier design
- [x] Closing (forward-ref Ch 11/12)

## Verify
- [x] `make examples` — all 18 transcripts generated; read each, tuned inputs
- [x] `make usage` — Hold + DownValues usage cards present
- [x] `make check-links` — OK, every \B{} resolves
- [x] `make pdf` — clean (271 pages), no undefined refs/errors; C excerpts + diagrams render
- [x] Facts spot-checked against src/ (expr.h/eval/attr/match/symtab)

## Wrap
- [x] `ROADMAP.md` Ch 10 scope rewritten + status → Verified
- [x] Changelog note in `docs/spec/changelog/2026-09-21.md`
- [x] Fixed Ch 7 DownValues "front" callout → "order of specificity" (+ xref §10.5)

## Review
- Comprehensive 10-section internals chapter written to the book's verified-example
  standard: every REPL transcript is real binary output; the five C excerpts are
  hand-quoted (new `ccode` style, clearly marked not build-verified).
- Nice emergent teaching point: the batch pipe leaves `%`/`Out[n]`/`$Line` symbolic,
  which the chapter uses as direct evidence that In/Out history is a REPL-loop layer
  on top of a stateless evaluator.
- Corrected a genuine inaccuracy: DownValues are ordered by specificity (not recency);
  fixed both the new chapter's claim and Ch 7's older callout so the book is consistent.
- Contributor-facing (book + docs only); no `src/` change, no `$VersionNumber` bump.
- Not committed (no commit requested). Untracked/modified: chapters/10-internals.tex,
  chapters/07-programming.tex, mathilda.sty, examples/10-internals/, ROADMAP.md,
  docs/spec/changelog/2026-09-21.md, tasks/.
