# Task: Implement NonCommutativeMultiply (`**`)

Plan: `~/.claude/plans/pasted-content-id-d470-let-s-implement-structured-cloud.md`

Associative (Flat) but non-commutative multiplication, operator `**`, attributes
`{Flat, OneIdentity, Protected}`, no automatic simplification beyond flattening.
No new C builtin — attribute-driven flattening in the evaluator does the work.

## Edits
- [x] `src/version.h` — bump 0.297 → 0.298
- [x] `src/parse.c` — `OP_NONCOMMULT` enum; `**` lexer arm (prec 5900); added NCM to parse-time n-ary flatten branch
- [x] `src/sym_names.{h,c}` — `SYM_NonCommutativeMultiply` decl/def/init
- [x] `src/print.c` — `get_expr_prec` 5900 (≥2 args); infix `" ** "` in the operator block
- [x] `src/attr.c` — `{"NonCommutativeMultiply", ATTR_FLAT|ATTR_ONEIDENTITY|ATTR_PROTECTED}`
- [x] `src/info.c` — docstring (after Dot)

## Tests
- [x] `tests/test_noncommutativemultiply.c` — 10 groups: parse/flatten/non-comm/one-arg/no-simp/print/precedence/attrs/pattern/memory — ALL PASS
- [x] `tests/CMakeLists.txt` — registered `noncommutativemultiply_tests`

## Docs / versioning
- [x] `docs/spec/builtins/arithmetic.md` — NonCommutativeMultiply section
- [x] `docs/spec/operators.md` — precedence row (5900)
- [x] `docs/spec/changelog/2026-10-05.md` — v0.298 entry (prepended)

## Verification
- [x] `make -j` build (GCC-16, clean); REPL smoke-test — every transcript line matches Mathematica exactly
- [x] built + ran `noncommutativemultiply_tests` (pass); regression set parse/print/evaluate/unevaluated/flatten_at all pass
- [x] valgrind: pure-NCM paths clean at core_init baseline (13,440 B/420 blocks, identical to parse_tests/evaluate_tests); +2 blocks in full test isolated to pre-existing Simplify-family leak
- [x] `make check-c99` PASS, `check-messages` PASS, `check-packed-aware` PASS
- [x] `check-array-exactness` PASS (0 MIXED); `check-nd-surfaces` PASS (NCM not flagged)
- [x] `check-fastpath-sweep` — exit 1 is PRE-EXISTING (stale OFF_BUFFER ~v0.156; 61 NEW heads all unrelated post-v0.156 builtins). NCM NOT in the list → no exempt needed. Gate is not in per-push CI.
- [x] rebuilt code-review-graph (incremental, ok)

## Review (v0.298, complete)

NonCommutativeMultiply (`**`) implemented end to end. Associative (Flat) but
non-commutative generalized multiplication, attributes `{Flat, OneIdentity,
Protected}`, no automatic simplification beyond flattening.

**Design:** no new C builtin — flattening is entirely attribute-driven in the
evaluator (`eval.c` Flat path, gated only on head-is-symbol + ATTR_FLAT). Nine
existing files touched additively + one new test file. Precedence 5900 placed in
the empty gap between Dot (5300) and Power (6500).

**Behavior verified** against the user's Mathematica transcript (REPL smoke +
unit tests): parse-time n-ary flatten (even under Hold), associativity,
non-commutativity (`a**b` ≠ `b**a`, Equal stays unevaluated), one-arg stays,
`{0**a,1**a}` stay, Expand/Simplify/FullSimplify no-ops, infix printing with
precedence-aware parens, and the full `Plus<Times<Dot<NCM<Power` chain. Flat +
OneIdentity pattern matching works (the `DOperator[L1__**L2_,…]` shape).

**Tests:** `tests/test_noncommutativemultiply.c` (10 groups) PASS; regression set
parse/print/evaluate/unevaluated/flatten_at PASS.

**Leaks:** pure-NCM paths valgrind clean at the core_init baseline (13,440 B/420
blocks — proven by isolation; identical to parse_tests/evaluate_tests). The +2
blocks in the full test are the pre-existing Simplify-family leak (the test
calls Simplify/FullSimplify only to assert NCM is a no-op under them).

**Gates:** check-c99 / check-messages / check-packed-aware / check-array-exactness
/ check-nd-surfaces all PASS. check-fastpath-sweep red is pre-existing backlog
(see above), NCM not implicated.

**Docs/version:** arithmetic.md section, operators.md row (5900), changelog
v0.298 prepended, version bumped 0.297 → 0.298. Graph rebuilt.

**Not committed (awaiting user):** per release-tagging rule this is a substantive
change → commit with `; v0.298` and tag `v0.298`. Pre-existing fastpath-sweep
backlog (61 heads, stale OFF_BUFFER) is a separate maintenance task, out of scope.
