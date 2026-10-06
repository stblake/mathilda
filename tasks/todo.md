# UpValues trial — todo

Feature-flagged (`#if UP_VALUES`, default ON via `-DUP_VALUES=1`). Plan: `UPVALUES_PLAN.md`.

## Phase 0 — build flag, symbols, struct field, teardown
- [x] makefile: `UP_VALUES ?= 1` → `-DUP_VALUES=1`
- [x] tests/CMakeLists.txt: `add_compile_definitions(UP_VALUES=1)`
- [x] symtab.h: `Rule* up_values;` field + decls + `symtab_up_value_count`
- [x] symtab.c: teardown (reset_node_payload, clear_symbol, count helper)
- [x] sym_names.{h,c}: SYM_UpValues/UpSet/UpSetDelayed/TagSet/TagSetDelayed/TagUnset/Definition/DownValues/OwnValues
- [x] Gate: `make UP_VALUES=1` and `make UP_VALUES=0` both build

## Phase 1 — symtab layer
- [x] symtab_add_up_value / get / apply_up_values
- [x] symtab_remove_matching_up_value
- [x] symtab_set_{up,down,own}_values

## Phase 2 — parser & printer
- [x] parse.c: OP_UPSET/UPSETDELAYED/TAGSET; lex `^=`,`^:=`,`/:`; TagSet led branch
- [x] print.c: prec + print_standard/print_tex + Definition render

## Phase 3 — evaluator
- [x] upvalue hook before apply_down_values_def
- [x] assignment dispatch for Up*/Tag*; build_assignment_target helper
- [x] handlers: apply_up_assignment / apply_tag_assignment / apply_tag_unset
- [x] {Up,Down,Own}Values[sym]=list intercept
- [x] classify_tag_position looks through Condition; element_dispatch_head handles `_tag`

## Phase 4 — builtins, attrs, docstrings, scoping
- [x] builtin_up_values (symbol + string)
- [x] register heads; attributes; docstrings
- [x] Clear/ClearAll/Remove/Unset up_values (via symtab teardown)
- [x] Block/Module save-restore

## Phase 5 — inspection
- [x] Definition builtin + print special-case (inert; FullForm = Definition[s])
- [x] ?name / Information extension

## Phase 6 — tests & valgrind
- [x] tests/test_upvalues.c (25 cases, all pass)
- [x] valgrind clean (no Mathilda-frame leaks, no invalid access)
- [x] check-messages; check-c99 pass
- [x] no regressions (eval/symtab/parse/regression/match/core_algebra/comparisons/trace green)
- [x] UP_VALUES=0 clean build, behaviour unchanged

## Phase 7 — docs & version
- [x] docs/spec/builtins/{assignment-and-rules,expression-information}.md
- [x] docs/spec/changelog/2026-10-05.md
- [x] version.h bump 0.294 → 0.295

---

## Review

Shipped the UpValues subsystem as a faithful Mathematica recreation, fully gated
behind `#if UP_VALUES` (default ON). All behaviours from the supplied Wolfram
reference verified via `-file` scripts and the 25-case unit suite.

**What works** (all matching the reference): UpSet/UpSetDelayed (single &
multi-symbol install, blank-head tags `a_mod`), TagSet/TagSetDelayed (up/down/own
classification, `tagnf`, Condition-wrapped LHS), TagUnset, modular arithmetic
end-to-end (`mod[0,5]`), upvalue-before-downvalue precedence, all-heads upvalue
with Hold firing / HoldComplete suppressing, UpValues[] reader (symbol + string +
`Names` map), all three `=list` setters (incl. reorder and copy idioms),
Definition (inert; displays own/down/up), ?name/Information display, and
Block/Clear/Remove scoping.

**Efficiency:** the per-call hook early-outs on a global `symtab_up_value_count`,
so a program with no up-values pays one load+branch before DownValue dispatch.
Candidate collection is level-one only, deduped, non-materializing, and reuses
the existing arity/first-arg-head dispatch pre-filter.

**Memory:** valgrind on the suite shows zero Mathilda-frame leaks and no invalid
reads/writes/frees/double-frees (the only "definitely lost" blocks are macOS
dyld/Objective-C startup, the standard valgrind-on-macOS false positives).

**Documented divergences** (pre-existing Mathilda conventions): rule order follows
specificity sort (insertion order the tie-break); immediate-vs-delayed is not
recorded, so Definition/?name render down/up values delayed. SubValues
(`f[x][y]` tags) are out of scope and declined with a message.

**Non-regressions:** the 4 `iter_tests` failures and `modular_tests` exit-1 are
PRE-EXISTING (float formatting `3.0` vs `3.`, Block dynamic-scope capture) —
confirmed identical on the `UP_VALUES=0` binary where this feature compiles out.

**Follow-ons (not done):** a `delayed` bit on `Rule` would let Definition render
`=` vs `:=` faithfully; book coverage deferred.

## Stress-test campaign (adversarial)

Behavioural corpus + independent diff bug-hunt + valgrind churn. Found & fixed 5
correctness bugs (4 in the feature, 1 pre-existing core):
1. Delayed `rhs /; test` not lifted onto the LHS for `^:=` / `TagSetDelayed`
   (guard never filtered) — `move_rhs_condition`.
2. `DownValues[s]=…`/`OwnValues[s]=…`/`UpValues[s]=…` bypassed Protected/Locked.
3. `element_dispatch_head` polluted `Blank`/pattern-wrapper heads; added `UpSet::nosym`
   and top-level `Condition` unwrap in `collect_level_one_symbols`.
4. Orderless upvalues fired order-dependently (first-arg-head dispatch filter is
   unsound for Orderless) and bound in non-WMA order — sort before match + skip the
   filter for Orderless.
5. **Pre-existing core bug:** `pattern_alpha_normalize` corrupted stored patterns on a
   non-matching `Unset`/`TagUnset` of a guarded rule (refcount-aliased in-place
   rename) — fixed with `deep_copy_private`; cures DownValues/OwnValues/UpValues.

Verified: 34-case unit suite green; behavioural + robustness corpora 0 FAILs;
valgrind clean under heavy install/fire/clear/Block/Unset churn (no Mathilda-frame
leaks, no invalid access); subsystems (Integrate/Simplify/Solve/Factor/Sum/Series/D)
byte-identical with and without an upvalue defined; no regressions across
eval/symtab/match/regression/cond/comparisons/etc.; `UP_VALUES=0` and `=1` both
build clean; check-messages + check-c99 pass.
