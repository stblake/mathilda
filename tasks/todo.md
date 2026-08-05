# O(1) Associations + Compile[]/auto-compile support

Full plan: `/Users/user/.claude/plans/i-would-like-to-twinkly-umbrella.md`

## Part A — Core interpreter: O(1) Associations (DONE)

- [x] `src/assoc_index.{c,h}`: persistent key→position hash index (`AssocIndex`)
- [x] `src/expr.h`: add `index` pointer to the `EXPR_FUNCTION` union arm (free slack; `sizeof(Expr)` unchanged at 48)
- [x] `src/expr.c`: init NULL in `expr_new_function` + `expr_unshare`; free in `expr_free`
- [x] `src/assoc.c`: `assoc_lookup_value` with lazy index build; route `Lookup`/`KeyExistsQ`/`KeyMemberQ`/`KeyFreeQ`
- [x] `src/part.c`: route `assoc_part_single`; fix `delete_path` in-place aliasing bug
- [x] `src/eval.c`: route `<|…|>[key]` and `Key[k][assoc]` accessors; **in-loop timestamp fixed-point exit** (node stability so the index persists)
- [x] `tests/bench_assoc.c`: C-level O(1) primitive gate + interpreter `Map[Lookup]` gate; `assoc_index.c` added to `COMMON_SRC`
- [x] docs: `data-structures.md` + changelog `2026-08-03.md`; lessons captured

### Review — Part A

- **Result:** single-key lookup is amortised **O(1)** (was O(n) scan). Primitive
  gate: hit ratio 1.04, miss 0.88 across a 2× size doubling. Interpreter
  `Total[Map[Lookup[a,#]&, keys]]`: **>100×** faster on 10⁴ entries, O(1) per probe.
- **Cost:** `sizeof(Expr)` unchanged; numeric hot loops perf-neutral (logistic map
  ~16 ns/iter, unchanged); no memory leaks (valgrind identical to baseline; 2000-cycle
  stress flat).
- **Tests:** association + eval-timestamp/eager-exit/sharing + list/patterns/match/iter/
  replace/sort/funcprog/packed/unevaluated suites all pass; `make check-c99` clean.
- ~~**Known limitation:** `Do`/`Table`/`Fold` re-evaluate a loop-invariant value
  O(n) per step (eval-clock churn).~~ **FIXED (Part A′, 2026-08-05)** — see below.

## Part A′ — Loop-invariant eval-clock churn fix (DONE)

The Part A caveat is resolved: `Do`/`Table`/`Fold` over a loop-invariant literal
value are now **O(1) per iteration**, general (lists/matrices/nested data, not just
associations).

- [x] `src/eval.c`: `g_last_rule_change_clock` finer epoch + `eval_rule_epoch_*`;
      GROUND flag in `last_evaluated_at` top bit (zero `sizeof(Expr)` cost);
      `ground_head` (six pure constructors) + `node_compute_ground`; short-circuits
      + stamp site; `eval_node_stamp`/`eval_node_is_ground` accessors.
- [x] `src/symtab.c` (down-value/clear/remove mark epoch; own-value soft unless
      Protected), `src/attr.c` + `src/core.c` (attribute/Protect sites → rule epoch).
- [x] `tests/bench_assoc.c` Do-loop O(1) gate (ratio 1.04); `tests/test_eval_timestamps.c`
      GROUND correctness (flag present/absent, survives soft bump, invalidated by
      rule change, delayed symbolic assoc never stale).
- [x] Verified: 20+ evaluator/pattern/rule suites pass; `make check-c99` clean;
      valgrind byte-identical to baseline (zero new allocations).
- **Soundness:** GROUND is restricted to six inert structural constructors whose
  canonical form is a pure function of args; straddle-safe consumption of child
  bits; a wrong answer would need redefining a constructor head without advancing
  the rule epoch, which no mutation path does. Details in changelog 2026-08-03.md.

## Part B — Compile[] & auto-compilation for Associations (NOT STARTED)

Built on Part A. Staged milestones; the compiler "cliff" means each is independently
shippable with no regression (an unsupported construct just bails to the interpreter).

- [x] **B1** Read-only parameter bag — DONE (2026-08-05). `CT_ASSOC` type + tightened
      `CT_IS_ARRAY`; `{p, _Association[, _valtype]}` argspec; borrowed handle in a
      scalar-bank register (index pre-built at marshalling); `ASSOC_LOOKUP`/`HASKEY`
      (K_KERN1, pure → CSE/LICM-hoisted), `ASSOC_LEN`, `ASSOC_VALUES`; program-owned
      `AssocSpec` pool (interned for CSE); const-fold of literal/global associations
      (auto-compile story); runtime-varying key bails (= B2). `tests/test_compile_assoc.c`
      (13 tests), valgrind = baseline, check-c99 + compile-coverage clean.
- [x] **B2** Runtime-varying (integer/real) keys — DONE (2026-08-05). `ASSOC_LOOKUP_DYN`
      (K_KERN2, pure); key emitted as a machine scalar in `R[b]`, `flags = result_ct|key_ct<<4`;
      malloc-free `assoc_lookup_value_i64`/`_real` probe Part A's index with a stack `Expr`
      (2 M probes / 20 k bag ≈ 130 ms). Array/string runtime keys still bail.
- [x] **B3** Associations as first-class VM values — DONE (2026-08-05). `KeyDrop`/`KeyTake`
      (`ASSOC_KEYSEL`) and `Counts` (`ASSOC_COUNTS`, array→assoc) produce OWNED associations
      in the array bank; `cf_unbox`/result-extraction return them; native cores
      (`assoc_key_select`, `assoc_counts_ndarray`, no evaluator). **Composition** works to
      any depth (`Lookup[KeyDrop[p,k],j]`, `KeyTake[KeyDrop[…]]`, `Total[Values[Counts[v]]]`)
      via `materialize_assoc_src` + free-source discipline (in-instruction for array-producers,
      free_if_tmp for scalar readers). **Deferred:** `KeyUnion` (→list), `PositionIndex`
      (list-valued) — don't fit the machine-scalar value model.
- [ ] **B4** Higher-order transforms (`Merge`/`GroupBy`/…) via a compiled callee.
- [ ] **B5** Functional key-update; symbol-rebinding `AssociateTo` bails.
