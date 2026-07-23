# Trace: produce Mathematica's nested evaluation trace  [DONE 2026-07-23]

Plan: `/Users/user/.claude/plans/trace-doesn-t-trace-an-reactive-hammock.md`

- [x] Replace flat depth-gated `TraceCollector` with a frame stack in `src/eval.c`
      (`TraceFrame`, `g_trace_top`/`g_trace_root_result`/`g_tracing`, push/pop/splice)
- [x] Record hooks in `evaluate()` loop into the current frame (dedup-by-last)
- [x] Snapshot the reassembled `f[evaluated_args]` form in `evaluate_step` (line 960)
      so `8 + 16 + 1` shows, not the pre-arg-eval `2^3 + 4^2 + 1`
- [x] Suppression counter (`g_trace_suppress`) around the builtin call and
      `Listable` threading so Range-internal / thread-element evals stay atomic
- [x] `eval_collect_trace` sets up/tears down the stack; resets suppression (Trace
      is itself a builtin, called under the caller's suppression)
- [x] `src/trace.c`: recursive HoldForm wrapping + recursive form-filter (flatten)
- [x] Both reported cases match exactly; DownValue/chain/nested-Trace verified
- [x] `trace_tests` pass (updated `x+1` multistep + `2*3+1,_Integer`→`{6,7}`;
      added nested test); eval/evaluate/core/iter/match/unevaluated/eager-exit/
      timestamps suites pass; no core-eval regression (change is `g_tracing`-gated)
- [x] valgrind: 0 leaks from trace code (11648B/146blk == macOS dyld baseline)
- [x] Docs: `docs/spec/builtins/expression-information.md` + changelog `2026-07-20.md`

Review: The evaluator's semantics are unchanged — the collector was rewritten to
*observe* the evaluator's existing recursion (frames per `evaluate()` call), plus
two observation points (reassembled-form snapshot; suppression of builtin/Listable
internal evals). Untraced hot path adds one predicted-false bool read per
`evaluate()` call.

---

# Implement Sow and Reap  [DONE 2026-07-23]

Plan: `/Users/user/.claude/plans/we-should-implement-sow-twinkly-spring.md`

- [x] Add SYM_Sow / SYM_Reap to sym_names.{h,c}
- [x] Write src/sowreap.h (prototypes)
- [x] Write src/sowreap.c (collector stack, builtin_sow, builtin_reap)
- [x] Register Sow/Reap in core.c (builtins, attributes, docstrings)
- [x] Write tests/test_sow_reap.c (11 test groups)
- [x] Add sowreap.c to COMMON_SRC + test block in tests/CMakeLists.txt
- [x] Build REPL clean (gcc -std=c99 -Wall -Wextra, EXIT 0), spot-checked all spec examples
- [x] Build + run sow_reap_tests -> all pass
- [x] valgrind: 0 leaks in Mathilda source (all definite-lost are macOS dyld/Accelerate baseline)
- [~] Docs: docs/spec tree does NOT exist in this checkout (not git-tracked) -> nothing to update;
      in-system docs are the docstrings (verified via Information[Sow]/Information[Reap]).

## Review

- Sow/Reap implemented in a self-contained `src/sowreap.c` with a file-global
  stack of `ReapFrame` (each frame on builtin_reap's C stack, linked via `prev`;
  balanced LIFO push/pop). Values/tag-groups are singly linked lists of Expr as
  requested -> O(1) ordered append, exact Sow order preserved.
- Full semantics: default/tagged/tag-list Sow, single-pattern and pattern-list
  Reap, f[tag,{vals}] handler form, tag-selective reaping, nested-Reap routing
  (innermost matching frame wins; non-match propagates outward), Throw-across-Reap
  unwind. All spec examples reproduce exactly.
- KNOWN DIVERGENCE (not a Sow bug): Mathilda's `Sum` evaluates its summand
  symbolically once and closed-forms it, so `Sow`/`Print` inside `Sum` fire once
  (WL iterates). `Do`/`Table`/`Module`+`Do` iterate faithfully; the test uses
  `Do` for the "reap across a loop" case. Fixing Sum is out of scope.
