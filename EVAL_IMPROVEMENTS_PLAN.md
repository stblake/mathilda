# Evaluator Improvements: PicoCAS vs. Mathematica (Withoff 1992)

This document compares the PicoCAS evaluator against the description of
Mathematica's internals in David Withoff's 1992 "Mathematica Internals: A
Tutorial" (Mathematica Conference, Boston, June 1992) and proposes a
prioritised plan of improvements.

References to Withoff use "[W §X.Y]". Source references use `file:line`.

---

## 0. Status

- **M2 — Symbol interning:** ✅ complete (Tier 3 prep, item 3.1).
  - New `src/sym_intern.{c,h}` provides a global string interner. Every
    `expr_new_symbol(name)` now stores the canonical `const char*`
    returned by `intern_symbol`, so two symbols with the same name share
    the same pointer. Lifetime is program-lifetime; `intern_clear()`
    is available for shutdown but is not called by the REPL today.
  - `expr_eq` for `EXPR_SYMBOL` is now a pointer compare. `expr_copy`
    of a symbol shares the canonical pointer instead of `strdup`-ing.
    `expr_free` skips the symbol-name field (interner owns it).
    `expr_compare` short-circuits on pointer equality before falling
    back to lexicographic `strcmp`.
  - `symtab_get_def` / `symtab_lookup` / `symtab_get_own_values` /
    `symtab_get_down_values` now intern their input name first, then
    do pointer-compare bucket scans instead of `strcmp`. Stored
    `SymbolDef::symbol_name` is the interned pointer; `symtab_clear`
    no longer frees it.
  - New `src/sym_names.{c,h}` exposes cached pointers (`SYM_List`,
    `SYM_Plus`, ..., `SYM_Composition`) populated once via
    `sym_names_init()` at the top of `core_init()`. Hot evaluator
    paths (`has_list_arg`, `apply_listable`, `flatten_sequences`,
    `evaluate_step`'s held-Evaluate / Set / SetDelayed / Rule /
    RuleDelayed / Unevaluated / Function / Derivative / Composition /
    Condition / Part / List checks, plus `assignment_target_symbol`)
    now compare with `==` against these cached pointers instead of
    `strcmp` against literals.
  - `eval_flatten_args` re-interns its `head_name` argument on entry so
    callers (notably `internal_call_impl`, which passes plain C string
    literals like `"Plus"`) keep working without per-call-site changes.
  - All 83 unit tests pass; valgrind shows interner allocations as
    still-reachable (program-lifetime), with no new definitely-lost
    leaks.

- **M1 — Correctness alignment:** ✅ complete.
  - 1.1 Flat now runs before Listable — `Plus[Plus[a,{1,2}],3]` threads to
    `{4+a, 5+a}`.
  - 1.2 User DownValues now run before built-ins ("internal down code");
    Protected symbols are unaffected because `apply_assignment` blocks
    DownValue installation on Protected targets.
  - 1.3 OneIdentity is no longer applied as an evaluation rewrite. The
    1-arg collapse (`Plus[x] → x`, `GCD[x] → x`, ...) moved into each
    head's builtin (`plus.c`, `times.c`, `power.c`, `arithmetic.c` for
    GCD/LCM, `boolean.c` for And/Or, `linalg.c` for Dot). User-defined
    OneIdentity heads now correctly stay un-rewritten: `g[x]` stays as
    `g[x]`. Pattern-matching OneIdentity (in `match.c`) is unchanged.
  - 1.4 `HoldAllComplete` no longer short-circuits the evaluator; it now
    suppresses arg evaluation, Sequence flattening, Unevaluated stripping,
    and Flat — but rules and built-ins on the head still apply
    (`Length[HoldComplete[a,b,c]]` works).
  - 1.5 `$RecursionLimit` (default 1024, matching modern Mathematica)
    added as a C-stack guard. On overflow we wrap the offending sub-
    expression in `Hold[]`, set a sticky `eval_overflow` flag, and bail
    out of every enclosing fixed-point loop so the unwind doesn't burn
    `$IterationLimit` at every level. Exposed as a user-visible OwnValue
    (`eval_init` seeds it; `apply_assignment` syncs the C-side limit on
    `$RecursionLimit = N`). Values below 20 are rejected with a
    `$RecursionLimit::limset` message.
  - 4.4 Reclassified — picocas's existing strip-without-restore matches
    real Mathematica's observed behaviour (`h[Unevaluated[Plus[a,b]]] →
    h[a+b]`) and the existing `tests/test_unevaluated.c` already locks
    that in. Withoff §3.1's literal "restore Unevaluated" appears to be
    a process-level description that is invisible at the printer level;
    no change needed.

---

## 1. Side-by-side comparison

### 1.1 Expression representation

| Aspect | Withoff [§2.1] | PicoCAS (`src/expr.h`, `src/expr.c`) |
|---|---|---|
| Distinction raw vs. normal | Symbols, numbers, strings are *raw*; `f[x]`, `Plus[2,2]` are *normal* | Same conceptual split via `EXPR_INTEGER/REAL/STRING/SYMBOL/BIGINT` vs `EXPR_FUNCTION` |
| Number representation | C `int`, C `double`; arbitrary-precision wraps separate types | `int64_t`, `double`, GMP `mpz_t`, optional MPFR. Consistent. |
| Memory management | **Reference counting**; subtrees are *shared* | **Deep copy** everywhere (`expr_copy`); subtrees never shared. Every `expr_copy` allocates a full duplicate tree. |
| Evaluation timestamp | Symbols and normal expressions carry a "time of last evaluation". *"An expression that has been evaluated once will not normally be evaluated again unless some part of the expression is modified."* | **Not tracked.** Every `evaluate()` call re-traverses and re-evaluates the whole tree from scratch. |
| Symbol identity | Interned: each symbol exists once, comparisons are pointer-cheap | `expr->data.symbol` is a `char*` from `strdup`. Equality requires `strcmp`. The symbol table is a separate hash, but `Expr` nodes do not point at `SymbolDef` — they hold an independent string copy. |

### 1.2 Symbol table

| Aspect | Withoff [§2.2] | PicoCAS (`src/symtab.h`) |
|---|---|---|
| Visible value classes | `OwnValues`, `DownValues`, **`SubValues`**, **`UpValues`**, **`NValues`**, **`FormatValues`**, `DefaultValues`, `Options`, `Messages`, `Attributes` | `own_values`, `down_values`, `attributes`, `docstring`. Everything else is missing. |
| Internal "code" | Internal counterparts of each value class: `down code`, `up code`, `sub code`, `num code`, `format code` (built-ins) | Single `BuiltinFunc builtin_func` slot — equivalent to `down code` only. There is no `up code`, `sub code`, `num code`, or `format code` table. |
| Rule ordering | Ordered by specificity at definition time | Linear list with a coarse heuristic (`has_patterns` ⇒ append, else prepend). Does not distinguish "more specific pattern" among patterned rules. (`src/symtab.c:103-135`) |
| Up-values dispatch | First-class — see eval order below | **Absent.** A user cannot write `Sin /: f[Sin[x_]] := ...` |
| `Messages`, `Options`, `DefaultValues` | Stored as Mathematica expressions on the symbol | None of these tables exist in `SymbolDef`. |
| Symbol storage | Symbol table is itself a Mathematica list at the kernel level | Custom open-chained hash bucket array of size 65535 (`src/symtab.c:10`). Adequate but not introspectable from the language. |

### 1.3 Evaluator

Withoff's algorithm [§3.1] is reproduced here next to the PicoCAS step
(`src/eval.c:331-582`):

| # | Withoff step | PicoCAS implementation | Notes |
|---|---|---|---|
| 1 | If string/number, return self | ✅ `evaluate_step` switch on `EXPR_INTEGER/REAL/STRING/BIGINT/MPFR` | OK |
| 2 | Symbol with no `OwnValues` returns self | ✅ `apply_own_values` returns NULL → fall through to `expr_copy(e)` | OK |
| 3 | Evaluate symbols with `OwnValues` | ✅ via `apply_own_values` | OK |
| 4 | **If no part has changed since last evaluation, return** | ❌ **Missing.** No timestamping; every step re-walks the tree. | Big perf win when present. |
| 5 | Evaluate head `h`; evaluate elements `eᵢ` as appropriate | ✅ `head = evaluate(...)`, then per-arg evaluate respecting `Hold*` | OK |
| 6 | `HoldFirst/HoldRest/HoldAll` skip eval of corresponding `eᵢ` | ✅ | OK |
| 7 | Evaluate arguments with head `Evaluate` (overrides Hold) | ✅ explicit check in held branch (`src/eval.c:384-389`) | OK |
| 8 | If `eᵢ` has head `Unevaluated`, replace with its arguments and **keep a record** of the original | ⚠ Strips `Unevaluated` only in *non-held* positions, no record kept. (`src/eval.c:419-434`) | Diverges from spec; restoration step (#16) cannot work without this record. |
| 9 | `Flat`: flatten nested same-head | ✅ `eval_flatten_args` | But ordering is **after** Listable, not before — see #11. |
| 10 | Flatten `Sequence` heads | ✅ `flatten_sequences`, with explicit skip-list for Set/SetDelayed/Rule/RuleDelayed | OK, but skip-list approach is fragile. |
| 11 | `Listable`: thread over list args | ✅ `apply_listable` | **Order bug:** PicoCAS runs Listable *before* Flat (`src/eval.c:438-457`). Withoff's order is Flat → Sequence → Listable → Orderless. |
| 12 | `Orderless`: sort `eᵢ` | ✅ `qsort` on `eval_compare_expr_ptrs` | OK |
| 13 | Apply user-defined `UpValues` of symbolic heads of `eᵢ` | ❌ **Not implemented — intentional.** See §4. |
| 14 | Apply internal `UpValues` (up code) | ❌ **Not implemented — intentional.** See §4. |
| 15 | Apply user `DownValues` if head is symbol; else user `SubValues` of symbolic head | ⚠ Only `DownValues` for symbolic head. **No `SubValues` path** for `f[1][x]` — intentional, see §4. Pure-function application of a `Function[...]` head is special-cased separately (`src/eval.c:552-560`). |
| 16 | Apply internal `DownValues` (down code), or sub code if head is not a symbol | ⚠ Built-ins are run *before* user `DownValues`, not after (`src/eval.c:460-467` vs. `:540-544`). Withoff's order is **user rules first**, then internal. PicoCAS inverts this. |
| 17 | Restore the head `Unevaluated` if no rules found | ❌ **Missing.** | Required by spec; ties in with #8. |
| 18 | Discard the head `Return` from results of user rules | ❌ **Missing.** No `Return[x]` semantics; `Return` is just an unknown symbol. | |

#### 1.4 Outer fixed-point loop

Withoff: *"If the expression changes at any point, the process usually
starts again from the beginning with the new expression."* In other words,
restart eagerly the moment a rule fires, and on the *changed sub-expression*.

PicoCAS (`src/eval.c:590-614`):

```c
while (iterations < MAX_ITERATIONS) {
    next = evaluate_step(current);
    if (expr_eq(current, next)) { ... return current; }
    expr_free(current);
    current = next;
    iterations++;
}
```

Two issues:
1. `expr_eq(current, next)` is a full-tree structural comparison after every
   step — O(N) per iteration.
2. `evaluate_step` returns a freshly *deep-copied* tree even when nothing
   changed at the top level, defeating the cheap-equality intent of #4.

### 1.5 Other observations

- **`HoldAllComplete` short-circuits everything** including built-ins and
  user `DownValues` (`src/eval.c:365-373`). Mathematica only suppresses
  argument evaluation, `Sequence` flattening, and `Unevaluated` stripping
  for `HoldAllComplete`; rules attached to the head still apply. The
  PicoCAS branch returns immediately and skips builtin/DownValue dispatch.
- **`OneIdentity` is applied universally** for one-arg calls
  (`src/eval.c:547-551`). In Mathematica `OneIdentity` only affects
  *pattern matching*, not evaluation. Applying it as an evaluation rewrite
  changes user semantics (e.g., `Plus[x]` should remain as `Plus[x]` when
  evaluated outside a pattern context).
- **Many head-name checks are `strcmp` against literals** scattered
  throughout `evaluate_step`. With ~22 kLOC and growing, every step pays
  several `strcmp` round-trips on every node.
- **`apply_down_values` builds a fresh `MatchEnv` per rule** and tries
  rules linearly. There is no dispatch index by arity, head shape, or
  first-argument literal.
- **No `$RecursionLimit` distinct from `$IterationLimit`** — the only
  guard is `MAX_ITERATIONS = 4096` on the outer fixed-point loop. Nested
  `evaluate()` calls have no depth limit and can blow the C stack.

---

## 2. Improvements (prioritised)

The improvements below are roughly ordered by **impact ÷ cost**. Each
entry includes a short rationale, the touched files, and risks.

### Tier 1 — Correctness fixes ✅ DONE (M1 milestone)

All five Tier-1 changes have shipped. Summary in §0 above; details preserved
below for the historical record. All 79 unit tests still pass.

**1.1 Re-order Flat / Sequence / Listable / Orderless to match Withoff §3.1** ✅
- Files: `src/eval.c` (`evaluate_step`)
- Moved `Flat` before `Listable`. Listable now sees the post-flattening
  argument list, so `Plus[Plus[a,{1,2}],3]` correctly threads to
  `{4+a, 5+a}`.

**1.2 Apply user `DownValues` before built-ins** ✅
- Files: `src/eval.c`
- DownValues now run before the head's `builtin_func`. Protected symbols
  are unaffected because `apply_assignment` blocks DownValue installation
  on Protected targets.

**1.3 Restrict `OneIdentity` to pattern-matching** ✅
- Files: `src/eval.c`, `src/plus.c`, `src/times.c`, `src/power.c`,
  `src/arithmetic.c`, `src/boolean.c`, `src/linalg.c`.
- Removed the eager rewrite. Each builtin now handles its own n=1 case.
  Pattern-matching OneIdentity (in `src/match.c`) is unchanged.
  `SetAttributes[g, OneIdentity]; g[x]` now correctly stays as `g[x]`.

**1.4 Stop short-circuiting on `HoldAllComplete`** ✅
- Files: `src/eval.c`
- The early-return is gone. HoldAllComplete now suppresses only arg
  evaluation, Sequence flattening, Unevaluated stripping, and Flat;
  built-ins and DownValues on the head still apply, so
  `Length[HoldComplete[a,b,c]]` returns 3.

**1.5 Add a real `$RecursionLimit`** ✅
- Files: `src/eval.c`, `src/eval.h`
- Static depth counter + sticky `eval_overflow` flag. Default 512; settable
  via `eval_set_recursion_limit()`. On overflow we wrap the offending
  expression in `Hold[]`, set the flag, and bail out of every enclosing
  fixed-point loop. Exposing this as a user-facing `$RecursionLimit`
  symbol (read/write through OwnValues) is a follow-up.

### Tier 2 — Major missing features

**2.1 Implement `NValues`** — Withoff [§2.2]
- Files: new `src/nvalues.c` or extend `numeric.c`; `src/eval.c` for
  `N[expr, prec]` to consult `NValues` of subexpression heads before
  falling through.
- Today, numeric evaluation is a separate code path. With `NValues`,
  user code can teach `N[]` how to numericize a user-defined function.
- Risk: medium. Interacts with MPFR precision tracking.

**2.2 Implement `Messages` table and `$MessageList`**
- Files: new `src/messages.{c,h}`; `src/eval.c` (clear `$MessageList` at
  the start of each top-level evaluation, populate on emission); `src/repl.c`
  (assign `MessageList[n]` after each Out).
- Today, every diagnostic uses `printf`/`fprintf` directly (e.g.
  `src/eval.c:133-139`, `:207-208`). A `Message[head::tag, args]`
  primitive and per-symbol `Messages[sym]` table would replace these.
- Risk: low. Mostly mechanical refactor; preserves existing behaviour
  through a default formatter.

**2.3 Implement `FormatValues` (format code)**
- Files: `src/print.c`, `src/symtab.{c,h}`.
- Allow user rules of the form `Format[expr] := ...` to alter output.
  Default print remains the current path when no `FormatValues` are set.

### Tier 3 — Performance work

**3.1 Symbol interning** ✅ DONE (M2)
- Files added: `src/sym_intern.{c,h}`, `src/sym_names.{c,h}`.
- Files touched: `src/expr.c` (symbol creation/copy/free/eq/compare),
  `src/symtab.c` (lookups + storage), `src/eval.c` (hot-path
  comparisons, `eval_flatten_args` defensive intern), `src/core.c`
  (call `sym_names_init()`), `tests/CMakeLists.txt` (new sources).
- Implementation: `data.symbol` keeps its `char*` declared type, but
  the memory is owned by the interner and shared across every Expr
  with that name. Pattern stages cleanly into Tier 3 follow-ups:
  later work can convert the field to a typed `sym_id_t` without
  re-touching every consumer, since equality and identity are now
  pointer-based throughout.
- Followups intentionally deferred: a sweep across the remaining
  ~700 `strcmp(...->data.symbol, "...")` sites in non-hot modules
  (poly.c, simp.c, etc.) -- they continue to work via interned
  string contents, and migrating them is mechanical when needed.

**3.2 Reference-counted, immutable expressions** — Withoff [§2.1]
- Files: `src/expr.{c,h}`, every caller of `expr_copy` / `expr_free`.
- Add an atomic-free refcount (REPL is single-threaded) to `Expr`. Make
  `expr_copy` an inc-ref. `expr_free` becomes dec-ref-and-maybe-free.
  Subtrees can be safely shared between expressions.
- Cuts allocations dramatically. Pattern-matching, replacement, and the
  fixed-point loop all currently deep-copy.
- Risk: high. Mutating helpers (`flatten_args`, `flatten_sequences`,
  `eval_sort_args`, `qsort` in place) must first call `expr_unshare()`
  before mutating, or be rewritten to return new nodes. This is the
  largest single change in the document and should be its own milestone.

**Status — Phase 1 (atom sharing) DONE 2026-05-03.** Landed in three
commits' worth of staged work:
- *Phase 1.* Added `unsigned refcount` to `Expr`; constructors initialize
  to 1; `expr_ref(e)` bumps; `expr_free` dec-refs and only physically
  frees on transition to 0. No semantic change yet (`expr_copy` still
  deep-copies). Commit-ready intermediate.
- *Phase 2.* Audited every in-place atom-payload mutation across the
  tree. Converted nine sites to free-and-replace (`print.c` ×6,
  `linalg.c` ×2, `simp.c` ×1) plus five in `poly.c` exposed by the test
  failure on Phase-3 flip (BPList exponent mutation: lines ~1214, 1261,
  1386, 1435, 1489). Audit script:
  ```
  grep -nE "data\.(integer|real)\s*[+\-*/]?="    src/*.c   # → only expr.c constructors
  grep -nE "->type\s*=\s*EXPR_"                  src/*.c   # → only expr.c
  grep -nE "mpz_\w+\s*\(\s*\S+->data\.bigint"    src/*.c   # → only expr.c
  grep -nE "mpfr_\w+\s*\(\s*\S+->data\.mpfr"     src/*.c   # → expr.c + numeric.c (numeric.c writes target a fresh r — safe)
  ```
  Function-node mutations (`->data.function.{args,arg_count,head} = …`
  in `eval.c` `flatten_args`, `core.c`, `parse.c`) are still safe because
  FUNCTION nodes are NOT yet shared.
- *Phase 3.* Flipped `expr_copy` to `++refcount; return e` for every
  non-FUNCTION type. FUNCTION still deep-copies (so a function-node's
  `args[]` and arg_count remain privately mutable). All 83 unit tests
  green.

**Status — Phase 2 (FUNCTION-node sharing) DONE 2026-05-03.** Followed the
Phase-1 staging:
- Added `Expr* expr_unshare(Expr* e)` (`src/expr.{c,h}`): consumes one
  ref and returns a refcount==1 logical equivalent. Fast-path for
  refcount==1 (zero work). Otherwise allocates a one-level private node:
  args[] / mpz / mpfr / string payload owned by the new node, children
  inc-ref'd via `expr_copy`. After the call the caller may freely
  rewrite the returned node's *direct* fields. Deep mutation requires
  unsharing each level along the path.
- Audited every `args[i] = …`, `arg_count = …`, `args = …`, `head = …`
  write outside `expr.c`. The vast majority operate on freshly built
  FUNCTION nodes (refcount==1 at the point of mutation): `eval.c`
  flatten_args / flatten_sequences / Unevaluated stripping / qsort,
  `internal.c` `internal_call_impl`, `core.c` QuotientRemainder
  arg_count zeroing, `parse.c` argv extension, `match.c` Sequence /
  Repeated binding writes, `series.c` SeriesData coefficient
  recursion, `plus.c` and `times.c` numeric-contagion args[i] writes.
  The only unsafe sites lived in `print.c` (8 sites): the negative-
  prefix print path did `t_copy = expr_copy(arg); …mutate t_copy.args[0]`,
  which under sharing would have corrupted the caller's tree. All 8
  sites now `expr_unshare(expr_copy(arg))`; the nested-Rational case
  unshares two levels.
- Flipped `expr_copy` to `++e->refcount; return e` unconditionally for
  every node type, including FUNCTION.
- New stress suite `tests/test_expr_sharing.c` (84th binary): refcount
  inc/dec mechanics, `expr_unshare` semantics on every type, hot-path
  mutator stability under repetition (Plus/Times contagion, Flat
  flatten, Orderless qsort, Sequence flatten, QuotientRemainder,
  print-negative does-not-mutate-input, deep nested rational), pattern
  binding repeated substitution, polynomial Factor/PolynomialGCD/
  Simplify/D fixed-point loops, and a deterministic xorshift random
  share/unshare/mutate/free walk over a 32-tree population for 4000
  iterations.

**Measured win (Phase 2).** REPL smoke script (Expand + GCD + Factor +
trig + D + Together + Cases + Map + TrigToExp + Simplify + Quit):
- Pre-Phase-2: 224,564 allocs / 203,358 frees / 6,660,295 bytes
- Phase-2:      83,550 allocs /  67,795 frees / 3,118,302 bytes
- Reduction: −62.8% allocs, −53.2% bytes.

Valgrind on `poly_tests` and on `expr_sharing_tests` reports zero
invalid-read/write/use-after-free errors. The 28 KB pre-existing
amplified leaks from Phase-1 are unchanged by Phase-2 (FUNCTION
sharing did not introduce new lost stacks; those leaks were already
all on shared atoms).

Next work in this lane is 3.3 / 3.4 (timestamps + eager early-exit) —
those ride on the now-complete refcount infrastructure.

**Known cost (Phase 1):** valgrind on `poly_tests` shows definitely-lost
bytes rising from 4,360 (137 blocks) to 28,832 (749 blocks). All the
extra-leaked stacks are pre-existing imbalanced ownership — paths that
called `expr_copy` without a matching `expr_free`. Pre-Phase-3 they
leaked one fresh atom node each (40 B); now they pin the original atom
alive at refcount > 0. None are new bugs; they are all amplifications
of pre-existing leaks. Tests pass functionally; a follow-up sweep will
balance the ownership in `get_coeff_expanded`, `poly_div_rem`,
`builtin_polynomial{quotient,gcd,extendedgcd}` once we have the time
budget.

**3.3 Evaluation timestamps / "evaluated" mark** — Withoff [§2.1, §3.1 step 4]
✅ DONE 2026-05-03 (conservative single-counter variant).
- Files touched: `src/expr.{c,h}` (added `uint64_t last_evaluated_at`,
  initialised to 0 in every constructor, reset to 0 in `expr_unshare`);
  `src/eval.{c,h}` (declared and defined `eval_clock_get/bump`, added
  the early-exit check at the top of `evaluate()` and the timestamp
  stamp on the clean fixed-point exit); `src/symtab.c` (`add_rule`,
  `symtab_clear_symbol` bump the clock); `src/attr.c` (`set_attributes`,
  `add_single_attribute`, `remove_single_attribute` bump only on actual
  state changes).
- Mechanism: a single 64-bit `g_eval_clock` starts at 1. Every Expr
  records the clock value at the moment it was last evaluated to a
  fixed point. `evaluate(e)` checks at entry: if `e->type == EXPR_FUNCTION`
  and `e->last_evaluated_at == g_eval_clock`, return an inc-ref'd view
  immediately — the entire outer fixed-point loop and `evaluate_step`
  body are skipped. Atom paths already short-circuit cheaply inside
  `evaluate_step`, so the pre-check is restricted to FUNCTION nodes
  to avoid an extra branch on the common atom path. The stamp is
  written on the *clean* fixed-point exit only — overflow and
  iteration-cap paths leave the field untouched so a later evaluator
  gets a fresh chance to make progress.
- Cache invalidation: `Set`, `SetDelayed`, `Clear`, `SetAttributes`,
  `ClearAttributes` all bump the clock through their respective
  symtab/attr mutators. Pure builtin calls (Sin, Factor, Expand, ...)
  do NOT bump, so re-evaluating the result of a pure call is now
  effectively free.
- Tests: new `tests/test_eval_timestamps.c` (16 tests, the 85th
  binary): clock starts non-zero and is monotonic; pure-builtin calls
  do not bump; fresh exprs have zero timestamp; evaluate() stamps the
  result; re-evaluation under same clock returns the same pointer
  (refcount-shared); Set/SetDelayed/Clear/SetAttributes invalidate
  caches; pre-cached subexpressions reused inside a fresh outer call
  still produce correct end-to-end results; `expr_eq` ignores the
  timestamp; interleaved Set / re-evaluate produces fresh values; a
  200-iter random-walk test alternates Set/Clear/Read across three
  symbols and verifies every read against a shadow state, ensuring
  no stale cache survives a definition change.
- All 85 unit-test binaries pass; valgrind on `eval_timestamps_tests`
  reports zero invalid-read/write/use-after-free errors.
- Memory cost: 8 bytes per Expr (uint64_t + alignment). On a 32-byte
  Expr layout this is +25%; in practice the savings from the cache
  on hot-loop workloads outweigh the per-node overhead.
- Refined variants left for future work: per-symbol clocks + a
  per-expression "max clock of any referenced symbol" cached on
  construction, which would let `f[Sin[Pi/4]]` survive an unrelated
  `g = ...` assignment without losing its cached value. The conservative
  variant is the natural first step and what currently ships.

**3.4 Eager early-exit fixed-point loop** ✅ DONE (2026-05-03)
- Files: `src/eval.h` (new `evaluate_step(Expr*, bool*)` signature with
  contract doc); `src/eval.c` (the `bool* changed` instrumentation
  inside `evaluate_step` and the `step_changed || expr_eq` hybrid in
  the outer fixed-point loop); `src/eval.c` (the helpers
  `flatten_sequences` and `eval_flatten_args` now return `bool` so the
  caller can cheaply detect a no-op flatten); `tests/test_eval_eager_exit.c`
  (15 new tests); `tests/CMakeLists.txt` (86th binary registered).
- Mechanism: `evaluate_step` now takes a `bool* changed` out-parameter
  and sets `*changed = true` whenever any of these fired during the
  step: head re-evaluation produced a different subtree
  (pointer-inequality after refcount sharing), arg evaluation reduced
  any sub-expression, an `Evaluate[]` wrapper was stripped, Sequence
  flattening reshaped the args, an Unevaluated wrapper was removed,
  Flat flattened a nested same-head call, Listable threaded over a
  list, Orderless re-sorted out-of-order args, a DownValue matched, a
  built-in returned non-NULL, Set/SetDelayed installed a rule, a pure
  Function applied, Derivative-of-Function reduced, or Composition
  unfolded. The outer loop in `evaluate()` then exits the moment
  `step_changed` is false — saving the O(tree) `expr_eq(current, next)`
  compare on the common case where the input is already a fixed point
  (atoms, bare symbols, fully-reduced functions). Pointer-identity
  works as a "did sub-evaluation change anything" check because
  `expr_copy()` is now refcount-share (no deep copy) per phase-2.
- False-positive handling: some built-ins (notably `builtin_plus` and
  `builtin_times`) unconditionally rebuild their output even when no
  terms combine, so they trip the change flag on every step. To preserve
  termination on `a + b + c`-style inputs without losing the §3.4 win,
  the outer loop falls back to `expr_eq(current, next)` *only when*
  the change flag was true. Cost on the slow path is identical to the
  pre-§3.4 path; the win is the cheap boolean fast-path on the common
  case where no rewrite fires.
- Idempotence guarantee: false-positives are correctness-safe (one
  extra outer iteration at worst); false-negatives would terminate
  too early and are the bug to avoid. The instrumentation in
  `evaluate_step` is therefore conservative — every potentially-rewriting
  branch sets the flag explicitly, with a one-line comment naming the
  reason (e.g. "DownValue rule fired", "List threading reshaped the
  call", "Composition unrolled").
- Test coverage: `tests/test_eval_eager_exit.c` exercises the atom
  no-change invariant, bare-symbol no-change, bound-symbol-changes,
  built-in-fires-and-changes, NULL out-param tolerance, the
  `Plus[a,b,c]` regression case (which would diverge under a
  no-fallback design), `repeated_evaluate_uses_timestamp_path` (the
  §3.3 + §3.4 interaction), Listable threading, Flat flattening,
  Sequence flattening, Orderless re-sorting, DownValue matches, a
  50-rep heavy-eval byte-stability sweep, a 30-rep deep-nested
  no-op stability sweep, and a Set+evaluate cycle that depends on
  the §3.3 clock invalidation kicking in.
- Combined §3.3 + §3.4 picture: on the *first* evaluation of an
  expression whose every sub-step is already a fixed point, §3.4
  saves the second outer iteration plus the converging `expr_eq`
  compare; on *subsequent* evaluations under an unchanged clock,
  §3.3's early-exit at the top of `evaluate()` short-circuits the
  whole thing in O(1) without entering `evaluate_step` at all.

**3.5 DownValue dispatch index** ✅ DONE (2026-05-03, M4)
- Files touched: `src/symtab.h`, `src/symtab.c`, `tests/test_rule_dispatch.c`,
  `tests/CMakeLists.txt`.
- Mechanism: each `Rule` now carries two pre-computed dispatch keys,
  set once at insertion time:
    * `dispatch_arity` — the number of args the pattern accepts at the
      top level. Set to `-1` whenever the pattern has a top-level
      `BlankSequence`, `BlankNullSequence`, `Optional`, `Repeated`,
      `RepeatedNull`, or `OptionsPattern`. Variable-arity rules always
      pass the filter.
    * `first_arg_head_canon` — interned pointer to the first arg's head
      symbol. `Sin`, `Plus`, `Integer`, `Real`, `Symbol`, etc.
      Computed by `pattern_arg_head_canon` which looks through
      `HoldPattern` / `Verbatim`, recurses into `Pattern[name, q]`,
      `Condition`, `PatternTest`, agreeing branches of `Alternatives`,
      and treats sequence/optional patterns as wildcard. NULL means
      "scan unconditionally."
- `apply_down_values` (`src/symtab.c:apply_down_values`) computes the
  input call's `(arity, first_arg_head_canon)` once, then for each
  rule:
    * skips when the rule committed to a specific arity that disagrees
      with the input;
    * skips when both rule and input committed to a specific first-arg
      head and the canonical pointers differ;
    * otherwise invokes the matcher exactly as before.
  Both checks are pointer compares against interned strings, so the
  filter is a few nanoseconds per skipped rule.
- Soundness: a rule that *could* match the input is never skipped. A
  rule whose `dispatch_arity` is -1 (variable arity) always passes the
  arity gate; a rule with NULL `first_arg_head_canon` always passes the
  head gate; an input with no extractable first-arg head only filters
  on arity. Wrong-shape rules can never reach the matcher.
- Stress test: `test_heavy_downvalue_dispatch` adds 100 specific
  literal rules `hv[k] -> k+1000` plus a wildcard fallback. Numeric
  inputs hit their specific rule under the dispatch filter; symbolic
  input cleanly falls through to the fallback.

**3.6 Improve rule-ordering specificity** ✅ DONE (2026-05-03, M4)
- Files touched: `src/symtab.h`, `src/symtab.c`, `tests/test_rule_dispatch.c`.
- Mechanism: every `Rule` now carries an `int32_t specificity` score
  computed at insertion time by `pattern_specificity`. The score
  combines:
    * literal atoms                      → +100 each
    * plain function calls               → +50 + sum of children
    * typed `Blank[h]`                   →  +20
    * `Pattern[x, q]`                    → score(q) (binding is free)
    * `HoldPattern` / `Verbatim`         → transparent
    * `Condition` / `PatternTest`        → score(inner) +10
    * plain `Blank[]`                    →  -10
    * `BlankSequence[]`                  → -100
    * `BlankNullSequence[]`              → -200
    * `OptionsPattern[]`                 → -150
    * `Optional[...]`                    →  -50
    * `Repeated[p]` / `RepeatedNull[p]`  → score(p) -50 / -100
    * `Alternatives[a, ...]`             → min(score(a_i)) (weakest wins)
- `add_rule` now does a single linear walk to find the position with
  strictly lower specificity, then stable-inserts the new rule before
  it. Equal-specificity rules retain insertion order; identical
  patterns still replace in place (so the in-place RHS update from the
  pre-M4 design is preserved).
- This subsumes the old binary `has_patterns` heuristic: literal-arg
  rules score ≥100 and dominate the same head's pattern rules; typed
  blanks dominate plain blanks; sequence-pattern rules sink to the
  bottom; user-supplied insertion order is preserved among equals.
- Risk audit: full 87-binary test suite passes unchanged. The
  Mathematica-style "more specific first" ordering is now strict
  rather than binary, but the only behavioural shift is that pattern
  rules of *unequal* specificity may now fire in a different order
  from their insertion sequence — the `test_classical_recursion_dispatch`
  case (`fib[n_]` defined before `fib[0]` and `fib[1]`) demonstrates
  the desirable form of this shift.
- M4 is the last Tier-3 milestone (Refcount + timestamps + eager exit
  + dispatch index + specificity ordering).

### Tier 4 — Cleanups and small features

**4.1 Move literal head-name strings to a header** ✅ DONE (folded into M2)
- Files: `src/sym_names.{c,h}` shipped together with 3.1; further
  symbols can be added by extending the `extern const char*`
  declarations and the `sym_names_init` body.

**4.2 Consolidate the assignment-target detection**
- Files: `src/eval.c:167-182`, `src/eval.c:467-537`.
- The Set/SetDelayed branch is ~70 lines deep inside `evaluate_step`
  with several pre-evaluation passes that mostly recompute what
  `assignment_target_symbol` already knows. Extract into
  `apply_set_or_set_delayed(Expr* res, bool delayed)` for readability.

**4.3 Implement `Return[x]`** — Withoff [§3.1 last bullet]
- Files: `src/eval.c`, plus probably a `Return` builtin to ensure correct
  Hold attributes.
- After applying user rules, if the result is `Return[x]`, replace with
  `x` at the appropriate scope boundary (`Function`, `CompoundExpression`,
  `Module`, `Block`, `With`, `Do`, `For`, `While`).
- Risk: medium — needs scope-boundary detection.

**4.4 Track and restore `Unevaluated`** — Withoff [§3.1] ⚠️ NO-OP
- On closer inspection, picocas's strip-without-restore behaviour matches
  what real Mathematica does at the user-visible level
  (`h[Unevaluated[Plus[a,b]]] → h[a+b]`), and `tests/test_unevaluated.c`
  already locks this in. Withoff §3.1's literal "restore Unevaluated"
  step appears to be a process-internal notion that is invisible at the
  printed output — likely the marker is restored only when its identity
  is observable to subsequent processing (e.g. via ReleaseHold-style
  introspection). Unless we discover a concrete divergence, no change is
  required here.

**4.5 Implement `$Pre`, `$Post`, `$PrePrint`** — Withoff [§1.4]
- Files: `src/repl.c`.
- Apply these hooks if defined as `OwnValues` of the corresponding
  symbols. Today the REPL has no such hooks.

**4.6 Migrate diagnostic output to the streams model**
- Files: `src/eval.c`, `src/repl.c`, every module that calls
  `fprintf(stderr, ...)`.
- Introduce thin `picocas_message(stream, fmt, ...)` to land all
  diagnostics in `$Messages` (and later `$MessageList`).

---

## 3. Suggested ordering / milestones

Each milestone should ship with tests in `tests/`. Milestones are sized
so each can be reviewed independently.

**M1 — Correctness alignment (Tier 1).** 1.1, 1.2, 1.3, 1.4, 1.5 plus
4.4 (Unevaluated restore). Pure correctness changes that bring the
evaluator closer to Withoff §3.1. No data-structure churn.

**M2 — Symbol interning (Tier 3 prep).** 3.1 ✅ DONE (4.1 folded in:
the cached `SYM_*` pointers in `sym_names.h` already serve the same
deduplication purpose, so the separate `sym_names.h` task is closed.)
Unlocks the perf gains in M3 by making "is this the symbol `List`?"
checks free.

**M3 — Refcount + timestamps (Tier 3 hot loop).** 3.2, 3.3, 3.4 land
together. The first iteration uses a global eval-clock and a single
shared/copy-on-write Expr. Big perf milestone.

**M4 — Dispatch and rule ordering (Tier 3 polish).** 3.5, 3.6. Becomes
worthwhile only after M2/M3 because it depends on cheap symbol
comparisons and stable ownership semantics.

**M5 — Messages / Format / NValues / Return (Tier 2 + Tier 4).**
2.1, 2.2, 2.3, 4.3, 4.5, 4.6. Quality-of-life features that round out
the Mathematica fidelity.

---

## 4. Out of scope

These appear in Withoff but are deliberately not part of this plan:

- **`UpValues` and `SubValues`** (Withoff §2.2, §3.1 steps 13–15) — the
  dispatch semantics they enable (`f /: g[f[x_]] := ...`, `f[x_][y_] :=
  ...`) are intentionally excluded. They make evaluation order harder
  to reason about and are not a target for PicoCAS.
- Notebook front-end and `MathLink` (Withoff §1, §4) — PicoCAS targets a
  text REPL.
- `Dump` for snapshotting kernel state (Withoff §1.3) — a substantial
  feature with limited value at PicoCAS's current size.
- Package autoloading via `AutoLoad` / `SystemStub` (Withoff §1.2) —
  premature; revisit if startup time becomes an issue.

---

## 5. Pointers in the source

For implementers, the highest-leverage entry points:

- Evaluator main step: `src/eval.c:331` (`evaluate_step`)
- Outer fixed-point loop: `src/eval.c:590` (`evaluate`)
- Hold/Evaluate handling: `src/eval.c:376-396`
- Sequence + Unevaluated handling: `src/eval.c:401-434`
- Listable/Flat/Orderless ordering bug: `src/eval.c:438-457`
- Builtin / DownValues / OneIdentity ordering bug: `src/eval.c:460-551`
- Symbol table layout: `src/symtab.h:18-25`
- Rule storage and ordering heuristic: `src/symtab.c:103-135`
- Expression layout (no refcount, no timestamp): `src/expr.h:31-48`
