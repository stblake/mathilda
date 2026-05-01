# Evaluator Improvements: PicoCAS vs. Mathematica (Withoff 1992)

This document compares the PicoCAS evaluator against the description of
Mathematica's internals in David Withoff's 1992 "Mathematica Internals: A
Tutorial" (Mathematica Conference, Boston, June 1992) and proposes a
prioritised plan of improvements.

References to Withoff use "[W §X.Y]". Source references use `file:line`.

---

## 0. Status

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
  - 1.5 `$RecursionLimit` (default 512) added as a C-stack guard. On
    overflow we wrap the offending sub-expression in `Hold[]`, set a
    sticky `eval_overflow` flag, and bail out of every enclosing
    fixed-point loop so the unwind doesn't burn `$IterationLimit` at
    every level. Settable via `eval_set_recursion_limit()`. Wiring this
    to a user-visible `$RecursionLimit` symbol is a pending tier-4 item.
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
| 13 | Apply user-defined `UpValues` of symbolic heads of `eᵢ` | ❌ **Missing.** | Major gap. |
| 14 | Apply internal `UpValues` (up code) | ❌ **Missing.** | Major gap. |
| 15 | Apply user `DownValues` if head is symbol; else user `SubValues` of symbolic head | ⚠ Only `DownValues` for symbolic head. **No `SubValues` path** for `f[1][x]`. Pure-function application of a `Function[...]` head is special-cased separately (`src/eval.c:552-560`), but user cannot attach rules to `f` for `f[x][y]` shape. |
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

**2.1 Implement `UpValues`** — Withoff [§2.2, §3.1]
- Files: `src/symtab.{c,h}` (add `up_values` field, `apply_up_values`),
  `src/eval.c` (insert step between Orderless and DownValues), `src/replace.c`
  (parse `f /: g[f[...]] := ...` via `TagSet` / `TagSetDelayed`).
- Adds a new step in `evaluate_step` that, for each symbolic head among
  the `eᵢ`, tries that symbol's `up_values` against the *parent*
  expression. First match wins. UpValues fire **before** DownValues per
  Withoff.
- Test idea: `Sin /: f[Sin[x_]] := g[x]` then `f[Sin[a]]` → `g[a]`.
- Risk: medium. Touches the parser too (`TagSet`, `TagSetDelayed`,
  `/:` operator), and `apply_assignment` to dispatch into `UpValues`
  storage when the LHS uses tag-set syntax.

**2.2 Implement `SubValues`** — Withoff [§2.2, §3.1]
- Files: `src/symtab.{c,h}`, `src/eval.c`.
- For an expression whose head is itself a normal expression (e.g.
  `f[a][b]`), look up `SubValues` on the symbolic head of the head
  (i.e. `f`). Currently PicoCAS handles this only through the
  pure-function special case and `Derivative[...]` (`src/eval.c:552-575`).
- Test idea: `f[x_][y_] := x + y` then `f[1][2]` → `3`. This currently
  installs as a `DownValues` rule on `f` and only matches because
  `apply_down_values` walks the whole expression — verify and adjust.

**2.3 Implement `NValues`** — Withoff [§2.2]
- Files: new `src/nvalues.c` or extend `numeric.c`; `src/eval.c` for
  `N[expr, prec]` to consult `NValues` of subexpression heads before
  falling through.
- Today, numeric evaluation is a separate code path. With `NValues`,
  user code can teach `N[]` how to numericize a user-defined function.
- Risk: medium. Interacts with MPFR precision tracking.

**2.4 Implement `Messages` table and `$MessageList`**
- Files: new `src/messages.{c,h}`; `src/eval.c` (clear `$MessageList` at
  the start of each top-level evaluation, populate on emission); `src/repl.c`
  (assign `MessageList[n]` after each Out).
- Today, every diagnostic uses `printf`/`fprintf` directly (e.g.
  `src/eval.c:133-139`, `:207-208`). A `Message[head::tag, args]`
  primitive and per-symbol `Messages[sym]` table would replace these.
- Risk: low. Mostly mechanical refactor; preserves existing behaviour
  through a default formatter.

**2.5 Implement `FormatValues` (format code)**
- Files: `src/print.c`, `src/symtab.{c,h}`.
- Allow user rules of the form `Format[expr] := ...` to alter output.
  Default print remains the current path when no `FormatValues` are set.

### Tier 3 — Performance work

**3.1 Symbol interning**
- Files: `src/expr.{c,h}`, `src/symtab.{c,h}`, plus a sweep over the
  tree (`grep "data.symbol" src/`).
- Replace `Expr->data.symbol: char*` with `Expr->data.symbol: SymbolDef*`
  (or a stable `const char*` returned by an interner). All "is the head
  named `List`" checks become pointer compares against well-known
  `SymbolDef*` constants cached at init time.
- Wins: every `strcmp(head->data.symbol, "List")` (there are dozens in
  `evaluate_step` alone) becomes a single pointer compare.
- Risk: high — touches every module. Stage by introducing a `sym_id_t`
  type alongside the existing string, migrate hot paths, then drop the
  string field.

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

**3.3 Evaluation timestamps / "evaluated" mark** — Withoff [§2.1, §3.1 step 4]
- Files: `src/expr.{c,h}`, `src/eval.c`, `src/symtab.c`.
- Maintain a global monotonic counter `eval_clock` incremented whenever
  *any* symbol is mutated (own/down/up/sub-value added or cleared,
  attribute changed). Each `Expr` carries `last_evaluated_at` (uint64).
- In `evaluate_step`, before doing any work, if `e->last_evaluated_at`
  is greater than or equal to the most recent change of any symbol the
  expression depends on, return `e` as-is.
- Conservative variant: a single global counter is checked; bumping any
  symbol invalidates everyone. Already a huge win because most evaluation
  cycles do not modify the symbol table at all.
- Refined variant: per-symbol clocks plus a per-expression "max clock of
  any referenced symbol" cached on construction.
- Cuts re-evaluation cost on the second and subsequent iterations of the
  outer fixed-point loop, and on REPL reuse via `In[n]`/`Out[n]`.
- Risk: medium-high once cache invalidation is dependency-aware. Start
  with the global-counter variant.

**3.4 Eager early-exit fixed-point loop**
- Files: `src/eval.c:590-614`.
- Replace `expr_eq(current, next)` with a "did anything fire" boolean
  threaded through `evaluate_step`. When nothing fired, we are done —
  no full structural compare needed.
- Companion: when `evaluate_step` does not produce a change, return the
  *same* pointer (refcount-bumped, or shared) rather than a deep copy.
- Risk: low if combined with 3.2; medium standalone (need to be careful
  not to leak the original).

**3.5 DownValue dispatch index**
- Files: `src/symtab.c`, `src/match.c`.
- Bucket rules on a cheap key — most commonly arity, then the symbolic
  head of the *first* argument when present. `apply_down_values` then
  scans only the relevant bucket.
- Risk: low. Buckets fall back to linear scan for rules whose key is
  ambiguous (sequence patterns).

**3.6 Improve rule-ordering specificity**
- Files: `src/symtab.c:103-135`.
- Replace the binary "has patterns or not" sort with a real pattern
  specificity score (number of literal sub-expressions, presence of
  sequence blanks, presence of `Condition`/`PatternTest`). Match
  Mathematica's `DownValues` ordering documentation.
- Risk: medium. Some test fixtures may rely on insertion order.

### Tier 4 — Cleanups and small features

**4.1 Move literal head-name strings to a header**
- Files: new `src/sym_names.h` with `extern const char* const SYM_LIST;`
  etc., plus a small init that interns them.
- Even before full interning (3.1), this de-duplicates literals.

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

**M2 — UpValues / SubValues (Tier 2 core).** 2.1 and 2.2. Together they
make PicoCAS able to host the dispatch patterns Mathematica users
expect for operator overloading and curried-call rules.

**M3 — Symbol interning (Tier 3 prep).** 3.1 and 4.1. Mechanical, but
unlocks the perf gains in M4 by making "is this the symbol `List`?"
checks free.

**M4 — Refcount + timestamps (Tier 3 hot loop).** 3.2, 3.3, 3.4 land
together. The first iteration uses a global eval-clock and a single
shared/copy-on-write Expr. Big perf milestone.

**M5 — Dispatch and rule ordering (Tier 3 polish).** 3.5, 3.6. Becomes
worthwhile only after M3/M4 because it depends on cheap symbol
comparisons and stable ownership semantics.

**M6 — Messages / Format / NValues / Return (Tier 2 trailing + Tier 4).**
2.3, 2.4, 2.5, 4.3, 4.5, 4.6. Quality-of-life features that round out
the Mathematica fidelity.

---

## 4. Out of scope

These appear in Withoff but are deliberately not part of this plan:

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
