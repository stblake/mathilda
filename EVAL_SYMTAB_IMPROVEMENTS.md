# Next-Generation Evaluator & Symbol-Table Architecture

A full implementation spec for the next evolution of Mathilda's evaluator
(`src/eval.c`) and symbol table (`src/symtab.c`) — detailed enough to hand to an
implementer.

---

## Context

The evaluator and symbol table were built when Mathilda was a small CAS. They have since
been hardened by the milestones documented in `plans/EVAL_IMPROVEMENTS_PLAN.md` (M1
correctness, M2 interning, M3 refcount+COW+eval-clock cache, M4 dispatch-index +
specificity ordering). Those milestones are **done and are not re-litigated here** — this
document is the *next* generation, targeting the structural gaps they left behind.

The core problem this addresses: **a symbol has no direct link to its own definition.**
`EXPR_SYMBOL` stores a bare interned `char*`. Every time the evaluator needs a symbol's
attributes, builtin pointer, or DownValues, it re-derives the definition through *two*
hash tables (the interner at 8191 buckets, then `symtab` at 65535 buckets), and
`get_attributes()` additionally does a **linear `strcmp` scan of a ~140-entry table** —
all on the hottest path in the system (`evaluate_step`, once per function node per pass).

**Intended outcome:** turn symbol→definition resolution into an O(1) pointer dereference
with zero hashing and zero string comparison on the eval hot path, via a **symbol-object
model** (the design Lisp and the Wolfram kernel both use). Layer targeted allocator and
rule-indexing wins on top. Target: measurable reduction in `evaluate()` wall-time and
allocations on the existing benchmark corpus, with zero behavioral change.

---

## Current Architecture (terse — see EVAL_IMPROVEMENTS_PLAN.md for M1–M5 detail)

- **Expr** (`expr.h:64`): tagged union with `refcount` (M3 COW) and `last_evaluated_at`
  (M3 eval-clock cache). `EXPR_SYMBOL` holds `char* symbol` — an interned pointer.
- **Interner** (`sym_intern.c`): 8191-bucket chained hash, djb2. `InternEntry{char* name;
  int is_system; next;}` — **no payload beyond the string**.
- **Symtab** (`symtab.c`): 65535-bucket chained hash, djb2 on the interned pointer's bytes.
  `SymbolDef{symbol_name; own_values; down_values; builtin_func; attributes; docstring;
  default_options;}`. Rules are specificity-sorted singly-linked lists with a precomputed
  `dispatch_arity` + `first_arg_head_canon` filter (M4).
- **Attributes** (`attr.c:143`): `get_attributes()` = linear `strcmp` over `builtin_attrs[]`
  (~140 static entries) **OR-merged** with `symtab_get_def()->attributes`.
- **Eval** (`eval.c`): fixed-point loop; per-step it calls `evaluate(head)`,
  `get_attributes(head->symbol)`, then `apply_down_values`, `symtab_get_def(head)->builtin_func`.

---

## Ranked Inefficiencies (impact × frequency)

| # | Issue | Site | Cost per eval-step | Freq |
|---|-------|------|--------------------|------|
| **1** | `get_attributes()` linear `strcmp` scan of ~140-entry table | `attr.c:149` | O(140·len) strcmp | every fn node |
| **2** | Symbol→def resolution goes through 2 hash tables, re-hashing an already-canonical pointer | `symtab.c:42`, `attr.c:157` | 1 intern hash + 1 symtab hash + chain walk | every fn node (≥2×: attrs + builtin) |
| **3** | Base attributes (static) merged with dynamic attributes at *runtime* every call, never cached | `attr.c:148–160` | table scan + symtab lookup + OR | every fn node |
| **4** | `apply_down_values` / builtin dispatch each re-resolve the def independently | `eval.c:875,883` | 2× `symtab_get_def` | every fn node |
| 5 | OwnValues have **no** dispatch filter — full list scan + `match()` per rule | `symtab.c:1040` | O(n·match) | every symbol node |
| 6 | DownValue list is linear even after M4 filter; no index by (arity, head) | `symtab.c:1002` | O(n) scan | fn nodes w/ rules |
| 7 | `expr_hash`/`expr_eq`/`expr_compare` recomputed O(tree), never cached on node | `expr.c:425,495`, `sort.c:233` | O(tree) | assoc lookups, Orderless sort, fixed-point fallback |
| 8 | Every `Rule`/`SymEntry`/`SymbolDef` individually `malloc`'d — fragmentation, poor locality | `symtab.c` | alloc churn | rule insert / lookup |
| 9 | O(n²) rule insertion (dup-scan + insertion-sort walk) | `symtab.c:651` | O(n²) over n adds | startup / def time |
| 10 | Eval-clock cache is global — any Set/Clear invalidates *all* cached evaluations | `eval.c:84` | coarse invalidation | after each assignment |

Issues **1–4** dominate: they fire on *every function node on every evaluation pass* and
are all consequences of the same root cause — the missing symbol→definition link.

---

## Centerpiece Redesign: The Symbol-Object Model

### Principle

A symbol *is* its definition cell. The canonical, interned representation of a symbol name
carries a direct pointer to its `SymbolDef`. `EXPR_SYMBOL` caches that `SymbolDef*`. The
evaluator then reaches attributes / builtin / rules by pointer dereference — no hashing, no
`strcmp`, no second table.

This unifies the interner and the symbol table into one structure and makes the SymbolDef
the single source of truth for *all* per-symbol metadata, including **base attributes
precomputed once** (killing inefficiency #1 and #3 outright).

### 1. Unify InternEntry and SymbolDef

Replace the two parallel tables with one. The interner becomes the symbol table; the
`SymbolDef` *is* the intern entry.

```c
/* symtab.h — the unified symbol cell. Interned + defined in one node. */
typedef struct SymbolDef {
    const char*  symbol_name;    /* owned here now (was owned by interner) */
    uint32_t     name_hash;      /* cached djb2 of name — computed once at creation */
    uint8_t      is_system;      /* migrated from InternEntry */

    /* definition payload (unchanged fields) */
    Rule*        own_values;
    Rule*        down_values;
    BuiltinFunc  builtin_func;
    uint32_t     attributes;     /* DYNAMIC attributes only (user SetAttributes) */
    uint32_t     base_attributes;/* NEW: static builtin attrs, set once at registration */
    uint32_t     eff_attributes; /* NEW: base | dynamic, recomputed only on mutation */
    char*        docstring;
    Expr*        default_options;

    struct SymbolDef* next;      /* hash chain (replaces SymEntry) */
} SymbolDef;
```

- `SymEntry` is deleted — `SymbolDef` chains itself (removes one malloc per symbol, #8).
- `intern_symbol(name)` becomes a thin wrapper: `return symtab_intern(name)->symbol_name;`
  so every existing caller of `intern_symbol` keeps compiling unchanged.
- `symtab_intern(name)` (new, the real work): hash → chain walk (pointer-compare then
  `strcmp` fallback for fresh literals) → create-if-absent. Returns `SymbolDef*`.
- The interner's `is_system` / `intern_mark_all_system` / `intern_is_system` /
  `intern_for_each` / `intern_clear` migrate onto the unified table (walk `SymbolDef` chains).
- **Single table, single hash.** Kill `src/sym_intern.c`'s separate bucket array; keep the
  file as a compatibility shim implementing `intern_symbol` over the unified table, OR fold
  it into `symtab.c` and have `sym_intern.h` declare the shim. (Shim keeps the ~1100
  `intern_symbol` / `SYM_*` call sites and `sym_names.c` untouched.)

### 2. Precompute base attributes at registration (kills #1, #3)

The `builtin_attrs[]` table in `attr.c` becomes a *one-time seeding* pass, not a per-call
scan.

- Add `attr_seed_base_attributes(void)`: iterate `builtin_attrs[]` once at init, calling
  `symtab_get_def(name)->base_attributes = attrs;` and set `eff_attributes = base | dynamic`.
  Call it from `attr_init()` / `core_init()` after `sym_names_init()`.
- `get_attributes(name)` collapses to:
  ```c
  uint32_t get_attributes(const char* name) {
      if (!name) return ATTR_NONE;
      SymbolDef* def = symtab_get_def(name);   /* still needed for the char* API */
      return def ? def->eff_attributes : ATTR_NONE;
  }
  ```
- **New fast overload for the hot path** (see §4): `get_attributes_def(SymbolDef* def)` →
  `return def->eff_attributes;` — pure field read, zero lookup.
- `set_attributes()` / `ClearAttributes` / `SetAttributes` update `def->attributes` then
  recompute `def->eff_attributes = def->base_attributes | def->attributes`, and bump the
  eval clock (already done). Base attrs are never lost by `Clear[f]` (only dynamic cleared).

### 3. Give EXPR_SYMBOL a cached SymbolDef*

```c
/* expr.h — EXPR_SYMBOL variant */
struct {
    const char*       name;   /* interned (== def->symbol_name); kept for API compat */
    struct SymbolDef* def;    /* NEW: resolved cell, or NULL until first resolve */
} symbol;
```

- `expr_new_symbol(name)` interns the name and, cheaply, stores the resolved `def`
  (since `symtab_intern` already returns it — no extra lookup).
- Equality/printing/hash still key on `name` (pointer identity preserved), so `expr_eq`,
  `expr_hash`, `expr_compare` are **unchanged** and the `char*`-based public API is intact.
- `def` is a *non-owning cache*. Because a `SymbolDef` lives for the whole program (like the
  interned string does today — see `intern_clear` only at shutdown), the pointer is stable.
  `Clear[f]`/`Remove[f]` must **not** free the `SymbolDef` node itself (only its rules /
  docstring / options), or must null the node's payload in place — see Migration risks.
- **This is a union-size-neutral change** on 64-bit (symbol was `char*` = 8 bytes; now two
  pointers = 16 bytes, still ≤ the FUNCTION variant's 3 pointers). Confirm `sizeof(Expr)`
  does not grow the union (it will not — FUNCTION is the largest variant).

### 4. Rewire the eval hot path to pointer derefs

In `evaluate_step` (`eval.c:724–891`), once `head` is an `EXPR_SYMBOL`, resolve the def
**once** and thread it through:

```c
SymbolDef* hdef = NULL;
if (head->type == EXPR_SYMBOL) {
    hdef = head->data.symbol.def;
    if (!hdef) hdef = head->data.symbol.def = symtab_get_def(head->data.symbol.name);
    attrs = hdef->eff_attributes;                 /* was get_attributes() scan+lookup */
}
...
Expr* down = apply_down_values_def(hdef, res);    /* pass def, no re-lookup */
...
if (hdef->builtin_func) { ret = hdef->builtin_func(res); ... }  /* was symtab_get_def again */
```

- Add `apply_down_values_def(SymbolDef*, Expr*)` and keep `apply_down_values(Expr*)` as a
  wrapper (resolves def then delegates) for non-hot callers.
- Eliminates issues #2 and #4 entirely on the hot path: **zero hash-table touches per
  function node** after the first resolve.

### Expected effect

Per function node per pass: from `{2 djb2 hashes + 2 chain walks + 140-entry strcmp scan}`
down to `{1 pointer load (cached def) + 1 field read (eff_attributes)}`. The first-ever
resolve of a given symbol node still costs one lookup; every subsequent pass is O(1).

---

## Layered Incremental Wins (on top of the symbol-object core)

### A. OwnValues dispatch fast path (#5)
Most symbols have exactly zero or one OwnValue, and the common case is a bare
`x = value` (a trivial `Blank`-free LHS). Add to `apply_own_values`:
- Fast path: if `def->own_values == NULL` → return NULL immediately (already cheap once
  `def` is in hand — no lookup).
- If the single rule's pattern is the symbol itself (no pattern vars), skip `env_new`/`match`
  and return `expr_copy(rule->replacement)` directly. Guard on a precomputed
  `rule->dispatch_arity == 0 && rule->specificity` "literal LHS" flag.

### B. DownValue (arity,head) bucket index (#6)
For symbols with many DownValues (integral tables have dozens under one head), replace the
single specificity-sorted list with a **small index**: hash-map from
`(dispatch_arity, first_arg_head_canon)` → sub-list, plus one "wildcard" list for
variable-arity/any-head rules. `apply_down_values_def` probes the matching bucket + wildcard
list only. Keep specificity ordering *within* each bucket. Only build the index lazily once a
symbol crosses a threshold (e.g. > 8 DownValues) to avoid overhead for ordinary user `f[x_]`.

### C. Arena/pool allocation for Rule and SymbolDef (#8, #9)
- Pool-allocate `Rule` and `SymbolDef` nodes from a growing arena (`src/symtab_arena.c`) —
  removes per-node `malloc`, improves chain-walk locality. Nodes are never individually
  freed during a session (Clear frees *rule payload* Exprs, not the node structs → free-list
  reuse). `symtab_clear` frees whole arena blocks.
- Optional: cap dup-scan in `add_rule` (#9) using the (arity,head) index from (B) so
  insertion dup-detection is O(bucket) not O(n).

### D. Cache expr_hash on the node (#7) — *optional, measure first*
Add `uint64_t hash_cache` (0 = unset) to Expr, populated lazily by `expr_hash`. Invalidate
on `expr_unshare` (already resets `last_evaluated_at`; reset hash there too). Benefits
Association lookups and Orderless sorting. **Gate on benchmark** — only land if
`bench_assoc` or eval corpus shows net win, since it grows `Expr` by 8 bytes.

### E. Per-symbol eval-clock (#10) — *deferred / research*
Global clock invalidation is coarse but correct and simple; WL itself uses coarse
invalidation. Documented as future work, not scheduled — do **not** implement speculatively
(prior guidance: no premature complexity). Revisit only if profiling shows cache-thrash
after assignments dominates a real workload.

---

## Phased Implementation Roadmap

Each phase compiles cleanly (`-std=c99 -Wall -Wextra`), passes the full test suite, and is
valgrind-clean before the next begins. Behavior is **identical** throughout — these are pure
performance/structure changes.

### Phase 0 — Baseline & guardrails  ✅ DONE (2026-07-13)
- **`tests/bench_eval.c`** built and wired into ctest (mirrors `bench_assoc`:
  median-of-trials wall time, machine-normalized against `Total[Range[40000]]`, gates at
  `SLOWDOWN_MAX = 2.5×` a committed baseline). Six workloads chosen to expose the dispatch
  tax rather than a heavy builtin body. Baseline recorded (Apple M-series, USE_ECM=OFF):

  | workload | expr | norm (before) |
  |---|---|---|
  | pure dispatch (8k distinct) | `Length[Table[hh[k,k],{k,8000}]]` | 0.60 |
  | downvalue rewrite (Nest 4k) | `Nest[g,0,4000]` (`g[x_]:=x+1`) | 0.19 |
  | listable thread (Cos 4k) | `Length[Cos[Range[4000]]]` | 0.15 |
  | orderless sort (Plus 1200) | `Length[Plus @@ Table[c[k] x,{k,1200}]]` | 0.92 |
  | replacerepeated fold (400) | `First[Range[400] //. {a_,b_,r___}:>{a+b,r}]` | 3.00 |
  | nested plus collapse (3k) | `Nest[#+a&,x,3000] /. a->0` | 0.39 |

  `norm` = workload µs ÷ calibration µs; stable to a few % run-to-run. Phases 1 & 3 should
  push these *down* — re-record after each intended speedup so the gate keeps catching
  accidental regressions on the new floor. The **pure-dispatch row is the primary Phase 1/3
  signal** (~2.7 µs/node over 8000 distinct nodes, almost entirely dispatch).
- **Workload-design gotchas found (important for interpreting results):**
  1. A `Nest[Function[u,f[u,u]],x,d]` "tree" is *not* a dispatch test — `replace_bindings`
     shares args by refcount, so it collapses to a ~13-node DAG that the eval-clock cache
     evaluates in O(1). Use `Table[hh[k,k],{k,n}]` with distinct integer args for genuinely
     distinct nodes.
  2. The `//.` fold is O(n²); keep it small (400) or it dominates total runtime.
- **Alloc counting: DEFERRED to Phase 5.** The plan assumed "existing M3 instrumentation" —
  there is none (`expr.c` calls `malloc` directly; no counter). Portable interception
  (malloc override) is fragile across Linux/macOS, and Phases 1 & 3 are CPU-dispatch wins,
  not allocation wins. A lightweight opt-in `expr_alloc_count()` counter will be added with
  the arena work (Phase 5), which is the phase whose success is actually measured in allocs.

### Phase 1 — Fold base attributes into the def  ✅ DONE (2026-07-13)
**Design changed from the original spec — simpler and more robust.** The spec proposed a
cached `eff_attributes = base | dynamic` recomputed at every mutation. Investigation found
`def->attributes` is written directly at **13+ runtime sites** (Block-style save/restore in
findmin, nsum, nlimit, nresidue, nseries, nderiv, nint, modular; Protect/Unprotect/ClearAll
in core.c), so a cached `eff` would impose a fragile 13-site sync obligation. Instead:

- Added **`base_attributes` (cached floor) + `base_seeded` flag** to `SymbolDef` — no `eff`
  field. `get_attributes()` folds the static `builtin_attrs[]` scan into `base_attributes`
  **once per symbol** (lazily on first call; eagerly in `attr_init` for table symbols) and
  returns `base_attributes | def->attributes` — an O(1) field read + OR, computed on read.
- **Zero sync obligation:** every mutation site still writes only `def->attributes`;
  `base_attributes` is immutable after seeding, so `base | dynamic` is always correct. No
  mutation site changed.
- **Bit-for-bit identical** to the old `scan | dynamic`, including the un-clearable base
  floor (Unprotect-resistance of builtins). `attributes` keeps its exact meaning, so
  match.c's direct `def->attributes` reads are untouched.
- **Files:** `src/symtab.h` (2 fields), `src/symtab.c` (init them), `src/attr.c`
  (`get_attributes` fast path + `base_attributes_of` helper + eager seed in `attr_init`).
- **Result:** `bench_eval` pure-dispatch row **0.60 → 0.40 norm (−33%)**; downvalue rewrite
  −16%, listable −13%; dispatch-light rows (orderless sort, replacerepeated fold) flat, as
  predicted. Full relevant suite green (clearall/protect, core, eval, match, symtab,
  eval_timestamps, regression, list). Struct grew → `make clean` required (done). New floor
  re-recorded in `bench_eval.c`.

### Phase 2 — Unify interner + symtab into one table  ✅ DONE (2026-07-13)
The 8191-bucket interner and 65535-bucket symbol table are now **one** 65535-bucket table
of self-chaining `SymbolDef` nodes (the old `SymEntry` indirection is gone). Resolving a
symbol to its definition no longer hashes the name twice across two tables.

- `SymbolDef` gained `is_system` (migrated from `InternEntry`), `has_def`, and `next`, and
  now **owns** its canonical name string. `symbol` union of `EXPR_SYMBOL` is unchanged
  (still a bare `char*`), so nothing else in the tree needed touching.
- **Key insight — the `has_def` flag.** After the merge, interning a name and *defining* it
  would collapse to the same node, which would change `Names[]` (enumerate every interned
  name) and `symtab_lookup` (report interned-only names as existing). `has_def` distinguishes
  the two: `intern_symbol` creates a node with `has_def=0`; `symtab_get_def` sets `has_def=1`.
  `symtab_for_each`/`Names[]` skip `has_def=0` nodes and `symtab_lookup` returns them as NULL,
  so the observable set is **bit-for-bit the pre-merge "was symtab_get_def ever called" set**.
- **Ownership / lifetime.** `symtab_clear` and `symtab_remove_symbol` reset a node to
  interned-only (`reset_node_payload`, `has_def=0`) but **keep the node and its name string**
  — so `SYM_*` and any live `Expr` holding the interned pointer stay valid (callers clear
  between subtests and keep evaluating; `SYM_*` must survive). `symtab_init` likewise resets
  payloads instead of `memset`-ing the array (a memset would strand the interned names).
  `intern_clear` is the only full teardown (frees nodes + names). The interning API
  (`intern_symbol`, `intern_is_system`, `intern_mark_all_system`, `intern_clear`) is
  implemented in `symtab.c` over the unified table; `sym_intern.c` is now an empty stub.
- **Files:** `src/symtab.h`, `src/symtab.c` (unified node + `node_intern`/`node_find`/
  `reset_node_payload` + migrated interning API), `src/sym_intern.c` (gutted).
- **Verified:** clean build; `names`, `symtab` (3× re-init), `context` (intern_is_system /
  symtab_lookup), `core`, `eval`, `match`, `clearall_remove_protect`, `eval_timestamps`,
  `regression`, `rule_dispatch` (calls `intern_clear`) all green. Valgrind: **zero** leaked
  bytes attributable to Mathilda source — a def/remove/clear session leaks byte-identically
  to a bare `Sin[1.0]` (both = fixed macOS dyld/Accelerate baseline). `bench_eval` passes;
  pure-dispatch nudged 0.42 → ~0.36 (the merge shaves the second hash off `symtab_get_def`,
  still called by `get_attributes`). Struct changed → `make clean` required.

### Phase 3a — Resolve the head's def once, thread it through  ✅ DONE (2026-07-13)
`evaluate_step` now resolves the head's `SymbolDef*` **once** and reuses it for attribute
lookup, DownValue dispatch, and builtin dispatch, instead of re-resolving (re-hashing) the
same head up to three times per node per pass.

- New def-threaded fast paths: `get_attributes_def(SymbolDef*)` (attr.c) and
  `apply_down_values_def(SymbolDef*, Expr*)` (symtab.c). The name-taking `get_attributes` /
  `apply_down_values` become thin wrappers (`resolve → delegate`) for non-hot callers.
- `evaluate_step` resolves `hdef = symtab_get_def(head)` once and uses `hdef->builtin_func`
  and `apply_down_values_def(hdef, …)` directly. The Phase 2 def node is stable (never
  freed/reallocated), so the cached pointer stays valid even if the symbol is redefined
  mid-step; each use reads its fields fresh.
- **Files:** `src/eval.c`, `src/attr.c/h`, `src/symtab.c/h`. No struct change (no `make
  clean` needed). **Verified:** green on eval, match, core, symtab, clearall_remove_protect,
  names, regression, rule_dispatch, eval_timestamps; no new allocations (no leak surface).

### Phase 3b — Cache SymbolDef* on EXPR_SYMBOL  ✅ DONE (2026-07-13; perf-neutral, kept for completeness)
The `symbol` union member is now `struct { char* name; struct SymbolDef* def; }`; `name`
stays at union offset 0 so the historical EXPR_STRING type-pun still holds and `sizeof(Expr)`
is unchanged (the struct is 16 bytes, still ≤ the 24-byte FUNCTION variant). The evaluator
resolves the head's def once, caches it on the node (`head->data.symbol.def`), and every
later pass over that node is a pointer load.

- **Mechanical rename:** 2864 `data.symbol` → `data.symbol.name` sites across 356 files
  (perl, portable `\b`). Only **3** were in-place writers, each of which now also maintains
  the cache: `expr_new_symbol` (`.def=NULL`, lazy), `expr_unshare` (copies `.def` — same
  symbol), `alpha_rename` in symtab.c (`.def=NULL` — the name changed).
- **Why it is safe (no staleness):** Phase 2 made SymbolDef nodes *stable* — reset in place,
  never reallocated. So a cached `def` pointer always points to the current node for that
  name; `Clear`/`Remove`/re-`core_init` mutate that same node's fields in place, which the
  cache reads fresh. The only full invalidation is `intern_clear` (frees all nodes), whose
  contract already declares every interned pointer invalid. Valgrind confirms: **0 invalid
  reads/writes**, leak total byte-identical to baseline.
- **Measured result: no change.** Pure-dispatch is flat at ~14,200 µs across Phases
  1/2/3a/3b — exactly as predicted. After Phase 1 removed the strcmp scan, the eliminated
  lookup is <2% of per-node cost, and for distinct nodes it does not even amortize. 3b is
  correct and clean but delivers no measurable speedup; it was implemented for architectural
  completeness at the maintainer's request.
- **Files:** `src/expr.h` (union), `src/expr.c` (constructor + unshare), `src/symtab.c`
  (alpha_rename), `src/eval.c` (cache read/fill), plus the mechanical rename across `src/` +
  `tests/`. Struct definition changed → `make clean` required.

### Phase 4 — OwnValue fast path + optional DownValue index
- Implement (A) OwnValue literal-LHS fast path.
- Implement (B) lazy (arity,head) DownValue index behind the > 8-rule threshold; keep
  specificity within buckets; verify integral-table-heavy heads still match identically.
- **Files:** `src/symtab.c`. **Test:** `test_downvalue*.c`, integral-table regressions
  (`Method->"CRCTable"`), `test_condition_downvalue.c`.

### Phase 5 — Arena allocation (optional, gated on Phase 0 alloc numbers)
- Pool `Rule`/`SymbolDef` via `src/symtab_arena.c`. **Files:** new arena module,
  `src/symtab.c`. **Test:** valgrind (no leaks; arena freed wholesale in `symtab_clear`).

### Phase 6 — expr_hash node cache (optional, gated on bench)
- Only if `bench_assoc` / `bench_eval` show a win. **Files:** `src/expr.c/h`. Reset in
  `expr_unshare`. **Test:** `bench_assoc` scaling gates must stay green.

---

## Migration Risks & Mitigations

- **SymbolDef lifetime vs. cached `def`:** `Clear[f]`/`Remove[f]` must not free the
  `SymbolDef` node (Exprs cache it). Change `symtab_remove_symbol` to clear payload
  (rules/docstring/options, reset `builtin_func`, `attributes`, keep `base_attributes` +
  name) rather than unlink+free the node. `symtab_for_each`/`Names[]` semantics: a
  payload-cleared node still exists — filter emptied nodes in enumeration if needed to
  preserve current `Names[]` output.
- **Init-order contract** (SYM_ pointers NULL until `sym_names_init`): a symbol Expr built
  before the table exists must resolve `def` lazily on first eval, never at construction.
  `expr_new_symbol` must tolerate a NULL table (store `def = NULL`, resolve later).
- **Incremental-build ABI skew** (NDArray lesson): changing `Expr`/`SymbolDef` layout →
  phantom segfaults from stale objects. Mandate `make clean` in the phase commits and note
  it in the changelog.
- **`sizeof(Expr)`:** confirm the symbol variant stays within the union's largest member;
  Phase 6's hash field is the only intentional growth (gated).
- **Ownership move of the name string** (Phase 2): the interner used to own it; now the
  SymbolDef does. Exactly one owner must free at shutdown — audit `intern_clear` becomes a
  no-op or delegates to `symtab_clear`.

---

## Verification

1. **Correctness (non-negotiable, every phase):**
   - `cd tests/build && cmake .. && make -j && for t in *_tests; do ./$t; done` — full suite
     green (grep output for `FAIL:` per the NDEBUG-soft-assert convention).
   - New builtins list-sync: no new files, so `COMMON_SRC` unaffected except Phase 5 arena
     (add `symtab_arena.c` to `tests/CMakeLists.txt` COMMON_SRC).
   - `valgrind --leak-check=full ./Mathilda` on a scripted session (define symbols, Clear,
     Remove, redefine, Names[]) — diff against a `Sin[1.0]` baseline (macOS dyld noise).
2. **Behavioral invariance:** run the integrate/simplify/sum regression corpora and diff
   output against `main` — must be byte-identical.
3. **Performance (the point):**
   - `tests/build/bench_assoc` scaling + absolute gates stay green.
   - New `tests/bench_eval.c`: report `evaluate()` wall-time and alloc-count deltas vs.
     Phase 0 baseline. Phases 1+3 should show the largest drop (attribute + def resolution).
   - Confirm the `get_attributes` strcmp scan and second hash lookup are gone via a profiler
     sample (perf/Instruments) on the symbol-heavy corpus.
4. **Docs:** update `Mathilda_spec.md` changelog table + this week's
   `docs/spec/changelog/<Monday>.md`; note the `make clean` requirement and the
   interner/symtab unification in `SPEC.md` §3.3.

---

## Non-Goals

- UpValues / SubValues (explicitly out of scope in EVAL_IMPROVEMENTS_PLAN.md; unchanged).
- Per-symbol fine-grained eval-clock (deferred; documented as future research only).
- Any change to evaluation *semantics*, pattern-matching results, or the public `char*`
  symbol API.
