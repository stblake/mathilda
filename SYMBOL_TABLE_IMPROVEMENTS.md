# Symbol table improvements

Status: every proposal has been worked through.
- **Done:** 1 (v0.6056); 2, 4 and 5 (v0.6057, in the form the evidence
  supported).
- **Declined, on measurements:** 3.

Each section records the outcome. Written 2026-10-27 while reducing the
per-call interpreter overhead of issue #112 (v0.6055). It records what the
current design costs, measured, and the redesigns worth doing. A smaller,
order-preserving change was already made and is described under
[What v0.6055 already did](#what-v06055-already-did).

## The current design

`src/symtab.c` has one table for interning and definitions.
- **Nodes.** A `SymbolDef` node exists as soon as a name is interned. It owns
  the name string and every rule list (own, down, sub, up, format and N
  values), plus attributes, the builtin pointer, options and kernel
  descriptors. `has_def` separates a real definition from an interned-only
  name.
- **Identity.** A symbol is identified by its interned `char*`. `SYM_*`
  globals are such pointers, and the code compares names by pointer
  throughout.
- **Lookup from an Expr.** An `EXPR_SYMBOL` caches its `SymbolDef*` lazily
  (`data.symbol.def`).
- **Nodes are never freed** except by `intern_clear`, the full teardown used
  between test runs.
  - `Remove`/`ClearAll` reset a node to interned-only.
  - `Module` mints a fresh `name$N` per local per call and never reclaims it.

## What it costs: measured

All figures are for the issue #112 skeleton: the trace's call structure with
integer arithmetic, n = 14, 190k `Module` calls.

| | Mathilda v0.6055 | reference |
|---|---|---|
| symbols left after one run | **+111,005** | **+11** |
| memory left after one run | **+44 MB** | +0.005 MB |
| run time | 0.66 s | 0.71 s |

- **The reference reclaims locals.** Its `Module` locals carry `Temporary`, and
  a local is removed once nothing references it. Only the 11 locals still
  referenced survive.
- **Mathilda keeps them all, for the whole session.** Every local stays as a
  `SymbolDef`, with its name string, its `own_values` rule and the rule's
  pattern and replacement trees.
  - A long-running session or a loop over `Module` grows without bound.
  - The interner's working set grows with it: every new name walks a longer
    chain.

Profile of the skeleton, before v0.6055:
- `intern_symbol` called from `builtin_module` (minting the local's name) was
  about 10% of samples, plus `snprintf` about 3%.
- That cost was almost all chain walking over cache-missing nodes. 2M locals
  in 65535 buckets averages about 35 nodes per chain, and each miss
  `strcmp`'d the chain.
- Rule binding (`rb_rec`, `substitute_scoping`) was another 15–20%. It copies
  the `Module` body twice per call: once to substitute the pattern variables,
  once to rename the locals.

## What v0.6055 already did

These changes are behaviour-preserving and order-preserving where it matters.
- **Stored hashes.** `SymbolDef.name_hash` holds the full 32-bit hash, so a
  chain walk compares integers and calls `strcmp` only on a hash hit.
- **A growable table.** Buckets are a power of two, doubled when the table
  holds more than twice as many nodes (`symtab_grow`). Nodes never move, so
  every `SymbolDef*` stays valid.
- **Deterministic disassembly names.** The only consumer that depended on
  bucket order was `compile/disasm.c`'s kernel-name lookup, which took the
  first symbol carrying a kernel. It now takes the alphabetically first.
  `Names` already sorted its output.

The result was 12% off the skeleton. It does nothing about the growth itself.

## Proposals, in order of value

### 1. Reclaim `Module` locals (Temporary semantics) — DONE in v0.6056

How it was implemented (`src/symtab.c`, "Module-local reclamation"):
- **Fresh names only.** `Module` flags a local `module_temp` only when its name
  was never interned, so no uncounted node can exist.
- **Counting.** Every node naming the local caches its `SymbolDef`;
  `temp_nodes` counts live nodes, adjusted in `expr_new_symbol`,
  `expr_unshare`, `alpha_rename` and `expr_free`.
- **The reclaim test.** A local is reclaimed when a walk of its own rule lists
  finds all of its nodes, each through refcount-1 parents. A shared parent, or
  a rule list detached by `Block` or an iterator, keeps it.
- **Timing.** A local is tried at its own `Module`'s exit. An escaped local is
  queued and retried only when the outermost evaluation ends, and only if
  something was queued or a queued local lost a node since the last sweep.
  That restriction comes from the audit: a sweep in a nested frame could free
  a `SymbolDef` that an outer rule loop is still reading.
- **Names are kept.** Name strings go to a graveyard freed by `intern_clear`,
  so no raw name pointer dangles or aliases a later symbol.
- **Kill switch.** `MATHILDA_NO_TEMP_RECLAIM=1`.
- **Result.** After the trace, +9 symbols and +2.7 MB, against +111,005 and
  +44 MB before (the reference: +11). `Module` costs 50 B per call (the kept
  name strings), down from 548. Clean under AddressSanitizer.

The original proposal follows.

This is the largest win: memory, interner working set, and fidelity to the
reference.

- **Reference behaviour.** A `Module` local is `Temporary`. When no expression
  refers to it any more, the symbol and its values are removed.
  `Module[{x}, x]` returns `x$123`, which survives while referenced.
- **Design.**
  1. Give `SymbolDef` a reference count of live `Expr` symbol nodes that point
     at it. `expr_new_symbol` / `expr_free` of an `EXPR_SYMBOL` would
     inc/dec it.
     - Only symbols flagged `Temporary` need it. A flag check keeps the cost off
       every other symbol.
     - Rule lists hold `Expr`s, so references from definitions are counted for
       free.
  2. When a Temporary symbol's count reaches 0 (or at `Module` exit, with a
     deferred sweep), unlink the node from its bucket, free its rules, and
     free the name.
  3. Name reuse needs care. A freed `name$N` must not be minted again while
     anything caches the old pointer: `SYM_*` never applies to locals, but
     `data.symbol.def` caches do. The counter is monotone, so names are not
     reused today, and freeing a node whose count is 0 means no `Expr` caches
     it.
- **Risks to audit.**
  - Every place that stores a bare `const char*` symbol name outside an `Expr`
    would dangle when the node is freed. Examples: `Names` collection,
    `MatchEnv.symbols`, memo tables keyed by name, `VarSet` in
    `flint_bridge.c`. Each needs to copy, hold a reference, or be proven
    short-lived.
  - `eval.c`'s `inert_observed` flags and the literal index (`lit_index`) live
    on the node, and they die with it, which is correct.
  - `Block` saves `old_own` and other rule lists by symbol *name*
    (`mathilda_strdup`). That is safe, but a Temporary symbol that is blocked
    must not be reclaimed mid-`Block`.
- **Verification.**
  - `symcount.wl` (this document's measurement) should go from +111,005 to
    about +10.
  - Under valgrind or ASan, a `Module`-heavy loop should be flat in memory.
  - The full ctest suite.
  - A new test: `Module[{x}, x]` survives while held and is reclaimed after
    `Clear`.

### 2. Lazy names for `Module` locals — DONE in v0.6057 (reduced form)

**Outcome.** Fully lazy names were not built: a local's name is needed by
canonical ordering (`symbol_cmp`), hashing (`expr_hash` hashes the name
bytes) and printing. So any local that enters an Orderless sum or product,
a hash-keyed structure, or output would materialise its name anyway, and the
laziness would add a deferred path to 4814 name reads for little saving.

What was removed instead is the redundant interning. Creating a local hashed
and walked its fresh name four times: a freshness check, `symtab_get_def`,
`expr_new_symbol(buf)` (a stack buffer never hits the pointer cache), and
`symtab_add_own_value(buf)`. `$ModuleNumber` was also re-interned from a
string literal on every call. Now:
- `symtab_intern_def` does one walk and reports "created";
- `expr_new_symbol_def` builds the node from the `SymbolDef`;
- `symtab_add_own_value_def` sets the value from the `SymbolDef`;
- `$ModuleNumber` is updated through `SYM_DollarModuleNumber`.

Measured: no change. Once v0.6055's growable table and v0.6056's
reclamation kept the table small, interning was no longer a measurable
cost. It is kept as the simpler code.

The original proposal follows.

Creating a local currently costs `snprintf` + intern + `calloc` + `strdup`. An
alternative is anonymous `SymbolDef`s that are not entered into the name table
at all.
- The printed name `name$N` is materialised only when something asks for it:
  printing, `ToString`, `Names`, `SymbolName`, or a symbol-table lookup by
  string.
- The local is reached only through `data.symbol.def`, which `Module` sets when
  it builds the renamed body.
- This removes the interner from the per-call path entirely.
- It needs (1)'s ownership model, and an audit of code that reads
  `data.symbol.name` expecting a table-registered string. Every pointer
  comparison against `SYM_*` stays valid, because locals are never `SYM_*`.

### 3. Integer symbol ids instead of `char*` identity — DECLINED (measured)

**Outcome.** Not done, because the precondition the proposal set ("only if
profiles show string-keyed side tables mattering") does not hold.
- **Symbol identity is already as cheap as ids.** It is an interned-pointer
  compare, exactly the cost of comparing ids.
- **The string costs that did show have been removed:**
  - interner chain walks, by v0.6055;
  - per-call re-interning, by item 2.
- **The matcher's binding environment (`env_get`/`env_set`) is
  pointer-first.** Its `strcmp` fallback, kept for the two literal keys
  `$OptionsPattern$` and `$OptionsPatternHead$` and for regex group names, was
  about 1.5% of a rule-application profile and usually exits at the first
  character.
- **The trade.** Converting 4814 `data.symbol.name` sites would add broad risk
  for no measurable gain.

**Revisit if** a profile shows a dense per-symbol side table (an array
indexed by symbol) being needed.

The original proposal follows.

Identity by `char*` makes every "is this the same symbol" a pointer compare,
which is cheap. But hashing needs the string, and any per-symbol side table
must be a hash map. Alternatives:
- A dense `uint32_t id` per `SymbolDef`, with side tables (e.g. a pattern
  matcher's binding slots, `Expand`'s generator sets) indexed by id.
- A `SymbolDef*` in `Expr` as the primary identity, with the name hanging off
  it.

This is a large mechanical change, touching 4814 `data.symbol.name`
sites (counted 2026-10-27). It is worth it only together with (2), and only if profiles show
string-keyed side tables mattering.

### 4. Copy the `Module` body once, not twice — DONE in v0.6057

**Outcome.**
- **No single pass.** Merging the rule substitution and the `Module` renaming
  into one pass would be unsound: `Module` must rename when it is evaluated
  (`$ModuleNumber`, nested scopes), not when the rule fires.
- **Sharing instead.** `substitute_scoping` now shares every subtree no local
  occurs in, including iterator specs and binding lists, as `replace_bindings`
  has done since v0.6052. Each pass therefore copies only the paths that
  change.
- **Faster lookups.** Interned names are compared by pointer instead of
  `strcmp`.
- **The larger cost was allocation, not copying.** `rb_rec` and
  `substitute_scoping` each made three, respectively one, heap allocations per
  node. They now use stack buffers for arity up to 8, since
  `expr_new_function` copies the array.
- **Measured.**
  - rule+`Module` loop: 2.60 → 1.81 µs per call (the reference: 2.27).
  - Bare `Module`: 1.58 → 1.50 µs.
  - The #112 skeleton: about 3.1 → 2.65 s for six runs.

The original proposal follows.

This is about the rule layer, not the table itself. Applying
`diracTr[a_, rest___] := Module[...]` costs two copies:
- `replace_bindings` copies the right-hand side to substitute `a` and `rest`;
- `builtin_module` copies the body again (`substitute_scoping`) to rename the
  locals.

`builtin_module` could rename into the tree `replace_bindings` just produced,
since that tree is uniquely owned (`refcount == 1`) when fresh. Or one pass
could do both substitutions, given a `MatchEnv` extended with the
local → `local$N` renames.

### 5. Insertion-ordered iteration — DONE in v0.6057

**Outcome.** Every `SymbolDef` sits on an intrusive doubly-linked
insertion-order list (`ins_prev`/`ins_next`). Unlinking on reclamation is
O(1). `symtab_for_each` and the other whole-table walks (init reset, clear,
`intern_clear`, `intern_mark_all_system`) follow definition order,
independent of the bucket array and its resizes. `disasm` keeps its
alphabetical tie-break, which is harmless and deterministic either way.

The original proposal follows.

Anything that enumerates the table (`symtab_for_each`) sees bucket order,
which v0.6055 changed and any future resize changes again.
- An intrusive insertion-order list (a `next_inserted` link) would give a
  stable, meaningful order (definition order, as the reference's `Names` is
  not) for free.
- That removes the need for consumers like `disasm` to impose their own order.

## Suggested plan for the session

1. Measure the baseline: `symcount.wl`, the #112 benchmark, and a
   `Module`-in-a-loop RSS curve.
2. Implement (1), reference counting for Temporary symbols only. Audit the
   bare-name holders listed above, and add the reclamation tests.
3. If the interner still shows in profiles, do (2).
4. Do (4) independently of the table work: it is local to `modular.c` and
   `match.c`.
5. Treat (3) and (5) as optional, decided on profiles after 1–4.
