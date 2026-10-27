# Symbol table improvements

Status: proposal, for a dedicated session. Written 2026-10-27 while reducing the
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

### 1. Reclaim `Module` locals (Temporary semantics)

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

### 2. Lazy names for `Module` locals

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

### 3. Integer symbol ids instead of `char*` identity

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

### 4. Copy the `Module` body once, not twice

This is about the rule layer, not the table itself. Applying
`diracTr[a_, rest___] := Module[...]` costs two copies:
- `replace_bindings` copies the right-hand side to substitute `a` and `rest`;
- `builtin_module` copies the body again (`substitute_scoping`) to rename the
  locals.

`builtin_module` could rename into the tree `replace_bindings` just produced,
since that tree is uniquely owned (`refcount == 1`) when fresh. Or one pass
could do both substitutions, given a `MatchEnv` extended with the
local → `local$N` renames.

### 5. Insertion-ordered iteration

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
