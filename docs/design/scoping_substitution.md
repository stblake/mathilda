# Scoping substitution: binder summaries and the Module clock

Status: implemented, v0.6653–v0.6657. This document replaces the proposal
`EVAL_CHANGE_PROPOSAL.md` (2026-11-21) as the reference. The proposal still holds the
measurements behind the decision (§2) and the rejected options (§5). This file is the
design as built: the invariants, the switches, and the three places where the build
departs from the proposal's prototype.

## 1. What problem this solves

Every substituting walk visited every node of the tree it was given:
- `substitute_scoping`, the Module/With renaming (`src/modular.c`);
- `rb_rec`, rule application (`src/match.c`);
- `substitute_names` and `substitute_slots`, pure-function application
  (`src/purefunc.c`).

v0.6057 made those walks share unchanged subtrees instead of copying them, but they
still *visited* them. On the slow `.m`-bound integration cases, 97% of the visits came
back unchanged. The worst was `goodPrime` in `TorsionModP.m`: a rule injects a large
context `Association` into a Module body four times, and the Module walk then traverses
all of it on every call, though capture avoidance guarantees it holds no local.

There was a second cost. Module ran two global evaluation-clock bumps per call: one to
rewrite `$ModuleNumber` and one for each fresh local's first binding. Each bump
invalidated every exact-clock fixed-point stamp in the session and failed every open
dependency frame.

## 2. Option D: Module does not move the clock

**The first binding of a fresh local** goes through `symtab_bind_fresh_local`
(`src/symtab.c`), an OwnValue store with no clock bump. "Fresh" means
`symtab_intern_def` just created the name. No node named it before Module built its
replacement node, so no stamp, INERT certificate or dependency frame can depend on its
binding. A pre-existing `x$n` (a user's own, or a reset `$ModuleNumber`) still bumps.
Later `Set`s in the body are ordinary bindings and still bump.

**`$ModuleNumber`.** The C counter `module_number` is authoritative. The symbol holds a
delayed OwnValue, `$ModuleNumber :> Internal`ModuleNumber[]`, and that builtin reads the
counter.
- `Internal`ModuleNumber` is not `eval_pure`, so a dependency frame that reads
  `$ModuleNumber` is tainted rather than certified against a value that now moves without
  a clock bump.
- A user `Set`, `SetDelayed` or `Clear` replaces the sentinel rule.
  `module_number_sync()` notices on the next Module or Unique call (the stored
  replacement is no longer the sentinel node), adopts an integer value, and reinstalls
  the sentinel.
- `Module`, `Unique` and `unique_make_batch` are the only readers and writers.

Known divergence: `OwnValues[$ModuleNumber]` shows the sentinel rule. A
`Block[{$ModuleNumber = n}, …]` that runs Modules leaves the counter where the Block's
Modules left it, instead of restoring the outer value.

`MATHILDA_MODULE_CLOCK_BUMP=1` restores both bumps, for A/B measurement.

## 3. Binder summaries (`src/symsum.{c,h}`)

### 3.1 The word

A 64-bit word in a FUNCTION node's `data.function.index` slot. The slot is shared with
the Association index, which is a malloc-aligned pointer, so the low bit tells the two
apart. `sizeof(Expr)` is unchanged.

| bits | meaning |
|---|---|
| 0 | tag: 1 = summary (0 = `AssocIndex*` or empty) |
| 1–16 | epoch the summary was computed in |
| 17–62 | 46 binder bits |
| 63 | SCOPE: the subtree has a `Module`/`With`/`Function` head |

Read the Association index only through `expr_fn_assoc_index()`, which returns NULL for
an odd word, and write it only through `expr_fn_set_assoc_index()`. An Association that
carries an index has no summary: walks descend into it, and its entries are summarised.

### 3.2 Which symbols get bits

- **Binder names.** These are names that a Module/With/Function binds (registered by the
  walker's wrapper before it walks), or that a rule binds as a pattern variable.
  Registration goes through `symsum_register_binder`, which assigns bits round-robin in
  first-registration order, so the assignment is deterministic for a given load order.
- **The SCOPE bit**, shared by the three scoping heads.
- **Every other symbol** contributes nothing.

Giving every symbol a hashed bit was simulated in the proposal: it saturates, because
every head is a symbol, and it removes only 38% of the visits. With binder-only bits, a
data value summarises to almost nothing and is skipped at its root.

### 3.3 When a summary is valid

A summary is valid when its tag is set, its epoch equals `g_symsum_epoch`, and the fast
paths are on. A symbol summarised *before* it became a binder is missing from that
summary. So registering a binder advances the epoch, which invalidates every summary,
but only when a summary of the current epoch saw that symbol as a bitless leaf. The
"seen" marker is `SymbolDef.ss_seen_epoch`, set in `symsum_leaf_bits` during a fill.

Most names are registered before any summary sees them, so registering them costs
nothing. If the 16-bit epoch is ever exhausted, the fast paths turn off for the rest of
the session.

**Departure 1.** The prototype stored an 8-bit epoch, compared it modulo 256, and bumped
it on every registration, 294 times in one Trager run. Summaries from 256 epochs earlier
therefore aliased the live epoch and were trusted. The prototype's measured gains came
from a build that could skip unsoundly.

### 3.4 Where a summary is filled

Summaries are filled only on **frozen** nodes, and are computed bottom-up from the
children's during the walk that is happening anyway (`symsum_fill`: O(arity), and only if
every FUNCTION child already has a valid summary). A node is frozen if it is either:
- reached through a node with `refcount ≥ 2`, which makes it immutable together with
  everything below it under copy-on-write; or
- part of a stored rule's right-hand side (`replace_bindings_stored`).

A fresh, rebuilt node never gets a summary.

A walk's *root* is frozen when its parent is shared. The body of a Function, Module or
With node is held only by that node, so the body's own refcount is 1 even when the
Function stored in an OwnValue is shared. `apply_pure_function` therefore passes
`head->refcount >= 2`, and Module and With pass `res->refcount >= 2`. Without this, a
stored Function's body was never summarised.

Invalidation:
- `expr_new_function` and pool recycling start the slot at NULL.
- `expr_invalidate_hash` clears a summary; every in-place rewrite already calls it.
- `expr_unshare` clears it in both branches. In the `refcount == 1` branch the node is
  returned for in-place mutation, and it may have been summarised while it was still
  shared.

**Departure 3.** The prototype treated the root of *every* `replace_bindings` call as
frozen, including a `ReplaceAll` rule owned by the expression being evaluated.

### 3.5 Consumers

Each consumer skips a FUNCTION subtree whose valid summary is disjoint from its mask, and
returns the node shared.

| walker | mask | note |
|---|---|---|
| `substitute_scoping` | bits of the env's old names | symbol visits pre-filtered by their own bit |
| `rb_rec` | pattern-variable bits \| SCOPE | see departure 2 |
| danger set (`rb_danger_build`) | — | value summaries ∩ the scope's local bits decide whether names must be collected |
| `substitute_names` / `substitute_slots` | parameter bits / `Slot` + `SlotSequence` | |

**Departure 2.** `scoping_capture_avoid` alpha-renames any binding scope whose local
collides with the danger set, *whether or not the scope contains a pattern variable*:
`f[x_] := {x, Module[{e}, e]}; f[e + 1]` renames the Module's `e`. A skip keyed on the
pattern variables alone would leave that Module unrenamed, and the result's local names
would change. Including SCOPE in the mask skips only subtrees that hold neither a pattern
variable nor a binding scope, which is exactly what the slow walk leaves untouched.

### 3.6 Measured (v0.6652 → v0.6657)

Instructions retired, from `tools/scope_bench/ab.sh`, with every answer hash identical:

| case | v0.6652 | v0.6657 | change |
|---|---|---|---|
| Trager `((1+x^2)(2+x^4)^(1/6))/x` | 75.49 G | 64.63 G | −14.4% |
| `cM-1-292` (Bronstein, `Log` tower) | 26.00 G | 22.67 G | −12.8% |
| the other nine corpus cases | | | −2.6% … +1.7% |

Steps on Trager: D −3.6%, 3a −9.9%, 3b −0.8%, 3c −0.5%, 3d −0.2%.

Micro rows (`micro.m`, µs per call, 10,003-leaf dead branch):

| row | v0.6652 | v0.6657 |
|---|---|---|
| rule only | 86 | 0.7 |
| Module | 184 | 2.3 |
| With | 181 | 1.0 |
| Function | 84 | 0.7 |
| Modules nested 8 deep | 100 | 8.5 |

Every row is now independent of the size of the dead branch.

## 4. Switches

| switch | effect |
|---|---|
| `MATHILDA_NO_SYMSUM=1` | every consumer off; the A/B baseline and the field fallback |
| `MATHILDA_SYMSUM_DIFF=1` | every fast walk is re-run slowly and compared with `expr_eq`; aborts on a difference (rule applications that renamed are not compared, since renames mint fresh names) |
| `-DMATHILDA_SYMSUM_VERIFY` | every summary a walker relies on is recomputed exactly; aborts if it is stale |
| `-DMTH_NO_SYMSUM` | consumers compiled out; forced on platforms without 64-bit pointers |
| `MATHILDA_MODULE_CLOCK_BUMP=1` | Option D off |

## 5. Invariants a future change must keep

- Never read or write `data.function.index` directly. Use `expr_fn_assoc_index()` and
  `expr_fn_set_assoc_index()`.
- A rewrite of a live node's head or arguments in place must call
  `expr_invalidate_hash`. This was already true for the structural hash, and it now also
  protects the summary.
- A new substituting walker, or a new scoping head, must take part. Either register the
  names it binds and give the head the SCOPE bit, or the summaries' meaning ("no binder
  name occurs below") no longer covers it.
- Copy-on-write: a node with `refcount ≥ 2` is immutable together with everything below
  it.

## 6. Verification

- `tests/test_scoping_semantics.c` holds the scoping probes, with expected values taken
  from the pre-change binary.
- `tests/test_symsum.c` covers the substrate.
- `tools/scope_bench/` has `micro.m`, `semantics.m`, and `ab.sh` (an instruction-count
  A/B over the 11-case corpus).
- Run the unit tier under `MATHILDA_SYMSUM_DIFF=1` and in a `-DMATHILDA_SYMSUM_VERIFY`
  build.
