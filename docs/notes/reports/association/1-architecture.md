## 1. Architecture & implementation

### 1.1 Representation: plain expression plus a side index

There is no `EXPR_ASSOCIATION` type. An association is an ordinary `EXPR_FUNCTION` node, `Association[Rule[k1,v1], ...]` (src/assoc.h:7-15). `is_association()` only checks the head (src/assoc.c:105-107). The one addition is an acceleration pointer placed in the union's spare space, so `sizeof(Expr)` does not change (src/expr.h:200-212):

```c
struct { struct Expr* head; struct Expr** args; size_t arg_count;
         struct AssocIndex* index; } function;
```

It is set to NULL on allocation (src/expr.c:267), reset on a physical copy (src/expr.c:531), and freed with the node (src/expr.c:581). `expr_eq`, `expr_hash` and `expr_compare` ignore it (src/expr.h:206-207). The parser builds the raw node: `<|` and `|>` are structural delimiters, not operators (src/parse.c:700-738, 871). Entries are parsed as ordinary expressions, and `builtin_association` does the validation.

### 1.2 The hash indexes

There are two indexes with the same table layout:

| | `KeyIndex` (transient) | `AssocIndex` (persistent) |
|---|---|---|
| Where | stack, src/assoc.c:34-77 | heap, on the node, src/assoc_index.c:17-22 |
| Used by | construction, bulk Lookup, KeyDrop/KeyTake, KeyUnion, Counts, GroupBy | single-key Lookup, KeyExistsQ, `a[k]`, `a[[k]]`, `Key[k][a]` |
| Table | open addressing, linear probe, `size_t` slots holding 1-based entry positions, 0 = empty | same |
| Capacity | smallest power of two above `2n+1`, minimum 8; load factor below 0.5 (src/assoc.c:42-47, src/assoc_index.c:26-31) | same |
| Growth | none: sized once, never rehashed | none: rebuilt from scratch instead |
| Key equality | `expr_eq` (SameQ-like) | `expr_eq` (src/assoc_index.c:73) |

The hash is `expr_hash`, an FNV-1a structural hash that is cached per node in `hash_cache` (src/expr.c:870-896, src/expr.h:170-182). Keys are compared exactly: `Lookup[<|1->"int", 1.->"real"|>, {1, 1.}]` gives `{int, real}`.

**Index lifecycle.** The index is built lazily, on the first single-key read (src/assoc.c:297-317). It is not built at canonicalisation, because `evaluate()` keeps the original node and discards the rebuilt one (src/assoc.c:225-232, src/assoc_index.h:20-27). It is never invalidated. It stays correct only because canonical associations are never modified in place: every update produces a new node, so a stale index can only be discarded, never read (src/assoc_index.h:29-33). The compiled VM builds the index ahead of time at its marshalling boundary, so worker threads never trigger the lazy build (src/compile/compile_vm.c:389, 433, 451).

### 1.3 Semantics

- **Order and duplicates.** `assoc_from_rules` keeps each key at the position of its first occurrence and gives it the value of its last occurrence (src/assoc.c:142-160). `<|x->1, y->2, x->3|>` gives `<|x -> 3, y -> 2|>`.
- **Fast path for canonical input.** If the input has no splicing and no duplicate keys, the builtin returns NULL and the node is left as it is (src/assoc.c:216-233). Splicing of lists and associations happens in `collect_entries` (src/assoc.c:182-199).
- **Rule vs RuleDelayed is inconsistent.** A canonical literal keeps `RuleDelayed`. Any rebuild goes through `make_rule` (src/assoc.c:151, 156), which turns delayed entries into `Rule`:

| Input | FullForm output |
|---|---|
| `<|p :> 1+1|>` | `Association[RuleDelayed[p, Plus[1, 1]]]` |
| `<|p :> 1+1, q -> 2, p :> 2+2|>` | `Association[Rule[p, 4], Rule[q, 2]]` |
| `Association[{p :> 1+1}]` | `Association[Rule[p, 2]]` |

  `KeyDrop` does keep `RuleDelayed`, because it copies whole entries (src/assoc.c:505). Point updates also keep it, because `assoc_entry_with_value` copies the head (src/assoc.c:118-121).

### 1.4 Memory and mutation

`expr_copy` only bumps a reference count, so entries are shared between associations. Every value rewrite must go through `assoc_entry_with_value` (src/assoc.h:134-145). The header records a real past bug: `u2 = u1; u2[["x"]] = 9` used to corrupt `u1`. The fix is in place (src/part.c:99-115), and a runtime check confirms it: `b=<|x->1|>; c=b; c[x]=9` leaves `b = <|x -> 1|>` and gives `c = <|x -> 9|>`.

`a[k] = v` is handled inside `Set`. When the symbol's OwnValue is an Association, the assignment is rewritten as `Part[a, Key[k]] = v` (src/eval.c:1098-1133). `expr_part_assign` then evaluates the symbol, copies all n arguments, finds the key with a **linear scan that ignores the index**, appends the key if it is missing, and stores a fresh association as the new OwnValue (src/part.c:218-293, 406-428). `AssociateTo` re-canonicalises the whole association on every call (src/assoc.c:1057-1089).

Hazards:
- The lazy build casts away `const` and writes into a node that may be shared (src/assoc.c:316). This is safe only while evaluation is single-threaded.
- `expr_unshare` returns a node unchanged when its refcount is 1 (src/expr.h:245-254). Code that unshares and then edits `args[]` in place would keep a stale index. The only protection is the documented rule that nobody does this.
- Unset does not remove a key. `u=<|p->1|>; u[p]=.` leaves `<|p -> 1|>`, and `KeyDropFrom` does not exist (`Names["KeyDropFrom"]` is `{}`).

### 1.5 Evaluation, attributes, printing

`Association` carries only `{Protected}` (src/assoc.c:1477-1478). It is not atomic and not `HoldAllComplete`: `AtomQ[<|x->1|>]` is `False`, and values are evaluated, so `<|p -> Print["side"]|>` prints `side` and stores `Null`. `ReplaceAll` reaches into keys, and a rename that creates a duplicate collapses it: `<|x->1, y->2|> /. x->y` gives `<|y -> 2|>`. Loop reuse depends on the evaluator's timestamp and "ground" fixed-point short-circuit, which lists `Association` as a pure constructor (src/eval.c:180-190, 2368-2390). `<|...|>[k]` and `Key[k][a]` are special cases inside `evaluate_step` (src/eval.c:2037-2075, 2076-2100). `KeyValuePattern` is implemented in the matcher (src/match.c:593-610). The printer uses `<|...|>` only when every argument is a two-argument Rule or RuleDelayed (src/print.c:34-45, 842-850); otherwise it falls back to `Association[...]`.

### 1.6 Complexity (measured on the built binary)

| Operation | Cost | Evidence |
|---|---|---|
| `a[k]`, `Lookup`, `KeyExistsQ` | O(1) after an O(n) first build | 10^4 lookups: 0.8 ms at n=10^4, 1.1 ms at n=8×10^4 |
| Construction, bulk `Lookup[a,{ks}]` | O(n + m) | src/assoc.c:133-168, 408-432 |
| `a[k] = v` (overwrite or insert) | **O(n)** per operation | 1000 overwrites: 0.63 s at n=10^4, 2.79 s at n=4×10^4 |
| n/10 inserts into an empty association | **O(n²)** overall | 0.014 / 0.064 / 0.30 / 1.64 s for n = 1, 2, 4, 8 ×10^4 |
| `KeyDrop`, `KeyTake` | O(n + m) | 100 calls: 68 ms at n=10^4, 273 ms at n=4×10^4 |
| `Keys`, `Values`, iteration | O(n), deep-copy by refcount bump | src/assoc.c:250-270 |

### 1.7 Design risks

- **The index and the scan disagree on malformed nodes.** `assoc_index_build` accepts any two-argument node (src/assoc_index.c:44-46). `assoc_scan` requires a Rule or RuleDelayed (src/assoc.c:286). As a result, `Lookup[Association[f[1,2]], 1]` returns `2` and `KeyExistsQ` returns `True`, while `Keys[...]` stays unevaluated.
- **Comments contradict the code.** src/assoc.c:291-296 says single-key reads never mutate or build lazily, and src/expr.h:205 says the index is "built at canonicalisation". Both are the opposite of what the code does (src/assoc.c:299-317).
- **Key set ignores `Unset`.** `a[k] =.` has no effect on the association.
- **Order is part of identity.** `<|p->1,q->2|> === <|q->2,p->1|>` is `False`, and the same comparison with `==` stays unevaluated.

### Key findings

- [strength] A hybrid design: a plain `Association[Rule...]` Expr plus a lazily built `AssocIndex` stored in unused union space. The generic expression toolchain works unchanged, and single-key reads are O(1) (src/expr.h:211, src/assoc.c:297-317).
- [strength] Hash-driven canonicalisation keeps first-occurrence order and last-value-wins in O(n). The refcount-aliasing bug on update is fixed and confirmed at runtime (src/assoc.c:133-168, src/part.c:99-115).
- [risk] Mutation is copy-everything. `a[k]=v` is O(n) and a linear scan, and filling an association incrementally is O(n²). Measured: 1.64 s for 8,000 inserts into an association that grows toward 8,000 keys.
- [risk] `RuleDelayed` survives only on the no-change fast path. Any duplicate or splice turns it into an evaluated `Rule` (src/assoc.c:151, 156).
- [risk] The index is correct only because nothing mutates in place, a rule no code enforces. The lazy build writes to a `const`, possibly shared node, and the malformed-entry checks differ between the index and the scan (src/assoc_index.c:44 vs src/assoc.c:286).
- [gap] Associations are not atomic and not `HoldAllComplete` (only `Protected`), so values evaluate and `ReplaceAll` rewrites keys.
- [gap] `a[k] =.` does not delete the key, and `KeyDropFrom` is missing.
- [gap] Header and comment documentation disagrees with the implemented index lifecycle (src/expr.h:205, src/assoc.c:291-296).
