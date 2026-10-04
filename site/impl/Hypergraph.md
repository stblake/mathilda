---
references:
  - "C. Berge, *Hypergraphs: Combinatorics of Finite Sets*, North-Holland Mathematical Library 45 (Elsevier, 1989)."
source: src/graph/hyp_util.c
---
**Algorithm.** `builtin_hypergraph` is a canonicalising constructor with three
accepted forms. `Hypergraph[{e1, ...}]` derives the vertices in first-appearance
order; `Hypergraph[{v1, ...}, {e1, ...}]` takes them explicitly (duplicates
dropped, first occurrence kept); `Hypergraph[g]` turns each edge `u<->v`/`u->v`
of a `Graph` into the hyperedge `{u, v}`. A hyperedge must be a `List`; anything
else leaves the call unevaluated. A canonical, already-valid `Hypergraph` is its
own fixed point — the builtin returns `NULL` for it. Vertices are resolved by an
**integer fast path** when every element is a machine integer in a range at most
`4n + 1024` wide (`Range[n]`, Wolfram-model states): first appearance by direct
addressing into a `seen` bitmap, which skips the hash index for the dominant
`Σ|e|` element-resolution cost; otherwise through the `GraphVIdx` open-addressing
hash. Construction then seeds the validated-hypergraph memo (`hyp_memo_from`),
resolving every hyperedge element to a vertex index.

**Data structures.** The object is a plain `Expr` tree
`Hypergraph[List, List]` — no new `EXPR_*` tag. The memo (`HYP_MEMO_SLOTS = 4`,
LRU-evicted, keyed by node pointer and holding an `expr_copy` reference so the
node stays alive and immutable) stores the vertex index, the hyperedges as a
vertex-index CSR in two forms — `eoff/ev` (raw, repeats kept) and `soff/sv`
(distinct vertices, aliasing the raw arrays when no hyperedge repeats a vertex)
— and, built lazily on first need, the vertex→hyperedge incidence CSR
(`voff/ve`). `hyp_int_list` returns a packed `int64` buffer above the packing
threshold.

**Complexity / limits.** Construction is `O(Σ|e|)`; validation by every later
head is `O(1)` on a memo hit. The integer fast path roughly halves construction.
Hyperedge counts and total incidence are capped at `INT32_MAX`. The memo is
module-static and lock-free: hypergraph builtins run only on the evaluator
thread.
