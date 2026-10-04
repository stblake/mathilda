---
references:
  - "E. Tomita, Y. Sutani, T. Higashi, S. Takahashi and M. Wakatsuki, *A simple and faster branch-and-bound algorithm for finding a maximum clique*, WALCOM 2010, LNCS 5942, 191-203."
  - "P. San Segundo, D. Rodríguez-Losada and A. Jiménez, *An exact bit-parallel algorithm for the maximum clique problem*, Comput. Oper. Res. **38** (2011) 571-581."
  - "C. Bron and J. Kerbosch, *Algorithm 457: finding all cliques of an undirected graph*, Comm. ACM **16** (1973) 575-577."
source: src/graph/galg_clique.c
---
**Algorithm.** `builtin_find_clique` first reduces `g` to the simple
mutual-adjacency graph `gc_mutual_graph` — undirected edges, and directed pairs
joined both ways, so a clique in a digraph needs the arcs in both directions.
With no size spec it calls `galg_max_clique`, a bitset branch-and-bound
(Tomita's MCS with San Segundo's BBMC bit-parallel colouring bound): vertices are
numbered by reverse degeneracy order (`gc_degeneracy`, an `O(n+m)` bucket peel),
greedy colouring supplies the pruning bound `depth + colour > incumbent`, and the
incumbent clique is returned sorted. A size spec (`k`, `{k}`, `{kmin,kmax}`, with
an optional count or `All`) switches to pivoted Bron-Kerbosch enumeration
(`gc_bk`), which lists **maximal** cliques in the size window largest-first; so
`FindClique[g, 2]` is `{}` when every maximal clique is larger. For `count == 1`
the enumeration tries sizes downward so the `kmin` prune bites hardest.

**Data structures.** Adjacency is a row-per-vertex bitset (`GcBits`, 64-bit
words), built induced on the working vertex set by `gc_bits_induced` (with an
optional complement for the independent-set sibling). Branch-and-bound keeps a
per-depth `P`-stack of candidate bitsets plus lazily allocated colouring scratch
(`O(ω·n)`, not `O(n²)`); Bron-Kerbosch concatenates each reported clique into a
growable `res`/`rstart` buffer and sorts the runs by `(size, lexicographic)` to
reproduce Wolfram's listing order. Above `GC_GLOBAL_MAX = 3000` vertices the
maximum-clique search is decomposed by degeneracy — every clique lies in
`{v} ∪ N⁺(v)` for its earliest vertex, a small local bitset problem skipped
outright when its core number cannot beat the incumbent.

**Complexity / limits.** Exponential in the worst case but pruned hard in
practice; capped at `GC_MAX_NODES = 5·10⁷` search nodes and `GC_BITSET_MAX = 8192`
vertices for the spec/complement paths, polled against the `TimeConstrained`
deadline every 4096 nodes. An exhausted budget leaves the call **unevaluated** —
never a non-maximum answer. `FindClique[g]` returns `{c}` (a one-element list
holding one maximum clique); the branch-and-bound proves it maximum but, like
Mathematica, does not specify which maximum clique.
