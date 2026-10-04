---
source: src/list/splitby.c
---
**Algorithm.** `SplitBy[list, f]` (`src/list/splitby.c`) walks `list` once, evaluating the
key `f[element]` for each element and starting a new run whenever the key differs from the
previous element's key. Only **adjacent** elements are ever grouped — this is the difference
from `GatherBy`, which collects equal keys from anywhere in the list, and from `Split`, which
compares adjacent elements (or a two-argument test) directly rather than through an evaluated
key. `SplitBy[list, {f1, f2, ...}]` splits by `f1`, then recursively splits each resulting
run by `f2`, nesting one level deeper per function.

**Data structures.** Ordinary `Expr` trees. The walk keeps the current run as a growing
`Expr**` vector and the previous key for the adjacency comparison; each completed run is
copied out under a fresh copy of the list's head (`splitby_make_run`). Keys are compared
structurally, so two adjacent elements whose keys evaluate to the same unevaluated form still
group together. The list's own head is preserved on the outer result and on each run.

**Complexity / limits.** `O(n)` key evaluations and comparisons for a single function; the
`{f1, ..., fk}` form multiplies by the nesting depth. `SplitBy` preserves element order
within every run (it never sorts), so the runs read left-to-right exactly as the input does.
