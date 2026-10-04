---
source: src/list/subsets.c
---
**Algorithm.** `Subsets[list]` (`src/list/subsets.c`) generates all subsets of `list` — the
power set — ordered first by increasing length and then lexicographically by original element
position within each length. The length can be restricted: `Subsets[list, n]` gives lengths
`0` through `n`, `Subsets[list, {n}]` exactly `n`, `Subsets[list, {nmin, nmax}]` the
inclusive range (a third element in the spec gives a length step). `Subsets[list, spec, s]`
returns only the first `s` subsets the spec would produce.

**Data structures.** Ordinary `Expr` trees. Each subset is emitted as a copy of the list's
elements at a chosen set of positions, carrying the list's own head (so `Subsets[f[a, b]]`
keeps head `f`). Subsets of a given length are produced by an index-combination generator
that advances position tuples in lexicographic order.

**Complexity / limits.** The full power set has `2^Length[list]` members, exponential in the
input, so the generator is **lazy**: it produces subsets on demand and the `, s` cap lets a
caller take just the first `s` without materialising the rest — the natural way to peek at
the start of an otherwise intractable enumeration. Order is deterministic (by length, then
lexicographic position), never sorted by value.
