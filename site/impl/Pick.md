---
source: src/list/pick.c
---
**Algorithm.** `Pick[expr, sel]` / `Pick[expr, sel, patt]` keeps the elements of
`expr` whose positionally-corresponding element of the selector `sel` matches `patt`
(literal `True` for the two-argument form). The selector mirrors `expr`'s structure,
so `pick_rec` descends both trees together: at each element, a selector that matches
`patt` keeps the whole `expr` element; a *compound* selector that did not match as a
whole recurses into the pair; an *atomic* non-matching selector drops the element.
The head at every level comes from `expr`, never `sel`. An association picks by
position over its entries — the value is compared and the entry (key included) kept;
an association used *as* a selector is atomic.

**Data structures.** The matcher (`match` via `pick_selects`) with a throwaway
`MatchEnv`; a per-level `kept` array rebuilt into a node with `expr`'s head.

**Complexity / limits.** O(size) with a match test per element. Any structural
disagreement — differing lengths, or a compound selector against an atomic
expression — is detected exactly (even arbitrarily deep) and aborts the whole call to
`NULL`, so `Pick[…]` is returned unevaluated, matching Mathematica's `Pick::incomp`.
