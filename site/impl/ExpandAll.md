---
source: src/expand.c
---
**Algorithm.** `builtin_expand_all` drives the recursive `expr_expand_all_impl`.
Where `Expand` distributes only at the top level, `ExpandAll` reaches every
subexpression — function heads and arguments, exponents, and the bases of
denominators — expanding bottom-up: it recurses into the head and each argument
first (an `Inequality`'s operator-symbol slots at odd indices are passed through
untouched), rebuilds the node, and then applies the top-level distributor
`expr_expand_impl` at that level. A denominator factor `Power[base, -k]`
(`k > 0`) is handled specially, mirroring `ExpandDenominator`: recurse into the
`base`, expand it to the positive power, and keep the reciprocal, so a
surrounding `Times` later distributes the numerator across the expanded
denominator. `ExpandAll[expr, patt]` leaves any part free of `patt` unchanged
(`expr_contains_patt` guards each descent). Like `Expand`, an expansion too large
to fit in memory returns `Overflow[]` rather than silently declining
(`overflow_mode` is set on the user-facing entry).

**Data structures.** Pure `Expr`-tree recursion; each level allocates a fresh
`Expr**` argument array, rebuilds through `eval_and_free`, and hands the result
to `expr_expand_impl` (which itself uses the FLINT polynomial multiplier for
polynomial-over-`Q` factors). No persistent state beyond the recursion stack.

**Complexity / limits.** Bounded by the size of the fully-expanded tree; the
`Overflow[]` guard uses the same Newton-box size estimate as `Expand`. A symbolic
structural head — no packed/NDArray or `Compile[]` path. Threads over equations,
inequalities, logic functions, and lists; `Protected`.
