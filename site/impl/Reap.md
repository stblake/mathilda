---
source: src/sowreap.c
---
**Algorithm.** `Reap[expr]` evaluates `expr` and returns `{value, {sown...}}`, collecting
every expression `Sow` deposits during that evaluation. `builtin_reap` (`src/sowreap.c`,
registered inline in `core_init`) mirrors `Catch`'s shell: it is `HoldFirst`, so it pushes a
fresh `ReapFrame` onto the file-global `g_reap_top` stack, calls `evaluate` on the held body
(the body is *borrowed*, not copied), pops the frame, and then — if the result is an
in-flight `Throw` (`eval_is_inflight_throw`) — frees the collectors and propagates. `Sow`
routes each sown value to the innermost frame whose tag patterns match, grouping values by
`expr_eq` on the tag. `Reap[expr, patt]` reaps only matching tags; `Reap[expr, {p1, ...}]`
makes one sublist per pattern; `Reap[expr, patt, f]` returns `f[tag, {e...}]` per tag group.

**Data structures.** Each `ReapFrame` lives on `builtin_reap`'s own C stack, linked through
a `prev` pointer in balanced LIFO order (no heap allocation for the frame). Per-tag values
and the tag groups are singly linked lists of `Expr`, giving `O(1)` ordered append that
preserves exact `Sow` order. The output list is assembled by *moving* the value pointers out
of those lists (nulling each node) and then freeing the frame — no deep copy.

**Complexity / limits.** `O(n)` in the number of sown values. One gotcha follows from
Mathilda's symbolic `Sum`: `Sum` evaluates its summand symbolically **once** and
closed-forms the result, so a `Sow` inside `Sum` fires once, not once per index. To reap one
value per iteration, use `Do`, `Table` or an explicit loop, which iterate faithfully — e.g.
`Reap[Do[Sow[i^2], {i, 4}]]`.
