# Reap

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Reap[expr]`**

Evaluates expr and returns {value, {sown...}}, collecting every expression sown by Sow during the evaluation. Reap\[expr, patt\] reaps only tags matching patt; Reap\[expr, {p1, ...}\] makes one sublist per pattern; Reap\[expr, patt, f\] returns f\[tag, {e...}\] per tag.

## Examples (4)

Every input below was run against the current Mathilda build and its output recorded.

### Applications (4)

{final value, {sown values in order}}

```mathematica
In[1]:= Reap[Sow[1]; Sow[2]; 99]
Out[1]= {99, {{1, 2}}}
```

Do iterates faithfully, so Sow fires each step

```mathematica
In[2]:= Reap[Do[Sow[i^2], {i, 4}]]
Out[2]= {Null, {{1, 4, 9, 16}}}
```

Grouped by tag, in first-seen order

```mathematica
In[3]:= Reap[Sow[1, x]; Sow[2, y]; Sow[3, x], _]
Out[3]= {3, {{1, 3}, {2}}}
```

Nothing sown -> an empty collection

```mathematica
In[4]:= Reap[42]
Out[4]= {42, {}}
```

## Algorithm

sowreap.c — Sow / Reap: dynamic-scope accumulation of intermediate results.

Mathematica semantics --------------------- Reap[expr] evaluates expr and returns {value, second}, where every value passed to Sow during that evaluation is collected. Sow[e] returns e and (as a side effect) records e in the nearest enclosing Reap that matches its tag. Because a single Reap may capture many Sows, this needs genuine dynamic state — a stack of active Reap "frames" — unlike Catch/Throw, whose first-throw-wins semantics let it work statelessly via an in-band sentinel.

Tags and patterns -----------------

```text
  Sow[e]                 -> tag None
  Sow[e, tag]            -> single tag
  Sow[e, {t1, t2, ...}]  -> one record per element (repeats allowed, so the
                            same value can appear several times)
```

A sown (tag, value) is routed to the innermost frame at least one of whose patterns matches tag; within that frame it is appended to every matching pattern-slot, grouped by structurally-identical tag in first-encounter order.

```text
Reap[expr]            == Reap[expr, _]      (single pattern _)
Reap[expr, patt]      single pattern; `second` is a flat list of entries
Reap[expr, {p1..pk}]  `second` has one slot per pattern (nesting one deeper)
Reap[expr, patt, f]   each entry is f[tag, {values}] instead of {values}
```

Data structures --------------- Per the design, values and tag-groups are singly linked lists of Expr, which give O(1) ordered append and preserve Sow order exactly. Each Reap frame lives on builtin_reap's C stack and is linked into a file-global stack via

```text
`prev`, so push/pop is a pointer write and the frame struct itself needs no
```

heap allocation. The system is single-threaded, so a plain global head is sufficient. Every return path of builtin_reap pops its frame, keeping the stack strictly balanced (LIFO).

## Implementation notes

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

**Attributes:** `HoldFirst`, `Protected`.

## References

- Source: [`src/sowreap.c`](https://github.com/stblake/mathilda/blob/main/src/sowreap.c)
- Specification index: [`Mathilda_spec.md`](https://github.com/stblake/mathilda/blob/main/Mathilda_spec.md)
- Tests: [`tests/test_assoc_atomicity.c`](https://github.com/stblake/mathilda/blob/main/tests/test_assoc_atomicity.c)
- Tests: [`tests/test_sow_reap.c`](https://github.com/stblake/mathilda/blob/main/tests/test_sow_reap.c)

## Notes & additional examples

### Notes

`Reap[expr]` evaluates `expr` and returns `{value, {sown...}}`, gathering everything `Sow`
deposits during the evaluation; `Reap[expr, patt]` keeps only matching tags, and the
third-argument form applies `f[tag, {e...}]` to each tag group. Values come back in exact
`Sow` order, grouped by tag.

One subtlety: Mathilda's `Sum` evaluates its summand **symbolically once** and closed-forms
the answer, so a `Sow` inside `Sum` fires a single time, not once per index. Use `Do`,
`Table` or an explicit loop when you want one sown value per iteration — hence the `Do`
example rather than a `Sum`.
