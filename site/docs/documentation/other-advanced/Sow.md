# Sow

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Sow[e]`**

Collects e into the nearest enclosing Reap and returns e. Sow\[e, tag\] collects e for the Reap whose pattern matches tag; Sow\[e, {tag1, tag2, ...}\] collects e once for each tag\_i.

## Examples (4)

Every input below was run against the current Mathilda build and its output recorded.

### Applications (4)

{lastValue, {{everything sown}}}

```mathematica
In[1]:= Reap[Sow[1]; Sow[2]; Sow[3]]
Out[1]= {3, {{1, 2, 3}}}
```

Collect squares across a loop

```mathematica
In[2]:= Reap[Do[Sow[i^2], {i, 1, 4}]]
Out[2]= {Null, {{1, 4, 9, 16}}}
```

Grouped by tag in first-encounter order

```mathematica
In[3]:= Reap[Sow[1, "a"]; Sow[2, "b"]; Sow[3, "a"], _]
Out[3]= {3, {{1, 3}, {2}}}
```

Returns its value; collection is a side effect

```mathematica
In[4]:= Sow[42]
Out[4]= 42
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

**Algorithm.** `builtin_sow` records its (already evaluated) first argument into the
nearest enclosing `Reap` whose pattern matches the tag, and returns a copy of the value
— `Sow` is transparent in an expression, its collection a side effect. The tag is set
by the arity: `Sow[e]` uses the tag `None`; `Sow[e, tag]` uses that single tag;
`Sow[e, {t1, t2, ...}]` records `e` once per element (repeats allowed, so the same
value can be routed to several tags). Each `(tag, value)` is handed to `sow_route`,
which walks the active `Reap` frames from innermost outward and deposits the value in
the first frame at least one of whose patterns matches the tag, appended to every
matching pattern-slot and grouped by structurally-identical tag in first-encounter
order.

**Data structures.** A file-global LIFO stack of `ReapFrame`s (each owned by a live
`builtin_reap` C-stack call, linked by a `prev` pointer — no heap frame needed). Within
a frame, values are singly linked `SowNode` lists grouped into `SowGroup`s per tag, with
an `O(1)` tail append that preserves `Sow` order exactly. The system is single-threaded,
so one global head suffices. `Sow` outside any `Reap` has no frame to route to, so it
simply returns its value.

**Complexity / limits.** Each `Sow` costs the match against the active frames' patterns
plus an `O(1)` append. Attributes: `Protected` (note: `Sow` is *not* `HoldFirst` — its
argument is evaluated normally; only `Reap` holds its first argument). The counterpart
`Reap[expr]` returns `{value, collected}` and empties the frame on every exit path,
keeping the stack strictly balanced.

**Attributes:** `Protected`.

## References

- Source: [`src/sowreap.c`](https://github.com/stblake/mathilda/blob/main/src/sowreap.c)
- Specification index: [`Mathilda_spec.md`](https://github.com/stblake/mathilda/blob/main/Mathilda_spec.md)
- Tests: [`tests/test_sow_reap.c`](https://github.com/stblake/mathilda/blob/main/tests/test_sow_reap.c)

## Notes & additional examples

### Notes

`Sow[e]` records `e` into the nearest enclosing `Reap` and returns `e` unchanged — it is
transparent in an expression, so you can drop it into existing code to harvest
intermediate values. `Sow[e, tag]` routes `e` to the `Reap` whose pattern matches `tag`;
`Sow[e, {t1, t2, ...}]` records `e` once per tag.

`Reap[expr]` returns `{value, collected}` where `value` is the result of `expr` and
`collected` holds everything sown, grouped by tag in the order each tag was first seen.
Outside any `Reap`, `Sow` simply returns its argument (nothing is collected). Note that
`Sow` evaluates its argument normally — it is `Reap` that is `HoldFirst`. `Sow` is
`Protected`.
