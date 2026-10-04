---
source: src/sowreap.c
---
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
