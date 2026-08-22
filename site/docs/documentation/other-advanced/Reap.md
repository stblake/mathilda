# Reap

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Reap[expr]`**

Evaluates expr and returns {value, {sown...}}, collecting every expression sown by Sow during the evaluation. Reap\[expr, patt\] reaps only tags matching patt; Reap\[expr, {p1, ...}\] makes one sublist per pattern; Reap\[expr, patt, f\] returns f\[tag, {e...}\] per tag.

## Examples

_No verified examples yet for this function._

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

**Attributes:** `HoldFirst`, `Protected`.

## References

- Source: [`src/core.c`](https://github.com/stblake/mathilda/blob/main/src/core.c)
- Specification index: [`Mathilda_spec.md`](https://github.com/stblake/mathilda/blob/main/Mathilda_spec.md)
- Tests: [`tests/test_sow_reap.c`](https://github.com/stblake/mathilda/blob/main/tests/test_sow_reap.c)
