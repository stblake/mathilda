---
source: src/parse.c
---
**Definition.** `List[e1, e2, ...]`, written `{e1, e2, ...}`, is the fundamental
ordered-container head. It is pure structure: there is no `List` builtin — the evaluator
simply evaluates each element and keeps it. Vectors are lists, matrices are lists of
lists, and the structural operators (`Part`, `Map`, `Take`, `Drop`, `Length`, ...) are
defined to act on `List`.

**Representation.** A bare `EXPR_SYMBOL` head (interned `SYM_List`) on an
`EXPR_FUNCTION` node. The parser (`src/parse.c`) rewrites the `{...}` surface syntax to a
`List[...]` node, and the printer (`src/print.c`) renders `List[...]` back as `{...}`;
`FullForm[{a, b, c}]` exposes the underlying `List[a, b, c]`. The packed-array
substrate (`src/pack.c`) can store a numeric `List` as a dense machine buffer
transparently, but the logical head is always `List`.

**Usage & limits.** Elements are evaluated normally and kept in the given order — `List`
has **no** `Orderless` attribute, so `{3, 1, 2}` is not sorted. In fact `Attributes[List]`
is empty (`{}`): it carries no attributes at all. Its uniformity is the point —
everything is an expression, so one set of generic tools works on every list.
