# Flatten

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
Flatten[list]
    flattens out nested lists, collapsing every level into a flat list
    with the same head as the top level.
Flatten[list, n]
    flattens only the top n levels.
Flatten[list, n, h]
    flattens only sublists whose head matches h, leaving other heads
    in place.
```

## Examples

_No verified examples yet for this function._

## Implementation notes

**Algorithm.** `builtin_flatten` (in `src/list.c`) accepts `Flatten[list]`, `Flatten[list, n]` (level cap), and `Flatten[list, n, h]` (custom head). It iterates over the top-level arguments calling the recursive worker `flatten_rec`, which splices the children of any subexpression whose head equals the flattening head `h` (default `List`) up into the output, descending up to `n` levels (n = -1 means unlimited). The collected arguments are gathered into a growable buffer and reassembled under the original head.

**Data structures.** A dynamically grown `Expr**` accumulator (`results`, with `count`/`cap`) holds the deep-copied leaf expressions before the final `expr_new_function` rebuild.

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/list.c`](https://github.com/stblake/mathilda/blob/main/src/list.c)
- Specification: [`docs/spec/builtins/structural-manipulation.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/structural-manipulation.md)
