# Subtract

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
x - y or Subtract[x, y] represents x - y; rewritten by the evaluator
to Plus[x, Times[-1, y]] so it inherits Plus's flattening and ordering.
```

## Examples

_No verified examples yet for this function._

## Implementation notes

`builtin_subtract` is a thin two-argument rewrite: `a - b` becomes `Plus[a, Times[-1, b]]`. It does no arithmetic itself — the returned `Plus`/`Times` tree is canonicalised and folded by the evaluator's `Plus`/`Times` machinery. Non-binary calls return `NULL`.

**Attributes:** `Listable`, `NumericFunction`, `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/arithmetic.c`](https://github.com/stblake/mathilda/blob/main/src/arithmetic.c)
- Specification index: [`Mathilda_spec.md`](https://github.com/stblake/mathilda/blob/main/Mathilda_spec.md)
