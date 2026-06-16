# List

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
{e1, e2, ...} or List[e1, e2, ...]
    represents an ordered list of the elements ei.
List is the fundamental container head: vectors are lists, matrices are
lists of lists, and the structural operators (Part, Map, Take, Drop,
Length, ...) act on List. Elements are evaluated normally and kept in the
given order (List has no Orderless attribute). The parser writes the {...}
syntax to List, and the printer renders List[...] back as {...}.
```

## Examples

_No verified examples yet for this function._

## Implementation notes

**Attributes:** none registered.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification index: [`Mathilda_spec.md`](https://github.com/stblake/mathilda/blob/main/Mathilda_spec.md)
