# Subsets

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
Subsets[list]
    Gives all subsets of list (the power set), ordered by
    increasing length and lexicographically by element position within
    each length. The head of list is kept on the subsets.
Subsets[list, n]
    Gives subsets of length 0 through n.
Subsets[list, {n}]
    Gives subsets of length exactly n.
Subsets[list, {nmin, nmax}]
    Gives subsets whose length lies in the
    inclusive range nmin to nmax; a third element gives a length step.
Subsets[list, spec, s]
    Gives only the first s subsets spec would
    produce, generated lazily.
```

## Examples

_No verified examples yet for this function._

## Implementation notes

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/list/list_init.c`](https://github.com/stblake/mathilda/blob/main/src/list/list_init.c)
- Specification index: [`Mathilda_spec.md`](https://github.com/stblake/mathilda/blob/main/Mathilda_spec.md)
