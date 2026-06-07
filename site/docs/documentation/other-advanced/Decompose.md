# Decompose

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
Decompose[poly, x]
    decomposes the univariate polynomial poly into the deepest possible
    composition {p1, p2, ..., pk} such that poly == p1[p2[...[pk[x]]...]],
    with each pi a polynomial of degree >= 2 in x.
    Returns {poly} if no nontrivial decomposition exists.
```

## Examples

_No verified examples yet for this function._

## Implementation notes

**Attributes:** `Listable`, `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification index: [`Mathilda_spec.md`](https://github.com/stblake/mathilda/blob/main/Mathilda_spec.md)
