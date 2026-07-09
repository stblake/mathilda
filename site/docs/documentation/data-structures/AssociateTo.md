# AssociateTo

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
AssociateTo[s, key -> val]  |  AssociateTo[s, {rules}]
    Adds or updates key-value pairs in the association held by symbol s,
    modifying s in place.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= asc = <|"a" -> 1|>; AssociateTo[asc, "b" -> 2]; asc
Out[1]= <|"a" -> 1, "b" -> 2|>
```

## Implementation notes

**Attributes:** `HoldFirst`, `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/assoc.c`](https://github.com/stblake/mathilda/blob/main/src/assoc.c)
- Specification: [`docs/spec/builtins/data-structures.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/data-structures.md)
