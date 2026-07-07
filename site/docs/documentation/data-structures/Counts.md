# Counts

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
Counts[list]
    Gives <|element -> count, ...|> tallying each distinct
    element. Hash-indexed: one O(n) pass.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= Counts[{1, 2, 2, 3, 3, 3}]
Out[1]= <|1 -> 1, 2 -> 2, 3 -> 3|>

In[2]:= Counts[<|"a" -> 1, "b" -> 1, "c" -> 2|>]
Out[2]= <|1 -> 2, 2 -> 1|>
```

## Implementation notes

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/assoc.c`](https://github.com/stblake/mathilda/blob/main/src/assoc.c)
- Specification: [`docs/spec/builtins/data-structures.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/data-structures.md)
