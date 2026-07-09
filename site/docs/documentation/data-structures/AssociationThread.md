# AssociationThread

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
AssociationThread[{k...}, {v...}]  |  AssociationThread[keys -> values]
    Builds <|k1 -> v1, ...|> from parallel key and value lists.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= AssociationThread[{"a", "b"}, {1, 2}]
Out[1]= <|"a" -> 1, "b" -> 2|>

In[2]:= AssociationThread[{"a", "b"} -> {1, 2}]
Out[2]= <|"a" -> 1, "b" -> 2|>
```

## Implementation notes

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/assoc.c`](https://github.com/stblake/mathilda/blob/main/src/assoc.c)
- Specification: [`docs/spec/builtins/data-structures.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/data-structures.md)
