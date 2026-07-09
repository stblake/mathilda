# KeyUnion

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
KeyUnion[{assoc1, assoc2, ...}]
    Gives the list of associations padded to the union of all their keys;
    a key absent from an association is filled with Missing["KeyAbsent", key].
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= KeyUnion[{<|"a" -> 1, "b" -> 2|>, <|"b" -> 3, "c" -> 4|>}]
Out[1]= {<|"a" -> 1, "b" -> 2, "c" -> Missing["KeyAbsent", "c"]|>, <|"a" -> Missing["KeyAbsent", "a"], "b" -> 3, "c" -> 4|>}
```

## Implementation notes

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/assoc.c`](https://github.com/stblake/mathilda/blob/main/src/assoc.c)
- Specification: [`docs/spec/builtins/data-structures.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/data-structures.md)
