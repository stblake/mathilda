# KeyTake

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
KeyTake[assoc, {k1, ...}]
    Gives the association of only the specified keys (order preserved).
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= KeyTake[<|"a" -> 1, "b" -> 2, "c" -> 3|>, {"c", "a"}]
Out[1]= <|"a" -> 1, "c" -> 3|>

In[2]:= KeyTake[{<|"a" -> 1, "b" -> 2|>, <|"a" -> 3, "b" -> 4|>}, {"a"}]
Out[2]= {<|"a" -> 1|>, <|"a" -> 3|>}
```

## Implementation notes

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/assoc.c`](https://github.com/stblake/mathilda/blob/main/src/assoc.c)
- Specification: [`docs/spec/builtins/data-structures.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/data-structures.md)
