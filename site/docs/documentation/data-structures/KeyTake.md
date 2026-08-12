# KeyTake

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`KeyTake[assoc, {k1, ...}]`**

Gives the association of only the specified keys (order preserved).

## Examples (2)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (2)

```mathematica
In[1]:= KeyTake[<|"a" -> 1, "b" -> 2, "c" -> 3|>, {"c", "a"}]
Out[1]= <|"a" -> 1, "c" -> 3|>

In[2]:= KeyTake[{<|"a" -> 1, "b" -> 2|>, <|"a" -> 3, "b" -> 4|>}, {"a"}]
Out[2]= {<|"a" -> 1|>, <|"a" -> 3|>}
```

## Implementation notes

**Attributes:** `Protected`.

## References

- Source: [`src/assoc.c`](https://github.com/stblake/mathilda/blob/main/src/assoc.c)
- Specification: [`docs/spec/builtins/data-structures.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/data-structures.md)
- Tests: [`tests/test_association.c`](https://github.com/stblake/mathilda/blob/main/tests/test_association.c)
- Tests: [`tests/test_compile_assoc.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compile_assoc.c)
