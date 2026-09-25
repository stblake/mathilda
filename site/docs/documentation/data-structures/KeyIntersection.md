# KeyIntersection

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`KeyIntersection[{assoc1, assoc2, ...}]`**

Gives the list of associations restricted to the keys common to all of them, each in the key order of assoc1. Elements may also be rules or lists of rules.

## Examples (2)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (2)

```mathematica
In[1]:= KeyIntersection[{<|"a" -> 1, "b" -> 2|>, <|"b" -> 3, "c" -> 4|>}]
Out[1]= {<|"b" -> 2|>, <|"b" -> 3|>}

In[2]:= KeyIntersection[{<|"a" -> 1, "b" -> 2|>, <|"b" -> 3, "a" -> 4|>}]
Out[2]= {<|"a" -> 1, "b" -> 2|>, <|"a" -> 4, "b" -> 3|>}
```

## Implementation notes

**Attributes:** `Protected`.

## References

- Source: [`src/assoc_ops.c`](https://github.com/stblake/mathilda/blob/main/src/assoc_ops.c)
- Specification: [`docs/spec/builtins/data-structures.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/data-structures.md)
- Tests: [`tests/test_assoc_ops.c`](https://github.com/stblake/mathilda/blob/main/tests/test_assoc_ops.c)
