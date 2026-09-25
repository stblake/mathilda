# KeyComplement

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`KeyComplement[{assoc1, assoc2, ...}]`**

Gives the entries of assoc1 whose keys occur in none of the other associations. Elements may also be rules or lists of rules.

## Examples (1)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (1)

```mathematica
In[1]:= KeyComplement[{<|"a" -> 1, "b" -> 2, "c" -> 3|>, <|"b" -> 0|>, <|"c" -> 0|>}]
Out[1]= <|"a" -> 1|>
```

## Implementation notes

**Attributes:** `Protected`.

## References

- Source: [`src/assoc_ops.c`](https://github.com/stblake/mathilda/blob/main/src/assoc_ops.c)
- Specification: [`docs/spec/builtins/data-structures.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/data-structures.md)
- Tests: [`tests/test_assoc_ops.c`](https://github.com/stblake/mathilda/blob/main/tests/test_assoc_ops.c)
