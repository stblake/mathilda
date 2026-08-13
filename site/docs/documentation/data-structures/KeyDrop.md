# KeyDrop

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`KeyDrop[assoc, key]  |  KeyDrop[assoc, {k1, ...}]`**

Gives assoc with the specified keys removed (order preserved).

## Examples (3)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (2)

```mathematica
In[1]:= KeyDrop[<|"a" -> 1, "b" -> 2, "c" -> 3|>, "b"]
Out[1]= <|"a" -> 1, "c" -> 3|>

In[2]:= KeyDrop[<|"a" -> 1, "b" -> 2, "c" -> 3|>, {"a", "c"}]
Out[2]= <|"b" -> 2|>
```

### Scope (1)

```mathematica
In[3]:= KeyDrop[{<|"a" -> 1, "b" -> 2|>, <|"a" -> 3, "b" -> 4|>}, "b"]
Out[3]= {<|"a" -> 1|>, <|"a" -> 3|>}
```

## Implementation notes

**Attributes:** `Protected`.

## References

**See also:** [KeyTake](../../data-structures/KeyTake/)

- Source: [`src/assoc.c`](https://github.com/stblake/mathilda/blob/main/src/assoc.c)
- Specification: [`docs/spec/builtins/data-structures.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/data-structures.md)
- Tests: [`tests/test_association.c`](https://github.com/stblake/mathilda/blob/main/tests/test_association.c)
- Tests: [`tests/test_compile_assoc.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compile_assoc.c)
