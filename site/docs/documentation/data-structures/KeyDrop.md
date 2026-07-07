# KeyDrop

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
KeyDrop[assoc, key]  |  KeyDrop[assoc, {k1, ...}]
    Gives assoc with the specified keys removed (order preserved).
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= KeyDrop[<|"a" -> 1, "b" -> 2, "c" -> 3|>, "b"]
Out[1]= <|"a" -> 1, "c" -> 3|>

In[2]:= KeyDrop[<|"a" -> 1, "b" -> 2, "c" -> 3|>, {"a", "c"}]
Out[2]= <|"b" -> 2|>
```

```mathematica
In[1]:= KeyDrop[{<|"a" -> 1, "b" -> 2|>, <|"a" -> 3, "b" -> 4|>}, "b"]
Out[1]= {<|"a" -> 1|>, <|"a" -> 3|>}
```

## Implementation notes

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/assoc.c`](https://github.com/stblake/mathilda/blob/main/src/assoc.c)
- Specification: [`docs/spec/builtins/data-structures.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/data-structures.md)
