# KeySort

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
KeySort[assoc]
    Sorts an association into canonical key order.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= KeySort[<|"c" -> 3, "a" -> 1, "b" -> 2|>]
Out[1]= <|"a" -> 1, "b" -> 2, "c" -> 3|>
```

## Implementation notes

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/assoc.c`](https://github.com/stblake/mathilda/blob/main/src/assoc.c)
- Specification: [`docs/spec/builtins/data-structures.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/data-structures.md)
