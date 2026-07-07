# CountsBy

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
CountsBy[list, f]
    Gives <|f[x] -> count, ...|> tallying elements by f[x].
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= CountsBy[Range[10], EvenQ]
Out[1]= <|False -> 5, True -> 5|>
```

## Implementation notes

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/assoc.c`](https://github.com/stblake/mathilda/blob/main/src/assoc.c)
- Specification: [`docs/spec/builtins/data-structures.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/data-structures.md)
