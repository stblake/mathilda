# KeySelect

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
KeySelect[assoc, pred]
    Keeps the entries whose key satisfies pred.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= KeySelect[<|1 -> 10, 2 -> 20, 3 -> 30|>, EvenQ]
Out[1]= <|2 -> 20|>
```

## Implementation notes

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/assoc.c`](https://github.com/stblake/mathilda/blob/main/src/assoc.c)
- Specification: [`docs/spec/builtins/data-structures.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/data-structures.md)
