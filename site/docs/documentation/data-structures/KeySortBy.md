# KeySortBy

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
KeySortBy[assoc, f]
    Sorts an association by f applied to each key (stable).
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= KeySortBy[<|"bbb" -> 1, "a" -> 2, "cc" -> 3|>, StringLength]
Out[1]= <|"a" -> 2, "cc" -> 3, "bbb" -> 1|>
```

## Implementation notes

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/assoc.c`](https://github.com/stblake/mathilda/blob/main/src/assoc.c)
- Specification: [`docs/spec/builtins/data-structures.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/data-structures.md)
