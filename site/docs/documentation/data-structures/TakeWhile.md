# TakeWhile

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
TakeWhile[list, crit]
    Gives the longest leading run of elements e for
    which crit[e] is True. Over an association, tests the values and keeps
    the matching leading entries (keys preserved).
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= TakeWhile[<|"a" -> 1, "b" -> 2, "c" -> 5, "d" -> 1|>, # < 3 &]
Out[1]= <|"a" -> 1, "b" -> 2|>

In[2]:= LengthWhile[<|"a" -> 1, "b" -> 2, "c" -> 5|>, # < 3 &]
Out[2]= 2
```

## Implementation notes

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/core.c`](https://github.com/stblake/mathilda/blob/main/src/core.c)
- Specification: [`docs/spec/builtins/data-structures.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/data-structures.md)
