# Lookup

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
Lookup[assoc, key]
    Gives the value for key, or Missing["KeyAbsent", key].
Lookup[assoc, key, default]
    Uses default when key is absent.
Lookup[assoc, {k1, k2, ...}]
    Looks up several keys at once (O(n+m)).
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= Lookup[<|"a" -> 1, "b" -> 2|>, "b"]
Out[1]= 2

In[2]:= Lookup[<|"a" -> 1|>, "z", 0]
Out[2]= 0

In[3]:= Lookup[{<|"a" -> 1, "b" -> 2|>, <|"a" -> 3|>}, "a", 0]
Out[3]= {1, 3}
```

## Implementation notes

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/assoc.c`](https://github.com/stblake/mathilda/blob/main/src/assoc.c)
- Specification: [`docs/spec/builtins/data-structures.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/data-structures.md)
