# FlattenAt

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`FlattenAt[list, n]`**

flattens out the sublist at position n of list, splicing its elements into list; a negative n counts from the end.

**`FlattenAt[expr, {i, j, ...}]`**

flattens out the part of expr at the position {i, j, ...}.

**`FlattenAt[expr, {{i1, ...}, {i2, ...}, ...}]`**

flattens out the parts of expr at several positions. The head of the spliced part is removed; FlattenAt works on any head, not just List.

## Examples (2)

Every input below was run against the current Mathilda build and its output recorded.

### Worked examples (2)

```mathematica
In[1]:= FlattenAt[{a,{b,c},{d,e},{f}}, 2]
Out[1]= {a, b, c, {d, e}, {f}}

In[2]:= FlattenAt[{a,{b,c},{d,e},{f}}, {{2},{4}}]
Out[2]= {a, b, c, {d, e}, f}
```

## Implementation notes

**Attributes:** `Protected`.

## References

**See also:** [List](../../other-advanced/List/), [MapAt](../../data-structures/MapAt/), [ReplaceAt](../../assignment-and-rules/ReplaceAt/)

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/structural-manipulation.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/structural-manipulation.md)
- Tests: [`tests/test_flatten_at.c`](https://github.com/stblake/mathilda/blob/main/tests/test_flatten_at.c)
