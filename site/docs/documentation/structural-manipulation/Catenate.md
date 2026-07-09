# Catenate

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
Catenate[{e1, e2, ...}]
    Concatenates the ei (which must share a head)
    into one, flattening a single level. A list of associations merges into
    one association (later keys win).
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= Catenate[{{1, 2}, {3, 4}, {5}}]
Out[1]= {1, 2, 3, 4, 5}

In[2]:= Catenate[{<|"a" -> 1, "b" -> 2|>, <|"b" -> 3, "c" -> 4|>}]
Out[2]= <|"a" -> 1, "b" -> 3, "c" -> 4|>

In[3]:= Catenate[Values[GroupBy[Range[6], EvenQ]]]
Out[3]= {1, 3, 5, 2, 4, 6}
```

## Implementation notes

- `Protected`.
- All elements must share a head; returns unevaluated for mixed heads.
- `Catenate[{}]` gives `{}`.

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/list/list_init.c`](https://github.com/stblake/mathilda/blob/main/src/list/list_init.c)
- Specification: [`docs/spec/builtins/structural-manipulation.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/structural-manipulation.md)
