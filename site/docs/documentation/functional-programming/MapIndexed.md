# MapIndexed

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
MapIndexed[f, list]
    Gives {f[e1, {1}], f[e2, {2}], ...}. Over an
    association, f[value, {Key[k]}] keeping keys.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= MapIndexed[f, {10, 20, 30}]
Out[1]= {f[10, {1}], f[20, {2}], f[30, {3}]}

In[2]:= MapIndexed[f, <|"a" -> 10, "b" -> 20|>]
Out[2]= <|"a" -> f[10, {Key["a"]}], "b" -> f[20, {Key["b"]}]|>
```

## Implementation notes

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/core.c`](https://github.com/stblake/mathilda/blob/main/src/core.c)
- Specification: [`docs/spec/builtins/functional-programming.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/functional-programming.md)
