# KeyValueMap

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`KeyValueMap[f, assoc]`**

Gives {f\[k1, v1\], f\[k2, v2\], ...}.

## Examples (2)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (2)

```mathematica
In[1]:= KeyValueMap[f, <|a -> 1, b -> 2|>]
Out[1]= {f[a, 1], f[b, 2]}

In[2]:= KeyValueMap[Plus, <|1 -> 10, 2 -> 20|>]
Out[2]= {11, 22}
```

## Implementation notes

**Attributes:** `Protected`.

## References

- Source: [`src/assoc.c`](https://github.com/stblake/mathilda/blob/main/src/assoc.c)
- Specification: [`docs/spec/builtins/data-structures.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/data-structures.md)
- Tests: [`tests/test_association.c`](https://github.com/stblake/mathilda/blob/main/tests/test_association.c)
