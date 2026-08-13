# AssociationMap

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`AssociationMap[f, {k1, k2, ...}]`**

Gives \<|k1 -\> f\[k1\], k2 -\> f\[k2\], ...|\>.

## Examples (1)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (1)

```mathematica
In[1]:= AssociationMap[#^2 &, {1, 2, 3, 4}]
Out[1]= <|1 -> 1, 2 -> 4, 3 -> 9, 4 -> 16|>
```

## Implementation notes

**Attributes:** `Protected`.

## References

- Source: [`src/assoc.c`](https://github.com/stblake/mathilda/blob/main/src/assoc.c)
- Specification: [`docs/spec/builtins/data-structures.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/data-structures.md)
- Tests: [`tests/test_association.c`](https://github.com/stblake/mathilda/blob/main/tests/test_association.c)
