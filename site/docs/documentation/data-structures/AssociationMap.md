# AssociationMap

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`AssociationMap[f, {k1, k2, ...}]`**

Gives \<|k1 -\> f\[k1\], k2 -\> f\[k2\], ...|\>.

**`AssociationMap[f, assoc]`**

Applies f to each rule k -\> v of assoc; the results (rules, lists of rules or associations) form the new association.

## Examples (2)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (2)

```mathematica
In[1]:= AssociationMap[#^2 &, {1, 2, 3, 4}]
Out[1]= <|1 -> 1, 2 -> 4, 3 -> 9, 4 -> 16|>

In[2]:= AssociationMap[Reverse, <|"a" -> 1, "b" -> 2|>]
Out[2]= <|1 -> "a", 2 -> "b"|>
```

## Implementation notes

**Attributes:** `Protected`.

## References

**See also:** [Nothing](../../lists-and-iteration/Nothing/)

- Source: [`src/assoc_ops.c`](https://github.com/stblake/mathilda/blob/main/src/assoc_ops.c)
- Specification: [`docs/spec/builtins/data-structures.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/data-structures.md)
- Tests: [`tests/test_assoc_forms.c`](https://github.com/stblake/mathilda/blob/main/tests/test_assoc_forms.c)
- Tests: [`tests/test_assoc_ops.c`](https://github.com/stblake/mathilda/blob/main/tests/test_assoc_ops.c)
- Tests: [`tests/test_association.c`](https://github.com/stblake/mathilda/blob/main/tests/test_association.c)
