# SubsetQ

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`SubsetQ[a, b]`**

Gives True if every element of b occurs in a (multiplicity ignored). Lists and associations (compared by value) may be mixed; other expressions must share a head.

## Examples (3)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (3)

```mathematica
In[1]:= SubsetQ[{1, 2, 3}, {3, 1}]
Out[1]= True

In[2]:= SubsetQ[{1, 2}, {1, 4}]
Out[2]= False

In[3]:= SubsetQ[<|"a" -> 1, "b" -> 2|>, {2}]
Out[3]= True
```

## Implementation notes

**Attributes:** `Protected`.

## References

- Source: [`src/assoc_ops.c`](https://github.com/stblake/mathilda/blob/main/src/assoc_ops.c)
- Specification: [`docs/spec/builtins/data-structures.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/data-structures.md)
- Tests: [`tests/test_assoc_ops.c`](https://github.com/stblake/mathilda/blob/main/tests/test_assoc_ops.c)
