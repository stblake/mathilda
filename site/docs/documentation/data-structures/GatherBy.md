# GatherBy

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`GatherBy[list, f]`**

Gathers elements with equal f\[element\] into sublists, in first-appearance order: {{group1}, {group2}, ...}.

## Examples (2)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (2)

```mathematica
In[1]:= GatherBy[{1, 2, 3, 4, 5, 6}, EvenQ]
Out[1]= {{1, 3, 5}, {2, 4, 6}}

In[2]:= GatherBy[<|"a" -> 1, "b" -> 2, "c" -> 3, "d" -> 4|>, EvenQ]
Out[2]= {<|"a" -> 1, "c" -> 3|>, <|"b" -> 2, "d" -> 4|>}
```

## Implementation notes

**Attributes:** `Protected`.

## References

**See also:** [GroupBy](../../data-structures/GroupBy/)

- Source: [`src/assoc.c`](https://github.com/stblake/mathilda/blob/main/src/assoc.c)
- Specification: [`docs/spec/builtins/data-structures.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/data-structures.md)
- Tests: [`tests/test_association.c`](https://github.com/stblake/mathilda/blob/main/tests/test_association.c)
- Tests: [`tests/test_list.c`](https://github.com/stblake/mathilda/blob/main/tests/test_list.c)
