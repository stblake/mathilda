# CountDistinct

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`CountDistinct[expr]`**

Gives the number of distinct elements of expr (of its values, for an association). One hash pass, O(n).

## Examples (2)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (2)

```mathematica
In[1]:= CountDistinct[{1, 2, 1, 3, 2}]
Out[1]= 3

In[2]:= CountDistinct[<|"a" -> 1, "b" -> 1, "c" -> 2|>]
Out[2]= 2
```

## Implementation notes

**Attributes:** `Protected`.

## References

**See also:** [SameQ](../../comparisons/SameQ/)

- Source: [`src/assoc_ops.c`](https://github.com/stblake/mathilda/blob/main/src/assoc_ops.c)
- Specification: [`docs/spec/builtins/data-structures.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/data-structures.md)
- Tests: [`tests/test_assoc_ops.c`](https://github.com/stblake/mathilda/blob/main/tests/test_assoc_ops.c)
