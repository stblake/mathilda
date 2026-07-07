# AllTrue

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
AllTrue[list, test]
    Gives True if test[e] is True for every element e
    (True for an empty list). Over an association, tests the values.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= AllTrue[{2, 4, 6}, EvenQ]
Out[1]= True

In[2]:= AnyTrue[{1, 3, 4}, EvenQ]
Out[2]= True

In[3]:= NoneTrue[{1, 3, 5}, EvenQ]
Out[3]= True
```

## Implementation notes

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/core.c`](https://github.com/stblake/mathilda/blob/main/src/core.c)
- Specification: [`docs/spec/builtins/functional-programming.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/functional-programming.md)
