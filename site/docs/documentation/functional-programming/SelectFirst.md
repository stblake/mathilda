# SelectFirst

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
SelectFirst[list, pred]
    Gives the first element e of list for which
    pred[e] is True, or Missing["NotFound"]. SelectFirst[list, pred, default]
    uses default. Over an association, tests values and returns the first match.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= SelectFirst[{1, 3, 4, 5, 6}, EvenQ]
Out[1]= 4

In[2]:= FirstCase[{1, 2, 3, 4}, _?EvenQ]
Out[2]= 2

In[3]:= SelectFirst[{1, 3, 5}, EvenQ, None]
Out[3]= None
```

## Implementation notes

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/core.c`](https://github.com/stblake/mathilda/blob/main/src/core.c)
- Specification: [`docs/spec/builtins/functional-programming.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/functional-programming.md)
