# UnitVector

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
UnitVector[k]
    gives the 2-D unit vector in the k-th direction.
UnitVector[n, k]
    gives the n-D unit vector: a length-n list with a 1 in position k
    and 0s elsewhere.
Components are exact integers unless WorkingPrecision requests
MachinePrecision or a digit count.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= UnitVector[1]
Out[1]= {1, 0}

In[2]:= UnitVector[3, 2]
Out[2]= {0, 1, 0}

In[3]:= UnitVector[2, WorkingPrecision -> MachinePrecision]
Out[3]= {0.0, 1.0}
```

## Implementation notes

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/lists-and-iteration.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/lists-and-iteration.md)
