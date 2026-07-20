# ConstantArray

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
ConstantArray[c, n]
    generates a list of n copies of the element c.
ConstantArray[c, {n1, n2, ...}]
    generates an n1 x n2 x ... nested array of copies of c.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= ConstantArray[c, 5]
Out[1]= {c, c, c, c, c}

In[2]:= ConstantArray[0, {2, 3}]
Out[2]= {{0, 0, 0}, {0, 0, 0}}
```

## Implementation notes

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/lists-and-iteration.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/lists-and-iteration.md)
