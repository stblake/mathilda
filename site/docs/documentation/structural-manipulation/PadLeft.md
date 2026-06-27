# PadLeft

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
PadLeft[list, n]
    makes a list of length n by padding list with zeros on the left.
PadLeft[list, n, x]
    pads by repeating the element x.
PadLeft[list, n, {x1, x2, ...}]
    pads by cyclically repeating the elements xi.
PadLeft[list, n, padding, m]
    leaves a margin of m elements of padding on the right.
PadLeft[list, {n1, n2, ...}]
    makes a nested list with length ni at level i.
PadLeft[list]
    pads a ragged array list with zeros to make it full.
A negative length pads on the right; a negative margin truncates trailing elements. The head of list need not be List.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= PadLeft[{a, b, c}, 10, {x, y, z}]
Out[1]= {z, x, y, z, x, y, z, a, b, c}

In[2]:= PadRight[{{a, b}, {c}}, {3, 5}]
Out[2]= {{a, b, 0, 0, 0}, {c, 0, 0, 0, 0}, {0, 0, 0, 0, 0}}
```

## Implementation notes

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/structural-manipulation.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/structural-manipulation.md)
