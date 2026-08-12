# PadRight

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`PadRight[list, n]`**

makes a list of length n by padding list with zeros on the right.

**`PadRight[list, n, x]`**

pads by repeating the element x.

**`PadRight[list, n, {x1, x2, ...}]`**

pads by cyclically repeating the elements xi.

**`PadRight[list, n, padding, m]`**

leaves a margin of m elements of padding on the left.

**`PadRight[list, {n1, n2, ...}]`**

makes a nested list with length ni at level i.

**`PadRight[list]`**

pads a ragged array list with zeros to make it full.

<details>
<summary>Notes</summary>

A negative length pads on the left; a negative margin truncates leading elements. The head of list need not be List.

</details>

## Examples (2)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (2)

```mathematica
In[1]:= PadLeft[{a, b, c}, 10, {x, y, z}]
Out[1]= {z, x, y, z, x, y, z, a, b, c}

In[2]:= PadRight[{{a, b}, {c}}, {3, 5}]
Out[2]= {{a, b, 0, 0, 0}, {c, 0, 0, 0, 0}, {0, 0, 0, 0, 0}}
```

## Implementation notes

**Attributes:** `Protected`.

## See also

[PadLeft](../../structural-manipulation/PadLeft/), [List](../../other-advanced/List/)

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/structural-manipulation.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/structural-manipulation.md)
- Tests: [`tests/test_packed_list.c`](https://github.com/stblake/mathilda/blob/main/tests/test_packed_list.c)
- Tests: [`tests/test_padright.c`](https://github.com/stblake/mathilda/blob/main/tests/test_padright.c)
