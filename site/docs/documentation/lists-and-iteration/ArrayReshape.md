# ArrayReshape

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`ArrayReshape[list, dims]`**

arranges the flattened elements of list into a rectangular array of dimensions dims, dropping extra elements or padding with 0 as needed.

**`ArrayReshape[list, dims, padding]`**

uses the given padding scheme (as in ArrayPad) when list is too short.

## Examples (2)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (2)

```mathematica
In[1]:= ArrayReshape[{a, b, c, d, e, f}, {2, 3}]
Out[1]= {{a, b, c}, {d, e, f}}

In[2]:= ArrayReshape[{1, 2, 3, 4, 5, 6, 7}, {5, 3}, x]
Out[2]= {{1, 2, 3}, {4, 5, 6}, {7, x, x}, {x, x, x}, {x, x, x}}
```

## Implementation notes

**Attributes:** `Protected`.

## References

**See also:** [List](../../other-advanced/List/), [ArrayPad](../../lists-and-iteration/ArrayPad/), [NDArray](../../linear-algebra/NDArray/)

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/lists-and-iteration.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/lists-and-iteration.md)
- Tests: [`tests/test_array_reshape.c`](https://github.com/stblake/mathilda/blob/main/tests/test_array_reshape.c)
