# Delete

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
Delete[expr, n] deletes the element at position n in expr.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= Insert[{a, b, c}, x, 2]
Out[1]= {a, x, b, c}

In[2]:= Delete[{a, b, c}, 2]
Out[2]= {a, c}
```

## Implementation notes

`builtin_delete` (in `src/part.c`) drives the recursive helper `delete_path`, which walks an integer position (or position path) into the expression tree and removes the targeted element by rebuilding the enclosing function with all arguments except that index. Negative indices count from the end, position `0` targets the head (replaced by a `Sequence[...]` of the remaining parts), and out-of-range indices leave the structure unchanged.

- `pos` can be a single index, a list (path), or a list of paths.
- `Delete[expr, 0]` replaces the head with `Sequence`.

**Attributes:** none registered.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/part.c`](https://github.com/stblake/mathilda/blob/main/src/part.c)
- Specification: [`docs/spec/builtins/structural-manipulation.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/structural-manipulation.md)
