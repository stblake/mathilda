# First

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
First[expr] gives the first element of expr.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= First[<|"a" -> 10, "b" -> 20|>]
Out[1]= 10

In[2]:= Rest[<|"a" -> 10, "b" -> 20, "c" -> 30|>]
Out[2]= <|"b" -> 20, "c" -> 30|>

In[3]:= Take[<|"a" -> 1, "b" -> 2, "c" -> 3|>, 2]
Out[3]= <|"a" -> 1, "b" -> 2|>
```

## Implementation notes

`builtin_first` (in `src/part.c`) takes a single argument and returns a deep copy of its first element (`args[0]`). It returns `NULL` (unevaluated) when the argument is atomic or has no elements.

**Attributes:** none registered.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/part.c`](https://github.com/stblake/mathilda/blob/main/src/part.c)
- Specification: [`docs/spec/builtins/data-structures.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/data-structures.md)

## Notes & additional examples

### Worked examples

```mathematica
In[1]:= First[{a,b,c}]
Out[1]= a
```

### Notes

`First[expr]` returns the first element (part 1) of any expression, not only lists.
