# Most

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
Most[expr] gives all but the last element of expr.
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

**Algorithm.** `builtin_most` returns a copy of the input with its last element dropped: it
copies args `0 .. n−2` into a new function node with the same head. Returns `NULL` (unevaluated)
for atoms or empty expressions.

**Attributes:** none registered.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/part.c`](https://github.com/stblake/mathilda/blob/main/src/part.c)
- Specification: [`docs/spec/builtins/data-structures.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/data-structures.md)

## Notes & additional examples

### Worked examples

```mathematica
In[1]:= Most[{a,b,c,d}]
Out[1]= {a, b, c}
```

### Notes

`Most[expr]` drops the last element; it is the complement of `Last`.
