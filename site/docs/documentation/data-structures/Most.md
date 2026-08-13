# Most

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Most[expr] gives all but the last element of expr.`**

## Examples (3)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (2)

```mathematica
In[1]:= Rest[<|"a" -> 10, "b" -> 20, "c" -> 30|>]
Out[1]= <|"b" -> 20, "c" -> 30|>

In[2]:= Take[<|"a" -> 1, "b" -> 2, "c" -> 3|>, 2]
Out[2]= <|"a" -> 1, "b" -> 2|>
```

### Applications (1)

```mathematica
In[3]:= Most[{a,b,c,d}]
Out[3]= {a, b, c}
```

## Implementation notes

**Algorithm.** `builtin_most` returns a copy of the input with its last element dropped: it
copies args `0 .. n−2` into a new function node with the same head. Returns `NULL` (unevaluated)
for atoms or empty expressions.

**Attributes:** none registered.

## References

**See also:** [First](../../data-structures/First/), [Last](../../data-structures/Last/), [Rest](../../data-structures/Rest/), [Take](../../data-structures/Take/), [Drop](../../data-structures/Drop/)

- Source: [`src/part.c`](https://github.com/stblake/mathilda/blob/main/src/part.c)
- Specification: [`docs/spec/builtins/data-structures.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/data-structures.md)
- Tests: [`tests/test_association.c`](https://github.com/stblake/mathilda/blob/main/tests/test_association.c)
- Tests: [`tests/test_compiledfunction.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compiledfunction.c)
- Tests: [`tests/test_eval.c`](https://github.com/stblake/mathilda/blob/main/tests/test_eval.c)
- Tests: [`tests/test_ndarray_selection.c`](https://github.com/stblake/mathilda/blob/main/tests/test_ndarray_selection.c)

## Notes & additional examples

### Notes

`Most[expr]` drops the last element; it is the complement of `Last`.
