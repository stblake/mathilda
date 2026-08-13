# Rest

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Rest[expr] gives all but the first element of expr.`**

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
In[3]:= Rest[{a,b,c,d}]
Out[3]= {b, c, d}
```

## Implementation notes

**Algorithm.** `builtin_rest` returns a copy of the input with its first element dropped: it
copies args `1 .. n−1` into a new function node with the same head. Returns `NULL`
(unevaluated) for atoms or empty expressions.

**Attributes:** none registered.

## References

**See also:** [First](../../data-structures/First/), [Last](../../data-structures/Last/), [Most](../../data-structures/Most/), [Take](../../data-structures/Take/), [Drop](../../data-structures/Drop/)

- Source: [`src/part.c`](https://github.com/stblake/mathilda/blob/main/src/part.c)
- Specification: [`docs/spec/builtins/data-structures.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/data-structures.md)
- Tests: [`tests/test_association.c`](https://github.com/stblake/mathilda/blob/main/tests/test_association.c)
- Tests: [`tests/test_compiledfunction.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compiledfunction.c)
- Tests: [`tests/test_eval.c`](https://github.com/stblake/mathilda/blob/main/tests/test_eval.c)
- Tests: [`tests/test_integer_exponent.c`](https://github.com/stblake/mathilda/blob/main/tests/test_integer_exponent.c)

## Notes & additional examples

### Notes

`Rest[expr]` drops the first element; it is the complement of `First`.
