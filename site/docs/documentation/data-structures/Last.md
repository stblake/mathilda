# Last

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Last[expr] gives the last element of expr.`**

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
In[1]:= Last[{a,b,c}]
Out[1]= c
```

## Implementation notes

`builtin_last` (in `src/part.c`) takes a single argument and returns a deep copy of its final element (`args[arg_count - 1]`). It returns `NULL` (unevaluated) when the argument is atomic or empty.

**Attributes:** none registered.

## See also

[First](../../data-structures/First/), [Rest](../../data-structures/Rest/), [Most](../../data-structures/Most/), [Take](../../data-structures/Take/), [Drop](../../data-structures/Drop/)

## References

- Source: [`src/part.c`](https://github.com/stblake/mathilda/blob/main/src/part.c)
- Specification: [`docs/spec/builtins/data-structures.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/data-structures.md)
- Tests: [`tests/test_association.c`](https://github.com/stblake/mathilda/blob/main/tests/test_association.c)
- Tests: [`tests/test_compile.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compile.c)
- Tests: [`tests/test_divisors.c`](https://github.com/stblake/mathilda/blob/main/tests/test_divisors.c)
- Tests: [`tests/test_eval.c`](https://github.com/stblake/mathilda/blob/main/tests/test_eval.c)

## Notes & additional examples

### Notes

`Last[expr]` returns the final element of any expression.
