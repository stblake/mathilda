# Length

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Length[expr]`**

gives the number of top-level elements in expr (the arity of its head).  Length of any atom is 0.

## Examples (1)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (1)

```mathematica
In[1]:= Length[{a, b, c}]
Out[1]= 3
```

## Performance

Against other systems, from the benchmark suite (same input, results cross-checked for agreement):

| case | Mathilda | Wolfram | Python |
|---|---:|---:|---:|
| GroebnerBasis cyclic-5 | 43.5 s | 18 s | 8.78e+03 s |
| NSolve degree 40 | 9.27 s | 0.779 s | 341 s |
| Factor x^120 - 1 | 7.01 s | 0.044 s | 3.85 s |
| Graph construction, 40000 edges | 5.81 s | 3.98 s | 14.7 s |
| Characters of 200k chars | 4.38 s | 2.54 s | 0.419 s |
| StringSplit on space, 200k chars | 3.12 s | 4.13 s | 0.723 s |

## Implementation notes

`builtin_length` (in `src/core.c`) returns an `EXPR_INTEGER` equal to the argument's `arg_count` when it is an `EXPR_FUNCTION`, and `0` for atoms (which have no parts).

- Returns the count of top-level arguments for functions.
- Returns `0` for all atomic expressions.

**Attributes:** none registered.

## References

- Source: [`src/core.c`](https://github.com/stblake/mathilda/blob/main/src/core.c)
- Specification: [`docs/spec/builtins/structural-manipulation.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/structural-manipulation.md)
- Tests: [`tests/test_association.c`](https://github.com/stblake/mathilda/blob/main/tests/test_association.c)
- Tests: [`tests/test_autocompile.c`](https://github.com/stblake/mathilda/blob/main/tests/test_autocompile.c)
- Tests: [`tests/test_boolean.c`](https://github.com/stblake/mathilda/blob/main/tests/test_boolean.c)
- Tests: [`tests/test_compile.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compile.c)
