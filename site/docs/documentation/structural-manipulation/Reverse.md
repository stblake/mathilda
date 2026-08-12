# Reverse

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Reverse[expr] reverses the order of elements in expr.`**

## Examples (2)

Every input below was run against the current Mathilda build and its output recorded.

### Applications (2)

```mathematica
In[1]:= Reverse[{1,2,3}]
Out[1]= {3, 2, 1}

In[2]:= Reverse[{a,b,c,d}]
Out[2]= {d, c, b, a}
```

## Performance

Against other systems, from the benchmark suite (same input, results cross-checked for agreement):

| case | Mathilda | Wolfram | Python |
|---|---:|---:|---:|
| Sort 4x10^6 | 42.2 s | 68.7 s | 111 s |
| gather v[[idx]], 4x10^6 | 16.8 s | 6.66 s | 7.18 s |
| Union of 4x10^6 integers | 12.4 s | 71.1 s | 376 s |
| Reverse 4x10^6 | 5.37 s | 0.297 s | 0.982 s |
| Join two 2x10^6 | 0.899 s | 0.6 s | 0.397 s |
| RotateLeft 4x10^6 by 1000 | 0.855 s | 0.307 s | 0.456 s |

## Implementation notes

**Algorithm.** `builtin_reverse` recurses through the expression with `reverse_rec`, reversing
the argument order at the levels selected by an optional level spec (`should_reverse_at_level`
matches an integer level, or any level listed in a `List` spec; default is level 1). At each
visited node it builds a new function with the same head, drawing children either forward or
reversed depending on whether the current level is selected, and recursing into each child.

**Attributes:** `Protected`.

## References

- Source: [`src/list.c`](https://github.com/stblake/mathilda/blob/main/src/list.c)
- Specification: [`docs/spec/builtins/structural-manipulation.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/structural-manipulation.md)
- Tests: [`tests/test_association.c`](https://github.com/stblake/mathilda/blob/main/tests/test_association.c)
- Tests: [`tests/test_compile.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compile.c)
- Tests: [`tests/test_convolutions.c`](https://github.com/stblake/mathilda/blob/main/tests/test_convolutions.c)
- Tests: [`tests/test_correlations.c`](https://github.com/stblake/mathilda/blob/main/tests/test_correlations.c)

## Notes & additional examples

### Notes

`Reverse[expr]` reverses the order of top-level elements.
