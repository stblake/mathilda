# First

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`First[expr] gives the first element of expr.`**

## Examples (4)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (3)

```mathematica
In[1]:= First[<|"a" -> 10, "b" -> 20|>] In[1b]:= First[<||>, 0] Out[1b]= 0

In[2]:= Rest[<|"a" -> 10, "b" -> 20, "c" -> 30|>]
Out[2]= <|"b" -> 20, "c" -> 30|>

In[3]:= Take[<|"a" -> 1, "b" -> 2, "c" -> 3|>, 2]
Out[3]= <|"a" -> 1, "b" -> 2|>
```

### Applications (1)

```mathematica
In[4]:= First[{a,b,c}]
Out[4]= a
```

## Performance

Against other systems, from the benchmark suite (same input, results cross-checked for agreement):

| case | Mathilda | Wolfram | Python |
|---|---:|---:|---:|
| Fourier 1200000 (mixed radix) | 12.7 s | 9.68 s | 8.53 s |
| ListConvolve 100000 x 2048 | 9.71 s | 1.1 s | 10.2 s |
| ListCorrelate 100000 x 2048 | 9.68 s | 1.1 s | 10.4 s |
| Fourier 262143 (awkward size) | 5.91 s | 5.13 s | 4.42 s |
| Fourier 2^18 (262144) | 4.46 s | 2.7 s | 2.35 s |
| NDS Van der Pol mu=1000 | 4.02 s | 0.234 s | 0.362 s |

## Implementation notes

`builtin_first` (in `src/part.c`) takes a single argument and returns a deep copy of its first element (`args[0]`). It returns `NULL` (unevaluated) when the argument is atomic or has no elements.

**Attributes:** none registered.

## References

**See also:** [Last](../../data-structures/Last/), [Rest](../../data-structures/Rest/), [Most](../../data-structures/Most/), [Take](../../data-structures/Take/), [Drop](../../data-structures/Drop/)

- Source: [`src/part.c`](https://github.com/stblake/mathilda/blob/main/src/part.c)
- Specification: [`docs/spec/builtins/data-structures.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/data-structures.md)
- Tests: [`tests/test_association.c`](https://github.com/stblake/mathilda/blob/main/tests/test_association.c)
- Tests: [`tests/test_basin_hopping.c`](https://github.com/stblake/mathilda/blob/main/tests/test_basin_hopping.c)
- Tests: [`tests/test_compile.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compile.c)
- Tests: [`tests/test_compiledfunction.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compiledfunction.c)

## Notes & additional examples

### Notes

`First[expr]` returns the first element (part 1) of any expression, not only lists.
