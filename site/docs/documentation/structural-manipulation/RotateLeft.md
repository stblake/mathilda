# RotateLeft

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`RotateLeft[expr, n] rotates the elements of expr n positions to the left.`**

## Examples (2)

Every input below was run against the current Mathilda build and its output recorded.

### Applications (2)

```mathematica
In[1]:= RotateLeft[{1,2,3,4},1]
Out[1]= {2, 3, 4, 1}

In[2]:= RotateLeft[{a,b,c,d,e},2]
Out[2]= {c, d, e, a, b}
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

**Algorithm.** `builtin_rotateleft` cyclically shifts elements toward the front by `n`
(default 1) using `rotate_rec`: at each level it computes the wrapped offset `((n mod len) +
len) mod len` and reads element `i` from source index `(i + offset) mod len`. A `List`-valued
`n` applies a per-level shift amount as the recursion descends into nested lists. `RotateRight`
is implemented by negating `n` and calling the same routine.

**Attributes:** `Protected`.

## See also

[RotateRight](../../structural-manipulation/RotateRight/)

## References

- Source: [`src/list.c`](https://github.com/stblake/mathilda/blob/main/src/list.c)
- Specification: [`docs/spec/builtins/structural-manipulation.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/structural-manipulation.md)
- Tests: [`tests/test_compiledfunction.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compiledfunction.c)
- Tests: [`tests/test_list.c`](https://github.com/stblake/mathilda/blob/main/tests/test_list.c)
- Tests: [`tests/test_ndarray_selection.c`](https://github.com/stblake/mathilda/blob/main/tests/test_ndarray_selection.c)
- Tests: [`tests/test_packed_list.c`](https://github.com/stblake/mathilda/blob/main/tests/test_packed_list.c)

## Notes & additional examples

### Notes

Cyclically shifts elements `n` positions to the left, wrapping the front elements to the back.
