# MapThread

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`MapThread[f, {{a1, a2, ...}, {b1, b2, ...}, ...}]`**

gives {f\[a1, b1, ...\], f\[a2, b2, ...\], ...}, applying f to corresponding elements of the lists.

**`MapThread[f, {e1, e2, ...}, n]`**

applies f to the parts of the ei at level n.

<details>
<summary>Notes</summary>

The ei must all have the same shape down through level n. MapThread is a generalization of Map to functions of several variables; it takes the function and its argument lists separately, unlike Thread. Lists of associations with identical keys thread over their values.

</details>

## Examples (3)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (3)

```mathematica
In[1]:= MapThread[f, {{a, b, c}, {x, y, z}}]
Out[1]= {f[a, x], f[b, y], f[c, z]}

In[2]:= MapThread[f, {{{a, b}, {c, d}}, {{u, v}, {s, t}}}, 2]
Out[2]= {{f[a, u], f[b, v]}, {f[c, s], f[d, t]}}

In[3]:= MapThread[Plus, {{a, b, c}, {u, v, w}, {x, y, z}}]
Out[3]= {a + u + x, b + v + y, c + w + z}
```

## Performance

Against other systems, from the benchmark suite (same input, results cross-checked for agreement):

| case | Mathilda | Wolfram | Python |
|---|---:|---:|---:|
| Clip to [0.25, 0.75] over 4x10^6 | 575 s | 1.95 s | 0.953 s |
| MapThread[Max] over 4x10^6 | 14.8 s | 692 s | 0.772 s |
| MapThread[Min] over 4x10^6 | 14.7 s | 687 s | 0.769 s |
| integer Mod over 4x10^6 | 3.88 s | 0.504 s | 3.28 s |
| a b + a over 4x10^6 | 0.754 s | 1.07 s | 1.41 s |
| a + b over 4x10^6 | 0.383 s | 0.516 s | 0.74 s |

## Implementation notes

**Attributes:** `Protected`.

## References

**See also:** [Map](../../data-structures/Map/), [Thread](../../functional-programming/Thread/), [NDArray](../../linear-algebra/NDArray/)

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/functional-programming.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/functional-programming.md)
- Tests: [`tests/test_image.c`](https://github.com/stblake/mathilda/blob/main/tests/test_image.c)
- Tests: [`tests/test_map_ndarray.c`](https://github.com/stblake/mathilda/blob/main/tests/test_map_ndarray.c)
- Tests: [`tests/test_mapthread.c`](https://github.com/stblake/mathilda/blob/main/tests/test_mapthread.c)
- Tests: [`tests/test_ndarray_selection.c`](https://github.com/stblake/mathilda/blob/main/tests/test_ndarray_selection.c)
