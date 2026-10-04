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

## Examples (7)

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

### Applications (4)

F is applied to corresponding elements

```mathematica
In[4]:= MapThread[f, {{a, b, c}, {x, y, z}}]
Out[4]= {f[a, x], f[b, y], f[c, z]}
```

Add two lists componentwise

```mathematica
In[5]:= MapThread[Plus, {{1, 2, 3}, {10, 20, 30}}]
Out[5]= {11, 22, 33}
```

A pure function of two arguments

```mathematica
In[6]:= MapThread[#1^#2 &, {{2, 3, 4}, {2, 2, 3}}]
Out[6]= {4, 9, 64}
```

Thread two levels deep

```mathematica
In[7]:= MapThread[f, {{{1, 2}, {3, 4}}, {{5, 6}, {7, 8}}}, 2]
Out[7]= {{f[1, 5], f[2, 6]}, {f[3, 7], f[4, 8]}}
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

**Algorithm.** `builtin_mapthread` applies `f` across *k* lists in parallel:
`MapThread[f, {{a1,a2,...}, {b1,b2,...}, ...}]` gives `{f[a1,b1,...],
f[a2,b2,...], ...}`. The outer container must be a `List` of the *k* expressions
to thread, and `f` takes *k* arguments. At the default level 1 the
`nd_mapthread2` fast path handles the common two-list case over packed buffers;
otherwise `mapthread_rec` descends `level` list levels in lock-step, building one
`f[...]` leaf per tuple of corresponding parts, then the whole result is
`evaluate`d once so `f`'s attributes fire and numeric leaves reduce. A structural
mismatch (the lists are not the same shape down through `level`) makes the
recursion return `NULL`, leaving the call unevaluated.

Any `NDArray` entry is materialised to a nested list before threading, so an
array threads exactly like the corresponding list; if an input was packed, the
result is repacked with that input's dtype (`map_try_repack`) — packed in, packed
out. `MapThread[f, {}]` is `{}`.

**Data structures.** The *k* thread entries are borrowed from the argument list
(copied into a scratch array only when an NDArray entry must be delisted).
`mapthread_rec` works on the borrowed sub-expression arrays, so no spine is copied
beyond the result being built.

**Complexity / limits.** `O(k · N)` leaf constructions for *N* tuples, plus one
final evaluation pass. The level argument must be a non-negative integer; the
entries must agree in shape down through that level or the call declines.

**Attributes:** `Protected`.

## References

**See also:** [Map](../../data-structures/Map/), [Thread](../../functional-programming/Thread/), [NDArray](../../linear-algebra/NDArray/)

- Source: [`src/funcprog.c`](https://github.com/stblake/mathilda/blob/main/src/funcprog.c)
- Specification: [`docs/spec/builtins/functional-programming.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/functional-programming.md)
- Tests: [`tests/test_image.c`](https://github.com/stblake/mathilda/blob/main/tests/test_image.c)
- Tests: [`tests/test_map_ndarray.c`](https://github.com/stblake/mathilda/blob/main/tests/test_map_ndarray.c)
- Tests: [`tests/test_mapthread.c`](https://github.com/stblake/mathilda/blob/main/tests/test_mapthread.c)
- Tests: [`tests/test_ndarray_selection.c`](https://github.com/stblake/mathilda/blob/main/tests/test_ndarray_selection.c)

## Notes & additional examples

### Notes

`MapThread[f, {l1, l2, ...}]` threads `f` over corresponding elements of the lists
`li`, so `f` receives one argument per list. It is the several-variable
generalisation of `Map`, and unlike `Thread` it takes the function and the argument
lists separately. The lists must have the same length; with the optional level `n`
they must agree in shape down through level `n`.

The pure-function form `#1^#2 &` shows the first and second threaded arguments as
`#1` and `#2`. Threading over packed arrays keeps the result packed.
