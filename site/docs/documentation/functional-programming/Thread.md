# Thread

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Thread[f[args]]`**

"threads" f over any lists that appear in args.

**`Thread[f[args], h]`**

threads f over any objects with head h that appear in args.

**`Thread[f[args], h, n]`**

threads f over objects with head h that appear in the first n args.

<details>
<summary>Notes</summary>

Functions with attribute Listable are automatically threaded over lists. All the elements in the specified args whose heads are h must be of the same length. Arguments that do not have head h are copied as many times as there are elements in the arguments that do have head h. Thread specifies argument positions using the standard sequence specification: All       all elements None      no elements n         elements 1 through n -n        last n elements {n}       element n only {m, n}    elements m through n inclusive {m, n, s} elements m through n in steps of s

</details>

## Examples (18)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (14)

```mathematica
In[1]:= Thread[f[{a, b, c}]]
Out[1]= {f[a], f[b], f[c]}

In[2]:= Thread[f[{a, b, c}, x]]
Out[2]= {f[a, x], f[b, x], f[c, x]}

In[3]:= Thread[f[{a, b, c}, {x, y, z}]]
Out[3]= {f[a, x], f[b, y], f[c, z]}

In[4]:= Thread[{a, b, c} == {x, y, z}]
Out[4]= {a == x, b == y, c == z}

In[5]:= Thread[Log[x == y], Equal]
Out[5]= Log[x] == Log[y]

In[6]:= Thread[f[{a, b}, {r, s}, {u, v}, {x, y}], List]
Out[6]= {f[a, r, u, x], f[b, s, v, y]}

In[7]:= Thread[f[{a, b}, {r, s}, {u, v}, {x, y}], List, None]
Out[7]= f[{a, b}, {r, s}, {u, v}, {x, y}]

In[8]:= Thread[f[{a, b}, {r, s}, {u, v}, {x, y}], List, 2]
Out[8]= {f[a, r, {u, v}, {x, y}], f[b, s, {u, v}, {x, y}]}

In[9]:= Thread[f[{a, b}, {r, s}, {u, v}, {x, y}], List, -2]
Out[9]= {f[{a, b}, {r, s}, u, x], f[{a, b}, {r, s}, v, y]}

In[10]:= Thread[f[{a, b}, {r, s}, {u, v}, {x, y}], List, {2, 4}]
Out[10]= {f[{a, b}, r, u, x], f[{a, b}, s, v, y]}

In[11]:= Thread[f[{a, b}, {r, s}, {u, v}, {x, y}], List, {1, -1, 2}]
Out[11]= {f[a, {r, s}, u, {x, y}], f[b, {r, s}, v, {x, y}]}

In[12]:= Thread[f[a + b, c + d]]
Out[12]= f[a + b, c + d]

In[13]:= Thread[f[a + b, c + d], Plus]
Out[13]= f[a, c] + f[b, d]

In[14]:= Thread[f[{a, b, c}, h, {x, y, z}]]
Out[14]= {f[a, h, x], f[b, h, y], f[c, h, z]}
```

### Applications (4)

```mathematica
In[15]:= Thread[f[{a, b, c}]]
Out[15]= {f[a], f[b], f[c]}

In[16]:= Thread[{x, y, z} -> {1, 2, 3}]
Out[16]= {x -> 1, y -> 2, z -> 3}

In[17]:= Thread[f[{a, b}, {c, d}, x]]
Out[17]= {f[a, c, x], f[b, d, x]}

In[18]:= Thread[Equal[{a, b, c}, {1, 2, 3}]]
Out[18]= {a == 1, b == 2, c == 3}
```

## Implementation notes

**Algorithm.** `builtin_thread` distributes a function over equal-length lists:
`Thread[f[{a,b},{c,d}]]` → `{f[a,c], f[b,d]}`. Arguments are `(expr, h, n)`: `h`
is the threading head (default `List`) and `n` is a position spec selecting which
of `f`'s `K` arguments take part. `thread_parse_spec` turns the spec into a
boolean mask of length `K`, handling `All`/`None`, an integer (first/last `|n|`
positions), and `{n}` / `{m,n}` / `{m,n,s}` index ranges with negative-from-end
indexing.

It then scans the masked, `h`-headed arguments to determine the common threading
length `L`; if two such arguments have different lengths the expression is
returned unchanged (a diagnostic is issued here in some systems; Mathilda elides it).
For each `k` in `0..L-1` it builds `f[...]` taking element `k` from every masked
threadable argument and copying all other arguments verbatim, then wraps the `L`
calls under `h[...]` and runs `evaluate()` so `f`'s attributes (`Listable`,
`OneIdentity`, ...) apply. Atoms and the no-threadable-argument case return a copy
of `expr`.

**Data structures.** `Expr`-tree only; a `bool* mask` of length `K`, plus a
`wrap_args` array of the `L` per-slice calls. The threading head uses the interned
`List` symbol (`expr_ref`) so `expr_eq` pointer comparisons work.

- `Protected`.
- Functions with attribute `Listable` are automatically threaded over lists.
- All the elements in the specified args whose heads are `h` must be of the same length; otherwise the expression is returned unchanged.
- Arguments that do not have head `h` are copied as many times as there are elements in the arguments that do have head `h`.
- The position specifier `n` uses the standard sequence specification:

  | Spec | Meaning |
  |------|---------|
  | `All` | all elements |
  | `None` | no elements |
  | `n` | elements 1 through `n` |
  | `-n` | last `n` elements |
  | `{n}` | element `n` only |
  | `{m, n}` | elements `m` through `n` inclusive |
  | `{m, n, s}` | elements `m` through `n` in steps of `s` |

**Attributes:** `Protected`.

## References

- Source: [`src/funcprog.c`](https://github.com/stblake/mathilda/blob/main/src/funcprog.c)
- Specification: [`docs/spec/builtins/functional-programming.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/functional-programming.md)
- Tests: [`tests/test_mapthread.c`](https://github.com/stblake/mathilda/blob/main/tests/test_mapthread.c)
- Tests: [`tests/test_reduce.c`](https://github.com/stblake/mathilda/blob/main/tests/test_reduce.c)
- Tests: [`tests/test_thread.c`](https://github.com/stblake/mathilda/blob/main/tests/test_thread.c)

## Notes & additional examples

### Notes

`Thread[f[args]]` distributes `f` over any lists in `args`, pairing them
positionally; non-list arguments are broadcast to every element. This is the
idiomatic way to build a rule list from parallel name/value lists
(`Thread[vars -> vals]`) or a system of equations from two vectors
(`Thread[Equal[lhs, rhs]]`). All list arguments must have the same length.
`Thread[f[args], h]` threads over a custom head `h` instead of `List`.
