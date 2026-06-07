# Thread

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
Thread[f[args]]
    "threads" f over any lists that appear in args.
Thread[f[args], h]
    threads f over any objects with head h that appear in args.
Thread[f[args], h, n]
    threads f over objects with head h that appear in the first n args.

Functions with attribute Listable are automatically threaded over
lists. All the elements in the specified args whose heads are h must
be of the same length. Arguments that do not have head h are copied
as many times as there are elements in the arguments that do have
head h.

Thread specifies argument positions using the standard sequence
specification:
    All       all elements
    None      no elements
    n         elements 1 through n
    -n        last n elements
    {n}       element n only
    {m, n}    elements m through n inclusive
    {m, n, s} elements m through n in steps of s
```

## Examples

All examples below are verified against the current Mathilda build.

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
```

## Implementation notes

- `Protected`.
- Functions with attribute `Listable` are automatically threaded over lists.
- All the elements in the specified args whose heads are `h` must be of the same length; otherwise the expression is returned unchanged.
- Arguments that do not have head `h` are copied as many times as there are elements in the arguments that do have head `h`.
- The position specifier `n` uses the standard sequence specification:

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/functional-programming.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/functional-programming.md)
