# Ordering

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Ordering[list] gives the positions in list at which each successive element of Sort[list] appears, so that list[[Ordering[list]]] is Sort[list].`**

**`Ordering[list, n] gives the positions of the n smallest elements; Ordering[list, -n] gives the positions of the n largest.`**

**`Ordering[list, seq] is equivalent to Take[Ordering[list], seq], where seq may be an integer n or -n, a {m, n} or {m, n, s} span, UpTo[k], or All.`**

**`Ordering[list, seq, p] orders using the ordering function p, as in Sort[list, p].`**

<details>
<summary>Notes</summary>

Ties are broken by original position (Ordering is stable). Ordering works on an expression with any head, and on an Association (ordering its values), always returning a list of integer positions. Ordering has a packed-array fast path and is compilable.

</details>

## Examples (5)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (5)

```mathematica
In[1]:= Ordering[{c, a, b}]
Out[1]= {2, 3, 1}

In[2]:= Ordering[{2, 6, 1, 9, 1, 2, 3}, 4]
Out[2]= {3, 5, 1, 6}

In[3]:= Ordering[{2, 6, 1, 9, 1, 2, 3}, -1]
Out[3]= {4}

In[4]:= Ordering[{2, 6, 1, 9, 1, 2, 3}, All, Greater]
Out[4]= {4, 2, 7, 6, 1, 5, 3}

In[5]:= Ordering[<|1 -> c, 2 -> a, 3 -> b|>]
Out[5]= {2, 3, 1}
```

## Implementation notes

- `Protected`.
- Uses the same internal canonical comparison (`expr_compare`) as `Sort`, and the same custom-ordering-function convention (`p` may return `1`, `0`, `-1`, `True`, or `False`).
- **Stable**: ties are broken by original position, so `Ordering[list, 1]` gives the position of the *first* minimum and `Ordering[{2, 2, 1}]` is `{3, 1, 2}`.
- The result is always a `List` of integer positions, regardless of `list`'s head — `Ordering[f[3, 1, 2]]` is `{2, 3, 1}`.
- Over an `Association`, orders by the **values** and returns their positions.
- Packed-array fast path: on a machine-number vector it argsorts the buffer directly (int64 argsort past `2^53` is exact), returning a packed int64 permutation.
- Compilable inside `Compile[]` and auto-compiled: `Ordering[vector]` lowers to a delegated buffer argsort whose result element type is always integer. A complex dtype, rank ≥ 2, or a custom comparator fall back to the interpreter.

**Attributes:** `Protected`.

## See also

[Sort](../../data-structures/Sort/), [List](../../other-advanced/List/), [Association](../../data-structures/Association/)

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/structural-manipulation.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/structural-manipulation.md)
- Tests: [`tests/test_compiledfunction.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compiledfunction.c)
- Tests: [`tests/test_packed_list.c`](https://github.com/stblake/mathilda/blob/main/tests/test_packed_list.c)
- Tests: [`tests/test_sort.c`](https://github.com/stblake/mathilda/blob/main/tests/test_sort.c)
