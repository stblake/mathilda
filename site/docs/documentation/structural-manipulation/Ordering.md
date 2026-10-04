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

Without p, ties are broken by original position (Ordering is stable). Ordering works on an expression with any head, and on an Association (ordering its values), always returning a list of integer positions. Ordering has a packed-array fast path and is compilable.

</details>

## Examples (9)

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

### Applications (4)

Positions that sort the list: list[[%]] is Sort[list]

```mathematica
In[6]:= Ordering[{c, a, b}]
Out[6]= {2, 3, 1}
```

Positions of the 4 smallest

```mathematica
In[7]:= Ordering[{2, 6, 1, 9, 1, 2, 3}, 4]
Out[7]= {3, 5, 1, 6}
```

Stable: ties keep input order

```mathematica
In[8]:= Ordering[{2, 2, 1}]
Out[8]= {3, 1, 2}
```

Order descending with a custom function

```mathematica
In[9]:= Ordering[{5, 3, 8, 1}, All, Greater]
Out[9]= {3, 1, 2, 4}
```

## Implementation notes

**Algorithm.** `builtin_ordering` gives the permutation of 1-based positions that
sorts `list`, so that `list[[Ordering[list]]] === Sort[list]`. It collects
borrowed subject pointers — for an `Association` each entry's *value* (the Rule's
second argument), otherwise the element itself — and sorts an index array `idx`.
Without a custom ordering function it runs libc `qsort` through
`ordering_index_compare`, which orders by `expr_compare` of the pointed-at
subjects and breaks ties by the original index ascending; that tiebreak is what
turns the unstable `qsort` into a **stable** argsort, so equal elements keep
their input order (`Ordering[{2, 2, 1}]` is `{3, 1, 2}`, and `Ordering[list, 1]`
names the first minimum). A custom `p` (3-arg form) uses the merge-sort
permutation `p_sort_perm`, placing ties exactly as `Sort[list, p]` does.
`Ordering[list, seq]` emits only a `Take`-style slice of the permutation
(`get_seq_spec_indices` handles an `Integer n`/`-n`, a `{m, n[, s]}` span,
`UpTo[k]`, or `All`). The result is always a `List` of `Integer` positions,
regardless of `list`'s head.

**Data structures.** A borrowed `Expr** subjects`, an `int64_t* idx` permutation,
and an optional `int64_t* sel` of selected ranks; the comparator reads a
file-static `(ordering_subjects, ordering_p)` context saved/restored around each
sort so a nested `Ordering` recurses correctly. The result is offered to the
packer (`pack_offer`), since a list of positions packs and is often used as
`Part` indices next.

**Complexity / limits.** `O(n log n)` comparisons. A packed or visible rank-1
`NDArray` takes the buffer fast path `ndstruct_ordering`, which argsorts the
machine words directly and returns a packed `int64` permutation — so, unlike
`Sort`, its result dtype is always `Integer` (it is on `pack.c`'s `AWARE` and
`INT64_OK` lists); rank ≥ 2, a complex dtype, or a custom comparator decline by
materialising and re-evaluating. It is also compilable:
`CompileDiagnostics[{{v, _Real, 1}}, Ordering[v]]` reports `Compiled -> True`
with `ResultType -> Array`, lowering to a delegated buffer argsort.
`Protected`.

- `Protected`.
- Uses the same internal canonical comparison (`expr_compare`) as `Sort`, and the same custom-ordering-function convention (`p` may return `1`, `0`, `-1`, `True`, or `False`).
- **Stable**: ties are broken by original position, so `Ordering[list, 1]` gives the position of the *first* minimum and `Ordering[{2, 2, 1}]` is `{3, 1, 2}`.
- The result is always a `List` of integer positions, regardless of `list`'s head — `Ordering[f[3, 1, 2]]` is `{2, 3, 1}`.
- Over an `Association`, orders by the **values** and returns their positions.
- Packed-array fast path: on a machine-number vector it argsorts the buffer directly (int64 argsort past `2^53` is exact), returning a packed int64 permutation.
- Compilable inside `Compile[]` and auto-compiled: `Ordering[vector]` lowers to a delegated buffer argsort whose result element type is always integer. A complex dtype, rank ≥ 2, or a custom comparator fall back to the interpreter.

**Attributes:** `Protected`.

## References

**See also:** [Sort](../../data-structures/Sort/), [List](../../other-advanced/List/), [Association](../../data-structures/Association/)

- Source: [`src/sort.c`](https://github.com/stblake/mathilda/blob/main/src/sort.c)
- Specification: [`docs/spec/builtins/structural-manipulation.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/structural-manipulation.md)
- Tests: [`tests/test_assoc_read.c`](https://github.com/stblake/mathilda/blob/main/tests/test_assoc_read.c)
- Tests: [`tests/test_bignum_rational_numeric.c`](https://github.com/stblake/mathilda/blob/main/tests/test_bignum_rational_numeric.c)
- Tests: [`tests/test_compiledfunction.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compiledfunction.c)
- Tests: [`tests/test_ml_classify.c`](https://github.com/stblake/mathilda/blob/main/tests/test_ml_classify.c)

## Notes & additional examples

### Notes

`Ordering[list]` gives the permutation of 1-based positions for which
`list[[Ordering[list]]]` is `Sort[list]`. The argsort is **stable** — ties are
broken by original position — so `Ordering[{2, 2, 1}]` is `{3, 1, 2}` and
`Ordering[list, 1]` names the *first* minimum. `Ordering[list, seq]` is
`Take[Ordering[list], seq]` (an integer `n`/`-n` for the `n` smallest/largest, a
`{m, n[, s]}` span, `UpTo[k]`, or `All`), and the 3-argument form orders by a
function `p` as in `Sort[list, p]`. The result is always a `List` of integers
regardless of `list`'s head; over an `Association` it orders by the values. A
machine-number vector takes a packed argsort (result dtype always integer), and
`Ordering[vector]` compiles.
