# OrderedQ

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`OrderedQ[h[e1, e2, ...]] gives True if the elements are in canonical order, and False otherwise.`**

**`OrderedQ[expr, p] uses the ordering function p to determine whether each pair of elements is in order.`**

## Examples (10)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (5)

```mathematica
In[1]:= OrderedQ[{1, 4, 2}]
Out[1]= False

In[2]:= OrderedQ[{"cat", "catfish", "fish"}]
Out[2]= True

In[3]:= OrderedQ[{1, Sqrt[2], 2, E, 3, Pi}, Less]
Out[3]= True

In[4]:= OrderedQ[{{a, 2}, {c, 1}, {d, 3}}, #1[[2]] < #2[[2]] &]
Out[4]= False

In[5]:= OrderedQ[f[b, a, c]]
Out[5]= False
```

### Applications (5)

```mathematica
In[6]:= OrderedQ[{1, 2, 3}]
Out[6]= True

In[7]:= OrderedQ[{3, 1, 2}]
Out[7]= False

In[8]:= OrderedQ[{5, 4, 3, 2, 1}, Greater]
Out[8]= True

In[9]:= Select[Permutations[{1, 2, 3}], OrderedQ]
Out[9]= {{1, 2, 3}}

In[10]:= OrderedQ[{x^2, x, 1}, (#1 > #2 &)]
Out[10]= True
```

## Implementation notes

**Algorithm.** `builtin_orderedq` scans adjacent element pairs of the list once and returns
`True` iff every pair is in canonical order. With no ordering function it uses `expr_compare`
(the canonical order described under `Sort`); with a second argument `p` it evaluates `p[a, b]`
and treats `True`/`1` as ordered. The first out-of-order pair short-circuits to `False`. Empty
and single-element lists are trivially `True`. This is the linear-scan predicate corresponding
to the same comparator that `Sort` uses.

- `Protected`.
- Uses the same internal canonical comparison logic as `Sort` by default.
- Custom ordering function `p` may return `1`, `0`, `-1`, `True`, or `False`.
- `OrderedQ` works with any expression head, not just `List`.
- Automatically handles 0- and 1-element lists.

**Attributes:** `Protected`.

## References

**See also:** [Sort](../../data-structures/Sort/), [List](../../other-advanced/List/)

- Source: [`src/sort.c`](https://github.com/stblake/mathilda/blob/main/src/sort.c)
- Specification: [`docs/spec/builtins/structural-manipulation.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/structural-manipulation.md)
- Tests: [`tests/test_ndarray_functions.c`](https://github.com/stblake/mathilda/blob/main/tests/test_ndarray_functions.c)
- Tests: [`tests/test_sort.c`](https://github.com/stblake/mathilda/blob/main/tests/test_sort.c)

## Notes & additional examples

### Notes

`OrderedQ[h[e1, e2, ...]]` returns `True` exactly when each adjacent pair is in
canonical order under `expr_compare` (the same total order the evaluator uses
for `Orderless` heads and `Sort`). Ties (equal elements) count as ordered, so a
list with repeats can still be `OrderedQ`. The two-argument form
`OrderedQ[expr, p]` replaces the comparator with any predicate `p[a, b]`,
letting you test for descending order, key-based order, or domain-specific
orderings.
