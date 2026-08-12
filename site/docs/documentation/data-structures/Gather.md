# Gather

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Gather[list]`**

Gathers identical elements of list into sublists, giving {{group1}, {group2}, ...}. Sublists appear in order of the first occurrence of their element, and elements keep their input order within a sublist. Equal elements are collected from anywhere in the list, not only from adjacent runs (unlike Split). Gather\[list\] is equivalent to GatherBy\[list, Identity\].

## Examples (3)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (3)

```mathematica
In[1]:= Gather[{1, 7, 3, 7, 2, 3, 9}]
Out[1]= {{1}, {7, 7}, {3, 3}, {2}, {9}}

In[2]:= Gather[{a, b, a}]
Out[2]= {{a, a}, {b}}

In[3]:= Gather[{}]
Out[3]= {}
```

## Algorithm

--------------------------------------------------------------------------- gather.c — Gather[list], the identity case of GatherBy.

Gather partitions a list into sublists of structurally identical elements: every element appears in exactly one sublist, two elements share a sublist iff expr_eq holds between them, sublists appear in order of the first occurrence of their element, and within a sublist elements keep their input order. Unlike Split, grouping is not restricted to adjacent runs:

```text
    Gather[{1, 7, 3, 7, 2, 3, 9}]  ->  {{1}, {7, 7}, {3, 3}, {2}, {9}}
    Gather[{a, b, a}]              ->  {{a, a}, {b}}
```

The grouping itself is not reimplemented here. assoc_gather_core (assoc.c) is the same hash-indexed O(n) engine that backs GatherBy; passing a NULL key function selects the identity key, so Gather[l] and GatherBy[l, Identity] agree by construction, and the identity path skips the n Identity[x] applications that spelling it as GatherBy[l, Identity] would evaluate. --------------------------------------------------------------------------

## Implementation notes

**Attributes:** `Protected`.

## See also

[Split](../../structural-manipulation/Split/)

## References

- Source: [`src/list/list_init.c`](https://github.com/stblake/mathilda/blob/main/src/list/list_init.c)
- Specification: [`docs/spec/builtins/data-structures.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/data-structures.md)
- Tests: [`tests/test_list.c`](https://github.com/stblake/mathilda/blob/main/tests/test_list.c)
