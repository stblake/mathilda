# Gather

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
Gather[list]
    Gathers identical elements of list into sublists, giving
    {{group1}, {group2}, ...}. Sublists appear in order of the first
    occurrence of their element, and elements keep their input order
    within a sublist. Equal elements are collected from anywhere in the
    list, not only from adjacent runs (unlike Split).
    Gather[list] is equivalent to GatherBy[list, Identity].
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= Gather[{1, 7, 3, 7, 2, 3, 9}]
Out[1]= {{1}, {7, 7}, {3, 3}, {2}, {9}}

In[2]:= Gather[{a, b, a}]
Out[2]= {{a, a}, {b}}

In[3]:= Gather[{}]
Out[3]= {}
```

## Implementation notes

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/list/list_init.c`](https://github.com/stblake/mathilda/blob/main/src/list/list_init.c)
- Specification: [`docs/spec/builtins/data-structures.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/data-structures.md)
