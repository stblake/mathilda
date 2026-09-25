# SortBy

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`SortBy[list, f]`**

Sorts the elements of list by the canonical order of f applied to each element.

**`SortBy[assoc, f]`**

Sorts an association by f applied to each value.

**`SortBy[list, f, p]`**

Compares the f values with the ordering function p.

**`SortBy[f]`**

Operator form: SortBy\[f\]\[expr\] is SortBy\[expr, f\].

## Examples (8)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (8)

```mathematica
In[1]:= Sort[<|"a" -> 3, "b" -> 1, "c" -> 2|>]
Out[1]= <|"b" -> 1, "c" -> 2, "a" -> 3|>

In[2]:= SortBy[<|"a" -> {9}, "b" -> {1}|>, First]
Out[2]= <|"b" -> {1}, "a" -> {9}|>

In[3]:= Total[<|"a" -> 3, "b" -> 1, "c" -> 2|>]
Out[3]= 6

In[4]:= Join[<|"a" -> 1, "b" -> 2|>, <|"b" -> 3, "c" -> 4|>]
Out[4]= <|"a" -> 1, "b" -> 3, "c" -> 4|>

In[5]:= Sort[<|"a" -> 2, "b" -> 3, "c" -> 1|>, Greater]
Out[5]= <|"b" -> 3, "a" -> 2, "c" -> 1|>

In[6]:= SortBy[<|"a" -> {1, 2}, "b" -> {0, 5}, "c" -> {1, 1}|>, First, Greater]
Out[6]= <|"a" -> {1, 2}, "c" -> {1, 1}, "b" -> {0, 5}|>

In[7]:= Ordering[<|"a" -> 2, "b" -> 2, "c" -> 1|>, All, Greater]
Out[7]= {2, 1, 3}

In[8]:= ReverseSort[<|"a" -> 2, "b" -> 2, "c" -> 1|>]
Out[8]= <|"a" -> 2, "b" -> 2, "c" -> 1|>
```

## Implementation notes

- With an ordering function, the sort is Mathematica 15's top-down merge sort
  (the merge keeps the left element unless `p[left, right]` is `False` or
  `-1`), so ties under a strict `p` such as `Greater` land exactly where
  Mathematica puts them. `Ordering[..., p]` and `KeySort[assoc, p]` share it.

**Attributes:** `Protected`.

## References

**See also:** [Sort](../../data-structures/Sort/), [Total](../../arithmetic/Total/), [Min](../../data-structures/Min/), [Max](../../data-structures/Max/), [Join](../../data-structures/Join/), [Greater](../../comparisons/Greater/)

- Source: [`src/core.c`](https://github.com/stblake/mathilda/blob/main/src/core.c)
- Specification: [`docs/spec/builtins/data-structures.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/data-structures.md)
- Tests: [`tests/test_assoc_forms.c`](https://github.com/stblake/mathilda/blob/main/tests/test_assoc_forms.c)
- Tests: [`tests/test_assoc_read.c`](https://github.com/stblake/mathilda/blob/main/tests/test_assoc_read.c)
- Tests: [`tests/test_association.c`](https://github.com/stblake/mathilda/blob/main/tests/test_association.c)
