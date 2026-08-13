# SortBy

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`SortBy[list, f]`**

Sorts the elements of list by the canonical order of f applied to each element.

**`SortBy[assoc, f]`**

Sorts an association by f applied to each value.

**`SortBy[f]`**

Operator form: SortBy\[f\]\[expr\] is SortBy\[expr, f\].

## Examples (4)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (4)

```mathematica
In[1]:= Sort[<|"a" -> 3, "b" -> 1, "c" -> 2|>]
Out[1]= <|"b" -> 1, "c" -> 2, "a" -> 3|>

In[2]:= SortBy[<|"a" -> {9}, "b" -> {1}|>, First]
Out[2]= <|"b" -> {1}, "a" -> {9}|>

In[3]:= Total[<|"a" -> 3, "b" -> 1, "c" -> 2|>]
Out[3]= 6

In[4]:= Join[<|"a" -> 1, "b" -> 2|>, <|"b" -> 3, "c" -> 4|>]
Out[4]= <|"a" -> 1, "b" -> 3, "c" -> 4|>
```

## Implementation notes

**Attributes:** `Protected`.

## References

**See also:** [Sort](../../data-structures/Sort/), [Total](../../arithmetic/Total/), [Min](../../data-structures/Min/), [Max](../../data-structures/Max/), [Join](../../data-structures/Join/)

- Source: [`src/core.c`](https://github.com/stblake/mathilda/blob/main/src/core.c)
- Specification: [`docs/spec/builtins/data-structures.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/data-structures.md)
- Tests: [`tests/test_association.c`](https://github.com/stblake/mathilda/blob/main/tests/test_association.c)
