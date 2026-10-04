# MinimalBy

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`MinimalBy[list, f]`**

Gives the element(s) of list for which f is minimal (all ties, in order). Over an association, gives the entries whose value minimises f. MinimalBy\[f\] is the operator form.

## Examples (5)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (2)

```mathematica
In[1]:= MaximalBy[{1, -5, 3, -5, 2}, Abs]
Out[1]= {-5, -5}

In[2]:= MinimalBy[<|"a" -> 1, "b" -> 3, "c" -> 2|>, Identity]
Out[2]= <|"a" -> 1|>
```

### Applications (3)

Both -2 and 2 tie for the smallest absolute value

```mathematica
In[3]:= MinimalBy[{5, -2, 3, 2}, Abs]
Out[3]= {-2, 2}
```

A unique minimum is still returned as a list

```mathematica
In[4]:= MinimalBy[{1, -5, 3, -2}, Abs]
Out[4]= {1}
```

Operator form MinimalBy[f][list]

```mathematica
In[5]:= MinimalBy[Abs][{5, -2, 3, 2}]
Out[5]= {-2, 2}
```

## Implementation notes

**Algorithm.** `builtin_minimal_by` is `maximal_minimal_by(res, 1)`, the mode-1
twin of `MaximalBy`. It evaluates the key `f[e]` for every element (for an
association, `f` of each *value*) into a parallel `keys[]` array, finds the least
key in one pass by `expr_compare`, then collects in a second pass every element
whose key equals that minimum, in original order — so all ties are returned. The
result keeps the input's head (a list stays a list; an association returns the
matching entries). `MinimalBy[f]` is the operator form.

The ranking is Mathilda's canonical order (`expr_compare`), so the keys may be any
expressions, not only numbers.

**Data structures.** An owned `keys[]` array of the *n* evaluated key `Expr`s and
an owned `out[]` array of the tied elements, both freed once the result is built.
The input argument array is borrowed.

**Complexity / limits.** `O(n)` key evaluations and `O(n)` comparisons; it selects
the minimum rather than sorting. The empty collection returns itself.

**Attributes:** `Protected`.

## References

**See also:** [MaximalBy](../../functional-programming/MaximalBy/)

- Source: [`src/sort.c`](https://github.com/stblake/mathilda/blob/main/src/sort.c)
- Specification: [`docs/spec/builtins/functional-programming.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/functional-programming.md)
- Tests: [`tests/test_assoc_forms.c`](https://github.com/stblake/mathilda/blob/main/tests/test_assoc_forms.c)
- Tests: [`tests/test_association.c`](https://github.com/stblake/mathilda/blob/main/tests/test_association.c)
- Tests: [`tests/test_list.c`](https://github.com/stblake/mathilda/blob/main/tests/test_list.c)

## Notes & additional examples

### Notes

`MinimalBy[list, f]` returns the element (or elements) of `list` for which `f` is
smallest, by Mathilda's canonical order — the dual of `MaximalBy`. All ties are
returned, in original order, so the result is always a list (`{-2, 2}` both have
absolute value 2). The operator form `MinimalBy[f]` defers the collection; over an
association the entries whose *value* minimises `f` are returned.
