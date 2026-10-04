# MaximalBy

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`MaximalBy[list, f]`**

Gives the element(s) of list for which f is maximal (all ties, in order). Over an association, gives the entries whose value maximises f. MaximalBy\[f\] is the operator form.

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

The element with the largest absolute value

```mathematica
In[3]:= MaximalBy[{1, -5, 3, -2}, Abs]
Out[3]= {-5}
```

All ties, in input order

```mathematica
In[4]:= MaximalBy[{"a", "bbb", "cc", "ddd"}, StringLength]
Out[4]= {"bbb", "ddd"}
```

Operator form MaximalBy[f][list]

```mathematica
In[5]:= MaximalBy[Abs][{1, -5, 3}]
Out[5]= {-5}
```

## Implementation notes

**Algorithm.** `builtin_maximal_by` is `maximal_minimal_by(res, 0)` (mode 1 is
`MinimalBy`). It evaluates the key `f[e]` for every element `e` — for an
association, `f` of each *value* — into a parallel `keys[]` array, finds the best
key in one pass by `expr_compare` (greatest for `MaximalBy`, least for
`MinimalBy`), then makes a second pass collecting **every** element whose key
equals the best, in original order. So all ties are returned, and the result keeps
the collection's head (a `List` stays a list, an association returns the matching
entries). `MaximalBy[f]` with one argument is the operator form `MaximalBy[#1,
f] &`.

The comparison is Mathilda's canonical order (`expr_compare`), the same total
order `Sort` uses, so keys need not be numeric — any expressions are ranked.

**Data structures.** An owned `keys[]` array of the *n* evaluated key `Expr`s, and
an owned `out[]` array of the tied elements; both freed after the result node is
built. The input argument array is borrowed.

**Complexity / limits.** `O(n)` key evaluations and `O(n)` comparisons — it finds
an extreme, so it does **not** sort. The empty collection returns itself.

**Attributes:** `Protected`.

## References

**See also:** [MinimalBy](../../functional-programming/MinimalBy/)

- Source: [`src/sort.c`](https://github.com/stblake/mathilda/blob/main/src/sort.c)
- Specification: [`docs/spec/builtins/functional-programming.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/functional-programming.md)
- Tests: [`tests/test_assoc_forms.c`](https://github.com/stblake/mathilda/blob/main/tests/test_assoc_forms.c)
- Tests: [`tests/test_association.c`](https://github.com/stblake/mathilda/blob/main/tests/test_association.c)

## Notes & additional examples

### Notes

`MaximalBy[list, f]` returns the element (or elements) of `list` for which `f` is
largest, measured by Mathilda's canonical order. All elements tying for the
maximum are returned together, in their original order, so the result is always a
list — here `{"bbb", "ddd"}` both have length 3. It selects an extreme rather than
sorting. The operator form `MaximalBy[f]` applies to a collection supplied later;
over an association the entries whose *value* maximises `f` are returned.
