# SubsetQ

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`SubsetQ[a, b]`**

Gives True if every element of b occurs in a (multiplicity ignored). Lists and associations (compared by value) may be mixed; other expressions must share a head.

## Examples (5)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (3)

```mathematica
In[1]:= SubsetQ[{1, 2, 3}, {3, 1}]
Out[1]= True

In[2]:= SubsetQ[{1, 2}, {1, 4}]
Out[2]= False

In[3]:= SubsetQ[<|"a" -> 1, "b" -> 2|>, {2}]
Out[3]= True
```

### Applications (2)

Is the second a subset of the first?

```mathematica
In[4]:= SubsetQ[{1, 2, 3, 4}, {2, 4}]
Out[4]= True
```

```mathematica
In[5]:= SubsetQ[{1, 2, 3}, {4}]
Out[5]= False
```

## Implementation notes

**Algorithm.** `builtin_subsetq` tests whether every element of `b` occurs in `a`,
ignoring multiplicity. It loads all of `a`'s elements into an `ExprSet` (a hash set),
then probes each element of `b`, short-circuiting on the first miss. Lists and
associations mix freely and are compared by element (by *value* for an association);
any other pair of expressions must share a head, or the `SubsetQ::heads` message
fires and the call is left unevaluated.

**Data structures.** An `ExprSet` hash set over `a`'s elements (membership by
`expr_eq`/`expr_hash`); `elem_at` abstracts list vs. association element access.

**Complexity / limits.** O(|a| + |b|) — one pass to build the set, one to probe.
Multiplicity is ignored, so `SubsetQ[{1, 1}, {1}]` and `SubsetQ[{1}, {1, 1}]` are both
`True`.

**Attributes:** `Protected`.

## References

- T. H. Cormen, C. E. Leiserson, R. L. Rivest and C. Stein, *Introduction to Algorithms*, 3rd ed. (MIT Press, 2009), ch. 11 (hash sets).
- Source: [`src/assoc_ops.c`](https://github.com/stblake/mathilda/blob/main/src/assoc_ops.c)
- Specification: [`docs/spec/builtins/data-structures.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/data-structures.md)
- Tests: [`tests/test_assoc_ops.c`](https://github.com/stblake/mathilda/blob/main/tests/test_assoc_ops.c)

## Notes & additional examples

### Notes

`SubsetQ[a, b]` is `True` when every element of `b` occurs in `a`; the test is built
on a hash set, so it is O(|a| + |b|). Multiplicity is ignored. Lists and associations
(compared by value) mix freely, but two other expressions must share a head.
