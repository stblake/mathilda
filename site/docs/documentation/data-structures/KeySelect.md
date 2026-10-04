# KeySelect

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`KeySelect[assoc, pred]`**

Keeps the entries whose key satisfies pred.

## Examples (3)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (1)

```mathematica
In[1]:= KeySelect[<|1 -> 10, 2 -> 20, 3 -> 30|>, EvenQ]
Out[1]= <|2 -> 20|>
```

### Applications (2)

Keep even-numbered keys

```mathematica
In[2]:= KeySelect[<|1 -> a, 2 -> b, 3 -> c, 4 -> d|>, EvenQ]
Out[2]= <|2 -> b, 4 -> d|>
```

```mathematica
In[3]:= KeySelect[<|1 -> a, 2 -> b, 3 -> c, 4 -> d|>, OddQ]
Out[3]= <|1 -> a, 3 -> c|>
```

## Implementation notes

**Algorithm.** `builtin_keyselect` keeps the entries whose *key* passes a predicate.
For each entry it forms `pred[key]` (`apply1`), evaluates it, and keeps the whole
entry only when the verdict is literally the symbol `True`; any other result drops
the entry. The surviving entries, in their original order, form a fresh
`Association`.

**Data structures.** A single `Expr**` output array sized to the entry count; the
predicate application is a throwaway evaluated per key.

**Complexity / limits.** O(n), one predicate evaluation per entry. The companion
`KeySelect` tests keys; `Select` (and `Discard`) over an association test the values
instead. A non-`True` verdict — including an unevaluated predicate call — counts as
rejection.

**Attributes:** `Protected`.

## References

- Source: [`src/assoc.c`](https://github.com/stblake/mathilda/blob/main/src/assoc.c)
- Specification: [`docs/spec/builtins/data-structures.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/data-structures.md)
- Tests: [`tests/test_assoc_forms.c`](https://github.com/stblake/mathilda/blob/main/tests/test_assoc_forms.c)
- Tests: [`tests/test_association.c`](https://github.com/stblake/mathilda/blob/main/tests/test_association.c)

## Notes & additional examples

### Notes

The predicate is applied to each key; an entry is kept only when the result is
exactly `True`. Order is preserved. To filter by the values rather than the keys, use
`Select` over the association.
