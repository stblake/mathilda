# CountsBy

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`CountsBy[list, f]`**

Gives \<|f\[x\] -\> count, ...|\> tallying elements by f\[x\].

## Examples (4)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (1)

```mathematica
In[1]:= CountsBy[Range[10], EvenQ]
Out[1]= <|False -> 5, True -> 5|>
```

### Applications (3)

```mathematica
In[2]:= CountsBy[{1, 2, 3, 4, 5}, EvenQ]
Out[2]= <|False -> 3, True -> 2|>

In[3]:= CountsBy[{-2, -1, 0, 1, 2}, Sign]
Out[3]= <|-1 -> 2, 0 -> 1, 1 -> 2|>

In[4]:= CountsBy[Range[10], Mod[#, 3] &]
Out[4]= <|1 -> 4, 2 -> 3, 0 -> 3|>
```

## Implementation notes

**Algorithm.** `builtin_countsby` returns `<|f[x] -> count|>`, tallying the
elements by the value of `f` applied to each. It evaluates `f[x]` once per
element (via the `apply1` helper) and keys the result into a `KeyIndex`: a new
`f`-value is recorded with count 1, a repeat increments the stored counter. The
distinct `f`-values and counts become the association's rules in first-appearance
order. An `Association` argument is tallied over its values by `f`
(`assoc_apply_over_values`).

**Data structures.** A `KeyIndex` open-addressing hash set over the *owned*
`f`-value keys, paired with an `int64_t` count array; the keys are adopted into
the result rules so no second copy is made.

**Complexity / limits.** `O(n)` evaluations of `f` plus `O(n)` hashing. Returns
`NULL` unless the argument is a `List` or association. There is no packed fast
path — `f` is an arbitrary symbolic function applied per element — so a packed
argument is materialised before counting.

**Attributes:** `Protected`.

## References

- T. H. Cormen, C. E. Leiserson, R. L. Rivest and C. Stein, *Introduction to Algorithms*, 3rd ed. (MIT Press, 2009), §11 (hash tables, open addressing).
- Source: [`src/assoc.c`](https://github.com/stblake/mathilda/blob/main/src/assoc.c)
- Specification: [`docs/spec/builtins/data-structures.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/data-structures.md)
- Tests: [`tests/test_assoc_forms.c`](https://github.com/stblake/mathilda/blob/main/tests/test_assoc_forms.c)
- Tests: [`tests/test_association.c`](https://github.com/stblake/mathilda/blob/main/tests/test_association.c)

## Notes & additional examples

### Notes

`CountsBy[list, f]` returns `<|f[x] -> count|>` — a histogram keyed by the value
of `f` applied to each element, in first-appearance order of the `f`-values. It
is `Counts` after bucketing by `f`: `CountsBy[list, Sign]` tallies how many
elements are negative, zero, and positive. Over an association the values are
tallied by `f`.
