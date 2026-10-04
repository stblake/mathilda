# PositionIndex

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`PositionIndex[list]`**

Gives \<|value -\> {positions}|\> mapping each distinct element to the list of 1-based positions where it occurs. O(n).

**`PositionIndex[assoc]`**

Gives \<|value -\> {keys}|\>, the keys at which each value occurs.

## Examples (3)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (2)

```mathematica
In[1]:= PositionIndex[{a, b, a, c, a, b}]
Out[1]= <|a -> {1, 3, 5}, b -> {2, 6}, c -> {4}|>

In[2]:= PositionIndex[<|"a" -> x, "b" -> y, "c" -> x|>]
Out[2]= <|x -> {"a", "c"}, y -> {"b"}|>
```

### Applications (1)

Each value -> the positions it occurs at

```mathematica
In[3]:= PositionIndex[{a, b, a, c, b, a}]
Out[3]= <|a -> {1, 3, 6}, b -> {2, 5}, c -> {4}|>
```

## Implementation notes

**Algorithm.** `builtin_positionindex` builds `<|value -> {positions}|>` in one hash
pass over the list. Each element is probed in a `KeyIndex`; a new element registers a
growable bucket, and every occurrence appends its 1-based position
(`expr_new_integer(i + 1)`). Distinct values keep first-appearance order. The
`assoc_ops_init` wrapper adds `PositionIndex[assoc]`, which maps each distinct value
to the list of *keys* at which it occurs.

**Data structures.** A single open-addressing `KeyIndex` over the distinct elements,
with a per-element growable `Expr**` bucket of integer positions (doubling on growth).

**Complexity / limits.** O(n), amortised O(1) per element. The result is a
re-canonicalised `Association`; because it indexes by value equality, unhashable or
symbolic elements are grouped by `expr_eq` just like concrete ones.

**Attributes:** `Protected`.

## References

- T. H. Cormen, C. E. Leiserson, R. L. Rivest and C. Stein, *Introduction to Algorithms*, 3rd ed. (MIT Press, 2009), ch. 11 (open-addressing hash tables).
- Source: [`src/assoc.c`](https://github.com/stblake/mathilda/blob/main/src/assoc.c)
- Specification: [`docs/spec/builtins/data-structures.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/data-structures.md)
- Tests: [`tests/test_assoc_ops.c`](https://github.com/stblake/mathilda/blob/main/tests/test_assoc_ops.c)
- Tests: [`tests/test_association.c`](https://github.com/stblake/mathilda/blob/main/tests/test_association.c)

## Notes & additional examples

### Notes

The result maps each distinct element to the sorted list of 1-based positions where
it appears, in first-appearance order of the values. It is the inverse view of a list
and runs in one O(n) hash pass. `PositionIndex[assoc]` instead maps each distinct
value to the list of keys that hold it.
