# Counts

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Counts[list]`**

Gives \<|element -\> count, ...|\> tallying each distinct element. Hash-indexed: one O(n) pass.

## Examples (5)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (2)

```mathematica
In[1]:= Counts[{1, 2, 2, 3, 3, 3}]
Out[1]= <|1 -> 1, 2 -> 2, 3 -> 3|>

In[2]:= Counts[<|"a" -> 1, "b" -> 1, "c" -> 2|>]
Out[2]= <|1 -> 2, 2 -> 1|>
```

### Applications (3)

```mathematica
In[3]:= Counts[{a, b, a, c, b, a}]
Out[3]= <|a -> 3, b -> 2, c -> 1|>

In[4]:= Counts[{1, 1, 2, 3, 3, 3}]
Out[4]= <|1 -> 2, 2 -> 1, 3 -> 3|>
```

Letter frequencies

```mathematica
In[5]:= Counts[Characters["mississippi"]]
Out[5]= <|"m" -> 1, "i" -> 4, "s" -> 4, "p" -> 2|>
```

## Implementation notes

**Algorithm.** `builtin_counts` returns `<|element -> count, ...|>` in order of
first appearance. For a boxed `List` it makes a single hash pass: each element is
looked up in a `KeyIndex`, a new distinct value is recorded with count 1 and
every repeat increments the stored counter, after which the distinct keys and
their counts are emitted as rules. An `Association` argument is counted over its
values (`Counts[Values[assoc]]`).

**Data structures.** The distinct-value table is a `KeyIndex` open-addressing
hash set over borrowed `Expr*` keys paired with an `int64_t` count array; the
result is a canonical `Association` built directly from the rule array.

**Complexity / limits.** `O(n)` for a list of `n` elements. A packed/`NDArray`
argument takes the fast path: `counts_from_ndarray` runs `ndred_tally`
(`ndreduce.c`) over the raw `int64`/`float64` words — direct-indexed when the
value range allows, hashed otherwise — and relabels its `{key, count}` pairs as
`key -> count` rules, so the per-element work happens on machine words rather
than boxed `Expr`s. Tally declines the dtypes and ranks it cannot key faithfully
(complex, rank > 1, non-finite floats) and returns the List-path answer, which is
rewritten the same way.

**Attributes:** `Protected`.

## References

**See also:** [NDArray](../../linear-algebra/NDArray/), [Tally](../../data-structures/Tally/), [Association](../../data-structures/Association/)

- T. H. Cormen, C. E. Leiserson, R. L. Rivest and C. Stein, *Introduction to Algorithms*, 3rd ed. (MIT Press, 2009), §11 (hash tables, open addressing).
- Source: [`src/assoc.c`](https://github.com/stblake/mathilda/blob/main/src/assoc.c)
- Specification: [`docs/spec/builtins/data-structures.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/data-structures.md)
- Tests: [`tests/test_assoc_ops.c`](https://github.com/stblake/mathilda/blob/main/tests/test_assoc_ops.c)
- Tests: [`tests/test_association.c`](https://github.com/stblake/mathilda/blob/main/tests/test_association.c)
- Tests: [`tests/test_compile_assoc.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compile_assoc.c)
- Tests: [`tests/test_ndarray_functions.c`](https://github.com/stblake/mathilda/blob/main/tests/test_ndarray_functions.c)

## Notes & additional examples

### Notes

`Counts[list]` returns `<|element -> count, ...|>` in order of first appearance —
the association form of `Tally`, which gives the same information as a list of
`{element, count}` pairs. It is the standard histogram primitive: counting
characters, residues, or category labels. A packed integer or real buffer is
counted on its machine words (through `Tally`'s direct-indexed or hashed count)
and relabelled as rules, so large numeric data stay on the buffer. Use `CountsBy`
to count by a function of each element.
