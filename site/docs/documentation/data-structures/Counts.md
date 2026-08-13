# Counts

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Counts[list]`**

Gives \<|element -\> count, ...|\> tallying each distinct element. Hash-indexed: one O(n) pass.

## Examples (2)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (2)

```mathematica
In[1]:= Counts[{1, 2, 2, 3, 3, 3}]
Out[1]= <|1 -> 1, 2 -> 2, 3 -> 3|>

In[2]:= Counts[<|"a" -> 1, "b" -> 1, "c" -> 2|>]
Out[2]= <|1 -> 2, 2 -> 1|>
```

## Implementation notes

**Attributes:** `Protected`.

## References

**See also:** [NDArray](../../linear-algebra/NDArray/), [Tally](../../data-structures/Tally/), [Association](../../data-structures/Association/)

- Source: [`src/assoc.c`](https://github.com/stblake/mathilda/blob/main/src/assoc.c)
- Specification: [`docs/spec/builtins/data-structures.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/data-structures.md)
- Tests: [`tests/test_association.c`](https://github.com/stblake/mathilda/blob/main/tests/test_association.c)
- Tests: [`tests/test_compile_assoc.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compile_assoc.c)
- Tests: [`tests/test_ndarray_functions.c`](https://github.com/stblake/mathilda/blob/main/tests/test_ndarray_functions.c)
