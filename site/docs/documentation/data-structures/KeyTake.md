# KeyTake

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`KeyTake[assoc, {k1, ...}]`**

Gives the association of only the specified keys, in the requested order.

## Examples (2)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (2)

```mathematica
In[1]:= KeyTake[<|"a" -> 1, "b" -> 2, "c" -> 3|>, {"c", "a"}]
Out[1]= <|"c" -> 3, "a" -> 1|>

In[2]:= KeyTake[{<|"a" -> 1, "b" -> 2|>, <|"a" -> 3, "b" -> 4|>}, {"a"}]
Out[2]= {<|"a" -> 1|>, <|"a" -> 3|>}
```

## Implementation notes

- Matches Mathematica 15: absent keys are skipped, and a key requested more
  than once is placed at its last occurrence (`{"c", "a", "c"}` gives
  `<|"a" -> .., "c" -> ..|>`).
- `O(n + m)`: the association's entries are hash-indexed once; `RuleDelayed`
  entries stay delayed.

**Attributes:** `Protected`.

## References

**See also:** [RuleDelayed](../../assignment-and-rules/RuleDelayed/)

- Source: [`src/assoc.c`](https://github.com/stblake/mathilda/blob/main/src/assoc.c)
- Specification: [`docs/spec/builtins/data-structures.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/data-structures.md)
- Tests: [`tests/test_assoc_forms.c`](https://github.com/stblake/mathilda/blob/main/tests/test_assoc_forms.c)
- Tests: [`tests/test_assoc_read.c`](https://github.com/stblake/mathilda/blob/main/tests/test_assoc_read.c)
- Tests: [`tests/test_association.c`](https://github.com/stblake/mathilda/blob/main/tests/test_association.c)
- Tests: [`tests/test_compile_assoc.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compile_assoc.c)
