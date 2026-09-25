# PositionIndex

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`PositionIndex[list]`**

Gives \<|value -\> {positions}|\> mapping each distinct element to the list of 1-based positions where it occurs. O(n).

**`PositionIndex[assoc]`**

Gives \<|value -\> {keys}|\>, the keys at which each value occurs.

## Examples (2)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (2)

```mathematica
In[1]:= PositionIndex[{a, b, a, c, a, b}]
Out[1]= <|a -> {1, 3, 5}, b -> {2, 6}, c -> {4}|>

In[2]:= PositionIndex[<|"a" -> x, "b" -> y, "c" -> x|>]
Out[2]= <|x -> {"a", "c"}, y -> {"b"}|>
```

## Implementation notes

**Attributes:** `Protected`.

## References

- Source: [`src/assoc_ops.c`](https://github.com/stblake/mathilda/blob/main/src/assoc_ops.c)
- Specification: [`docs/spec/builtins/data-structures.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/data-structures.md)
- Tests: [`tests/test_assoc_ops.c`](https://github.com/stblake/mathilda/blob/main/tests/test_assoc_ops.c)
- Tests: [`tests/test_association.c`](https://github.com/stblake/mathilda/blob/main/tests/test_association.c)
