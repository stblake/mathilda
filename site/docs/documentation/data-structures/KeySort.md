# KeySort

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`KeySort[assoc]`**

Sorts an association into canonical key order.

**`KeySort[assoc, p]`**

Sorts the entries by key using the ordering function p.

## Examples (2)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (2)

```mathematica
In[1]:= KeySort[<|"c" -> 3, "a" -> 1, "b" -> 2|>]
Out[1]= <|"a" -> 1, "b" -> 2, "c" -> 3|>

In[2]:= KeySort[<|"c" -> 1, "a" -> 2, "b" -> 3|>, -Order[#1, #2] &]
Out[2]= <|"c" -> 1, "b" -> 3, "a" -> 2|>
```

## Implementation notes

- `KeySort[assoc, p]` orders the keys exactly as `Sort[Keys[assoc], p]` would
  (the shared merge sort), so a `p` that does not evaluate to `False`/`-1`
  leaves the order alone, as in Mathematica 15.

**Attributes:** `Protected`.

## References

- Source: [`src/assoc.c`](https://github.com/stblake/mathilda/blob/main/src/assoc.c)
- Specification: [`docs/spec/builtins/data-structures.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/data-structures.md)
- Tests: [`tests/test_assoc_read.c`](https://github.com/stblake/mathilda/blob/main/tests/test_assoc_read.c)
- Tests: [`tests/test_association.c`](https://github.com/stblake/mathilda/blob/main/tests/test_association.c)
