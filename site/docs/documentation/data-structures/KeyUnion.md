# KeyUnion

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`KeyUnion[{assoc1, assoc2, ...}]`**

Gives the list of associations padded to the union of all their keys; a key absent from an association is filled with Missing\["KeyAbsent", key\].

**`KeyUnion[{assoc1, assoc2, ...}, f]`**

Fills each absent key k with f\[k\] instead.

## Examples (1)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (1)

```mathematica
In[1]:= KeyUnion[{<|"a" -> 1|>, <|"b" -> 2|>}, 0 &]
Out[1]= {<|"a" -> 1, "b" -> 0|>, <|"a" -> 0, "b" -> 2|>}
```

## Implementation notes

**Attributes:** `Protected`.

## References

**See also:** [Missing](../../data-structures/Missing/)

- Source: [`src/assoc_ops.c`](https://github.com/stblake/mathilda/blob/main/src/assoc_ops.c)
- Specification: [`docs/spec/builtins/data-structures.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/data-structures.md)
- Tests: [`tests/test_assoc_ops.c`](https://github.com/stblake/mathilda/blob/main/tests/test_assoc_ops.c)
- Tests: [`tests/test_association.c`](https://github.com/stblake/mathilda/blob/main/tests/test_association.c)
