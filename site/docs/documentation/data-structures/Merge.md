# Merge

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Merge[{assoc1, assoc2, ...}, f]`**

Combines associations, applying f to the list of values collected for each key (e.g. Merge\[{...}, Total\]). The list may also hold rules and lists of rules, each rule contributing one value.

## Examples (3)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (3)

```mathematica
In[1]:= Merge[{<|"a" -> 1|>, <|"a" -> 2, "b" -> 3|>}, Total]
Out[1]= <|"a" -> 3, "b" -> 3|>

In[2]:= Merge[<|"g1" -> <|"a" -> 1|>, "g2" -> <|"a" -> 2, "b" -> 3|>|>, Total]
Out[2]= <|"a" -> 3, "b" -> 3|>

In[3]:= Merge[{"a" -> 1, "b" -> 2, "a" -> 3}, Total]
Out[3]= <|"a" -> 4, "b" -> 2|>
```

## Implementation notes

**Attributes:** `Protected`.

## References

- Source: [`src/assoc_ops.c`](https://github.com/stblake/mathilda/blob/main/src/assoc_ops.c)
- Specification: [`docs/spec/builtins/data-structures.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/data-structures.md)
- Tests: [`tests/test_assoc_forms.c`](https://github.com/stblake/mathilda/blob/main/tests/test_assoc_forms.c)
- Tests: [`tests/test_assoc_ops.c`](https://github.com/stblake/mathilda/blob/main/tests/test_assoc_ops.c)
- Tests: [`tests/test_association.c`](https://github.com/stblake/mathilda/blob/main/tests/test_association.c)
