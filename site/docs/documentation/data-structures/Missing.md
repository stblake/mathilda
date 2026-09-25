# Missing

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Missing[]  |  Missing["reason"]  |  Missing["reason", data]`**

Represents missing data. Lookup and key access give Missing\["KeyAbsent", k\] for an absent key; KeyUnion and JoinAcross fill gaps with Missing\[...\]. Test with MissingQ, remove with DeleteMissing.

## Examples (3)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (3)

```mathematica
In[1]:= <|"a" -> 1|>["q"]
Out[1]= Missing["KeyAbsent", "q"]

In[2]:= Missing["NotAvailable"]
Out[2]= Missing["NotAvailable"]

In[3]:= DeleteMissing[{1, Missing["NotAvailable"], 3}]
Out[3]= {1, 3}
```

## Implementation notes

**Attributes:** `Protected`.

## References

**See also:** [KeyUnion](../../data-structures/KeyUnion/), [JoinAcross](../../data-structures/JoinAcross/), [MissingQ](../../data-structures/MissingQ/), [DeleteMissing](../../data-structures/DeleteMissing/)

- Source: [`src/assoc_ops.c`](https://github.com/stblake/mathilda/blob/main/src/assoc_ops.c)
- Specification: [`docs/spec/builtins/data-structures.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/data-structures.md)
- Tests: [`tests/test_assoc_forms.c`](https://github.com/stblake/mathilda/blob/main/tests/test_assoc_forms.c)
- Tests: [`tests/test_assoc_ops.c`](https://github.com/stblake/mathilda/blob/main/tests/test_assoc_ops.c)
- Tests: [`tests/test_assoc_read.c`](https://github.com/stblake/mathilda/blob/main/tests/test_assoc_read.c)
- Tests: [`tests/test_association.c`](https://github.com/stblake/mathilda/blob/main/tests/test_association.c)
