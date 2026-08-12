# KeyExistsQ

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`KeyExistsQ[assoc, key]`**

Gives True if key is present in assoc, else False.

## Examples (2)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (2)

```mathematica
In[1]:= KeyExistsQ[<|"a" -> 1|>, "a"]
Out[1]= True

In[2]:= KeyFreeQ[<|"a" -> 1|>, "b"]
Out[2]= True
```

## Implementation notes

**Attributes:** `Protected`.

## See also

[KeyMemberQ](../../data-structures/KeyMemberQ/), [KeyFreeQ](../../data-structures/KeyFreeQ/)

## References

- Source: [`src/assoc.c`](https://github.com/stblake/mathilda/blob/main/src/assoc.c)
- Specification: [`docs/spec/builtins/data-structures.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/data-structures.md)
- Tests: [`tests/test_association.c`](https://github.com/stblake/mathilda/blob/main/tests/test_association.c)
- Tests: [`tests/test_compile_assoc.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compile_assoc.c)
