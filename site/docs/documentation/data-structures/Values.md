# Values

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Values[assoc]`**

Gives a list of the values of an association (or rules).

## Examples (1)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (1)

```mathematica
In[1]:= Values[<|"a" -> 1, "b" -> 2|>]
Out[1]= {1, 2}
```

## Implementation notes

**Attributes:** `Protected`.

## References

- Source: [`src/assoc.c`](https://github.com/stblake/mathilda/blob/main/src/assoc.c)
- Specification: [`docs/spec/builtins/data-structures.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/data-structures.md)
- Tests: [`tests/test_association.c`](https://github.com/stblake/mathilda/blob/main/tests/test_association.c)
- Tests: [`tests/test_compile_assoc.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compile_assoc.c)
