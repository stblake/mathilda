# DeleteMissing

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`DeleteMissing[expr]`**

Removes all Missing\[...\] elements (equivalent to DeleteCases\[expr, \_Missing\]). Over an association, drops entries whose value is Missing\[...\].

## Examples (1)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (1)

```mathematica
In[1]:= DeleteMissing[Lookup[<|"a" -> 1, "b" -> 2|>, {"a", "z", "b"}]]
Out[1]= {1, 2}
```

## Implementation notes

**Attributes:** `Protected`.

## See also

[Lookup](../../data-structures/Lookup/)

## References

- Source: [`src/patterns.c`](https://github.com/stblake/mathilda/blob/main/src/patterns.c)
- Specification: [`docs/spec/builtins/data-structures.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/data-structures.md)
- Tests: [`tests/test_association.c`](https://github.com/stblake/mathilda/blob/main/tests/test_association.c)
