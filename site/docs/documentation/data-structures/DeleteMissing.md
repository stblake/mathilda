# DeleteMissing

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`DeleteMissing[expr]`**

Removes all Missing\[...\] elements (equivalent to DeleteCases\[expr, \_Missing\]). Over an association, drops entries whose value is Missing\[...\].

**`DeleteMissing[expr, n]`**

Removes Missing\[...\] elements at levels 1 through n (n may be Infinity); association values count as one level down.

**`DeleteMissing[expr, n, d]`**

Removes the elements at levels 1..n that contain a Missing\[...\] at depth d or less.

## Examples (3)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (3)

```mathematica
In[1]:= DeleteMissing[Lookup[<|"a" -> 1, "b" -> 2|>, {"a", "z", "b"}]]
Out[1]= {1, 2}

In[2]:= DeleteMissing[{1, {Missing[], 2}}, 2]
Out[2]= {1, {2}}

In[3]:= DeleteMissing[{{1, Missing[]}, {2}}, 1, 1]
Out[3]= {{2}}
```

## Implementation notes

**Attributes:** `Protected`.

## References

**See also:** [Lookup](../../data-structures/Lookup/)

- Source: [`src/assoc_ops.c`](https://github.com/stblake/mathilda/blob/main/src/assoc_ops.c)
- Specification: [`docs/spec/builtins/data-structures.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/data-structures.md)
- Tests: [`tests/test_assoc_ops.c`](https://github.com/stblake/mathilda/blob/main/tests/test_assoc_ops.c)
- Tests: [`tests/test_association.c`](https://github.com/stblake/mathilda/blob/main/tests/test_association.c)
