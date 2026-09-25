# Values

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Values[assoc]`**

Gives a list of the values of an association (or of a rule or list of rules).

**`Values[{assoc1, assoc2, ...}]`**

Threads over lists of associations and lists of rules.

**`Values[assoc, f]`**

Wraps each value: {f\[v1\], f\[v2\], ...}.

## Examples (3)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (3)

```mathematica
In[1]:= Values[<|"a" -> 1, "b" -> 2|>]
Out[1]= {1, 2}

In[2]:= Values[<|"a" -> 1, "b" -> 2|>, f]
Out[2]= {f[1], f[2]}

In[3]:= Values[{<|"a" -> 1|>, {"b" -> 2}}]
Out[3]= {{1}, {2}}
```

## Implementation notes

**Attributes:** `Protected`.

## References

**See also:** [Keys](../../data-structures/Keys/)

- Source: [`src/assoc_ops.c`](https://github.com/stblake/mathilda/blob/main/src/assoc_ops.c)
- Specification: [`docs/spec/builtins/data-structures.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/data-structures.md)
- Tests: [`tests/test_assoc_ops.c`](https://github.com/stblake/mathilda/blob/main/tests/test_assoc_ops.c)
- Tests: [`tests/test_association.c`](https://github.com/stblake/mathilda/blob/main/tests/test_association.c)
- Tests: [`tests/test_compile_assoc.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compile_assoc.c)
