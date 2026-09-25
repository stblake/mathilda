# Keys

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Keys[assoc]`**

Gives a list of the keys of an association (or of a rule or list of rules).

**`Keys[{assoc1, assoc2, ...}]`**

Threads over lists of associations and lists of rules.

**`Keys[assoc, f]`**

Wraps each key: {f\[k1\], f\[k2\], ...}.

## Examples (3)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (3)

```mathematica
In[1]:= Keys[<|"a" -> 1, "b" -> 2|>]
Out[1]= {"a", "b"}

In[2]:= Keys[<|"a" -> 1, "b" -> 2|>, f]
Out[2]= {f["a"], f["b"]}

In[3]:= Keys[{<|"a" -> 1|>, <|"b" -> 2, "c" -> 3|>}]
Out[3]= {{"a"}, {"b", "c"}}
```

## Implementation notes

**Attributes:** `Protected`.

## References

- Source: [`src/assoc_ops.c`](https://github.com/stblake/mathilda/blob/main/src/assoc_ops.c)
- Specification: [`docs/spec/builtins/data-structures.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/data-structures.md)
- Tests: [`tests/test_assoc_atomicity.c`](https://github.com/stblake/mathilda/blob/main/tests/test_assoc_atomicity.c)
- Tests: [`tests/test_assoc_ops.c`](https://github.com/stblake/mathilda/blob/main/tests/test_assoc_ops.c)
- Tests: [`tests/test_association.c`](https://github.com/stblake/mathilda/blob/main/tests/test_association.c)
