# Splice

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Splice[{e1, e2, ...}]`**

Is replaced by the sequence e1, e2, ... when it appears inside a List or an Association.

**`Splice[{e1, e2, ...}, h]`**

Splices into any head matching the pattern h.

## Examples (5)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (5)

```mathematica
In[1]:= {1, Splice[{2, 3}], 4}
Out[1]= {1, 2, 3, 4}

In[2]:= <|"a" -> 1, Splice[{"b" -> 2, "c" -> 3}]|>
Out[2]= <|"a" -> 1, "b" -> 2, "c" -> 3|>

In[3]:= Table[Splice[{i, -i}], {i, 3}]
Out[3]= {1, -1, 2, -2, 3, -3}

In[4]:= f[1, Splice[{2, 3}]]
Out[4]= f[1, Splice[{2, 3}]]

In[5]:= f[1, Splice[{2, 3}, _]]
Out[5]= f[1, 2, 3]
```

## Implementation notes

**Attributes:** `Protected`.

## References

**See also:** [List](../../other-advanced/List/), [Association](../../data-structures/Association/), [SequenceHold](../../expression-information/SequenceHold/), [HoldAllComplete](../../expression-information/HoldAllComplete/)

- Source: [`src/assoc_ops.c`](https://github.com/stblake/mathilda/blob/main/src/assoc_ops.c)
- Specification: [`docs/spec/builtins/data-structures.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/data-structures.md)
- Tests: [`tests/test_assoc_ops.c`](https://github.com/stblake/mathilda/blob/main/tests/test_assoc_ops.c)
