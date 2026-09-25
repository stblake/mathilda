# MissingQ

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`MissingQ[expr]`**

Gives True if expr has head Missing (Missing\[\], Missing\["reason"\], Missing\["KeyAbsent", k\], ...), and False otherwise.

## Examples (3)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (3)

```mathematica
In[1]:= MissingQ[Missing["KeyAbsent", "z"]]
Out[1]= True

In[2]:= MissingQ[<|"a" -> 1|>["z"]]
Out[2]= True

In[3]:= Select[{1, Missing[], 3}, Not @* MissingQ]
Out[3]= {1, 3}
```

## Implementation notes

**Attributes:** `Protected`.

## References

**See also:** [Missing](../../data-structures/Missing/)

- Source: [`src/assoc_ops.c`](https://github.com/stblake/mathilda/blob/main/src/assoc_ops.c)
- Specification: [`docs/spec/builtins/data-structures.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/data-structures.md)
- Tests: [`tests/test_assoc_ops.c`](https://github.com/stblake/mathilda/blob/main/tests/test_assoc_ops.c)
