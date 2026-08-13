# MinimalBy

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`MinimalBy[list, f]`**

Gives the element(s) of list for which f is minimal (all ties, in order). Over an association, gives the entries whose value minimises f. MinimalBy\[f\] is the operator form.

## Examples (2)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (2)

```mathematica
In[1]:= MaximalBy[{1, -5, 3, -5, 2}, Abs]
Out[1]= {-5, -5}

In[2]:= MinimalBy[<|"a" -> 1, "b" -> 3, "c" -> 2|>, Identity]
Out[2]= <|"a" -> 1|>
```

## Implementation notes

**Attributes:** `Protected`.

## References

**See also:** [MaximalBy](../../functional-programming/MaximalBy/)

- Source: [`src/core.c`](https://github.com/stblake/mathilda/blob/main/src/core.c)
- Specification: [`docs/spec/builtins/functional-programming.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/functional-programming.md)
- Tests: [`tests/test_association.c`](https://github.com/stblake/mathilda/blob/main/tests/test_association.c)
- Tests: [`tests/test_list.c`](https://github.com/stblake/mathilda/blob/main/tests/test_list.c)
