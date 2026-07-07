# MaximalBy

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
MaximalBy[list, f]
    Gives the element(s) of list for which f is maximal
    (all ties, in order). Over an association, gives the entries whose
    value maximises f. MaximalBy[f] is the operator form.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= MaximalBy[{1, -5, 3, -5, 2}, Abs]
Out[1]= {-5, -5}

In[2]:= MinimalBy[<|"a" -> 1, "b" -> 3, "c" -> 2|>, Identity]
Out[2]= <|"a" -> 1|>
```

## Implementation notes

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/core.c`](https://github.com/stblake/mathilda/blob/main/src/core.c)
- Specification: [`docs/spec/builtins/functional-programming.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/functional-programming.md)
