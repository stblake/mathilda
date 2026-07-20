# HoldAll

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
HoldAll
    is an attribute that specifies that all arguments to a function are to be maintained in an unevaluated form.
You can use Evaluate to evaluate the arguments of a HoldAll function in a controlled way.
Even when a function has attribute HoldAll, Sequence objects that appear in its arguments are still by default flattened; Unevaluated wrappers on a held argument are, however, left intact (use HoldComplete/HoldAllComplete to also suppress Sequence flattening).
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= Hold[1+1, 2+2]
Out[1]= Hold[1 + 1, 2 + 2]

In[2]:= Hold[Sequence[a, b], c]
Out[2]= Hold[a, b, c]

In[3]:= Hold[Unevaluated[1+2]]
Out[3]= Hold[Unevaluated[1 + 2]]

In[4]:= Hold[Evaluate[1+2], 3+4]
Out[4]= Hold[3, 3 + 4]
```

## Implementation notes

**Attributes:** none registered.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/expression-information.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/expression-information.md)
