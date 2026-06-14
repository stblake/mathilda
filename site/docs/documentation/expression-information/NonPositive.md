# NonPositive

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
NonPositive[x]
    gives True if x is a real number that is negative or zero, and False
if x is a manifestly positive real number or a non-real complex number.
For non-numeric x the expression is left unevaluated. NonPositive is
Listable, so it threads over lists element by element.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= NonPositive[{1.6, 3/4, Pi, 0, -5, 1 + I, Sin[10^5]}]
Out[1]= {False, False, False, True, True, False, False}

In[2]:= NonPositive[{x, Sin[y]}]
Out[2]= {NonPositive[x], NonPositive[Sin[y]]}

In[3]:= NonPositive[1 - Pi]
Out[3]= True
```

## Implementation notes

**Attributes:** `Listable`, `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/expression-information.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/expression-information.md)
