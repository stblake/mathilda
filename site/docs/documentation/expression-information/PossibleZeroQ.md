# PossibleZeroQ

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
PossibleZeroQ[expr] gives True if symbolic and numerical methods suggest that expr has value zero, and False otherwise.
The general problem of deciding whether an expression is zero is undecidable; PossibleZeroQ is a quick but not always accurate test.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= PossibleZeroQ[E^(I Pi/4) - (-1)^(1/4)]
Out[1]= True

In[2]:= PossibleZeroQ[(x + 1)(x - 1) - x^2 + 1]
Out[2]= True

In[3]:= PossibleZeroQ[(E + Pi)^2 - E^2 - Pi^2 - 2 E Pi]
Out[3]= True

In[4]:= PossibleZeroQ[E^Pi - Pi^E]
Out[4]= False

In[5]:= PossibleZeroQ[2^(2 I) - 2^(-2 I) - 2 I Sin[Log[4]]]
Out[5]= True

In[6]:= PossibleZeroQ[Sqrt[x^2] - x]
Out[6]= False

In[7]:= PossibleZeroQ[Sin[x]^2 + Cos[x]^2 - 1]
Out[7]= True
```

## Implementation notes

**Attributes:** `Listable`, `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/expression-information.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/expression-information.md)
