# TrigFactor

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
TrigFactor[expr]
    factors trigonometric functions in expr.
    TrigFactor operates on both circular and hyperbolic functions.
    TrigFactor factors polynomials in trigonometric functions and collapses
    Pythagorean, angle-addition, and double-angle identities where possible,
    broadly acting as the inverse of TrigExpand.
    TrigFactor automatically threads over lists, as well as equations,
    inequalities, and logic functions.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= TrigFactor[Sin[x]^2 + Cos[x]^2]
Out[1]= 1

In[2]:= TrigFactor[Cosh[x]^2 - Sinh[x]^2]
Out[2]= 1

In[3]:= TrigFactor[2 Sin[x] Cos[x]]
Out[3]= Sin[2 x]

In[4]:= TrigFactor[Cos[x]^2 - Sin[x]^2]
Out[4]= Cos[2 x]

In[5]:= TrigFactor[Sin[a] Cos[b] + Cos[a] Sin[b]]
Out[5]= Sin[a + b]

In[6]:= TrigFactor[Cos[a] Cos[b] + Sin[a] Sin[b]]
Out[6]= Cos[a - b]

In[7]:= TrigFactor[Sin[x]^2 + Tan[x]^2]
Out[7]= (1 + Cos[x]^2) Tan[x]^2

In[8]:= TrigFactor[Cosh[x]^2 - Cosh[x]^4]
Out[8]= -Cosh[x]^2 Sinh[x]^2
```

## Implementation notes

- `Listable`, `Protected`.
- Operates on both circular (`Sin`, `Cos`, `Tan`, `Cot`, `Sec`, `Csc`) and

**Attributes:** `Listable`, `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/elementary-functions.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/elementary-functions.md)
