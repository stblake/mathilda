# TrigExpand

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
TrigExpand[expr]
    expands out trigonometric functions in expr.
    TrigExpand operates on both circular and hyperbolic functions.
    TrigExpand splits up sums and integer multiples that appear in arguments of
    trigonometric functions, and then expands out products of trigonometric
    functions into sums of powers, using trigonometric identities when possible.
    TrigExpand automatically threads over lists, as well as equations,
    inequalities, and logic functions.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= TrigExpand[Sin[2 x]]
Out[1]= 2 Cos[x] Sin[x]

In[2]:= TrigExpand[Sin[x + y]]
Out[2]= Cos[x] Sin[y] + Sin[x] Cos[y]

In[3]:= TrigExpand[Sin[3 x]]
Out[3]= -Sin[x]^3 + 3 Cos[x]^2 Sin[x]

In[4]:= TrigExpand[Cos[x + y + z]]
Out[4]= Cos[x] Cos[y] Cos[z] - Cos[x] Sin[y] Sin[z] - Sin[x] Cos[y] Sin[z] - Sin[x] Sin[y] Cos[z]

In[5]:= TrigExpand[Sin[x]^2 + Cos[x]^2]
Out[5]= 1

In[6]:= TrigExpand[Sinh[4 x]]
Out[6]= 4 Cosh[x] Sinh[x]^3 + 4 Cosh[x]^3 Sinh[x]

In[7]:= TrigExpand[Cosh[x - y]]
Out[7]= Cosh[x] Cosh[y] - Sinh[x] Sinh[y]

In[8]:= TrigExpand[Tanh[2 t]]
Out[8]= 2 Cosh[t] Sinh[t] Sech[2 t]
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
