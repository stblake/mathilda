# TrigReduce

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
TrigReduce[expr]
    rewrites products and powers of trigonometric functions in expr in
    terms of trigonometric functions with combined arguments.
TrigReduce operates on both circular and hyperbolic functions; given a
trigonometric polynomial it typically yields a linear expression
involving trigonometric functions with more complicated arguments
(broadly the inverse of TrigExpand).
TrigReduce automatically threads over lists, equations, inequalities,
and logic functions.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= TrigReduce[2 Cos[x]^2]
Out[1]= 1 + Cos[2 x]

In[2]:= TrigReduce[2 Sin[x] Cos[y]]
Out[2]= Sin[x + y] + Sin[x - y]

In[3]:= TrigReduce[2 Cosh[x] Cosh[y]]
Out[3]= Cosh[x + y] + Cosh[x - y]

In[4]:= TrigReduce[Sin[a] (Cos[b] - Sin[b]) + Cos[a] (Sin[b] + Cos[b])]
Out[4]= Cos[a + b] + Sin[a + b]

In[5]:= TrigReduce[Tan[x] + Tan[y]]
Out[5]= Sec[x] Sec[y] Sin[x + y]

In[6]:= TrigReduce[Coth[x] + Coth[y]]
Out[6]= Csch[x] Csch[y] Sinh[x + y]

In[7]:= TrigReduce[Sin[x]^4]
Out[7]= 1/8 (3 + Cos[4 x] - 4 Cos[2 x])

In[8]:= TrigReduce[2 Sin[x + y] Cos[x - y]]
Out[8]= Sin[2 x] + Sin[2 y]
```

## Implementation notes

- Applies the classical product-to-sum identities (Sin·Cos, Sin·Sin,

**Attributes:** `Listable`, `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/elementary-functions.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/elementary-functions.md)
