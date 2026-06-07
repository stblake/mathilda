# ExpToTrig

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
ExpToTrig[expr]
    rewrites exponentials and logarithms in expr in terms of circular and
    hyperbolic trigonometric functions when possible.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= ExpToTrig[Exp[I x]]
Out[1]= Cos[x] + I Sin[x]

In[2]:= ExpToTrig[Log[1 + I x] - Log[1 - I x]]
Out[2]= (2*I) ArcTan[x]

In[3]:= ExpToTrig[Exp[I x] == -1]
Out[3]= Cos[x] + I Sin[x] == -1
```

## Implementation notes

- `Listable`, `Protected`.
- Tries when possible to give results that do not involve explicit complex numbers.
- ExpToTrig is natively the inverse of `TrigToExp`.
- Automatically threads over lists, equations, inequalities, and logic functions.

**Attributes:** `Listable`, `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/elementary-functions.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/elementary-functions.md)
