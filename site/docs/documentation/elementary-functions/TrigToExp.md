# TrigToExp

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
TrigToExp[expr]
    rewrites circular and hyperbolic trigonometric functions (and their
    inverses) in expr in terms of exponentials and logarithms.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= TrigToExp[Cos[x]]
Out[1]= 1/2 E^(-I x) + 1/2 E^(I x)

In[2]:= TrigToExp[ArcSin[x]]
Out[2]= -I Log[I x + Sqrt[1 - x^2]]

In[3]:= TrigToExp[{Sin[x], Cos[x], Tan[x]}]
Out[3]= {(-1/2*I) E^(I x) + (1/2*I) E^(-I x), 1/2 E^(-I x) + 1/2 E^(I x), -I E^(I x)/(E^(-I x) + E^(I x)) + I E^(-I x)/(E^(-I x) + E^(I x))}
```

## Implementation notes

- `Listable`, `Protected`.
- Operates on both circular and hyperbolic functions, and their inverses.
- Automatically threads over lists, equations, inequalities, and logic functions.

**Attributes:** `Listable`, `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/elementary-functions.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/elementary-functions.md)
