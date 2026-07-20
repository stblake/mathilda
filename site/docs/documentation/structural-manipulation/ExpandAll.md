# ExpandAll

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
ExpandAll[expr]
    expands out all products and integer powers in any part of expr.
ExpandAll[expr, patt]
    avoids expanding parts of expr that do not contain terms matching patt.
ExpandAll effectively maps Expand and ExpandDenominator onto every part of expr,
    including function heads, arguments, exponents, and denominators.
ExpandAll automatically threads over lists, as well as equations,
    inequalities, and logic functions.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= ExpandAll[1/(1+x)^3 + Sin[(1+x)^3]]
Out[1]= 1/(1 + 3 x + 3 x^2 + x^3) + Sin[1 + 3 x + 3 x^2 + x^3]

In[2]:= ExpandAll[(x+z)^2/(x+y)^2]
Out[2]= x^2/(x^2 + 2 x y + y^2) + 2 (x z)/(x^2 + 2 x y + y^2) + z^2/(x^2 + 2 x y + y^2)

In[3]:= ExpandAll[E^(I a (t-b))]
Out[3]= E^(-I a b + I a t)

In[4]:= ExpandAll[((1+a) (1+b))[x]]
Out[4]= (1 + a + b + a b)[x]

In[5]:= ExpandAll[(f[(x+y)^2] + g[(y+z)^2])^2, x]
Out[5]= f[x^2 + 2 x y + y^2]^2 + g[(y + z)^2]^2 + 2 f[x^2 + 2 x y + y^2] g[(y + z)^2]
```

## Implementation notes

- `Protected`.
- A thin recursive driver over the accelerated `Expand`: it descends into every

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/structural-manipulation.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/structural-manipulation.md)
