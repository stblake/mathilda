---
status: Stable
references:
---
### Worked examples

```mathematica
In[1]:= Simplify[(x^2 - 1)/(x - 1)]
Out[1]= 1 + x
```

```mathematica
In[1]:= Simplify[Sin[x]^2 + Cos[x]^2]
Out[1]= 1
```

```mathematica
In[1]:= Simplify[x + x + x]
Out[1]= 3 x
```

```mathematica
In[1]:= Simplify[Sqrt[x^2], x > 0]
Out[1]= x

In[2]:= Simplify[Sqrt[x^2], Element[x, Reals]]
Out[2]= Abs[x]
```

```mathematica
In[1]:= Simplify[Cosh[x]^2 - Sinh[x]^2]
Out[1]= 1

In[2]:= Simplify[Log[a b] - Log[a] - Log[b], {a > 0, b > 0}]
Out[2]= 0

In[3]:= Simplify[Cos[3 x]/Cos[x] - (2 Cos[2 x] - 1)]
Out[3]= 0
```

### Notes

`Simplify` tries a collection of transformations — `Together`, `Cancel`,
`Expand`, `Factor`, `Apart`, `TrigExpand`, `TrigFactor`, and a `TrigToExp` round
trip — and keeps the smallest result, so it can cancel `(x^2-1)/(x-1)` and reduce
the Pythagorean identity to `1`. A second argument supplies assumptions:
`Simplify[Sqrt[x^2], x > 0]` uses `x > 0` to drop the absolute value and return
`x`. Assumptions may be equations, inequalities, or domain statements like
`Element[x, Integers]`. `Simplify` threads automatically over lists, equations,
and logical combinations.
