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
