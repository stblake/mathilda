### Worked examples

With the specialist enabled (the default), a transcendental equation is solved by inverting its head; disabling it returns the equation unevaluated:

```mathematica
In[1]:= Solve[Sin[x] == 1/2, x]
Out[1]= {{x -> ConditionalExpression[2 C[1] Pi + 5/6 Pi, Element[C[1], Integers]]}, {x -> ConditionalExpression[2 C[1] Pi + 1/6 Pi, Element[C[1], Integers]]}}

In[2]:= Solve[Sin[x] == 1/2, x, InverseFunctions -> False]
Out[2]= Solve[Sin[x] == 1/2, x, InverseFunctions -> False]
```

Inverting `Exp` yields the full set of complex logarithm branches:

```mathematica
In[1]:= Solve[Exp[x] == 5, x]
Out[1]= {{x -> ConditionalExpression[Log[5] + (2*I) C[1] Pi, Element[C[1], Integers]]}}
```

### Notes

`InverseFunctions` is a `Solve` option (default `Automatic`, i.e. enabled) that
controls the inverse-function specialist for invertible heads such as `Sin`,
`Log`, `Exp`, and integer `Power`. With it on, `Sin[x] == 1/2` is inverted to
the full periodic solution set; setting `InverseFunctions -> False` disables the
specialist, so the equation is returned unevaluated.
