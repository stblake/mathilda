### Worked examples

```mathematica
In[1]:= Solve[Exp[x] == 2, x]
Out[1]= {{x -> ConditionalExpression[Log[2] + (2*I) C[1] Pi, Element[C[1], Integers]]}}

In[2]:= Solve[Exp[x] == 2, x, GeneratedParameters -> K]
Out[2]= {{x -> ConditionalExpression[Log[2] + (2*I) K[1] Pi, Element[K[1], Integers]]}}
```

```mathematica
In[1]:= Solve[Sin[x] == 0, x, GeneratedParameters -> n]
Out[1]= {{x -> ConditionalExpression[Pi + 2 Pi n[1], Element[n[1], Integers]]}, {x -> ConditionalExpression[2 Pi n[1], Element[n[1], Integers]]}}
```

### Notes

`GeneratedParameters` is a `Solve` option that names the head used for the
fresh integer-parameter symbols the inverse-function solver introduces for
equations with infinitely many solutions. The default `C` produces
`C[1], C[2], ...`; here `GeneratedParameters -> K` switches them to `K[1]`. The
same indexed parameter labels every branch of a periodic solution set:
`Solve[Sin[x] == 0, x, GeneratedParameters -> n]` returns the two families
`2 Pi n[1]` and `Pi + 2 Pi n[1]`, each carrying its own `Element[n[1], Integers]`
constraint.
