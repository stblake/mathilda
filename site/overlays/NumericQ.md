### Worked examples

```mathematica
In[1]:= NumericQ[Pi]
Out[1]= True

In[2]:= NumericQ[x]
Out[2]= False
```

A deep tree of constants and transcendental functions is recognised as numeric
without any value being computed:

```mathematica
In[1]:= NumericQ[Gamma[1/2] + Zeta[3]]
Out[1]= True
```

One non-numeric leaf is enough to spoil the whole expression:

```mathematica
In[1]:= NumericQ[x + 1]
Out[1]= False
```

The classification looks through `NumericFunction` heads recursively, so mixed
elementary and special functions of numeric arguments still qualify:

```mathematica
In[1]:= NumericQ[Sin[2] + Log[3]]
Out[1]= True
```

### Notes

An expression is numeric if it is an explicit number, a constant such as `Pi`, or
a `NumericFunction` whose arguments are all numeric. The test is structural and
recursive — it never evaluates the expression to a number — so `Gamma[1/2] +
Zeta[3]` is reported numeric while a single symbolic leaf such as `x` makes the
whole expression non-numeric.
