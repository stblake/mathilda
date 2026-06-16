### Worked examples

```mathematica
In[1]:= ToExpression["1 + 2*3"]
Out[1]= 7
```

```mathematica
In[1]:= ToExpression["D[ArcTan[x], x]"]
Out[1]= 1/(1 + x^2)
```

```mathematica
In[1]:= ToExpression["Series[Exp[x], {x, 0, 4}]"]
Out[1]= 1 + x + 1/2 x^2 + 1/6 x^3 + 1/24 x^4 + O[x]^5
```

```mathematica
In[1]:= ToExpression["x^2 + 1", InputForm, Hold]
Out[1]= Hold[x^2 + 1]
```

```mathematica
In[1]:= ToExpression["bad syntax ]["]
Out[1]= $Failed
```

### Notes

`ToExpression[input]` parses a string as Mathilda input and evaluates the result
— so the full symbolic engine is reachable from text, including `D`, `Series`,
`Solve`, and friends. The three-argument form wraps a head around the parsed tree
before evaluation; `ToExpression[s, InputForm, Hold]` is the standard way to get
the *unevaluated* parsed form. A syntax error yields `$Failed` rather than an
abort, making `ToExpression` safe to use on untrusted input or in `Map`.
