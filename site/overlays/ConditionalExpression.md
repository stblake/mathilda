### Worked examples

```mathematica
In[1]:= ConditionalExpression[x, True]
Out[1]= x
```

```mathematica
In[1]:= ConditionalExpression[1/x, x != 0] /. x -> 0
Out[1]= Undefined
```

```mathematica
In[1]:= ConditionalExpression[Sqrt[x^2], x > 0]
Out[1]= ConditionalExpression[Sqrt[x^2], x > 0]
```

```mathematica
In[1]:= ConditionalExpression[ConditionalExpression[e, a > 0], b > 0]
Out[1]= ConditionalExpression[e, a > 0 && b > 0]
```

### Notes

`ConditionalExpression[expr, cond]` carries a value together with the assumption under which it is valid. It evaluates to `expr` when the condition is `True` and to `Undefined` when it is `False`; nested forms merge their conditions with `And`.
