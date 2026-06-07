### Worked examples

```mathematica
In[1]:= 10/4
Out[1]= 5/2

In[2]:= Divide[10, 4]
Out[2]= 5/2

In[3]:= a/b
Out[3]= a/b
```

### Notes

`x / y` is rewritten by the evaluator to `Times[x, Power[y, -1]]`, so it inherits Times's flattening and ordering; integer quotients auto-reduce to an exact `Rational`.
