### Worked examples

```mathematica
In[1]:= NumericQ[Pi]
Out[1]= True

In[2]:= NumericQ[x]
Out[2]= False

In[3]:= NumericQ[Sqrt[2] + 1]
Out[3]= True
```

### Notes

An expression is numeric if it is an explicit number, a constant such as `Pi`, or a `NumericFunction` whose arguments are all numeric. A bare symbol like `x` is not numeric.
