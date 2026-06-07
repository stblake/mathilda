### Worked examples

```mathematica
In[1]:= 3 > 2
Out[1]= True

In[2]:= 5 > 5
Out[2]= False

In[3]:= Greater[5, 3]
Out[3]= True

In[4]:= x > 2
Out[4]= x > 2
```

### Notes

`>` is the operator form of `Greater`. On purely numeric arguments it decides to `True` or `False`; if an argument is symbolic and the relation cannot be settled, the expression is returned unevaluated.
