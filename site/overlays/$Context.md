### Worked examples

```mathematica
In[1]:= $Context
Out[1]= "Global`"

In[2]:= Begin["foo`"]
Out[2]= "foo`"

In[3]:= $Context
Out[3]= "foo`"
```

### Notes

`$Context` is the string giving the current context; new unqualified symbols are created here. It is read-only for direct assignment — change it through `Begin`, `BeginPackage`, `End`, and `EndPackage`.
