### Worked examples

```mathematica
In[1]:= $Context
Out[1]= "Global`"

In[2]:= Begin["foo`"]
Out[2]= "foo`"

In[3]:= End[]
Out[3]= "foo`"

In[4]:= $Context
Out[4]= "Global`"
```

### Notes

`End[]` closes the context opened by the matching `Begin[]`, restoring the previously active `$Context` and returning the context it just closed as a string.
