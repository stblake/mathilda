### Worked examples

```mathematica
In[1]:= Put[x^2 + 1, "/tmp/mathilda_demo.m"]
Out[1]= Null

In[2]:= Get["/tmp/mathilda_demo.m"]
Out[2]= 1 + x^2
```

### Notes

`Get["file"]` reads the expressions in a file, evaluates them in order, and returns the last result.
