### Worked examples

```mathematica
In[1]:= Put[x^2 + 1, "/tmp/mathilda_demo.m"]
Out[1]= Null

In[2]:= FilePrint["/tmp/mathilda_demo.m"]
1 + x^2
Out[2]= Null
```

### Notes

`Put[expr, "file"]` (equivalently `expr >> "file"`) writes `expr` to the file, replacing any prior contents. Read it back with `Get`.
