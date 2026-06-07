### Worked examples

```mathematica
In[1]:= Put[x^2 + 1, "/tmp/mathilda_demo.m"]
Out[1]= Null

In[2]:= FilePrint["/tmp/mathilda_demo.m"]
1 + x^2
Out[2]= Null
```

### Notes

`FilePrint["file"]` writes the raw textual contents of the file to stdout and returns `Null`; it does not parse or evaluate the file (use `Get` for that).
