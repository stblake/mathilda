### Worked examples

```mathematica
In[1]:= Put[x^2 + 1, "/tmp/mathilda_demo.m"]
Out[1]= Null

In[2]:= PutAppend[y, "/tmp/mathilda_demo.m"]
Out[2]= Null

In[3]:= FilePrint["/tmp/mathilda_demo.m"]
1 + x^2
y
Out[3]= Null
```

### Notes

`PutAppend[expr, "file"]` (equivalently `expr >>> "file"`) adds `expr` on a new line at the end of the file, creating the file if it does not exist.
