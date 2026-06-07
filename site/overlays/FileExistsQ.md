### Worked examples

```mathematica
In[1]:= Put[x^2 + 1, "/tmp/mathilda_demo.m"]
Out[1]= Null

In[2]:= FileExistsQ["/tmp/mathilda_demo.m"]
Out[2]= True

In[3]:= FileExistsQ["/tmp/does_not_exist.m"]
Out[3]= False
```

### Notes

`FileExistsQ["name"]` returns `True` if a file with the given name exists and `False` otherwise.
