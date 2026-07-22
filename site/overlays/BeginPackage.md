### Worked examples

```mathematica
In[1]:= $Context
Out[1]= "Global`"

In[2]:= BeginPackage["MyPkg`"]
Out[2]= "MyPkg`"

In[3]:= $Context
Out[3]= "MyPkg`"

In[4]:= $ContextPath
Out[4]= {"MyPkg`", "System`"}
```

### Notes

`BeginPackage["ctx`"]` sets `$Context` to `ctx`` and restricts `$ContextPath` to `{"ctx`", "System`"}`, matching the standard package prologue so only system symbols and the package's own symbols resolve under short names. Use the matching `EndPackage[]` to close it.
