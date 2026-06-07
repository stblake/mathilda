### Worked examples

```mathematica
In[1]:= BeginPackage["MyPkg`"]
Out[1]= "MyPkg`"

In[2]:= EndPackage[]
Out[2]= Null

In[3]:= $Context
Out[3]= "Global`"

In[4]:= $ContextPath
Out[4]= {"MyPkg`", "Global`", "System`"}
```

### Notes

`EndPackage[]` restores the context state saved by `BeginPackage` and prepends the just-closed package context to `$ContextPath`, so the package's exported symbols remain visible under their short names. It returns `Null`.
