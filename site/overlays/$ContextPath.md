### Worked examples

```mathematica
In[1]:= $ContextPath
Out[1]= {"Global`", "System`"}

In[2]:= BeginPackage["MyPkg`"]
Out[2]= "MyPkg`"

In[3]:= $ContextPath
Out[3]= {"MyPkg`", "System`"}
```

### Notes

`$ContextPath` is the ordered list of contexts searched to resolve bare identifiers to existing qualified symbols. It is read-only for direct assignment and is modified by `BeginPackage` and `EndPackage`.
