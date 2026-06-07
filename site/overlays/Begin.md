### Worked examples

```mathematica
In[1]:= $Context
Out[1]= "Global`"

In[2]:= Begin["foo`"]
Out[2]= "foo`"

In[3]:= $Context
Out[3]= "foo`"

In[4]:= End[]
Out[4]= "foo`"

In[5]:= $Context
Out[5]= "Global`"
```

### Notes

`Begin["ctx`"]` switches `$Context` to `ctx`` and saves the previous context so the matching `End[]` can restore it. An argument that starts with a backtick is taken relative to the current context (e.g. `Begin["`Private`"]` inside `MyPkg`` yields `MyPkg`Private``).
