### Worked examples

```mathematica
In[1]:= Array[f,5]
Out[1]= {f[1], f[2], f[3], f[4], f[5]}

In[2]:= Array[#^2&,4]
Out[2]= {1, 4, 9, 16}

In[3]:= Array[f,3,0]
Out[3]= {f[0], f[1], f[2]}
```

### Notes

`Array[f, n]` builds a length-`n` list by applying `f` to indices `1..n`; an optional third argument sets the starting index.
