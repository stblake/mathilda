### Worked examples

```mathematica
In[1]:= x = 5
Out[1]= 5

In[2]:= x + 1
Out[2]= 6

In[3]:= Clear[x]
Out[3]= Null

In[4]:= x + 1
Out[4]= 1 + x
```

### Notes

`Clear[s]` removes all OwnValues and DownValues attached to `s`, so the symbol becomes undefined again (attributes and the symbol itself are left intact).
