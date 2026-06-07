### Worked examples

```mathematica
In[1]:= Print["Hello, Mathilda!"]
"Hello, Mathilda!"
Out[1]= Null

In[2]:= Print[2 + 3]
5
Out[2]= Null

In[3]:= Print["x = ", 2^10]
"x = "1024
Out[3]= Null
```

### Notes

`Print` writes its evaluated arguments to stdout, concatenated with no separator, then returns `Null`. String arguments are printed with their surrounding quotes.
