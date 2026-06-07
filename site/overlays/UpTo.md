### Worked examples

```mathematica
In[1]:= Take[{a, b, c, d}, UpTo[2]]
Out[1]= {a, b}

In[2]:= Take[{a, b}, UpTo[5]]
Out[2]= {a, b}
```

### Notes

`UpTo[n]` is a count specification meaning "as many as `n`, but no error if fewer are available." With `Take`, requesting more elements than exist (Out[2]) returns all of them rather than failing.
