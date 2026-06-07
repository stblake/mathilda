### Worked examples

```mathematica
In[1]:= (s = 0; Do[s = s + i, {i, 5}]; s)
Out[1]= 15

In[2]:= (p = 1; Do[p = p*2, {3}]; p)
Out[2]= 8
```

### Notes

`Do[expr, {i, imax}]` runs `expr` with `i` taking values `1` through `imax`; `Do[expr, n]` simply repeats `expr` `n` times. `Do` returns `Null`, so it is used for side effects — read the accumulated value out of a variable afterwards.
