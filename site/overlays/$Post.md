### Worked examples

```mathematica
In[1]:= $Post = Framed
Out[1]= Framed[Framed]

In[2]:= 2 + 2
Out[2]= Framed[4]
```

### Notes

`$Post`, if set, is applied to every output expression after evaluation. Here
`$Post = Framed` wraps each result in `Framed[...]` (note the assignment's own
echo is also wrapped). Unset by default, in which case results are shown
unmodified.
