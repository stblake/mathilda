### Worked examples

```mathematica
In[1]:= Median[{5, 1, 3, 2, 4}]
Out[1]= 3
```

```mathematica
In[1]:= Median[{1, 2, 3, 4}]
Out[1]= 5/2
```

```mathematica
In[1]:= Median[Table[k^2, {k, 1, 10}]]
Out[1]= 61/2
```

### Notes

`Median[data]` returns the middle value of the sorted data. For an odd number of
elements it is the single central element (the first example sorts to
`{1, 2, 3, 4, 5}`, giving `3`); for an even number it is the exact average of the
two central elements, so `Median[{1, 2, 3, 4}]` is `5/2`. The result is kept in
exact arithmetic — the median of the first ten squares is `61/2`, the average of
the 5th and 6th sorted values `25` and `36`. `Median` expects numeric data; an
even-length list of unresolved symbols cannot be averaged and is left
unevaluated with a `Median::rectn` message.
