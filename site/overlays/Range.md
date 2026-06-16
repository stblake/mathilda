---
status: Stable
references:
---
### Worked examples

```mathematica
In[1]:= Range[5]
Out[1]= {1, 2, 3, 4, 5}
```

```mathematica
In[1]:= Range[2, 10, 2]
Out[1]= {2, 4, 6, 8, 10}
```

```mathematica
In[1]:= Range[0, 1, 1/4]
Out[1]= {0, 1/4, 1/2, 3/4, 1}
```

A negative step counts down, and `Range` chains naturally with the functional
operators it is built to feed — here the exact triangular numbers and the sum of
the first hundred integers:

```mathematica
In[1]:= Range[10, 1, -1]
Out[1]= {10, 9, 8, 7, 6, 5, 4, 3, 2, 1}

In[2]:= Map[#^2 &, Range[5]]
Out[2]= {1, 4, 9, 16, 25}

In[3]:= Total[Range[100]]
Out[3]= 5050
```

### Notes

`Range[n]` produces the integers `1` through `n`. `Range[n, m]` runs from `n` to
`m` in unit steps, and `Range[n, m, d]` uses step `d`. The step may be an exact
rational, in which case the result stays exact rather than being converted to
floating point. `Range` is the most direct way to build the index list that
`Map`, `Select`, or `Fold` then consume.
