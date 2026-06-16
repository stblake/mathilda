### Worked examples

```mathematica
In[1]:= Range[10][[2 ;; 8]]
Out[1]= {2, 3, 4, 5, 6, 7, 8}
```

```mathematica
In[1]:= Range[10][[2 ;; 8 ;; 2]]
Out[1]= {2, 4, 6, 8}

In[2]:= Range[10][[;; ;; 3]]
Out[2]= {1, 4, 7, 10}
```

### Notes

`i ;; j` is the span of elements `i` through `j`, and `i ;; j ;; k` steps by `k`. Inside `Part` (`[[...]]`) it extracts a contiguous or strided slice: `;; ;; 3` keeps every third element starting from the first, and omitting an endpoint defaults to the start or end of the list. Evaluated on its own a span stays inert as `Span[i, j]`, since its meaning is supplied by the indexing context it appears in.
