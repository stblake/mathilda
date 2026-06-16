### Worked examples

```mathematica
In[1]:= Take[{a, b, c, d}, UpTo[2]]
Out[1]= {a, b}

In[2]:= Take[{a, b}, UpTo[5]]
Out[2]= {a, b}
```

The defining behavior is graceful saturation: asking for far more elements than exist returns everything available instead of raising an error, which makes it safe inside generic pipelines where the length is unknown:

```mathematica
In[1]:= Take[Range[10], UpTo[100]]
Out[1]= {1, 2, 3, 4, 5, 6, 7, 8, 9, 10}
```

### Notes

`UpTo[n]` is a count specification meaning "as many as `n`, but no error if fewer are available." With `Take`, requesting more elements than exist (Out[2], and the saturating example above) returns all of them rather than failing.
