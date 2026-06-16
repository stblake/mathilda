---
status: Stable
references:
---
### Worked examples

```mathematica
In[1]:= NestList[f, x, 3]
Out[1]= {x, f[x], f[f[x]], f[f[f[x]]]}
```

```mathematica
In[1]:= NestList[2 # &, 1, 5]
Out[1]= {1, 2, 4, 8, 16, 32}
```

The successive convergents of Newton's iteration for `Sqrt[2]`, kept exact —
each entry roughly doubles the number of correct digits:

```mathematica
In[1]:= NestList[(# + 2/#)/2 &, 1, 4]
Out[1]= {1, 3/2, 17/12, 577/408, 665857/470832}
```

Building the first few convergents of a continued fraction symbolically:

```mathematica
In[1]:= NestList[1/(1 + #) &, x, 3]
Out[1]= {x, 1/(1 + x), 1/(1 + 1/(1 + x)), 1/(1 + 1/(1 + 1/(1 + x)))}
```

### Notes

`NestList[f, expr, n]` returns a list of length `n + 1`: the first element is
`expr` and each subsequent element applies one more `f`. It is the
value-collecting companion of `Nest`, useful for capturing an entire orbit or
sequence (powers of two, a counter, the iterates of a map). `n` must be a
non-negative integer, and each intermediate value is evaluated before being
appended.
