### Worked examples

```mathematica
In[1]:= Join[{1, 2}, {3, 4}, {5}]
```

```mathematica
In[1]:= Join[f[a], f[b, c]]  (* any shared head, not just List *)
```

```mathematica
In[1]:= Join[<|a -> 1|>, <|b -> 2|>]  (* associations merge, later values winning *)
```

### Notes

`Join[e1, e2, ...]` concatenates its arguments, which must share a head, into a
single expression under that head — the standard list-concatenation operator, and
it works equally on `f[...]` expressions. Associations are merged key-wise, with
later associations overriding earlier values on shared keys. A trailing integer
argument gives a level specification (default 1) so that `Join[..., n]`
concatenates at depth `n`. Use `Catenate` to flatten one list of parts rather
than several arguments.
