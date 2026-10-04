### Worked examples

```mathematica
In[1]:= MissingQ[Missing[]]
In[2]:= MissingQ[Lookup[<|a -> 1|>, b]]  (* the absent-key result is Missing *)
In[3]:= MissingQ[5]
```

### Notes

`MissingQ` is a plain head test: `True` for any `Missing[…]` expression and `False`
otherwise. It is the standard guard after a `Lookup` or key access that may not find
its key.
