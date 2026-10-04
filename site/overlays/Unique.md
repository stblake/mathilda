### Worked examples

```mathematica
In[1]:= Unique[]
In[2]:= Unique[x]
In[3]:= Unique["c"]
In[4]:= Unique[{p, q}]
```

### Notes

`Unique[]` generates a brand-new symbol, `Unique[x]` or `Unique["x"]` uses the
given name as a prefix, and `Unique[{a, b, ...}]` returns a list of fresh symbols
sharing one numeric suffix. Each name is formed from a prefix and a number drawn
from the `$ModuleNumber` counter — the same counter `Module` uses for its
locals.

Freshness is guaranteed by construction: the counter is advanced until the
candidate name is unused for every prefix in the batch, so a returned symbol has
never named anything before. Each generated symbol is `Temporary`.
