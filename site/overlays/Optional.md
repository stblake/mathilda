### Worked examples

A default-valued parameter: the second argument may be omitted, in which case
the default is supplied.

```mathematica
In[1]:= f[x_, y_ : 1] := x + y

In[2]:= f[3]
Out[2]= 4

In[3]:= f[3, 10]
Out[3]= 13

In[4]:= g[x_, y_ : 0, z_ : 0] := {x, y, z}

In[5]:= g[1, 2]
Out[5]= {1, 2, 0}
```

Several defaulted parameters give "optional trailing arguments" with sensible
fallbacks — here a quadratic whose linear and constant coefficients default to
`0` and `1`:

```mathematica
In[1]:= lin[a_, b_ : 0, c_ : 1] := a x^2 + b x + c

In[2]:= lin[2]
Out[2]= 1 + 2 x^2

In[3]:= lin[2, 3, 4]
Out[3]= 4 + 3 x + 2 x^2
```

The `_.` sugar draws the default from `Default[f]` at the call site, so a
pattern like `x_ + y_.` matches a bare term by treating the missing summand as
its additive identity `0` — the mechanism Mathematica uses to make rules robust
against absent structure:

```mathematica
In[1]:= p[x_ + y_.] := {x, y}

In[2]:= p[a]
Out[2]= {a, 0}
```

### Notes

`Optional[p, def]` (surface syntax `p : def`) lets a pattern argument be omitted, supplying `def` in its place — the standard way to give function definitions default-valued parameters (Out[2], Out[5]). The shorthand `patt_.` is sugar for `Optional[patt_, Default[f]]`, drawing the default from `Default[f]` at the call site so a single rule can match expressions with or without a given term.
