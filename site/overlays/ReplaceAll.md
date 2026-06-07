---
status: Stable
references:
---
### Worked examples

```mathematica
In[1]:= x + y /. x -> 2
Out[1]= 2 + y
```

```mathematica
In[1]:= {x, y, z} /. {x -> 1, z -> 3}
Out[1]= {1, y, 3}
```

```mathematica
In[1]:= x^2 + x /. x -> a + 1
Out[1]= 1 + a + (1 + a)^2
```

```mathematica
In[1]:= f[1, 2] /. f[a_, b_] -> a + b
Out[1]= 3
```

### Notes

`expr /. rules` is the shorthand for `ReplaceAll[expr, rules]`. It traverses
`expr` top-down and rewrites each subexpression that matches a rule's left-hand
side. A single rule or a list of rules may be supplied; with a list, the first
matching rule wins at each position. Replacement is not re-applied to the result
of a match (that is `ReplaceRepeated`, `//.`). Patterns with named blanks
(`f[a_, b_] -> a + b`) bind variables for use on the right-hand side.
