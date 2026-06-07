---
status: Stable
references:
  - "Harold Abelson and Gerald Jay Sussman, *Structure and Interpretation of Computer Programs*, 2nd ed., §3.1 (local state and lexical scoping)."
---
### Worked examples

```mathematica
In[1]:= Module[{x = 5}, x^2 + 1]
Out[1]= 26
```

```mathematica
In[1]:= Module[{a = 2, b = 3}, a*b + a + b]
Out[1]= 11
```

```mathematica
In[1]:= f[n_] := Module[{s = 0}, s = n^2 + n; s]; f[4]
Out[1]= 20
```

### Notes

`Module[{vars}, body]` introduces lexically scoped local variables, optionally
with initial values (`{x = 5}`). The locals are renamed to unique symbols so they
never collide with global definitions of the same name. Inside the body, locals
can be reassigned (`s = n^2 + n`) and a sequence of statements separated by `;`
is evaluated left to right, with the last expression returned. This makes
`Module` the standard tool for writing multi-step function definitions.
