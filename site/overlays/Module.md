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

```mathematica
In[1]:= g[n_] := Module[{f}, f[0] = 1; f[k_] := k*f[k - 1]; f[n]]; g[6]
Out[1]= 720
```

```mathematica
In[1]:= Module[{x = 1}, Do[x = x + 1/x, {5}]; x]
Out[1]= 969581/272890
```

### Notes

`Module[{vars}, body]` introduces lexically scoped local variables, optionally
with initial values (`{x = 5}`). The locals are renamed to unique symbols so they
never collide with global definitions of the same name. Inside the body, locals
can be reassigned (`s = n^2 + n`) and a sequence of statements separated by `;`
is evaluated left to right, with the last expression returned. This makes
`Module` the standard tool for writing multi-step function definitions. The local
symbols can even carry their own `DownValues`: in the factorial example a *local*
`f` is given a base case and a recursive rule, so the recursion runs entirely
inside the module without polluting any global `f`. The continued-fraction
example iterates `x -> x + 1/x` five times with `Do`, returning the exact rational
result thanks to arbitrary-precision arithmetic.
