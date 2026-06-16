---
status: Stable
references:
  - "Harold Abelson and Gerald Jay Sussman, *Structure and Interpretation of Computer Programs*, 2nd ed., §1.3.2 (lambda; constructing procedures)."
---
### Worked examples

```mathematica
In[1]:= (# + 1 &)[10]
Out[1]= 11
```

```mathematica
In[1]:= Function[x, x^2][5]
Out[1]= 25
```

```mathematica
In[1]:= Function[{x, y}, x + y][3, 4]
Out[1]= 7
```

```mathematica
In[1]:= f = #1 - #2 &; f[10, 3]
Out[1]= 7
```

```mathematica
In[1]:= Nest[# ^ 2 + 1 &, x, 3]
Out[1]= 1 + (1 + (1 + x^2)^2)^2
```

```mathematica
In[1]:= NestList[1/(1 + #) &, x, 3]
Out[1]= {x, 1/(1 + x), 1/(1 + 1/(1 + x)), 1/(1 + 1/(1 + 1/(1 + x)))}
```

```mathematica
In[1]:= Fold[#1 * 10 + #2 &, 0, {1, 2, 3, 4}]
Out[1]= 1234
```

```mathematica
In[1]:= (## &)[a, b, c]
Out[1]= Sequence[a, b, c]
```

### Notes

`body &` is the anonymous (pure) function shorthand, where `#`/`#1` denotes the
first argument, `#2` the second, and so on. The named forms `Function[x, body]`
and `Function[{x1, x2, ...}, body]` bind explicit parameter symbols. Parameter
binding is lexical: arguments are substituted into the body before evaluation,
and nested `Function`s shadow their own parameters. A pure function can be stored
in a symbol (`f = #1 - #2 &`) and called like any other function head. Combined
with `Nest`/`NestList`/`Fold` a single anonymous function drives an entire
iteration — `Nest[#^2 + 1 &, x, 3]` builds the symbolic 3-fold composition,
`NestList[1/(1 + #) &, x, 3]` unfolds the convergents of a continued fraction, and
`Fold[#1*10 + #2 &, 0, ...]` reassembles digits via Horner's rule. The slot
sequence `##` (and `##1`, `##2`, ...) splices all arguments at once, returning a
`Sequence`.
