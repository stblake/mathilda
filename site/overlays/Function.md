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

### Notes

`body &` is the anonymous (pure) function shorthand, where `#`/`#1` denotes the
first argument, `#2` the second, and so on. The named forms `Function[x, body]`
and `Function[{x1, x2, ...}, body]` bind explicit parameter symbols. Parameter
binding is lexical: arguments are substituted into the body before evaluation,
and nested `Function`s shadow their own parameters. A pure function can be stored
in a symbol (`f = #1 - #2 &`) and called like any other function head.
