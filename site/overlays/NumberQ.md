### Worked examples

```mathematica
In[1]:= NumberQ[3]
```

```mathematica
In[1]:= NumberQ[2/3]
```

```mathematica
In[1]:= NumberQ[1 + 2 I]
```

```mathematica
In[1]:= NumberQ[Pi]  (* a symbolic constant is not an explicit number *)
```

```mathematica
In[1]:= NumberQ[x]
```

### Notes

`NumberQ[expr]` is `True` for an **explicit** number — an integer, bigint, machine
or arbitrary-precision real, rational, or complex. It draws the line exactly where
`NumericQ` does not: `NumberQ[Pi]` is `False` because `Pi` is a symbol that merely
*has* a numeric value, whereas `NumericQ[Pi]` is `True`. Use `NumberQ` when you
need an already-evaluated literal number, and `NumericQ` when a symbolic constant
or a numeric-function call should also qualify.
