### Worked examples

```mathematica
In[1]:= Length[{a, b, c}]  (* the number of elements of a list *)
```

```mathematica
In[2]:= Length[f[x, y]]  (* the argument count of any function *)
```

```mathematica
In[3]:= Length[a + b + c]  (* Plus[a, b, c] has three arguments *)
```

```mathematica
In[4]:= Length[x]  (* an atom has no parts *)
```

### Notes

`Length[expr]` returns the number of top-level arguments (`arg_count`) when
`expr` is a function, and `0` for every atom — a symbol, number, or string has no
parts. Because a sum is `Plus[...]` and a list is `List[...]` internally, the
same count works on both; `Length` makes no distinction between a `List` and any
other head.
