### Worked examples

```mathematica
In[1]:= Head[f[x]]  (* the symbol applied as the head *)
```

```mathematica
In[2]:= Head[{1, 2, 3}]  (* a list is List[...] internally *)
```

```mathematica
In[3]:= Head[3/4]  (* an atom reports its type name *)
```

```mathematica
In[4]:= Head[a + b, f]  (* the two-argument form wraps the result in f *)
```

### Notes

`Head[expr]` returns the top-level wrapper via `expr_head`: the applied symbol or
expression for a function (`f` for `f[x]`, `Plus` for `a + b`, `List` for a
list), and the type symbol for an atom — `Integer`, `Real`, `Rational`,
`Complex`, `Symbol`, or `String`. The two-argument form `Head[expr, h]` returns
`h[Head[expr]]`, leaving the outer application for the evaluator to reduce. Head
extraction is also what `expr[[0]]` and `Extract[expr, {0}]` perform.
