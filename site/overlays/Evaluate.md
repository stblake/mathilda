### Worked examples

```mathematica
In[1]:= Hold[1 + 1]  (* Hold keeps its argument unevaluated *)
```

```mathematica
In[1]:= Hold[Evaluate[1 + 1]]  (* Evaluate forces evaluation even inside the held position *)
```

### Notes

`Evaluate[expr]` forces `expr` to be evaluated even when it sits in an argument position
that the enclosing head would otherwise hold. In an ordinary (unheld) position it is a
no-op, since the argument has already been evaluated by the standard pipeline — the C
builtin then just unwraps it.

With more than one argument it splices through `Sequence`, so `Evaluate[a, b]` behaves as
`Sequence[a, b]`.
