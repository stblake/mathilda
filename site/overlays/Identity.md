### Worked examples

```mathematica
In[1]:= Identity[x]  (* returns its argument unchanged *)
```

```mathematica
In[1]:= Identity[1 + 1]  (* its argument is evaluated normally first *)
```

```mathematica
In[1]:= Map[Identity, {1, 2, 3}]  (* the identity function is handy as a default callback *)
```

### Notes

`Identity[expr]` returns `expr` unchanged. It takes exactly one argument; any other arity
is left unevaluated. Having no held attributes, it evaluates its argument through the
normal pipeline before returning it, so `Identity[1 + 1]` is `2`.

It is most useful as a neutral function argument — a do-nothing callback where some
transformation is expected.
