### Worked examples

```mathematica
In[1]:= StringQ["AbC"]
```

```mathematica
In[1]:= StringQ[""]  (* the empty string still counts *)
```

```mathematica
In[1]:= StringQ[123]
```

```mathematica
In[1]:= StringQ[{"a", "b"}]  (* not Listable: a list is not a string *)
```

### Notes

`StringQ[expr]` is `True` exactly when `expr` is a string, including the empty
string `""`. It is deliberately **not** `Listable`, so a list of strings is tested
as a single object and gives `False` rather than threading element-wise — if you
want the element-wise test, map it. Any arity other than one is a malformed call:
`StringQ[]` emits `StringQ::argx` and stays unevaluated.
