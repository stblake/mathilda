### Worked examples

```mathematica
In[1]:= DeleteMissing[{1, Missing[], 2, Missing["x"], 3}]
```

```mathematica
In[1]:= DeleteMissing[{a, Missing[], b, c, Missing["NotAvailable"]}]
```

```mathematica
In[1]:= DeleteMissing[<|a -> 1, b -> Missing[], c -> 3|>]  (* drops entries whose value is Missing *)
```

### Notes

`DeleteMissing[expr]` removes every `Missing[...]` element — equivalent to
`DeleteCases[expr, _Missing]` — which is the standard way to clean a dataset of
absent-value markers. Over an association it drops the entries whose *value* is
missing, keeping the rest. The two- and three-argument forms
`DeleteMissing[expr, n]` and `DeleteMissing[expr, n, d]` restrict the work to
levels `1..n` and, with `d`, to elements that contain a missing value no deeper
than `d`.
