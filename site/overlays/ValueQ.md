### Worked examples

```mathematica
In[1]:= ValueQ[x]
```

```mathematica
In[1]:= x = 5; ValueQ[x]  (* now x carries an OwnValue *)
```

```mathematica
In[1]:= f[a_] := a^2; ValueQ[f[2]]  (* true as soon as the head f has a DownValue *)
```

```mathematica
In[1]:= ValueQ[f]  (* the bare symbol f has only DownValues, so no OwnValue *)
```

### Notes

`ValueQ[expr]` reports whether a value has been defined for `expr` **without
evaluating it** — it is `HoldAll`, so it inspects the symbol itself rather than the
value the symbol would evaluate to. A bare symbol is valued iff it carries an
`OwnValue` (an immediate `x = 5` or a delayed `y := RandomReal[]`); a compound
`f[...]` is valued iff its head `f` carries any `DownValue`, regardless of whether
the particular arguments match a rule — which is why `ValueQ[f[2]]` is `True` while
`ValueQ[f]` is `False`. Numbers, strings, and undefined symbols give `False`.
