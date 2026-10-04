### Worked examples

```mathematica
In[1]:= Equivalent[True, True]  (* all arguments share one truth value *)
```

```mathematica
In[1]:= Equivalent[True, False]  (* a True and a False together *)
```

```mathematica
In[1]:= Equivalent[True, a, b]  (* a literal True forces the remaining atoms *)
```

```mathematica
In[1]:= Equivalent[False, a]  (* a literal False negates the remaining atom *)
```

```mathematica
In[1]:= Equivalent[p, q]  (* distinct symbolic atoms stay symbolic *)
```

### Notes

`Equivalent[e1, e2, …]` is `True` when all of the `ei` share one truth value — all
`True` or all `False`. It is `Flat`, `Orderless` and `OneIdentity`; `Equivalent[]`
and `Equivalent[e]` are `True`.

Evaluation folds the literals and cancels duplicates: a `True` together with a
`False` gives `False`, a single literal forces each remaining atom (`True` leaves
it, `False` negates it) and the results are joined with `And`. `LogicalExpand`,
`Reduce` and `FindInstance` expand `Equivalent` to the cyclic conjunction
`Implies[a1, a2] && … && Implies[an, a1]`.
