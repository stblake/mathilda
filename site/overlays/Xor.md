### Worked examples

```mathematica
In[1]:= Xor[True, False]  (* an odd number of True arguments gives True *)
```

```mathematica
In[1]:= Xor[True, True, True]  (* three Trues: still odd *)
```

```mathematica
In[1]:= Xor[p, q, p]  (* duplicate arguments cancel in pairs *)
```

```mathematica
In[1]:= Xor[True, a]  (* a single True negates the remaining argument *)
```

```mathematica
In[1]:= Xor[p, q]  (* distinct symbolic arguments stay symbolic *)
```

### Notes

`Xor[e1, e2, …]` is `True` when an odd number of the `ei` are `True`. It is
`Flat`, `Orderless` and `OneIdentity`, so nested `Xor` flattens, arguments are
sorted canonically, and `Xor[e]` collapses to `e`.

Evaluation folds the literal Booleans and cancels duplicate arguments (`a` Xor `a`
is `False`): `Xor[]` is `False`, `Xor[True, False]` is `True`, `Xor[True, True]`
is `False`, and an odd count of consumed `True`s negates the surviving core, so
`Xor[True, a]` is `Not[a]`.
