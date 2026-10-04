### Worked examples

```mathematica
In[1]:= AssociationMap[f, {a, b, c}]
```

```mathematica
In[1]:= AssociationMap[#^2 &, {1, 2, 3}]  (* keys are the inputs, values the results *)
```

```mathematica
In[1]:= AssociationMap[Reverse, <|a -> 1, b -> 2|>]  (* over an association, f acts on each key -> value rule *)
```

### Notes

`AssociationMap[f, {k1, ...}]` builds `<|k1 -> f[k1], ...|>` — the keys are the
list elements and the values are `f` applied to them, which is the natural way to
tabulate a function over a set of inputs. Given an association instead of a key
list, `f` is applied to each entry *as a rule* `k -> v`, and its result (a rule,
a list of rules, an association, or `Nothing`) is spliced into the output — so
`AssociationMap[Reverse, assoc]` swaps keys and values.
