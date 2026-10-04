### Worked examples

```mathematica
In[1]:= AssociationComap[{f, g}, x]
```

```mathematica
In[1]:= AssociationComap[{Total, Max, Min}, {3, 1, 2}]  (* several summaries of one dataset, labelled by function *)
```

### Notes

`AssociationComap[{f1, f2, ...}, x]` is the reverse of `AssociationMap`: it
applies each function to the *same* value `x` and labels the results by the
functions, giving `<|f1 -> f1[x], f2 -> f2[x], ...|>`. It is a compact way to
compute several named summaries of one object — for instance `Total`, `Max`, and
`Min` of a list — in a single keyed record. The first argument must be a list of
functions.
