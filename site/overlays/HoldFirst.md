### Worked examples

```mathematica
In[1]:= Attributes[HoldFirst]  (* an attribute symbol carries itself, plus Protected *)
```

```mathematica
In[1]:= SetAttributes[holdfn, HoldFirst]  (* give a test function the attribute *)
```

```mathematica
In[1]:= holdfn[1 + 1, 2 + 2]  (* the first argument is held, the rest evaluate *)
```

```mathematica
In[1]:= Attributes[holdfn]  (* and the function now reports it *)
```

### Notes

`HoldFirst` is an evaluation attribute. Setting it on a symbol `f` makes the
evaluator keep `f`'s first argument unevaluated while evaluating the others, which
is exactly what assignment heads such as `Set`, `AppendTo` and `SetAttributes` rely
on. It is one bit (`ATTR_HOLDFIRST`) in a symbol's attribute mask; `SetAttributes`
/ `ClearAttributes` flip it and `Attributes` reads it back. `Evaluate` can override
a hold locally. It is the first-argument member of the `HoldFirst` / `HoldRest` /
`HoldAll` family.
