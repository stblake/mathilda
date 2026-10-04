### Worked examples

```mathematica
In[1]:= SetAttributes[f, HoldAll]  (* give f the HoldAll attribute *)
In[1]:= Attributes[f]
In[1]:= f[2 + 3]  (* the argument is held, so the sum is not computed *)
In[1]:= g[2 + 3]  (* an ordinary head evaluates its argument first *)
```

```mathematica
In[1]:= Attributes[SetDelayed]  (* := is defined to hold all of its arguments *)
```

### Notes

`HoldAll` is an attribute symbol, not a function. A head carrying it has all of its
arguments held unevaluated before the head runs — the mechanism behind `Hold`,
`SetDelayed` (`:=`), `Function` and the control-flow heads. It appears only inside
`Attributes[...]`, `SetAttributes[sym, HoldAll]` and `ClearAttributes[sym, HoldAll]`.

Internally `HoldAll` is the pair `HoldFirst | HoldRest`; `Attributes` reports that pair as
the single token `HoldAll`. Unlike `HoldAllComplete` it still allows upvalue lookup and
`Sequence`/`Unevaluated` processing on the held arguments.
