### Worked examples

```mathematica
In[1]:= {1, Splice[{2, 3}], 4}  (* splices into the surrounding list *)
In[2]:= f[a, Splice[{b, c}, _], d]  (* the head pattern _ lets it splice into any head *)
In[3]:= Splice[{1, 2}]  (* inert on its own *)
```

### Notes

The one-argument `Splice[{e1, e2, …}]` expands into the enclosing `List` or
`Association`; to splice into another head you give a head pattern, as in
`Splice[list, _]`. On its own, or in a head it does not apply to, `Splice` stays
unevaluated. The work happens in the evaluator's sequence-flattening pass.
