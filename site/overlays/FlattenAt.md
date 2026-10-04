### Worked examples

```mathematica
In[1]:= FlattenAt[{a, {b, c}, {d, e}, {f}}, 2]  (* splice the sublist at position 2 into its parent *)
```

```mathematica
In[2]:= FlattenAt[{a, {b, c}, {d, e}, {f}}, {{2}, {4}}]  (* several positions, each against the original *)
```

```mathematica
In[3]:= FlattenAt[f[g[1, 2], g[3, 4]], 1]  (* any head works, not just List *)
```

```mathematica
In[4]:= FlattenAt[{a, {b, c}, {d, e}, {f}}, -1]  (* a negative position counts from the end *)
```

### Notes

`FlattenAt` removes the head of the subexpression at each position and splices
its arguments into the surrounding expression. Position resolution is the shared
`MapAt`/`ReplaceAt` walker, so `{2}` (one deep path) and `{{2}, {4}}` (two
separate paths) mean exactly what they do there, and an out-of-range position
leaves the call unevaluated. Because the head at the position is *removed*,
`FlattenAt[{1, {{2}, {3}}, 4}, 2]` gives `{1, {2}, {3}, 4}` — to flatten *within*
a part instead, use `MapAt[Flatten, ...]`. Several positions are resolved against
the original expression with no index bookkeeping, since each targeted slot
briefly holds a single `Sequence` the evaluator splices on the next pass.
