### Worked examples

```mathematica
In[1]:= Catch[Throw[5]]  (* the thrown value is handed to the nearest Catch *)
```

```mathematica
In[1]:= g[x_] := If[x > 10, Throw[overflow], x!]; Catch[g[2] + g[11]]  (* the throw unwinds out of the Plus *)
```

```mathematica
In[1]:= Throw[orphan]  (* uncaught at top level: reported and wrapped in Hold *)
```

### Notes

`Throw[value]` stops evaluation and returns `value` to the nearest enclosing
`Catch`. Its arguments are evaluated first (`Throw` is not held), and the
`Throw[...]` node then propagates up through every intervening expression — the
second example shows it escaping a partially-evaluated `Plus`.

`Throw[value, tag]` is caught only by a `Catch[expr, form]` whose `form` matches
`tag`. An uncaught `Throw[value]` or `Throw[value, tag]` prints `Throw::nocatch`
and comes back as `Hold[Throw[...]]`; an uncaught `Throw[value, tag, f]` instead
returns `f[value, tag]`.
