### Worked examples

```mathematica
In[1]:= Catch[Do[If[i^2 > 20, Throw[i]], {i, 1, 10}]]  (* the first i whose square exceeds 20 *)
```

```mathematica
In[1]:= Catch[a; b; Throw[c]; d]  (* evaluation stops at the Throw, so d is never reached *)
```

```mathematica
In[1]:= Catch[Throw[x, u], u]  (* a tagged throw caught by a matching form *)
```

```mathematica
In[1]:= Catch[Throw[v, tg], tg, f]  (* the three-argument form returns f[value, tag] *)
```

```mathematica
In[1]:= Catch[Throw[val, mytag], othertag]  (* a non-matching tag is left uncaught *)
```

```mathematica
In[1]:= Catch[1 + 1]  (* with no throw, Catch is just the value of its body *)
```

### Notes

`Catch[expr]` catches the first `Throw` generated anywhere while `expr` is
evaluated and returns its value, or returns `expr`'s own value if nothing is
thrown. Because the throw propagates through *any* enclosing head — `Plus`,
`Map`, `Table`, a function call — `Catch` is the tool for a non-local exit, where
`Return` only escapes a scope boundary.

`Catch[expr, form]` catches only a `Throw[value, tag]` whose `tag` matches `form`
(the tag is re-evaluated before each comparison); a throw with a non-matching tag,
and any tagless `Throw[value]`, propagates to an outer `Catch`. `Catch[expr, form,
f]` returns `f[value, tag]`. The construct uses sentinel propagation rather than
`setjmp`/`longjmp`, so it is leak-free.
