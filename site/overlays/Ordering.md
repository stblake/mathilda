### Worked examples

```mathematica
In[1]:= Ordering[{c, a, b}]  (* positions that sort the list: list[[%]] is Sort[list] *)
```

```mathematica
In[2]:= Ordering[{2, 6, 1, 9, 1, 2, 3}, 4]  (* positions of the 4 smallest *)
```

```mathematica
In[3]:= Ordering[{2, 2, 1}]  (* stable: ties keep input order *)
```

```mathematica
In[4]:= Ordering[{5, 3, 8, 1}, All, Greater]  (* order descending with a custom function *)
```

### Notes

`Ordering[list]` gives the permutation of 1-based positions for which
`list[[Ordering[list]]]` is `Sort[list]`. The argsort is **stable** — ties are
broken by original position — so `Ordering[{2, 2, 1}]` is `{3, 1, 2}` and
`Ordering[list, 1]` names the *first* minimum. `Ordering[list, seq]` is
`Take[Ordering[list], seq]` (an integer `n`/`-n` for the `n` smallest/largest, a
`{m, n[, s]}` span, `UpTo[k]`, or `All`), and the 3-argument form orders by a
function `p` as in `Sort[list, p]`. The result is always a `List` of integers
regardless of `list`'s head; over an `Association` it orders by the values. A
machine-number vector takes a packed argsort (result dtype always integer), and
`Ordering[vector]` compiles.
