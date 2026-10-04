### Worked examples

```mathematica
In[1]:= Order[a, b]  (* a sorts before b: +1 *)
```

```mathematica
In[2]:= Order[b, a]  (* reversed: -1 *)
```

```mathematica
In[3]:= Order[a, a]  (* identical: 0 *)
```

```mathematica
In[4]:= Order[6, Pi]  (* structural: the Integer 6 sorts before the symbol Pi *)
```

```mathematica
In[5]:= Order[6, N[Pi]]  (* two numeric atoms, now compared by value *)
```

### Notes

`Order[e1, e2]` is the user-facing surface of the canonical comparator every
sorting routine (`Sort`, `SortBy`, `OrderedQ`, `Ordering`) is built on: `1` if
`e1` is before `e2`, `-1` if after, `0` if identical. The comparison is
**structural**, not by numerical value — `Order[6, Pi]` is `1` because an
`Integer` sorts before a symbol, whereas `Order[6, N[Pi]]` is `-1` because two
numeric atoms compare by value. It needs exactly two arguments (otherwise it
stays unevaluated) and is compilable: over machine numbers it lowers to
`Sign[e2 - e1]`, returning the same `{1, 0, -1}` integer.
