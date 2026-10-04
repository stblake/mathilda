### Worked examples

```mathematica
In[1]:= ChartStyle -> {Red, Blue}  (* an inert option keyword: it just names the left of a rule *)
```

```mathematica
In[1]:= MemberQ[Attributes[ChartStyle], Protected]  (* protected, with no value of its own *)
```

```mathematica
In[1]:= Head[BarChart[{3, 1, 2}, ChartStyle -> {Red, Green, Blue}]]  (* accepted by its owning plot builtin *)
```

### Notes

`ChartStyle` is an inert option keyword for `BarChart` and `Histogram` — a colour or list of colours cycled through the bars. It has no builtin and no
own-value; it is `Protected` and stays symbolic, read only as the left-hand side of
a `ChartStyle -> value` option rule while the owning builtin builds its `Graphics[...]`
result. The right-hand value evaluates normally before the builtin sees it.
