### Worked examples

```mathematica
In[1]:= Lighting -> None  (* an inert option keyword: it just names the left of a rule *)
```

```mathematica
In[1]:= MemberQ[Attributes[Lighting], Protected]  (* protected, with no value of its own *)
```

```mathematica
In[1]:= Head[Plot3D[x + y, {x, 0, 1}, {y, 0, 1}, Lighting -> None]]  (* accepted by its owning plot builtin *)
```

### Notes

`Lighting` is an inert option keyword for `Graphics3D`, `Plot3D` and `ParametricPlot3D` — surface shading for 3D graphics (Automatic = per-face Lambertian, None/False = flat raw colour). It has no builtin and no
own-value; it is `Protected` and stays symbolic, read only as the left-hand side of
a `Lighting -> value` option rule while the owning builtin builds its `Graphics3D[...]`
result. The right-hand value evaluates normally before the builtin sees it.
