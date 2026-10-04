### Worked examples

```mathematica
In[1]:= LightGreen  (* a pale pastel -- its own constant, not a tint of Green *)
```

```mathematica
In[1]:= FullForm[LightGreen]  (* an OwnValue, not an inert head -- already resolved *)
```

```mathematica
In[1]:= Head[Graphics[{LightGreen, Disk[]}]]  (* usable wherever a directive is expected *)
```

### Notes

`LightGreen` is the named colour `RGBColor[0.88, 1, 0.88]`. It is bound as a `Protected`
OwnValue, so it resolves to a real colour literal wherever a bare argument is evaluated:
`FullForm[LightGreen]` is `RGBColor[0.88, 1, 0.88]` and `Head[LightGreen]` is `RGBColor`.

The exact `1` green component prints as an integer while the `0.88` reds print as reals.
It is a separate constant from `Green`, one of the pale `Light*` pastel family. It carries
no opacity; use `Opacity[a]` or the four-argument `RGBColor[r, g, b, a]`.
