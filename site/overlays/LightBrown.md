### Worked examples

```mathematica
In[1]:= LightBrown  (* a pale pastel -- its own constant, not a tint of Brown *)
```

```mathematica
In[1]:= FullForm[LightBrown]  (* an OwnValue, not an inert head -- already resolved *)
```

```mathematica
In[1]:= Head[Graphics[{LightBrown, Disk[]}]]  (* usable wherever a directive is expected *)
```

### Notes

`LightBrown` is the named colour `RGBColor[0.94, 0.91, 0.88]`. It is bound as a `Protected`
OwnValue, so it resolves to a real colour literal wherever a bare argument is evaluated:
`FullForm[LightBrown]` is `RGBColor[0.94, 0.91, 0.88]` and `Head[LightBrown]` is `RGBColor`.

It is a *separate* constant from `Brown`, not a computed tint — the `Light*` family is its
own set of pale pastel entries (`LightRed`, `LightGreen`, `LightBlue`, …). It carries no
opacity; use `Opacity[a]` or the four-argument `RGBColor[r, g, b, a]`.
