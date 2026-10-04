### Worked examples

```mathematica
In[1]:= LightPink  (* a pale pastel -- distinct from Pink *)
```

```mathematica
In[1]:= FullForm[LightPink]  (* an OwnValue, not an inert head -- already resolved *)
```

```mathematica
In[1]:= Head[Graphics[{LightPink, Disk[]}]]  (* usable wherever a directive is expected *)
```

### Notes

`LightPink` is the named colour `RGBColor[1, 0.925, 0.925]`. It is bound as a `Protected`
OwnValue, so it resolves to a real colour literal wherever a bare argument is evaluated:
`FullForm[LightPink]` is `RGBColor[1, 0.925, 0.925]` and `Head[LightPink]` is `RGBColor`.

The exact `1` red component prints as an integer and the `0.925` green/blue as reals. It
is a distinct constant from `Pink` (`RGBColor[1, 0.5, 0.5]`), one of the pale `Light*`
pastel family. It carries no opacity; use `Opacity[a]` or the four-argument
`RGBColor[r, g, b, a]`.
