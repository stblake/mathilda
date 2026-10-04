### Worked examples

```mathematica
In[1]:= LightYellow  (* a pale pastel -- distinct from Yellow *)
```

```mathematica
In[1]:= FullForm[LightYellow]  (* an OwnValue, not an inert head -- already resolved *)
```

```mathematica
In[1]:= Head[Graphics[{LightYellow, Disk[]}]]  (* usable wherever a directive is expected *)
```

### Notes

`LightYellow` is the named colour `RGBColor[1, 1, 0.85]`. It is bound as a `Protected`
OwnValue, so it resolves to a real colour literal wherever a bare argument is evaluated:
`FullForm[LightYellow]` is `RGBColor[1, 1, 0.85]` and `Head[LightYellow]` is `RGBColor`.

The exact `1` red/green components print as integers and the `0.85` blue as a real. It is
a distinct constant from `Yellow` (`RGBColor[1, 1, 0]`), one of the pale `Light*` pastel
family. It carries no opacity; use `Opacity[a]` or the four-argument `RGBColor[r, g, b, a]`.
