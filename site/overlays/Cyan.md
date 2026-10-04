### Worked examples

```mathematica
In[1]:= Cyan  (* the named secondary evaluates straight to its RGBColor literal *)
```

```mathematica
In[1]:= FullForm[Cyan]  (* an OwnValue, not an inert head -- already resolved *)
```

```mathematica
In[1]:= Head[Graphics[{Cyan, Disk[]}]]  (* usable wherever a directive is expected *)
```

### Notes

`Cyan` is the named colour `RGBColor[0, 1, 1]`. It is bound as a `Protected` OwnValue, so
it resolves to a real colour literal wherever a bare argument is evaluated: `FullForm[Cyan]`
is `RGBColor[0, 1, 1]` and `Head[Cyan]` is `RGBColor`.

`Cyan`, `Magenta` (`RGBColor[1, 0, 1]`) and `Yellow` (`RGBColor[1, 1, 0]`) are the
secondary (subtractive) primaries — the colours a single-ink `CMYKColor` directive maps
to. The exact-`0`/`1` components print as integers. `Cyan` carries no opacity; use
`Opacity[a]` or the four-argument `RGBColor[r, g, b, a]`.
