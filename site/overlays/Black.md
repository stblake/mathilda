### Worked examples

```mathematica
In[1]:= Black  (* a greyscale name resolves to GrayLevel, not RGBColor *)
```

```mathematica
In[1]:= FullForm[Black]  (* an OwnValue, not an inert head -- already resolved *)
```

```mathematica
In[1]:= Head[Graphics[{Black, Line[{{0, 0}, {1, 1}}]}]]  (* usable wherever a directive is *)
```

### Notes

`Black` is the named colour `GrayLevel[0]`. Unlike the saturated primaries it uses
`GrayLevel`, not `RGBColor`, matching Mathematica's `InputForm`: `FullForm[Black]` is
`GrayLevel[0]` and `Head[Black]` is `GrayLevel`.

It is bound as a `Protected` OwnValue, so it resolves to a real colour literal wherever a
bare argument is evaluated, and can be dropped anywhere a style directive is accepted —
`Graphics`, `PlotStyle`, `FrameStyle`, `Directive`. Its companions are `White`
(`GrayLevel[1]`), `Gray` (`GrayLevel[0.5]`) and `LightGray` (`GrayLevel[0.85]`). It carries
no opacity; use `Opacity[a]` or `GrayLevel[g, a]`.
