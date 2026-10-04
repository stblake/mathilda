### Worked examples

```mathematica
In[1]:= Red  (* the named primary evaluates straight to its RGBColor literal *)
```

```mathematica
In[1]:= FullForm[Red]  (* an OwnValue, not an inert head -- it has already resolved *)
```

```mathematica
In[1]:= Head[Graphics[{Red, Disk[]}]]  (* usable wherever a colour directive is expected *)
```

```mathematica
In[1]:= FullForm[FrameStyle -> Red]  (* an option value is a resolved colour too *)
```

### Notes

`Red` is the named colour `RGBColor[1, 0, 0]`. It is bound as a `Protected` OwnValue, not
as an inert head, so it resolves to a real colour literal wherever a bare argument is
evaluated: `FullForm[Red]` is `RGBColor[1, 0, 0]` and `Head[Red]` is `RGBColor`.

That is what lets it be dropped anywhere a style directive is accepted — inside
`Graphics`, as a `PlotStyle`/`Epilog`/`FrameStyle` value, or in a `Directive`. The three
exact-`0`/`1` components print as integers, matching Mathematica's `InputForm`. `Red`
carries no opacity; use `Opacity[a]` or the four-argument `RGBColor[r, g, b, a]` for a
translucent colour.
