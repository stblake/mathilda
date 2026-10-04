### Worked examples

```mathematica
In[1]:= LightGray  (* a named colour evaluates to its colour directive *)
```

```mathematica
In[1]:= LightGray === GrayLevel[0.85]  (* same thing, confirmed *)
```

```mathematica
In[1]:= List @@ LightGray  (* its components *)
```

### Notes

`LightGray` is one of Mathilda's named colour constants. Rather than being an inert head, it
is a symbol whose OwnValue is the directive `GrayLevel[0.85]`, so it resolves to a real colour
wherever one is expected — as a `PlotStyle` / `VectorStyle` directive, inside
`Graphics` primitives, or in `ColorRules`.

A grey-scale constant, so it expands to a `GrayLevel` directive rather than `RGBColor`. The exact components `0` and `1` print as plain integers (matching
Mathematica's `InputForm`). `LightGray` is `Protected`.
