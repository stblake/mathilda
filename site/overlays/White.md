### Worked examples

```mathematica
In[1]:= White  (* a named colour evaluates to its colour directive *)
```

```mathematica
In[1]:= White === GrayLevel[1]  (* same thing, confirmed *)
```

```mathematica
In[1]:= List @@ White  (* its components *)
```

### Notes

`White` is one of Mathilda's named colour constants. Rather than being an inert head, it
is a symbol whose OwnValue is the directive `GrayLevel[1]`, so it resolves to a real colour
wherever one is expected — as a `PlotStyle` / `VectorStyle` directive, inside
`Graphics` primitives, or in `ColorRules`.

A grey-scale constant, so it expands to a `GrayLevel` directive rather than `RGBColor`. The exact components `0` and `1` print as plain integers (matching
Mathematica's `InputForm`). `White` is `Protected`.
