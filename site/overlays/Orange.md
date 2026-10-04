### Worked examples

```mathematica
In[1]:= Orange  (* a named colour evaluates to its colour directive *)
```

```mathematica
In[1]:= Orange === RGBColor[1, 0.5, 0]  (* same thing, confirmed *)
```

```mathematica
In[1]:= List @@ Orange  (* its components *)
```

### Notes

`Orange` is one of Mathilda's named colour constants. Rather than being an inert head, it
is a symbol whose OwnValue is the directive `RGBColor[1, 0.5, 0]`, so it resolves to a real colour
wherever one is expected — as a `PlotStyle` / `VectorStyle` directive, inside
`Graphics` primitives, or in `ColorRules`.

It expands to an `RGBColor` directive, like most of the named colours. The exact components `0` and `1` print as plain integers (matching
Mathematica's `InputForm`). `Orange` is `Protected`.
