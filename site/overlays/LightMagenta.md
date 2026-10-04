### Worked examples

```mathematica
In[1]:= LightMagenta  (* the named constant evaluates to its RGBColor literal *)
```

```mathematica
In[1]:= FullForm[LightMagenta]  (* its full form is a plain RGBColor triple *)
```

```mathematica
In[1]:= MemberQ[Attributes[LightMagenta], Protected]  (* a protected style constant *)
```

```mathematica
In[1]:= Graphics[{LightMagenta, Disk[]}]  (* used as a directive, styling the primitive that follows it *)
```

### Notes

`LightMagenta` is a named colour constant for `RGBColor[1, 0.9, 1]` (a pale magenta tint). Unlike `CMYKColor`, which is
an inert head, `LightMagenta` is an **OwnValue**: it *evaluates* to its `RGBColor[...]`
literal, which is why it resolves anywhere a real colour would — as a primitive-list
directive, a `PlotStyle`/`Background`/`FrameStyle` value, or a `ColorFunction`
result. A channel of exactly `0` or `1` prints as an integer; other channels print
as reals. It is `Protected`, and is one of the named colours defined together in
`graphics_init.c`.
