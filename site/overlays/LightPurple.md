### Worked examples

```mathematica
In[1]:= LightPurple  (* the named constant evaluates to its RGBColor literal *)
```

```mathematica
In[1]:= FullForm[LightPurple]  (* its full form is a plain RGBColor triple *)
```

```mathematica
In[1]:= MemberQ[Attributes[LightPurple], Protected]  (* a protected style constant *)
```

```mathematica
In[1]:= Graphics[{LightPurple, Disk[]}]  (* used as a directive, styling the primitive that follows it *)
```

### Notes

`LightPurple` is a named colour constant for `RGBColor[0.94, 0.88, 0.94]` (a pale purple tint). Unlike `CMYKColor`, which is
an inert head, `LightPurple` is an **OwnValue**: it *evaluates* to its `RGBColor[...]`
literal, which is why it resolves anywhere a real colour would — as a primitive-list
directive, a `PlotStyle`/`Background`/`FrameStyle` value, or a `ColorFunction`
result. A channel of exactly `0` or `1` prints as an integer; other channels print
as reals. It is `Protected`, and is one of the named colours defined together in
`graphics_init.c`.
