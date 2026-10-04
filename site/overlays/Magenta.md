### Worked examples

```mathematica
In[1]:= Magenta  (* the named constant evaluates to its RGBColor literal *)
```

```mathematica
In[1]:= FullForm[Magenta]  (* its full form is a plain RGBColor triple *)
```

```mathematica
In[1]:= MemberQ[Attributes[Magenta], Protected]  (* a protected style constant *)
```

```mathematica
In[1]:= Graphics[{Magenta, Disk[]}]  (* used as a directive, styling the primitive that follows it *)
```

### Notes

`Magenta` is a named colour constant for `RGBColor[1, 0, 1]` (pure magenta (full red and blue)). Unlike `CMYKColor`, which is
an inert head, `Magenta` is an **OwnValue**: it *evaluates* to its `RGBColor[...]`
literal, which is why it resolves anywhere a real colour would — as a primitive-list
directive, a `PlotStyle`/`Background`/`FrameStyle` value, or a `ColorFunction`
result. A channel of exactly `0` or `1` prints as an integer; other channels print
as reals. It is `Protected`, and is one of the named colours defined together in
`graphics_init.c`.
