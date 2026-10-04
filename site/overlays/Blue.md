### Worked examples

```mathematica
In[1]:= Blue  (* the named constant evaluates to its RGBColor literal *)
```

```mathematica
In[1]:= FullForm[Blue]  (* its full form is a plain RGBColor triple *)
```

```mathematica
In[1]:= MemberQ[Attributes[Blue], Protected]  (* a protected style constant *)
```

```mathematica
In[1]:= Graphics[{Blue, Disk[]}]  (* used as a directive, styling the primitive that follows it *)
```

### Notes

`Blue` is a named colour constant for `RGBColor[0, 0, 1]` (pure blue). Unlike `CMYKColor`, which is
an inert head, `Blue` is an **OwnValue**: it *evaluates* to its `RGBColor[...]`
literal, which is why it resolves anywhere a real colour would — as a primitive-list
directive, a `PlotStyle`/`Background`/`FrameStyle` value, or a `ColorFunction`
result. A channel of exactly `0` or `1` prints as an integer; other channels print
as reals. It is `Protected`, and is one of the named colours defined together in
`graphics_init.c`.
