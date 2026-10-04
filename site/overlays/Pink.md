### Worked examples

```mathematica
In[1]:= Pink  (* the named constant evaluates to its RGBColor literal *)
```

```mathematica
In[1]:= FullForm[Pink]  (* its full form is a plain RGBColor triple *)
```

```mathematica
In[1]:= MemberQ[Attributes[Pink], Protected]  (* a protected style constant *)
```

```mathematica
In[1]:= Graphics[{Pink, Disk[]}]  (* used as a directive, styling the primitive that follows it *)
```

### Notes

`Pink` is a named colour constant for `RGBColor[1, 0.5, 0.5]` (a soft pink). Unlike `CMYKColor`, which is
an inert head, `Pink` is an **OwnValue**: it *evaluates* to its `RGBColor[...]`
literal, which is why it resolves anywhere a real colour would — as a primitive-list
directive, a `PlotStyle`/`Background`/`FrameStyle` value, or a `ColorFunction`
result. A channel of exactly `0` or `1` prints as an integer; other channels print
as reals. It is `Protected`, and is one of the named colours defined together in
`graphics_init.c`.
