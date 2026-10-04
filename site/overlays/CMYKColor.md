### Worked examples

```mathematica
(* inert: CMYKColor stays symbolic, converted to RGB only at render time *)
In[1]:= CMYKColor[1, 0, 0, 0]
```

```mathematica
(* a cyan disk -- the directive styles the primitive that follows it *)
In[1]:= Graphics[{CMYKColor[1, 0, 0, 0], Disk[]}]
```

```mathematica
(* it is a Protected style directive, like RGBColor/GrayLevel/Hue *)
In[1]:= MemberQ[Attributes[CMYKColor], Protected]
```

```mathematica
(* its head is preserved -- there is no builtin that rewrites it *)
In[1]:= Head[CMYKColor[0.2, 0.4, 0.1, 0.3]]
```

### Notes

`CMYKColor` is an inert, `Protected` style directive for the subtractive (print)
colour model. It has no builtin and no own-values, so it never evaluates to an
`RGBColor` at the language level — it stays symbolic and is converted only at render
time, by `r = (1 - c)(1 - k)`, `g = (1 - m)(1 - k)`, `b = (1 - y)(1 - k)`. The forms
are `CMYKColor[c, m, y, k]`, `CMYKColor[c, m, y]` (`k = 0`), `CMYKColor[c, m, y, k,
a]` (with opacity), and the list forms `CMYKColor[{c, m, y, k}]` /
`CMYKColor[{c, m, y, k, a}]`; components and opacity outside `[0, 1]` are clipped.

It is recognised as a colour wherever `RGBColor`/`GrayLevel`/`Hue` are — as a
primitive-list directive, a `PlotStyle`/`Background`/`FrameStyle` value, or a
`ColorFunction` return value — and the same conversion is implemented in both the
CLI Raylib renderer and the notebook JSON frontend.
