### Worked examples

```mathematica
In[1]:= edge = Image[{{0., 0., 1., 1.}, {0., 0., 1., 1.}, {0., 0., 1., 1.}}];  (* a vertical step edge *)
```

```mathematica
In[1]:= Part[ImageData[GradientFilter[edge]], 2, 2]  (* the response peaks across the step *)
```

```mathematica
In[1]:= ImageDimensions[GradientFilter[edge]]  (* same size as the input, single channel *)
```

### Notes

`GradientFilter[image]` gives the gradient magnitude `Sqrt[dx^2 + dy^2]`, using the normalised
Sobel derivatives of `DerivativeFilter`. The magnitude rather than `|dx| + |dy|` because it is
**rotation invariant**: an edge at 45° reports the same strength as one at 0°, where the
absolute sum would report it `√2` times stronger and so bias every downstream threshold by
orientation.

A colour image is reduced to luminance first and differentiated once, rather than
differentiated per channel and combined by some arbitrary rule. The result is a single-channel
`"Real"` image.
