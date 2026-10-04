### Worked examples

```mathematica
In[1]:= ImageData[ColorReplace[Image[{{0., 1.}}], 1. -> 0.]]  (* the bright pixel is recoloured to 0 *)
```

```mathematica
In[1]:= ImageChannels[ColorReplace[Image[{{1.}}], 1. -> RGBColor[1., 0., 0.]]]  (* grey -> colour when the target is non-grey *)
```

```mathematica
In[1]:= ImageData[ColorReplace[Image[{{0., 0.5, 1.}}], 0.5 -> 0., 0.1]]  (* only pixels within 0.1 of 0.5 change *)
```

### Notes

`ColorReplace[image, old -> new]` replaces every pixel within a tolerance (default `0.02`) of
`old` by `new`; a list of rules applies several at once. Where rules overlap the **nearest**
target wins rather than the first, so the result does not depend on the order the rules were
written. Colours may be `RGBColor[r, g, b]`, `GrayLevel[v]`, a number, or `{r, g, b}`.

Replacing a grey image's colour with a non-grey one produces a three-channel image — flattening
`red` to its luminance would give grey when the caller asked for red. An alpha channel always
passes through: transparency is not a colour.
