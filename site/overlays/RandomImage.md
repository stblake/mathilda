### Worked examples

```mathematica
In[1]:= ImageDimensions[RandomImage[]]  (* the default is a 150x150 grey noise field *)
```

```mathematica
In[1]:= ImageChannels[RandomImage[1, {4, 4}, ColorSpace -> "RGB"]]  (* three independent channels *)
```

```mathematica
In[1]:= SeedRandom[1]; ImageData[RandomImage[1, {1, 2}]]  (* drawn from RandomReal's stream, so SeedRandom pins it *)
```

### Notes

`RandomImage[]` gives a 150 × 150 grey image of uniform noise on `[0, 1]`; `RandomImage[max]`
scales the range to `[0, max]`, `RandomImage[max, {w, h}]` sets the size (a single `n` means
`{n, n}`), and `ColorSpace -> "RGB"` gives three independent channels.

Samples are drawn from the **same** stream as `RandomReal`, so `SeedRandom` makes the result
reproducible. The samples are scaled by `max` but not clamped — a `"Real"` image may hold
values above 1, and clamping belongs in `Export`.
