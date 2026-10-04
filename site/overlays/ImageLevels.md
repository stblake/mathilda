### Worked examples

```mathematica
(* a "Bit" image uses its 2 natural levels: one 0 and three 1s *)
In[1]:= ImageLevels[Image[{{0, 1}, {1, 1}}, "Bit"]]
```

```mathematica
(* forcing the bin count with a second argument *)
In[1]:= ImageLevels[Image[{{0., 0.5}, {0.5, 1.}}], 2]
```

```mathematica
(* the counts sum to the pixel count exactly *)
In[1]:= Total[Last /@ ImageLevels[Image[{{0, 1}, {1, 1}}, "Bit"]]]
```

### Notes

`ImageLevels[image]` gives `{{level, count}, ...}`: the histogram as data, not a
plot — use `Histogram` over the result for a picture. `ImageLevels[image, n]`
uses `n` bins.

Levels are on the same unit scale as `ImageData`, so a level can be compared
against a pixel value without rescaling. A `"Bit"` image uses its 2 natural
levels and a `"Byte"` its 256, because those *are* the distinct values; a
`"Real"` image has no natural set and is binned into 256 over `[0, 1]`. The
counts sum to the pixel count exactly, every pixel landing in one bin. It accepts
volumes as well as planes.
