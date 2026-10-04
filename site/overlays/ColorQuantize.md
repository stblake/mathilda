### Worked examples

```mathematica
(* quantising three distinct colours down to two preserves the image dimensions *)
In[1]:= ImageDimensions[ColorQuantize[Image[{{{1., 0, 0}, {0, 1., 0}, {0, 0, 1.}}}], 2]]
```

```mathematica
(* the channel count is preserved *)
In[1]:= ImageChannels[ColorQuantize[Image[{{{1., 0, 0}, {0, 1., 0}}}], 2]]
```

```mathematica
(* deterministic: the same image quantised twice gives the same result *)
In[1]:= ImageData[ColorQuantize[Image[{{{1., 0, 0}, {0, 1., 0}}}], 2]] === ImageData[ColorQuantize[Image[{{{1., 0, 0}, {0, 1., 0}}}], 2]]
```

### Notes

`ColorQuantize[image, n]` reduces the image to at most `n` colours by **median
cut**: the box with the widest single-channel spread is split at its median until
`n` boxes remain, and each box collapses to its mean colour.

Widest spread, not most pixels, because a large box of nearly identical colours
does not need splitting and a small one spanning half the spectrum does. Median
cut, not k-means, because it is **deterministic** — a palette that depended on a
random stream could not be tested or documented. The channel count is preserved
and an alpha channel passes through.
