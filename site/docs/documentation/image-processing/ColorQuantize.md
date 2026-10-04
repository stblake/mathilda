# ColorQuantize

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`ColorQuantize[image, n] reduces the image to at most n colours by MEDIAN CUT: the box with the widest single-channel spread is split at its median until n boxes remain, and each collapses to its mean colour. Widest spread rather than most pixels, since a large box of nearly identical colours does not need splitting and a small one spanning half the spectrum does. Median cut rather than k-means because it is DETERMINISTIC -- a palette that depended on the random stream could not be tested or documented. The channel count is preserved and alpha passes through.`**

## Examples (11)

Every input below was run against the current Mathilda build and its output recorded.

### Basic Examples (4)

```mathematica
In[1]:= ramp = Image[Table[N[(i + j)/32], {i, 1, 16}, {j, 1, 16}], "Real"];

In[2]:= ColorQuantize[ramp, 4]
Out[2]= -Image-

In[3]:= Table[Length[Union[Flatten[ImageData[ColorQuantize[ramp, n]]]]], {n, 1, 4}]
Out[3]= {1, 2, 3, 4}

In[4]:= ColorQuantize[Image[Table[{N[i/16], N[j/16], 0.5}, {i, 1, 16}, {j, 1, 16}], "Real"], 6]
Out[4]= -Image-
```

### Properties & Relations (4)

```mathematica
In[5]:= ramp = Image[Table[N[(i + j)/32], {i, 1, 16}, {j, 1, 16}], "Real"];
```

The same input gives the same palette, every time

```mathematica
In[6]:= ImageData[ColorQuantize[ramp, 5]] === ImageData[ColorQuantize[ramp, 5]]
Out[6]= True
```

```mathematica
In[7]:= ImageDimensions[ColorQuantize[ramp, 4]] === ImageDimensions[ramp]
Out[7]= True
```

More colours than the image holds cannot invent any

```mathematica
In[8]:= Length[Union[Flatten[ImageData[ColorQuantize[Image[{{0., 1.}, {0., 1.}}, "Real"], 8]]]]] <= 2
Out[8]= True
```

### Applications (3)

```mathematica
In[9]:= ImageDimensions[ColorQuantize[Image[{{{1., 0, 0}, {0, 1., 0}, {0, 0, 1.}}}], 2]]
Out[9]= {3, 1}

In[10]:= ImageChannels[ColorQuantize[Image[{{{1., 0, 0}, {0, 1., 0}}}], 2]]
Out[10]= 3

In[11]:= ImageData[ColorQuantize[Image[{{{1., 0, 0}, {0, 1., 0}}}], 2]] === ImageData[ColorQuantize[Image[{{{1., 0, 0}, {0, 1., 0}}}], 2]]
Out[11]= True
```

## Algorithm

imagecolor.c -- ColorReplace, ColorQuantize and HistogramTransform.

Three heads that act on an image's COLOURS rather than its geometry, and they share the one thing that makes such operations awkward: a decision made per pixel needs a global view first. Replacing a colour needs a distance rule, quantising needs a palette derived from every pixel, and equalising needs the whole distribution. So each of these makes a pass to gather, then a pass to write — which is why none of them fits the filter machinery in imagefilter.c.

## Implementation notes

**Algorithm.** `builtin_colorquantize` reduces an image to at most `n` colours by
**median cut**. It loads the pixels to a flat buffer, copies each as an RGB
triple into `qbuf` (a grey pixel is replicated across the three channels), and
starts with one box covering all pixel indices. While there are fewer than `n`
boxes it repeatedly splits the box with the **widest single-channel spread** —
widest spread rather than most pixels, since a large box of nearly identical
colours does not need splitting and a small one spanning half the spectrum does.
To split, it `qsort`s that box's indices on the widest channel (ties broken by
index, so the palette does not depend on qsort's internal choices) and cuts at
the median. When every box is a single colour the loop stops early. Finally each
box collapses to its **mean** colour, written back to every pixel it owns; a grey
output takes the mean's luminance, and an alpha channel passes through unchanged.
Median cut is chosen over k-means precisely because it is **deterministic** — a
palette depending on a random seed could not be tested or documented.

**Data structures.** `qbuf` is `n_px × 3` doubles (the colours), `qidx` an index
permutation sorted per box, `boxes` an array of `{start, len}` ranges into
`qidx`, and `out` the `n_px × channels` result buffer. File-scope `qbuf`/`qidx`/
`qchan` are the comparator's context for `qsort`. The result is a `"Real"` image.

**Complexity / limits.** Each of the up-to-`n` splits sorts a box:
`O(n · n_px log n_px)` worst case, with the per-split scan for the widest channel
`O(n · n_px)`. `n` is `1..4096`. The channel count is preserved and alpha passes
through.

- `Protected`. **Median cut**: the box with the widest single-channel spread is split at its
  median until `n` boxes remain, and each collapses to its mean colour. Widest spread rather
  than most pixels — a large box of nearly identical colours does not need splitting, and a
  small one spanning half the spectrum does.
- Median cut rather than k-means because it is **deterministic**: a palette that depended on the
  random stream could be neither tested nor documented.
- Channel count and dimensions are preserved; alpha passes through.
- Asking for more colours than the image holds cannot invent any.

**Attributes:** `Protected`.

## References

- P. Heckbert, *Color image quantization for frame buffer display*, Computer Graphics (SIGGRAPH) **16** (1982) 297-307.
- Source: [`src/imagecolor.c`](https://github.com/stblake/mathilda/blob/main/src/imagecolor.c)
- Specification: [`docs/spec/builtins/image-processing.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/image-processing.md)
- Tests: [`tests/test_image.c`](https://github.com/stblake/mathilda/blob/main/tests/test_image.c)

## Notes & additional examples

### Notes

`ColorQuantize[image, n]` reduces the image to at most `n` colours by **median
cut**: the box with the widest single-channel spread is split at its median until
`n` boxes remain, and each box collapses to its mean colour.

Widest spread, not most pixels, because a large box of nearly identical colours
does not need splitting and a small one spanning half the spectrum does. Median
cut, not k-means, because it is **deterministic** — a palette that depended on a
random stream could not be tested or documented. The channel count is preserved
and an alpha channel passes through.
