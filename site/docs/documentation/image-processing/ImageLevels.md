# ImageLevels

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`ImageLevels[image] gives {{level, count}, ...}: the histogram as DATA, not a plot -- use Histogram over the result for a picture. ImageLevels[image, n] uses n bins. Levels are on the same unit scale as ImageData, so a level can be compared against a pixel value without rescaling. A "Bit" image uses its 2 natural levels and "Byte" its 256, because those ARE the distinct values; a "Real" image has no natural set and is binned into 256 over [0, 1]. The counts sum to the pixel count exactly, every pixel landing in one bin. Accepts volumes as well as planes.`**

## Examples (6)

Every input below was run against the current Mathilda build and its output recorded.

### Basic Examples (2)

```mathematica
In[1]:= bit = Image[Table[Boole[Mod[i + j, 2] == 0], {i, 1, 8}, {j, 1, 8}]];

In[2]:= ImageLevels[bit]
Out[2]= {{0.0, 32}, {1.0, 32}}
```

### Scope (1)

```mathematica
In[3]:= ImageLevels[Image[{{0, 1}, {1, 0}}]]
Out[3]= {{0.0, 2}, {1.0, 2}}
```

### Applications (3)

```mathematica
In[4]:= ImageLevels[Image[{{0, 1}, {1, 1}}, "Bit"]]
Out[4]= {{0.0, 1}, {1.0, 3}}

In[5]:= ImageLevels[Image[{{0., 0.5}, {0.5, 1.}}], 2]
Out[5]= {{0.0, 1}, {1.0, 3}}

In[6]:= Total[Last /@ ImageLevels[Image[{{0, 1}, {1, 1}}, "Bit"]]]
Out[6]= 4
```

## Implementation notes

**Algorithm.** `builtin_imagelevels` returns the pixel histogram as **data**, a
list of `{level, count}` pairs — not a plot (`Histogram` over the result gives
the picture). `ImageLevels[image]` chooses the bin count from the type: a `"Bit"`
image uses its 2 natural levels, a `"Byte"` image its 256, because those *are*
the distinct values; only a `"Real"` image, which has no natural set, is binned
into 256 over `[0, 1]` (`IMG_LEVEL_BINS`). `ImageLevels[image, n]` forces `n`
bins. It loads the whole buffer (plane or volume), clamps each value to `[0, 1]`,
maps it to a bin with `floor(v · (nbins − 1) + 0.5)`, and accumulates a count.
The returned level for bin `i` is `i/(nbins − 1)`, on the **same unit scale**
`ImageData` uses, so a level is directly comparable to a pixel value without
rescaling. Every pixel lands in exactly one bin, so the counts sum to the pixel
count exactly.

**Data structures.** A `calloc`'d `double` count array of `nbins`, filled by one
pass over the flat unit-scale buffer (`image_load` / `image3d_load`), then
marshalled into a `List` of two-element `{real level, integer count}` pairs.

**Complexity / limits.** `O(n + nbins)`, `n` the pixel-channel count. `n` bins
must be between `2` and `65536`. Accepts volumes as well as planes; a colour
image's channels are binned together into the same histogram.

**Attributes:** `Protected`.

## References

**See also:** [Histogram](../../graphics/Histogram/), [ImageData](../../image-processing/ImageData/)

- Source: [`src/imagefilter.c`](https://github.com/stblake/mathilda/blob/main/src/imagefilter.c)
- Specification: [`docs/spec/builtins/image-processing.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/image-processing.md)
- Tests: [`tests/test_image.c`](https://github.com/stblake/mathilda/blob/main/tests/test_image.c)

## Notes & additional examples

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
