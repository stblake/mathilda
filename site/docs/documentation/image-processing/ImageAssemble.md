# ImageAssemble

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`ImageAssemble[{{a, b}, {c, d}}] tiles a grid of images into one; ImageAssemble[{a, b}] makes a single row. Each tile keeps its natural size -- a row is as tall as its tallest tile and a column as wide as its widest, and any gap is left blank rather than stretched, since stretching would resample an image the caller did not ask to resize. Alpha survives if any tile had it.`**

## Examples (10)

Every input below was run against the current Mathilda build and its output recorded.

### Basic Examples (4)

```mathematica
In[1]:= a = Image[Table[N[(i + j)/32], {i, 1, 16}, {j, 1, 16}], "Real"];

In[2]:= ImageAssemble[{a, a}]
Out[2]= -Image-

In[3]:= ImageDimensions[ImageAssemble[{a, a}]]
Out[3]= {32, 16}

In[4]:= ImageDimensions[ImageAssemble[{{a, a}, {a, a}}]]
Out[4]= {32, 32}
```

### Applications (3)

```mathematica
In[5]:= a = Image[Table[N[Boole[(i - 16)^2 + (j - 16)^2 <= 100]], {i, 1, 32}, {j, 1, 32}], "Real"];
```

A contact sheet comparing one filter at four radii

```mathematica
In[6]:= ImageAssemble[{Table[GaussianFilter[a, r], {r, 1, 2}], Table[GaussianFilter[a, r], {r, 3, 4}]}]
Out[6]= -Image-
```

The same image before and after, side by side

```mathematica
In[7]:= ImageAssemble[{a, EdgeDetect[a]}]
Out[7]= -Image-
```

### Applications (3)

Two 2x1 tiles side by side: width 4

```mathematica
In[8]:= ImageDimensions[ImageAssemble[{Image[{{0., 1.}}], Image[{{1., 0.}}]}]]
Out[8]= {4, 1}
```

The tiles laid out left to right

```mathematica
In[9]:= ImageData[ImageAssemble[{Image[{{0., 1.}}], Image[{{1., 0.}}]}]]
Out[9]= {{0.0, 1.0, 1.0, 0.0}}
```

A 2x2 grid

```mathematica
In[10]:= ImageDimensions[ImageAssemble[{{Image[{{0.}}], Image[{{1.}}]}, {Image[{{1.}}], Image[{{0.}}]}}]]
Out[10]= {2, 2}
```

## Implementation notes

**Algorithm.** `builtin_imageassemble` tiles a grid of images into one.
`ImageAssemble[{{a, b}, {c, d}}]` is a 2-D grid; `ImageAssemble[{a, b}]` a single row (a flat
list is one row, the reading order a list already implies — distinguished by probing whether
the first element is itself an image). Each tile is decoded with `image_load`. Row heights and
column widths are taken as the maximum over the row/column (`row_h[i]`, `col_w[j]`), so tiles
of different sizes line up on a grid; the total canvas is the sum of those. Each tile is
blitted at its natural size into its cell — nothing is stretched, because stretching would
resample an image the caller did not ask to resize, and any leftover gap is left blank.

The combined colour-channel count is the maximum over all tiles (`promote`), with grey
replicated into colour via `sample`. An alpha channel is added to the sheet iff **any** tile
had one, and where a sheet has alpha its gaps are transparent rather than opaque black — the
one property a sprite sheet must keep.

**Data structures.** An array of decoded `Img` tiles plus parallel `row_h`/`col_w` arrays; a
single `calloc`'d `total_w · total_h · out_channels` canvas, wrapped by `image_build_real` as a
packed `"Real"` image (every image head returns a packed buffer — `make check-image-packing`).

**Complexity / limits.** `O(total pixels)` to clear the canvas plus `O(sum of tile pixels)` to
copy; a ragged grid is accepted (short rows simply leave cells blank). A non-image element
makes the whole call decline.

- `Protected`.
- Each tile keeps its **natural size**: a row is as tall as its tallest tile, a column as wide as
  its widest, and any gap is left blank rather than stretched — stretching would resample an image
  the caller did not ask to resize.
- Channels are promoted as in `ImageCompose`, and alpha survives if any tile had it.

**Attributes:** `Protected`.

## References

**See also:** [ImageCompose](../../image-processing/ImageCompose/)

- Source: [`src/imagecompose.c`](https://github.com/stblake/mathilda/blob/main/src/imagecompose.c)
- Specification: [`docs/spec/builtins/image-processing.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/image-processing.md)
- Tests: [`tests/test_image.c`](https://github.com/stblake/mathilda/blob/main/tests/test_image.c)

## Notes & additional examples

### Notes

`ImageAssemble[{{a, b}, {c, d}}]` tiles a grid; `ImageAssemble[{a, b}]` makes a single row. A
row is as tall as its tallest tile and a column as wide as its widest; each tile keeps its
natural size and any gap is left blank rather than stretched, since stretching would resample
an image the caller did not ask to resize.

A grey tile beside a colour one is promoted to colour, and the assembled sheet carries an alpha
channel if any tile did — with its gaps transparent, the property a sheet of sprites needs to
keep.
