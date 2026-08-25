# ImageAssemble

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`ImageAssemble[{{a, b}, {c, d}}] tiles a grid of images into one; ImageAssemble[{a, b}] makes a single row. Each tile keeps its natural size -- a row is as tall as its tallest tile and a column as wide as its widest, and any gap is left blank rather than stretched, since stretching would resample an image the caller did not ask to resize. Alpha survives if any tile had it.`**

## Examples (7)

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

## Implementation notes

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
