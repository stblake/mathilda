# AlphaChannel

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`AlphaChannel[image] gives the image's opacity as a one-channel image. An image with no alpha channel answers with an all-opaque one rather than declining: "how transparent is this?" has an answer for every image, and it is "not at all". Two channels are read as grey+alpha and four as RGB+alpha.`**

## Examples (7)

Every input below was run against the current Mathilda build and its output recorded.

### Basic Examples (4)

```mathematica
In[1]:= a = Image[Table[N[(i + j)/32], {i, 1, 16}, {j, 1, 16}], "Real"];

In[2]:= Union[Flatten[ImageData[AlphaChannel[a]]]]
Out[2]= {1.0}

In[3]:= Union[Flatten[ImageData[AlphaChannel[SetAlphaChannel[a, 0.25]]]]]
Out[3]= {0.25}

In[4]:= ImageChannels[AlphaChannel[Image[Table[{0.2, 0.4, 0.6}, {i, 1, 8}, {j, 1, 8}], "Real"]]]
Out[4]= 1
```

### Applications (3)

An RGBA pixel's alpha is channel 4

```mathematica
In[5]:= ImageData[AlphaChannel[Image[{{{1., 0., 0., 0.5}}}]]]
Out[5]= {{0.5}}
```

A grey image has no alpha: the result is all-opaque

```mathematica
In[6]:= AlphaChannel[Image[{{0., 1.}, {1., 0.}}]]
Out[6]= -Image-
```

Which reads back as 1 everywhere

```mathematica
In[7]:= ImageData[AlphaChannel[Image[{{0., 1.}, {1., 0.}}]]]
Out[7]= {{1.0, 1.0}, {1.0, 1.0}}
```

## Implementation notes

**Algorithm.** `builtin_alphachannel` returns an image's opacity as a one-channel image. It
decodes the input to a unit-interval buffer with `image_load`, then reads the opacity per
pixel: a two-channel image is grey+alpha (channel 1), a four-channel image is RGB+alpha
(channel 3), and any other channel count has no alpha, so it answers `1.0` — fully opaque.
That last case is a decision, not a decline: "how transparent is this?" has an answer for
every image, and for one without an alpha channel the answer is "not at all". The opacity
buffer is wrapped back up with `image_build_real`.

**Data structures.** One decoded `Img` (`{w, h, c, buf}` of unit-interval reals, row-major,
channels innermost) in, a fresh `width · height` buffer out. The result is a `"Real"`
single-channel image built as a packed `NDArray` by `image_build_real`, which is what keeps
every image head's output a packed buffer (enforced by `make check-image-packing`).

**Complexity / limits.** `O(width · height)` — one pass, picking one sample per pixel. The
grey/alpha and RGB/alpha channel conventions are the only ones recognised; a 3-channel image
is treated as alpha-free and reads back all-opaque.

- `Protected`.
- An image with **no** alpha channel answers with an all-opaque one rather than declining: "how
  transparent is this?" has an answer for every image, and it is "not at all".
- Two channels are read as grey+alpha, four as RGB+alpha.

**Attributes:** `Protected`.

## References

- Source: [`src/imagecompose.c`](https://github.com/stblake/mathilda/blob/main/src/imagecompose.c)
- Specification: [`docs/spec/builtins/image-processing.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/image-processing.md)
- Tests: [`tests/test_image.c`](https://github.com/stblake/mathilda/blob/main/tests/test_image.c)

## Notes & additional examples

### Notes

`AlphaChannel` extracts the opacity as a one-channel image. Two channels are read as
grey+alpha and four as RGB+alpha; any other channel count — a grey or an RGB image — is
treated as fully opaque, and the result is an all-`1` image rather than a decline. "How
transparent is this?" has an answer for every image.

The result is always a `"Real"` single-channel image handed back as a packed buffer, so it
feeds straight into another filter.
