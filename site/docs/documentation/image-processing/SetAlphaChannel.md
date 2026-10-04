# SetAlphaChannel

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`SetAlphaChannel[image] attaches a fully opaque alpha channel. SetAlphaChannel[image, a] sets one opacity everywhere when a is a number in [0, 1], or per pixel when a is an image of the same dimensions (read as grey, so a colour mask is not taken as its red channel alone). A mask of the wrong size is declined rather than resampled.`**

## Examples (12)

Every input below was run against the current Mathilda build and its output recorded.

### Basic Examples (5)

```mathematica
In[1]:= a = Image[Table[N[(i + j)/32], {i, 1, 16}, {j, 1, 16}], "Real"];

In[2]:= ImageChannels[SetAlphaChannel[a, 0.5]]
Out[2]= 2

In[3]:= ImageChannels[SetAlphaChannel[Image[Table[{0.5, 0.2, 0.9}, {i, 1, 8}, {j, 1, 8}], "Real"], 0.5]]
Out[3]= 4

In[4]:= SetAlphaChannel[a, Image[Table[N[j/16], {i, 1, 16}, {j, 1, 16}], "Real"]]
Out[4]= -Image-

In[5]:= Head[SetAlphaChannel[a, Image[{{0.5}}, "Real"]]]
Out[5]= SetAlphaChannel
```

### Properties & Relations (3)

```mathematica
In[6]:= a = Image[Table[N[(i + j)/32], {i, 1, 16}, {j, 1, 16}], "Real"];
```

A colour mask is averaged, not read as its red channel: {0.2, 0.4, 0.6} gives 0.4

```mathematica
In[7]:= Module[{m = Image[Table[{0.2, 0.4, 0.6}, {i, 1, 16}, {j, 1, 16}], "Real"], u}, u = Union[Flatten[ImageData[AlphaChannel[SetAlphaChannel[a, m]]]]]; Round[First[u], 0.0001]]
Out[7]= 0.4
```

Setting then removing gets back the channel count it started with

```mathematica
In[8]:= ImageChannels[RemoveAlphaChannel[SetAlphaChannel[a, 0.5]]] === ImageChannels[a]
Out[8]= True
```

### Applications (4)

A grey image gains an alpha channel: grey+alpha = 2

```mathematica
In[9]:= ImageChannels[SetAlphaChannel[Image[{{0., 1.}}]]]
Out[9]= 2
```

An RGB image gains alpha: 4 channels

```mathematica
In[10]:= ImageChannels[SetAlphaChannel[Image[{{{1., 0., 0.}}}]]]
Out[10]= 4
```

One opacity everywhere, written into the new channel

```mathematica
In[11]:= ImageData[SetAlphaChannel[Image[{{0., 1.}}], 0.5]]
Out[11]= {{{0.0, 0.5}, {1.0, 0.5}}}
```

A same-size mask sets per-pixel opacity

```mathematica
In[12]:= ImageData[SetAlphaChannel[Image[{{0., 1.}}], Image[{{0.25, 0.75}}]]]
Out[12]= {{{0.0, 0.25}, {1.0, 0.75}}}
```

## Implementation notes

**Algorithm.** `builtin_setalphachannel` attaches or replaces an image's opacity. The result
always carries one more channel than the source's colour-channel count (`colour_count` strips
any existing alpha first, so the output is colour+alpha, never doubled). The new opacity is
either:

- a constant — `SetAlphaChannel[image, a]` with `a` a number in `[0, 1]`, written into the
  alpha slot of every pixel; or
- a mask — `SetAlphaChannel[image, maskimage]`, which must match the source in size (a
  mismatched mask is declined, not silently resampled) and is read as **grey**: its own colour
  channels are averaged, so a colour mask is not quietly taken as its red channel alone; or
- fully opaque `1.0` when no second argument is given.

Colour samples are copied through `sample`, which replicates a grey source across colour
channels so the operation has one code path for grey and colour inputs.

**Data structures.** One or two decoded `Img` buffers in; a fresh
`width · height · (colour_count + 1)` buffer out, wrapped by `image_build_real` as a packed
`"Real"` image (every image head returns a packed buffer — `make check-image-packing`).

**Complexity / limits.** `O(width · height · channels)`, one pass. A mask of the wrong
dimensions is refused; the opacity constant must lie in `[0, 1]`.

- `Protected`.
- A mask is read as **grey** (its channels averaged), so a colour mask is not silently taken as
  its red channel alone.
- A mask of the wrong size is declined rather than resampled: a mismatch is a mistake, not a
  request to interpolate.

**Attributes:** `Protected`.

## References

- Source: [`src/imagecompose.c`](https://github.com/stblake/mathilda/blob/main/src/imagecompose.c)
- Specification: [`docs/spec/builtins/image-processing.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/image-processing.md)
- Tests: [`tests/test_image.c`](https://github.com/stblake/mathilda/blob/main/tests/test_image.c)

## Notes & additional examples

### Notes

`SetAlphaChannel[image]` attaches a fully opaque alpha; `SetAlphaChannel[image, a]` sets one
opacity everywhere (a number in `[0, 1]`); `SetAlphaChannel[image, mask]` sets it per pixel
from a same-size image read as grey — a colour mask is averaged rather than taken as its red
channel alone. A mask of the wrong size is declined, not resampled silently.

The result always carries exactly one alpha channel: an existing alpha is replaced, not
duplicated, so grey → grey+alpha (2 channels) and RGB → RGBA (4).
