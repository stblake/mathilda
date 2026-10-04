# RemoveAlphaChannel

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`RemoveAlphaChannel[image] drops the alpha channel. RemoveAlphaChannel[image, b] instead COMPOSITES over a background of brightness b, which is the difference between forgetting the transparency and resolving it: a half-transparent white pixel over black is grey, where dropping alpha would leave it white.`**

## Examples (7)

Every input below was run against the current Mathilda build and its output recorded.

### Basic Examples (4)

```mathematica
In[1]:= t = SetAlphaChannel[Image[{{1.0, 1.0}, {1.0, 1.0}}, "Real"], 0.5];

In[2]:= ImageChannels[RemoveAlphaChannel[t]]
Out[2]= 1

In[3]:= Round[Max[Flatten[ImageData[RemoveAlphaChannel[t]]]], 0.001]
Out[3]= 1.0

In[4]:= Round[Max[Flatten[ImageData[RemoveAlphaChannel[t, 0.]]]], 0.001]
Out[4]= 0.5
```

### Applications (3)

```mathematica
In[5]:= ImageChannels[RemoveAlphaChannel[Image[{{{1., 0, 0, 0.5}}}]]]
Out[5]= 3

In[6]:= ImageData[RemoveAlphaChannel[Image[{{{1., 1., 1., 0.5}}}], 0]]
Out[6]= {{{0.5, 0.5, 0.5}}}

In[7]:= ImageData[RemoveAlphaChannel[Image[{{{1., 0., 0., 0.5}}}]]]
Out[7]= {{{1.0, 0.0, 0.0}}}
```

## Implementation notes

**Algorithm.** `builtin_removealphachannel` drops an image's alpha channel.
`RemoveAlphaChannel[image]` simply discards it, keeping each colour channel's
stored value. `RemoveAlphaChannel[image, b]` instead **composites over** a
background of brightness `b`: for each pixel it reads the alpha `al` (1.0 for an
image with no alpha) and writes `v · al + b · (1 − al)` per colour channel —
which is the difference between *forgetting* the transparency and *resolving*
it. A half-transparent white pixel over black (`b = 0`) becomes grey, where
plain dropping would leave it white. The output channel count is
`colour_count(im.c)` — the colour channels without the alpha — so a 4-channel
RGBA image becomes 3-channel RGB and a 2-channel grey-plus-alpha becomes
1-channel grey.

**Data structures.** The image is read into an `Img` struct via `img_read`;
`sample(&im, x, y, k, cc)` reads channel `k`, and `has_alpha`/`colour_count`
derive the alpha presence and colour-channel count from the channel count. The
output is a flat `double` buffer of `width · height · colour_count`, built into a
`"Real"` image (`image_build_real`).

**Complexity / limits.** `O(width · height · channels)` — a single per-pixel
pass. The composite form needs a scalar background brightness; an image with no
alpha channel is returned with its colour channels intact (the composite reads
`al = 1`, leaving values unchanged).

- `Protected`.
- The two forms are genuinely different: a half-transparent white pixel over black is **grey**,
  where dropping alpha leaves it white. One forgets the transparency; the other resolves it.

**Attributes:** `Protected`.

## References

- Source: [`src/imagecompose.c`](https://github.com/stblake/mathilda/blob/main/src/imagecompose.c)
- Specification: [`docs/spec/builtins/image-processing.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/image-processing.md)
- Tests: [`tests/test_image.c`](https://github.com/stblake/mathilda/blob/main/tests/test_image.c)

## Notes & additional examples

### Notes

`RemoveAlphaChannel[image]` drops the alpha channel.
`RemoveAlphaChannel[image, b]` instead **composites** over a background of
brightness `b`, which is the difference between forgetting the transparency and
resolving it: a half-transparent white pixel over black is grey, where dropping
alpha would leave it white.

The composite writes `v · a + b · (1 − a)` per colour channel. The output keeps
only the colour channels — a 4-channel RGBA image becomes 3-channel RGB, a
2-channel grey-plus-alpha becomes 1-channel grey.
