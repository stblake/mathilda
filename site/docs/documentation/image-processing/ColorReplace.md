# ColorReplace

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`ColorReplace[image, old -> new] replaces every pixel within a tolerance of `old` by `new`; ColorReplace[image, {r1, r2, ...}] applies several rules and ColorReplace[image, rules, tol] sets the tolerance (default 0.02 -- at 0 only bit-identical colours match, which after any filtering is nothing at all). Distance is Euclidean in RGB, and where rules overlap the NEAREST wins rather than the first, so the answer does not depend on the order they were written. Colours may be RGBColor[r, g, b], GrayLevel[v], a number or {r, g, b}. Replacing a grey image's colour with a non-grey one produces a three-channel image, since flattening the new colour to its luminance would give grey when the caller asked for red. An alpha channel passes through: transparency is not a colour.`**

## Examples (13)

Every input below was run against the current Mathilda build and its output recorded.

### Basic Examples (4)

```mathematica
In[1]:= img = Image[Table[{If[j <= 4, 1., 0.], If[j <= 4, 0., 1.], 0.}, {i, 1, 8}, {j, 1, 8}], "Real"];

In[2]:= ColorReplace[img, RGBColor[1, 0, 0] -> RGBColor[0, 0, 1]]
Out[2]= -Image-

In[3]:= Union[Flatten[ImageData[ColorReplace[img, RGBColor[1, 0, 0] -> RGBColor[0, 0, 1]]]]]
Out[3]= {0.0, 1.0}

In[4]:= ColorReplace[img, {RGBColor[1, 0, 0] -> RGBColor[0, 1, 0], RGBColor[0, 1, 0] -> RGBColor[0, 0, 1]}]
Out[4]= -Image-
```

### Properties & Relations (6)

```mathematica
In[5]:= img = Image[Table[{If[j <= 4, 1., 0.], If[j <= 4, 0., 1.], 0.}, {i, 1, 8}, {j, 1, 8}], "Real"];

In[6]:= grey = Image[Table[N[j/8], {i, 1, 8}, {j, 1, 8}], "Real"];
```

A colour nothing matches leaves the image exactly alone

```mathematica
In[7]:= ImageData[ColorReplace[img, RGBColor[0.5, 0.5, 0.5] -> RGBColor[0, 0, 0]]] === ImageData[img]
Out[7]= True
```

Grey replaced by colour promotes to three channels; grey by grey does not

```mathematica
In[8]:= {ImageChannels[ColorReplace[grey, GrayLevel[0.5] -> RGBColor[1, 0, 0], 0.1]], ImageChannels[ColorReplace[grey, GrayLevel[0.5] -> GrayLevel[0.], 0.1]]}
Out[8]= {3, 1}
```

The colour forms are interchangeable

```mathematica
In[9]:= ImageData[ColorReplace[grey, 0.5 -> 0., 0.1]] === ImageData[ColorReplace[grey, GrayLevel[0.5] -> GrayLevel[0.], 0.1]]
Out[9]= True
```

A tolerance wide enough to reach every colour collapses the image to one

```mathematica
In[10]:= Length[Union[Flatten[ImageData[ColorReplace[img, RGBColor[1, 0, 0] -> RGBColor[0.25, 0.25, 0.25], 2.0]]]]]
Out[10]= 1
```

### Applications (3)

The bright pixel is recoloured to 0

```mathematica
In[11]:= ImageData[ColorReplace[Image[{{0., 1.}}], 1. -> 0.]]
Out[11]= {{0.0, 0.0}}
```

Grey -> colour when the target is non-grey

```mathematica
In[12]:= ImageChannels[ColorReplace[Image[{{1.}}], 1. -> RGBColor[1., 0., 0.]]]
Out[12]= 3
```

Only pixels within 0.1 of 0.5 change

```mathematica
In[13]:= ImageData[ColorReplace[Image[{{0., 0.5, 1.}}], 0.5 -> 0., 0.1]]
Out[13]= {{0.0, 0.0, 1.0}}
```

## Algorithm

imagecolor.c -- ColorReplace, ColorQuantize and HistogramTransform.

Three heads that act on an image's COLOURS rather than its geometry, and they share the one thing that makes such operations awkward: a decision made per pixel needs a global view first. Replacing a colour needs a distance rule, quantising needs a palette derived from every pixel, and equalising needs the whole distribution. So each of these makes a pass to gather, then a pass to write — which is why none of them fits the filter machinery in imagefilter.c.

## Implementation notes

**Algorithm.** `builtin_colorreplace` recolours pixels within a tolerance of a target colour.
The rule(s) — `old -> new`, or a list of them — are parsed by `read_rules`, which accepts
`RGBColor[r, g, b]`, `GrayLevel[v]`, a bare number, or `{r, g, b}` and expands grey specs to
three equal components so grey and colour pixels compare in one path. The image is decoded to a
unit buffer; the default tolerance is `0.02` (a stated constant, because at 0 only a
bit-identical colour matches, which after any filtering is nothing). For each pixel the
Euclidean RGB distance to every rule's source is measured, and the **nearest** rule within
tolerance wins — not the first, so overlapping rules do not make the answer depend on the order
they were written. The matched pixel takes the rule's target; an unmatched pixel is copied
through, and an alpha channel always passes through untouched (transparency is not a colour).

If any target is non-grey while the input is grey, the output is promoted to three channels —
flattening the new colour to its luminance would give grey where the caller asked for red.

**Data structures.** Fixed-capacity rule tables `from[64][3]`/`to[64][3]`; one decoded unit
buffer in, a fresh `width · height · out_channels` buffer out, wrapped by `image_build_real` as
a packed `"Real"` image (every image head returns a packed buffer — `make check-image-packing`).

**Complexity / limits.** `O(width · height · rules)`. At most 64 rules; the tolerance must be
non-negative. Distance is Euclidean in RGB only — there is no perceptual colour-space metric.

- `Protected`. Colours may be `RGBColor[r, g, b]`, `GrayLevel[v]`, a bare number, or `{r, g, b}`.
- Distance is Euclidean in RGB; the default tolerance is `0.02`. At `0` only bit-identical
  colours match, which after any filtering is nothing at all.
- Where rules overlap the **nearest** wins, not the first, so the result does not depend on the
  order they were written in.
- Replacing a grey image's colour with a non-grey one produces a **three-channel** image:
  flattening the new colour to its luminance would hand back grey when the caller asked for red.
  Grey-for-grey does not promote.
- An alpha channel passes through untouched — transparency is not a colour.

**Attributes:** `Protected`.

## References

- Source: [`src/imagecolor.c`](https://github.com/stblake/mathilda/blob/main/src/imagecolor.c)
- Specification: [`docs/spec/builtins/image-processing.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/image-processing.md)
- Tests: [`tests/test_image.c`](https://github.com/stblake/mathilda/blob/main/tests/test_image.c)

## Notes & additional examples

### Notes

`ColorReplace[image, old -> new]` replaces every pixel within a tolerance (default `0.02`) of
`old` by `new`; a list of rules applies several at once. Where rules overlap the **nearest**
target wins rather than the first, so the result does not depend on the order the rules were
written. Colours may be `RGBColor[r, g, b]`, `GrayLevel[v]`, a number, or `{r, g, b}`.

Replacing a grey image's colour with a non-grey one produces a three-channel image — flattening
`red` to its luminance would give grey when the caller asked for red. An alpha channel always
passes through: transparency is not a colour.
