# Image

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Image[data] is a raster image, normalising to the canonical Image[data, type]. The data is a rectangular height x width array of pixel values, or height x width x channels for a colour image, so it is indexed data[[y, x]] with rows running down the image -- note that ImageDimensions reports {width, height}, transposed relative to this. The type is inferred from the values: all-integer data in {0, 1} is "Bit", all-integer in 0..255 is "Byte", anything else is "Real". Image[data, type] states the type instead -- "Bit" (0 or 1), "Byte" (0..255), "Bit16" (0..65535) or "Real" ("Real32" and "Real64" are accepted as synonyms) -- and, as in Mathematica, COERCES the data to it: an integer type rounds each value to the nearest integer and clips it to the type's range, so Image[{{0, 300}}, "Byte"] stores {0, 255}; "Real" keeps any real value. Image[image, type] converts an image between types, preserving brightness. Data that is not a rectangular array of real numbers (ragged, non-numeric or complex) is left unevaluated.`**

## Examples (43)

Every input below was run against the current Mathilda build and its output recorded.

### Basic Examples (6)

```mathematica
In[1]:= Image[{{0., 1.}, {1., 0.}}]
Out[1]= -Image-

In[2]:= Image[{{0., 0.25, 0.5, 0.75, 1.}}]
Out[2]= -Image-

In[3]:= ImageDimensions[Image[{{0., 1., 0.5}, {1., 0., 0.25}}]]
Out[3]= {3, 2}

In[4]:= ImageChannels[Image[{{0., 1.}, {1., 0.}}]]
Out[4]= 1

In[5]:= ImageType[Image[{{0., 1.}, {1., 0.}}]]
Out[5]= "Real"

In[6]:= ImageData[Image[{{0., 1.}, {1., 0.}}]]
Out[6]= {{0.0, 1.0}, {1.0, 0.0}}
```

### Scope (32)

The type is inferred from the values

```mathematica
In[7]:= ImageType[Image[{{0, 1}, {1, 0}}]]
Out[7]= "Bit"
```

```mathematica
In[8]:= ImageType[Image[{{0, 128}, {255, 7}}]]
Out[8]= "Byte"

In[9]:= ImageType[Image[{{0., 0.5}}]]
Out[9]= "Real"
```

A stated type coerces the data: 300 clips to 255

```mathematica
In[10]:= ImageData[Image[{{0, 300}}, "Byte"], "Byte"]
Out[10]= {{0, 255}}
```

```mathematica
In[11]:= ImageData[Image[{{0.5, 2, 0.4}}, "Bit"], "Bit"]
Out[11]= {{1, 1, 0}}

In[12]:= ImageType[Image[{{0, 1}, {1, 0}}, "Byte"]]
Out[12]= "Byte"
```

Ragged data declines rather than being padded

```mathematica
In[13]:= Head[Image[{{1., 2.}, {3.}}]]
Out[13]= Image
```

```mathematica
In[14]:= Head[Image[{}]]
Out[14]= Image
```

Data is indexed [[y, x]] while ImageDimensions reports {width, height}

```mathematica
In[15]:= Dimensions[ImageData[Image[{{1., 2., 3.}, {4., 5., 6.}}]]]
Out[15]= {2, 3}
```

```mathematica
In[16]:= ImageDimensions[Image[{{1., 2., 3.}, {4., 5., 6.}}]]
Out[16]= {3, 2}

In[17]:= ImageData[Image[{{1., 2., 3.}, {4., 5., 6.}}]][[1, 3]]
Out[17]= 3.0
```

Three channels make a colour image

```mathematica
In[18]:= ImageChannels[Image[{{{1., 0., 0.}, {0., 1., 0.}}}]]
Out[18]= 3
```

```mathematica
In[19]:= Image[{{{1., 0., 0.}, {0., 1., 0.}, {0., 0., 1.}}}]
Out[19]= -Image-

In[20]:= ImageDimensions[Image[{{{1., 0., 0.}, {0., 1., 0.}}}]]
Out[20]= {2, 1}
```

Two or four channels carry alpha

```mathematica
In[21]:= ImageChannels[Image[{{{1., 0.5}, {0., 1.}}}]]
Out[21]= 2
```

```mathematica
In[22]:= ImageChannels[Image[{{{1., 0., 0., 0.5}}}]]
Out[22]= 4
```

ImageData reports STORED values with an explicit type

```mathematica
In[23]:= ImageData[Image[{{0, 255}}, "Byte"], "Byte"]
Out[23]= {{0, 255}}
```

```mathematica
In[24]:= ImageData[Image[{{0, 255}}, "Byte"]]
Out[24]= {{0.0, 1.0}}

In[25]:= ImageData[Image[{{1, 0}, {0, 1}}, "Bit"], "Bit"]
Out[25]= {{1, 0}, {0, 1}}
```

An already-canonical image is left alone, so evaluation reaches a fixed point

```mathematica
In[26]:= Image[Image[{{0., 1.}}]] === Image[{{0., 1.}}]
Out[26]= True
```

Pixels survive a round trip through ImageData exactly

```mathematica
In[27]:= Module[{d = {{0.1, 0.2}, {0.3, 0.4}}}, ImageData[Image[d]] === d]
Out[27]= True
```

```mathematica
In[28]:= Image[Table[N[i j]/9, {i, 3}, {j, 3}]]
Out[28]= -Image-

In[29]:= Image[Table[N[Mod[i + j, 2]], {i, 8}, {j, 8}]]
Out[29]= -Image-

In[30]:= Image[Table[N[(i - 1)/7], {i, 8}, {j, 8}]]
Out[30]= -Image-
```

A colour ramp

```mathematica
In[31]:= Image[Table[{N[(j - 1)/7], N[(i - 1)/7], 0.5}, {i, 8}, {j, 8}]]
Out[31]= -Image-
```

```mathematica
In[32]:= Image[Table[If[Mod[Quotient[i - 1, 2] + Quotient[j - 1, 2], 2] == 0, 0., 1.], {i, 8}, {j, 8}]]
Out[32]= -Image-

In[33]:= Image[Table[N[Boole[(i - 4.5)^2 + (j - 4.5)^2 <= 9]], {i, 8}, {j, 8}]]
Out[33]= -Image-

In[34]:= ImageQ[Image[{{0., 1.}}]]
Out[34]= True

In[35]:= ImageQ[{{0., 1.}}]
Out[35]= False

In[36]:= Image3DQ[Image[{{0., 1.}}]]
Out[36]= False
```

The storage is a packed buffer, not a tree of Expr nodes

```mathematica
In[37]:= Head[Part[Image[Table[N[i j]/64, {i, 8}, {j, 8}]], 1]]
Out[37]= NDArray
```

```mathematica
In[38]:= Part[Image[{{0., 1.}}], 2]
Out[38]= "Real"
```

### Applications (5)

A small real-typed checkerboard, stored as a packed buffer

```mathematica
In[39]:= Image[{{0., 1.}, {1., 0.}}]
Out[39]= -Image-
```

All values in {0,1}, so the type infers to Bit

```mathematica
In[40]:= ImageType[Image[{{0, 1}, {1, 0}}]]
Out[40]= "Bit"
```

Integers to 255 infer to Byte

```mathematica
In[41]:= ImageType[Image[{{0, 128, 255}}]]
Out[41]= "Byte"
```

A stated type rounds and clips: 300 -> 255, scaled out

```mathematica
In[42]:= ImageData[Image[{{0, 300}}, "Byte"]]
Out[42]= {{0.0, 1.0}}
```

Width x height, transposed from the data's rows

```mathematica
In[43]:= ImageDimensions[Image[{{0., 1., 0.}}]]
Out[43]= {3, 1}
```

## Implementation notes

**Algorithm.** `builtin_image` calls the shared constructor `img_construct`, which turns
every spelling into the canonical two-argument form `Image[data, type]`:

1. `Image[data]` — infer the type from the values (`img_shape` collects the range and an
   all-integer flag in one pass): all-integer in `{0, 1}` is `"Bit"`, all-integer in `0..255`
   is `"Byte"`, anything else is `"Real"`. Inference looks only at the data, so it is
   predictable — `Image[{{0,1}}]` is a bit image, `Image[{{0.,1.}}]` a real one.
2. `Image[data, type]` — coerce the data to the stated type. `img_coerce` rounds each value
   to the nearest integer (halves up) and clips to `[0, max]` for the integer types
   (`"Bit"` 1, `"Byte"` 255, `"Bit16"` 65535), and keeps any real for `"Real"`; `"Real32"`
   and `"Real64"` are accepted as synonyms of `"Real"`.
3. `Image[image, type]` — a conversion: the source is read in unit scale and re-quantised to
   the target type, so brightness is preserved.

The shape walk (`img_shape` / `img_shape_fast`) rejects a ragged array rather than padding
it — a ragged array is not an image, and every downstream filter indexes it as rectangular,
so a clear refusal here replaces an out-of-bounds read far away. Complex-valued data is
refused. `ImageData` reports stored values scaled back to the unit interval (`img_to_unit`):
a `"Byte"` 255 comes back as exactly `1.0`.

**Data structures.** The canonical node is `Image[data, "type"]`. `data` is a VISIBLE
`NDArray` (`present_as = NDA_HEAD_NDARRAY`) whenever `ndbuild_open` can pack it, otherwise the
equivalent nested `List`s; it is row-major, height × width (× channels, channels innermost).
The dtype follows the pixel type — `int64` for `"Bit"`/`"Byte"`/`"Bit16"`, `float64` for
`"Real"` — because `ImageData` reports stored values and a byte 200 must print as `200`, not
`200.`. Storing on the visible surface is deliberate: the evaluator's post-gate materialises a
resting *packed List* but never touches a visible `NDArray`, so a container that comes to rest
holding its pixels keeps its buffer (this is why `make check-image-packing` audits that every
image head hands back a packed buffer).

**Complexity / limits.** Construction is `O(pixels)` — every pixel is validated and coerced
once. A query such as `ImageDimensions` then uses `img_shape_fast`, which re-checks only
rectangularity (`O(height)`) and reads a packed shape in `O(1)`, so repeated size queries in a
filter pipeline do not re-walk the pixels. Data that is not a rectangular array of real
numbers is left unevaluated, which is what makes validity decidable by `ImageQ`.

**Attributes:** `Protected`.

## References

**See also:** [ImageQ](../../image-processing/ImageQ/), [ImageDimensions](../../image-processing/ImageDimensions/), [ImageData](../../image-processing/ImageData/), [ImageType](../../image-processing/ImageType/)

- Source: [`src/image.c`](https://github.com/stblake/mathilda/blob/main/src/image.c)
- Specification: [`docs/spec/builtins/image-processing.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/image-processing.md)
- Tests: [`tests/test_image.c`](https://github.com/stblake/mathilda/blob/main/tests/test_image.c)

## Notes & additional examples

### Notes

`Image[data]` normalises to the canonical `Image[data, type]`, which is also real Wolfram
syntax; `ImageQ` tests for exactly that form. The type is inferred from the values alone, so
`Image[{{0,1}}]` is a `"Bit"` image and `Image[{{0.,1.}}]` a `"Real"` one — a distinction a
caller can rely on.

The data is stored height × width (× channels), which is **transposed** relative to
`ImageDimensions`'s `{width, height}` — the single most common source of silently-wrong image
code, and Mathematica's convention. `ImageData` scales stored values back into `[0, 1]` using
the type's range, so the type is not decoration: it is what makes that scaling well defined.

Computed and constructed images share one representation — a visible packed `NDArray` — so
`===` compares two images with identical pixels as equal regardless of how each was built.
