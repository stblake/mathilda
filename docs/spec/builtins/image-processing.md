# Image Processing

The raster image subsystem. This page currently covers the **representation and its
accessors**; filtering (`ImageConvolve`, `GaussianFilter`, `Binarize`) and `Image3D` build on
these and land separately.

## Image

`Image[data]` is a raster image, normalising to the canonical `Image[data, type]`.
Attributes: `Protected`.

- `Image[data]` — type inferred from the values
- `Image[data, "Bit" | "Byte" | "Real"]` — type stated, and validated against the data

**The canonical form is `Image[data, type]`, which is also real Wolfram syntax.** Normalising to
it is what makes validity *decidable*: a builtin that returns `NULL` leaves its expression alone,
so a valid `Image[data]` and a nonsensical `Image["hello"]` would otherwise be
indistinguishable — both merely unevaluated. A canonical two-argument form is exactly what passed
validation, and `ImageQ` tests for it.

**Data layout, and the transposition that is the whole trap.** `data` is a rectangular
`height × width` array, indexed `data[[y, x]]` with rows running *down* the image; a colour image
is `height × width × channels`, interleaved. `ImageDimensions` reports `{width, height}` —
**transposed relative to this**. That is Mathematica's convention and the single most common
source of silently-wrong image code, so the tests use a non-square image throughout: a square one
cannot tell the two apart.

**Type inference looks at the values and nothing else**, so it is predictable from the data
alone:

| data | inferred type |
|---|---|
| all integers, values within `{0, 1}` | `"Bit"` |
| all integers, values within `0..255` | `"Byte"` |
| anything else | `"Real"` |

So `Image[{{0, 1}, {1, 0}}]` is a bit image while `Image[{{0., 1.}, {1., 0.}}]` is a real one —
the same numbers written differently, and a distinction a caller can rely on.

A **stated** type is refused where the data does not fit it: `Image[{{0, 300}}, "Byte"]` stays
unevaluated rather than reinterpreting 300, which would corrupt every later scaling.

**Ragged data declines** rather than being padded or truncated. A ragged array is not an image,
and every filter downstream indexes it as rectangular — so a clear refusal here beats an
out-of-bounds read somewhere far away.

## ImageQ

`ImageQ[expr]` gives `True` for a valid canonical image, `False` otherwise. Attributes:
`Protected`.

This is *how validity is tested*, because `Head` cannot do it: malformed input to `Image` stays
unevaluated, so both a valid image and a refused one have head `Image`.

## ImageDimensions

`ImageDimensions[image]` gives `{width, height}`. Attributes: `Protected`.

Transposed relative to `ImageData`, which is `height × width`.

## ImageChannels

`ImageChannels[image]` gives the number of colour channels — `1` for grey, otherwise the length
of each pixel's value list (`3` for RGB, `4` with alpha). Attributes: `Protected`.

A pixel list of a different length in one place is not a channel count, so that declines rather
than reporting the first pixel's length as if it were the image's.

## ImageType

`ImageType[image]` gives `"Bit"`, `"Byte"` or `"Real"`. Attributes: `Protected`.

The type is not decoration: it fixes the **range** of a stored value, which is what makes
`ImageData`'s scaling well defined.

## ImageData

`ImageData[image]` gives the pixels as reals in `[0, 1]`, scaling out the image's type.
`ImageData[image, type]` gives the stored values unscaled. Attributes: `Protected`.

A `"Byte"` 255 comes back as **exactly** `1.0` and a 0 as exactly `0.0`; the tests assert these
as equalities rather than tolerances, so a wrong divisor (254, or 256) fails outright.

```
In[1]:= ImageData[Image[{{0, 128, 255}}]]
Out[1]= {{0.0, 0.501961, 1.0}}

In[2]:= ImageData[Image[{{0, 128, 255}}], "Byte"]
Out[2]= {{0, 128, 255}}
```

Storing the original values and scaling **on read** keeps `Image[…]` round-trippable — what went
in is what `FullForm` shows — and makes the scaling an assertion about one function rather than a
property spread across two. Real data passes through untouched, *including* values outside
`[0, 1]`: storing faithfully beats clamping silently, which would destroy data the caller may
want back.

`ImageData[image, type]` accepts only the image's **own** type. Converting between types has its
own rounding decisions, and an accessor should not make them silently.

## Performance

Accessors are **O(height)**, not O(pixels). The first implementation routed every accessor
through the full validator, which made `ImageDimensions` cost **0.59 ms** on a 512×512 image —
asking how wide an image is touched all 262,144 pixels. Since the canonical form can only be
produced by `Image[]`, which validates every pixel once, the accessors re-check only
**rectangularity** (O(height) row-length comparisons). That is ~1000× faster:

| | before | after |
|---|---|---|
| `ImageDimensions` (512×512) | 0.59 ms | 0.58 µs |
| `ImageQ` (512×512) | 0.62 ms | 0.64 µs |

Rectangularity is still checked, and a test pins it: a hand-typed
`Image[{{1,2},{3}}, "Bit"]` must not pass as an image, because every filter downstream will index
it as rectangular.

`ImageData` on a 512×512 byte image is ~4.8 ms, dominated by building 262,144 `Expr` reals — the
cost of materialising an array as expressions, and the reason filters will operate on the stored
data rather than round-tripping through `ImageData`.
