# ImageData

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`ImageData[image] gives the pixel array as reals in [0, 1], scaling out the image's type -- a "Byte" 255 comes back as exactly 1.0. The array is height x width, or height x width x channels for a colour image, interleaved. ImageData[image, type] gives the stored values unscaled instead, where type must be the image's own type; converting between types is a separate operation with its own rounding, not something this does silently.`**

## Examples (5)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (2)

```mathematica
In[1]:= ImageData[Image[{{0, 128, 255}}]]
Out[1]= {{0.0, 0.501961, 1.0}}

In[2]:= ImageData[Image[{{0, 128, 255}}], "Byte"]
Out[2]= {{0, 128, 255}}
```

### Applications (3)

```mathematica
In[3]:= ImageData[Image[{{0, 255}, {255, 0}}, "Byte"]]
Out[3]= {{0.0, 1.0}, {1.0, 0.0}}

In[4]:= ImageData[Image[{{0, 255}}, "Byte"], "Byte"]
Out[4]= {{0, 255}}

In[5]:= ImageData[Image[{{0.25, 0.5}, {0.75, 1.}}]]
Out[5]= {{0.25, 0.5}, {0.75, 1.0}}
```

## Implementation notes

**Algorithm.** `builtin_imagedata` returns the pixel array as a nested list.
`ImageData[image]` gives **unit-interval** reals, scaling out the image's type:
the type's maximum (`255` for `"Byte"`, `65535` for `"Bit16"`, `1` for `"Bit"`
and `"Real"`) is divided out by `img_to_unit`, so a `"Byte"` 255 comes back as
exactly `1.0`. `ImageData[image, type]` instead returns the **stored values
unscaled** — but `type` must equal the image's own type, since converting between
types is a separate operation with its own rounding and this routine refuses to
do it silently (a mismatched `type` returns unevaluated). The nested result is
rebuilt by `img_scale_tree` / `nd_nest`, one array axis per recursive call with
stride the product of the remaining dimensions, so the same code serves a grey
plane (height × width), a colour plane (height × width × channels, interleaved)
and a colour volume (rank 4) without unrolling each rank by hand. The scaling
walk is shape-agnostic, so a volume needs only its type.

**Data structures.** Operates on the image node's argument 0, normally a packed
NDArray buffer (`is_ndarray`); leaves are read with `ndt_get` and rebuilt into
`List` trees. The unscaled `"Real"` float64 case needs no per-element conversion
at all — `img_to_unit` is the identity there. `ImageData` is on `pack.c`'s
`AWARE` list.

**Complexity / limits.** `O(height · width · channels)` (× depth for a volume) —
one scaled leaf per pixel-channel, plus the `List` nodes of the rebuilt tree.
`ImageData[image, type]` with `type` not equal to the image's own type is
declined rather than performing a hidden conversion.

**Attributes:** `Protected`.

## References

**See also:** [FullForm](../../expression-information/FullForm/)

- Source: [`src/image.c`](https://github.com/stblake/mathilda/blob/main/src/image.c)
- Specification: [`docs/spec/builtins/image-processing.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/image-processing.md)
- Tests: [`tests/test_image.c`](https://github.com/stblake/mathilda/blob/main/tests/test_image.c)

## Notes & additional examples

### Notes

`ImageData[image]` gives the pixel array as reals in `[0, 1]`, scaling out the
image's type — the array is height × width, or height × width × channels
(interleaved) for a colour image, the same shape the image was built from.

`ImageData[image, type]` gives the stored values **unscaled** instead, but only
when `type` is the image's own type: converting between types is a separate
operation with its own rounding, and this never does it silently — a mismatched
type returns unevaluated. The result shape is exactly the one `Image[...]` was
handed, which is what makes `Image[ImageData[img] ..., type]` a faithful
round-trip. It accepts a volume as well as a plane.
