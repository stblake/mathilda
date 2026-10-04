---
source: src/image.c
---
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
