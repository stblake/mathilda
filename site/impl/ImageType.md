---
source: src/image.c
---
**Algorithm.** `builtin_imagetype` returns the pixel type as a string. It tries
`image3d_info` first (so a volume is accepted) and then `image_info`, reading out
the stored `ImgType` enum and mapping it through `img_type_name` to one of
`"Bit"`, `"Byte"`, `"Bit16"` or `"Real"`. The type is not recomputed from the
data — it is the second argument of the canonical `Image[data, type]` /
`Image3D[data, type]` node, fixed at construction. That construction is where the
inference happens: all-integer data in `{0,1}` becomes `"Bit"`, all-integer in
`0..255` becomes `"Byte"`, anything else `"Real"`. The type's whole purpose is to
fix the **range** of a stored value — `"Bit"` is `{0,1}`, `"Byte"` is `0..255`,
`"Bit16"` is `0..65535`, `"Real"` is already the unit interval — which is what
makes `ImageData`'s scaling to `[0,1]` well defined.

**Data structures.** Reads the type string directly from argument 1 of the image
node; no pixel buffer is touched. `ImageType` is on `pack.c`'s `AWARE` list, so a
packed-buffer image answers without being unpacked.

**Complexity / limits.** `O(1)`. Returns unevaluated (`NULL`) for a non-image.
`"Real32"`/`"Real64"` are accepted as synonyms at construction but normalise to
`"Real"`, so that is what is reported back.
