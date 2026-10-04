---
source: src/image.c
---
**Algorithm.** `builtin_imagedimensions` reports an image's size. It tries `image3d_info`
first: a volume answers `{width, height, depth}` — three elements, **fully reversed** from the
`depth × height × width` storage order. A plane answers `{width, height}` via `image_info` —
**transposed** relative to `ImageData`, which returns a `height × width` array. Both are
Mathematica's conventions and the subsystem's most common trap, which is why every size test
uses a non-square image (a square one cannot tell the two axes apart).

The size comes from `img_shape_fast`, not the full validator, and that is a measured
decision: routing the query through the per-pixel validator made `ImageDimensions` cost
0.59 ms on a 512 × 512 image because asking how wide an image is touched all 262144 pixels. A
filter pipeline queries dimensions constantly, so `img_shape_fast` instead reads a packed
buffer's dims in `O(1)` and, for nested data, checks only that rows are the same length
(`O(height)`) — ~500× cheaper, with the per-pixel numeric check left where it is paid once, in
construction.

**Data structures.** Reads the stored `NDArray` dims (or the nested-`List` row lengths); the
result is a 2- or 3-element `List` of integers.

**Complexity / limits.** `O(1)` for a packed image, `O(height)` (or `O(depth · height)`) for a
nested one. Registered packed-aware, so a packed image's dims are read directly rather than the
buffer being materialised first.
