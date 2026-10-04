---
source: src/image.c
---
**Algorithm.** `builtin_image3d` is `img_construct(res, "Image3D", vol = true)` —
the volumetric twin of `Image`, sharing one constructor so the two ranks cannot
drift apart. `Image3D[data]` normalises to the canonical `Image3D[data, type]`:
the data is a **depth × height × width** array of voxels (or
depth × height × width × channels for colour), so it is indexed `data[[z, y,
x]]` with slices outermost. Type inference matches `Image`: all-integer voxels in
`{0, 1}` give `"Bit"`, all-integer in `0..255` give `"Byte"`, anything else
`"Real"`; a stated type **coerces** the data to it (round and clip for an integer
type). `Image3D[volume, type]` converts an existing volume between types,
preserving brightness. The data is canonicalised onto the visible NDArray surface
with stored values preserved exactly (`img_typed_storage`), and malformed input
(ragged, non-numeric, complex) is left unevaluated. Validity is tested by
`Image3DQ`, not `ImageQ` (which is `False` for a volume). The one trap the header
flags: `ImageDimensions` reports `{width, height, depth}` — **fully reversed**
from the `z, y, x` storage order, which is Mathematica's convention.

**Data structures.** A canonical `Image3D[data, type]` node whose `data` is a
packed NDArray (a flat depth · height · width · channels buffer, row-major,
channels innermost) with a dtype that follows the pixel type (`"Bit"`/`"Byte"`
hold integers, `"Real"` holds float64). `image3d_info` reads depth, height,
width, channels and type without copying; `image3d_load` materialises a flat
unit-scale buffer for the filters. `Image3D` is on `pack.c`'s `AWARE` list.

**Complexity / limits.** `O(depth · height · width · channels)` to canonicalise
and infer. An already-canonical volume on the NDArray surface is returned
unchanged (`NULL`) so the evaluator reaches a fixed point. `ImageQ`,
`ImageDimensions`, `ImageChannels`, `ImageType` and `ImageData` all accept either
rank; only `ImageQ`/`Image3DQ` are rank-specific.
