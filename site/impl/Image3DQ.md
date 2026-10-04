---
source: src/image.c
---
**Algorithm.** `builtin_image3dq` is a predicate: it calls `image3d_info` on its argument and
returns `True` or `False`. `image3d_info` accepts only the canonical volumetric form
`Image3D[data, "type"]` — the head must be `Image3D`, the second argument one of the canonical
type names, and the data a rectangular `depth × height × width` array (or
`depth × height × width × channels` for colour). `img3_shape_fast` checks rectangularity by
walking only slice and row lengths, and `img_data_storable` rejects complex storage; a packed
buffer answers its shape in `O(1)` from its dims. A malformed `Image3D[...]` returns `NULL`
from its constructor and so stays unevaluated, which is exactly why a separate predicate is
needed: an unevaluated `Image3D["hello"]` and a valid one would otherwise be
indistinguishable.

`Image3DQ` is the volumetric twin of `ImageQ`; `ImageQ` is `False` for a volume and
`Image3DQ` `False` for a plane, so the two partition the image heads cleanly.

**Data structures.** Reads the `Expr` tree without copying it. For a packed `Image3D` the
voxels live in a visible rank-3 (grey) or rank-4 (colour) `NDArray`, slices outermost, indexed
`data[[z, y, x]]`.

**Complexity / limits.** `O(1)` for a packed volume's shape plus the storable scan; a nested
volume pays `O(depth · height)` for the rectangularity check and `O(voxels)` for the storable
check. Returns a boolean, so it is registered packed-aware (the probe never materialises the
buffer it inspects).
