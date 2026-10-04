---
source: src/imagegeom.c
---
**Algorithm.** `builtin_imagereflect` flips an image about an axis by a **pure index
permutation** — no interpolation, so it is exact and self-inverse. `ImageReflect[image]`
reflects top-to-bottom (the default); a side argument selects the axis via `reflect_side_axis`:
`Top`/`Bottom` the height axis, `Left`/`Right` the width axis. Either name of a pair selects
the same axis, since reflecting "to the top" and "to the bottom" are the same operation. The
destination pixel `(x, y)` reads source `(w−1−x or x, h−1−y or y)` depending on the mode.

For an `Image3D` (`reflect3_run`), `Front`/`Back` additionally name the depth axis — the pair
Mathematica uses for volumes. Those two **decline on a plane**, which has no depth axis:
silently reinterpreting them as another axis would turn a caller's mistake into a wrong
picture.

Because it is a permutation, the algebra is exact: reflecting twice about the same axis is the
identity **bit for bit**, and reflections about different axes commute exactly.

**Data structures.** One decoded unit buffer in; a fresh same-size buffer out, wrapped by
`image_build_real` / `image3d_build_real` as a packed `"Real"` image (every image head returns
a packed buffer — `make check-image-packing`).

**Complexity / limits.** `O(pixels)`, one pass. The only inputs refused are `Front`/`Back` on a
plane and a non-side symbol.
