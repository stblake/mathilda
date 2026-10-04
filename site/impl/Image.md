---
source: src/image.c
---
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
