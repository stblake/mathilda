---
source: src/image.c
---
**Algorithm.** `builtin_imageq` is a one-argument predicate that returns
`True`/`False` and nothing symbolic: it calls `image_info(arg, NULL, NULL, NULL,
NULL)` and wraps the boolean. `image_info` is the single validity gate the whole
subsystem shares — it checks that the argument is a two-argument `Image[data,
type]` whose second argument is a type string (`"Bit"`, `"Byte"`, `"Bit16"`,
`"Real"`), that `img_shape_fast` finds a rectangular height × width (× channels)
array, and that `img_data_storable` confirms every leaf is in the range the type
fixes. `ImageQ` is the companion to the fact that malformed input to `Image[...]`
is left **unevaluated** rather than erroring — the constructor returns `NULL` for
ragged, non-numeric or complex data, so the head stays `Image[...]` and `ImageQ`
is how a caller tests whether that happened. A volume is deliberately `False`
here (it is `Image3DQ` that accepts one), since the two ranks are distinct
objects.

**Data structures.** Reads the canonical `Image` node — an `EXPR_FUNCTION` with
head `Image`, argument 0 the pixel array (normally a packed NDArray buffer),
argument 1 the type string. No buffer is loaded or copied; validation walks only
the shape and leaf scalars. `ImageQ` is on `pack.c`'s `AWARE` list, so a packed
pixel buffer is inspected in place rather than being unpacked into boxed `Expr`
nodes.

**Complexity / limits.** `O(1)` for the structural checks plus one pass over the
leaves for `img_data_storable`; no allocation. Returns `False`, never
unevaluated, for every non-image — it is a total predicate.
