---
source: src/imagecompose.c
---
**Algorithm.** `builtin_alphachannel` returns an image's opacity as a one-channel image. It
decodes the input to a unit-interval buffer with `image_load`, then reads the opacity per
pixel: a two-channel image is grey+alpha (channel 1), a four-channel image is RGB+alpha
(channel 3), and any other channel count has no alpha, so it answers `1.0` — fully opaque.
That last case is a decision, not a decline: "how transparent is this?" has an answer for
every image, and for one without an alpha channel the answer is "not at all". The opacity
buffer is wrapped back up with `image_build_real`.

**Data structures.** One decoded `Img` (`{w, h, c, buf}` of unit-interval reals, row-major,
channels innermost) in, a fresh `width · height` buffer out. The result is a `"Real"`
single-channel image built as a packed `NDArray` by `image_build_real`, which is what keeps
every image head's output a packed buffer (enforced by `make check-image-packing`).

**Complexity / limits.** `O(width · height)` — one pass, picking one sample per pixel. The
grey/alpha and RGB/alpha channel conventions are the only ones recognised; a 3-channel image
is treated as alpha-free and reads back all-opaque.
