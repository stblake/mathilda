---
source: src/imagefilter.c
---
**Algorithm.** `builtin_gradientfilter` returns the gradient **magnitude**
`Sqrt[dx² + dy²]`. Colour is reduced to luminance **first** (`img_grey_plane`), then
differentiated once — not differentiated per channel and combined, which would need an
arbitrary rule (max? sum? norm?). The two first derivatives come from the same normalised Sobel
stencils `DerivativeFilter` uses (`deriv_kernel(0,1)` and `deriv_kernel(1,0)`, each through
`convolve_dispatch`), and are combined pixelwise as the Euclidean norm.

The magnitude rather than `|dx| + |dy|` is chosen for **rotation invariance**: the length of a
vector does not depend on which way the axes point, so an edge at 45° reports the same strength
as one at 0°, where the absolute-sum alternative would report it `√2` times stronger and bias
every downstream threshold by orientation. A volume (`gradient3_run`) takes three separable
derivatives and combines them with one square root at the end — `sqrt` is not additive, so
taking it per component would give the sum of absolute derivatives, a different and larger
quantity.

**Data structures.** A grey buffer plus two derivative buffers `dx`, `dy`; the result through
`image_build_real` as a packed single-channel `"Real"` image (every image head returns a packed
buffer — `make check-image-packing`).

**Complexity / limits.** `O(pixels)` — two separable 3-tap convolutions plus a per-pixel
`hypot`. The border is replicate, as for all convolutions; the output is single-channel.
