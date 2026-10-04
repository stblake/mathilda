---
source: src/imagefilter.c
---
**Algorithm.** `builtin_derivativefilter` applies a separable derivative kernel.
`DerivativeFilter[image, {n, m}]` is the `n`-th derivative down the rows and the `m`-th across
the columns, each order 0–2. `deriv_kernel` builds the full 2-D kernel as the outer product of
two 1-D stencils from `deriv_stencil`:

- order 0: `{1, 2, 1}/4` — smoothing, normalised to sum 1 (so it preserves a constant, and
  being symmetric, a linear ramp exactly);
- order 1: `{+1/2, 0, −1/2}` — central difference. On `f(x) = c x` it gives exactly `c`;
- order 2: `{1, −2, 1}` — second difference. On `f(x) = c x²` it gives exactly `2c`.

So `{0, 1}` is Sobel-x and `{1, 0}` Sobel-y. Two subtleties the code pins: the order-1 stencil
is written `{+1/2, 0, −1/2}`, **pre-flipped**, because `ImageConvolve` reflects its kernel and
a reading-order stencil would compute the *negated* derivative — caught only by asserting an
exact signed value, since a gradient magnitude squares the sign away. And the stencils are
**normalised** (/4, /2), unlike the raw integer Sobel kernels that report eight times the true
slope. The full matrix is handed to the same `convolve_dispatch` every filter uses, which
re-derives the separable factorisation rather than trusting it. A volume takes `deriv3_run`
with a third stencil.

**Data structures.** One decoded unit buffer; the result through `image_build_real` as a packed
`"Real"` image (every image head returns a packed buffer — `make check-image-packing`).

**Complexity / limits.** The kernel factors to rank 1, so it costs `O(pixels · (kw + kh))`.
Orders are restricted to 0–2; the padding is replicate, as for all convolutions.
