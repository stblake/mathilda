---
source: src/imagefilter.c
---
**Algorithm.** `builtin_gaussianfilter` blurs an image with a Gaussian of radius `r`. It is
defined as exactly `ImageConvolve[image, GaussianMatrix[r]]` and **implemented** that way —
`builtin_gaussianmatrix` builds the `(2r+1) × (2r+1)` matrix (`exp(−(dx²+dy²)/(2σ²))`, σ = r/2
by default, normalised by the *realised* sum rather than the analytic `2πσ²` so a truncated
kernel still sums to 1 and does not darken the image), and the same `convolve_dispatch` every
filter uses runs it. Two independent implementations of one identity is how the identity
quietly stops holding, so there is only one. A volume takes `gauss3_kernel` + `convolve3_run`.

The convolution is true convolution (kernel reflected), with replicate ("Fixed") padding — a
constant image convolved with a kernel summing to 1 comes back unchanged everywhere, edges
included, where zero padding would darken them. The Gaussian is symmetric, so reflection is
invisible here; it matters only on asymmetric kernels. The result is always `"Real"` — a
Gaussian of bytes is not a byte.

**Data structures.** A decoded unit buffer and a dense kernel matrix; the result through
`image_build_real` as a packed `"Real"` image (every image head returns a packed buffer —
`make check-image-packing`). `convolve_dispatch` re-derives the kernel's separability, so a
rank-1 Gaussian costs `kw + kh` taps per pixel, not `kw · kh`.

**Complexity / limits.** `O(pixels · (kw + kh))` via the separable path. The radius must be a
non-negative integer (≤ 512 planar, ≤ 32 volumetric).
