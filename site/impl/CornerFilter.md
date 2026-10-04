---
references:
  - "C. Harris and M. Stephens, *A combined corner and edge detector*, Proc. 4th Alvey Vision Conf. (1988) 147-151."
  - "J. Shi and C. Tomasi, *Good features to track*, Proc. CVPR (1994) 593-600."
source: src/imagefilter.c
---
**Algorithm.** `builtin_cornerfilter` gives the corner strength at every pixel
from the eigenvalues of the Gaussian-weighted **structure tensor** (the
second-moment matrix of the gradient). `CornerFilter[image]` uses window radius
`2`; `CornerFilter[image, r]` sets it; `CornerFilter[image, r, method]` or a
`Method ->` option selects `"MinimumEigenvalue"` (the default, Shi–Tomasi's
`lambda_min`) or `"Harris"` (`det − 0.04 · trace²`). Options are stripped first
(`options_extract`), and an explicit positional method beats the registered
default, so `CornerFilter[img, 2, "Harris"]` is honoured rather than silently
overridden. Colour is reduced to luminance first, since a corner is a property of
brightness. The response (`corner_response` for a plane, `corner3_response` for a
volume) computes per-pixel gradients, forms the products `gx², gy², gx·gy`, and
Gaussian-smooths each over the window to build the tensor; then per pixel it
either takes the Harris combination or the smaller eigenvalue, clamping a
rounding-negative eigenvalue to `0` (the tensor is positive semidefinite). Both
measures score a **straight edge as zero**: every gradient in the window is
parallel, so the tensor has rank 1 and its determinant and smaller eigenvalue
vanish; both large eigenvalues mean a corner.

**Data structures.** Flat `double` buffers for the luminance plane, the two (or
three) gradient components, the tensor entries (`sxx`, `syy`, `sxy`, and for a
volume `szz`, `sxz`, `syz`), a separable Gaussian line, and the response. The
result is a single-channel `"Real"` image (`image_build_real`).

**Complexity / limits.** `O(width · height)` for the gradients and the
pointwise eigenvalue step, plus the separable Gaussian smoothing (`O(r)` per
pixel per axis). Radius `1..32`. `ImageCorners` builds on this response with
thresholding, non-maximum suppression and minimum-separation to return discrete
positions.
