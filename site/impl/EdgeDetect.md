---
references:
  - "J. Canny, *A computational approach to edge detection*, IEEE TPAMI **8** (1986) 679-698."
  - "N. Otsu, *A threshold selection method from gray-level histograms*, IEEE Trans. Systems, Man, and Cybernetics **9** (1979) 62-66."
source: src/imagefilter.c
---
**Algorithm.** `builtin_edgedetect` implements the **Canny** edge detector,
returning a `"Bit"` image. `EdgeDetect[image]` uses a default Gaussian smoothing
radius of 2; `EdgeDetect[image, r]` sets it (`r = 0` means no smoothing);
`EdgeDetect[image, r, t]` sets the high threshold explicitly. The image is reduced
to a luminance plane (`img_grey_plane`) and run through four stages:

1. **Smooth** — convolve with `GaussianMatrix[r]` (reusing that builder), because
   a derivative amplifies noise. `r = 0` skips it.
2. **Gradient** — the normalised Sobel pair (`deriv_kernel`) gives `dx`, `dy`;
   the magnitude is `sqrt(dx² + dy²)`.
3. **Non-maximum suppression** (`canny_nms`) along the gradient direction, which
   is what makes an edge one pixel wide rather than a thick band.
4. **Hysteresis** — keep any pixel above the high threshold, plus any pixel above
   `0.4 ×` high that is 8-connected (flood-filled via an explicit stack) to a
   strong one, so a real edge survives its faint stretches while isolated weak
   responses do not.

The high threshold defaults to **Otsu on the suppressed magnitude** (`img_otsu`
over `nms`), where the two classes really are edge-against-non-edge; on the raw
magnitude it would be dominated by the ridge flanks. The `0.4` low/high ratio is
the conventional choice, stated because it is a choice.

**Data structures.** Six flat `double`/`size_t` scratch arrays of `width ·
height` (`sm`, `dx`, `dy`, `mag`, `nms`, and the hysteresis stack `stk`), plus
the two Sobel kernels. The result is a 0/1 mask packed into a `"Bit"` image.

**Complexity / limits.** `O(width · height)` for the per-pixel stages plus the
separable Gaussian and Sobel convolutions; hysteresis is a linear flood fill over
the connected weak-edge pixels. Radius is capped at 64; a blank or uniform
gradient field (no two Otsu classes) honestly yields no edges.
