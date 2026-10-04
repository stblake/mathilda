---
references:
  - "C. Harris and M. Stephens, *A Combined Corner and Edge Detector*, Proc. 4th Alvey Vision Conf. (1988) 147-151."
  - "J. Shi and C. Tomasi, *Good Features to Track*, Proc. CVPR (1994) 593-600."
source: src/imagefilter.c
---
**Algorithm.** `builtin_imagecorners` finds corners via the **structure tensor**. A corner is
where the image gradient points in two independent directions, which is a statement about the
second-moment matrix `M = [[⟨Ix Ix⟩, ⟨Ix Iy⟩], [⟨Ix Iy⟩, ⟨Iy Iy⟩]]` of the gradient over a
Gaussian-weighted window. `structure_tensor` builds the three smoothed gradient products
(reusing the normalised Sobel stencils and a separable Gaussian through `convolve_dispatch`);
`corner_response` reads them as `MinimumEigenvalue` (Shi-Tomasi) `λ_min = ½(tr − √((Sxx−Syy)² +
4 Sxy²))`, written as a sum of squares under the root so cancellation cannot produce a NaN.
Along a straight edge `M` has rank 1, so `λ_min = 0` exactly — the property that separates a
corner detector from an edge detector.

`corner_peaks_list` then applies, in order: a threshold at `frac` of the largest response
(default 0.05); 3×3 non-maximum suppression; a minimum-separation filter, greedy in
**descending response** so a cluster's survivor is its strongest member; and finally the
feature cap (`MaxFeatures` option, or the 5th positional argument). The list is sorted
strongest first, ties broken by position for determinism.

**Data structures.** A grey plane, three tensor buffers and a response buffer; the result is a
`List` of `{row, column}` positions (1-based, indexing `ImageData` directly) — **not** an
image, so the `make check-image-packing` audit does not apply. Note this is *not* Mathematica's
`{x, y}` from the bottom-left.

**Complexity / limits.** `O(pixels)` for the response (separable convolutions) plus the greedy
separation pass over detected peaks. Window radius ≤ 32.
