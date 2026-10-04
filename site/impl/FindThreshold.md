---
references:
  - "N. Otsu, *A threshold selection method from gray-level histograms*, IEEE Trans. Systems, Man, and Cybernetics **9** (1979) 62-66."
source: src/imagefilter.c
---
**Algorithm.** `builtin_findthreshold` returns the single level that best
separates an image into two classes by **Otsu's method**. It reduces the image to
a Rec. 601 luminance plane (`img_grey_plane`) and calls `img_otsu`, which bins
the pixels into a 256-bin histogram (`IMG_OTSU_BINS`), with out-of-range values
clamped **into** the extreme bins rather than dropped so they still vote. It then
sweeps the threshold `t` over the bins, maintaining incrementally the class
weight `w0` and first moment `sum0`, and maximises the **between-class variance**
`sb = w0 · w1 · (mu0 − mu1)²` — algebraically equivalent to minimising the
weighted within-class variance but computable in one incremental pass. A NaN
pixel aborts (it has no bin); a strictly-greater comparison keeps the **lowest**
winning bin, so ties are deterministic. The returned threshold sits at the
**upper edge** of the winning bin, `(best_bin + 0.5)/255`, so a pixel in that bin
falls on the low side. An image whose pixels are all identical has only one
occupied cluster and no split exists, so the routine returns `false` and the head
is left unevaluated.

**Data structures.** A fixed 256-element `double` histogram on the stack; a
single flat luminance `double` plane loaded from the image. No per-threshold
allocation — the sweep is `O(bins)`.

**Complexity / limits.** `O(width · height)` to build the histogram plus
`O(bins)` for the sweep; the threshold is a `double` on the same unit scale as
`ImageData`. Returns unevaluated for a degenerate (single-value) image. The same
`img_otsu` backs `Binarize`'s default threshold and `EdgeDetect`'s high
threshold.
