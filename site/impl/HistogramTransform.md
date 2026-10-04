---
references:
  - "R. C. Gonzalez and R. E. Woods, *Digital Image Processing*, 3rd ed. (Pearson, 2008), §3.3 (histogram equalization)."
source: src/imagecolor.c
---
**Algorithm.** `builtin_histogramtransform` equalises the histogram, spreading
the brightness distribution toward uniform over 256 bins (`HT_BINS`). It builds
the histogram of the **luminance** (Rec. 601 `luma` for a colour image, the
single channel otherwise), accumulates it into a cumulative distribution `cdf`
normalised to `[0, 1]`, and remaps each pixel's luminance `l` to `nl = cdf[bin]`.
The new luminance is applied to every channel **as a ratio** `f = nl / l`, so hue
survives — equalising each channel independently is the other obvious choice and
it shifts colour, because it removes exactly the imbalance that makes an image
warm or cool. A black pixel (`l <= 1e-9`) has no ratio to scale and takes the new
luminance `nl` in every channel, which is grey — the only hue-free answer
available — and each channel is clamped to `[0, 1]`. An alpha channel passes
through unchanged.

**Data structures.** A fixed 256-element `size_t` histogram and `double` CDF on
the stack; one flat unit-scale pixel buffer (`image_load`) and one output buffer
of the same shape. Two passes over the pixels — one to build the histogram, one
to remap. The result is a `"Real"` image.

**Complexity / limits.** `O(n + bins)`, `n = width · height · channels`, bins
fixed at 256. Operates on a plane; the mapping is derived once from the global
luminance distribution and applied to every pixel.
