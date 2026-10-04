---
source: src/imagefilter.c
---
**Algorithm.** `builtin_gaussianmatrix` builds a `(2r+1) × (2r+1)` Gaussian
matrix normalised to sum **exactly 1**. `GaussianMatrix[r]` uses the default
standard deviation `sigma = r/2` (which puts the kernel's edge at two standard
deviations, where the Gaussian has fallen to ~13% of its peak);
`GaussianMatrix[{r, sigma}]` states sigma explicitly. The radius must be a
non-negative integer (a fractional radius has no matrix size); `r = 0` uses
`sigma = 1`. Each entry is `exp(-(dx² + dy²) / (2 sigma²))` with `dx, dy`
measured from the centre, accumulated into a running `sum`, and then every entry
is divided by that **realised** sum. Dividing by the realised sum rather than the
analytic `2 pi sigma²` is deliberate: the analytic constant is correct only for
an infinite kernel, and using it on a truncated one leaves the sum slightly under
1, which darkens an image a little on every pass. The matrix is returned as a
plain nested `List`, so it can be handed to `ImageConvolve`, `MeanFilter`'s path,
or used in arithmetic.

**Data structures.** A `n × n` flat `double` scratch buffer (`n = 2r+1`) filled
and normalised, then marshalled into nested `List` rows of `expr_new_real`
leaves. No image is involved — this is a kernel generator.

**Complexity / limits.** `O(n²)` to fill and normalise, `n = 2r+1`. The radius is
capped at 512. `GaussianFilter[image, r]` is defined as
`ImageConvolve[image, GaussianMatrix[r]]` through the same matrix, and
`EdgeDetect` reuses this builder for its smoothing stage.
