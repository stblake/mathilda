# GaussianMatrix

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`GaussianMatrix[r] gives a (2r+1) x (2r+1) Gaussian matrix normalised to sum 1. GaussianMatrix[{r, sigma}] states the standard deviation; it defaults to r/2, which puts the kernel's edge at two standard deviations. Normalisation divides by the realised sum rather than the analytic 2 pi sigma^2, because the analytic constant is correct only for an infinite kernel and using it on a truncated one leaves the sum under 1 -- which darkens an image slightly on every pass.`**

## Examples (6)

Every input below was run against the current Mathilda build and its output recorded.

### Basic Examples (2)

```mathematica
In[1]:= GaussianMatrix[1]
Out[1]= {{0.0113437, 0.0838195, 0.0113437}, {0.0838195, 0.619347, 0.0838195}, {0.0113437, 0.0838195, 0.0113437}}

In[2]:= Total[Flatten[GaussianMatrix[2]]]
Out[2]= 1.0
```

### Applications (4)

```mathematica
In[3]:= GaussianMatrix[1]
Out[3]= {{0.0113437, 0.0838195, 0.0113437}, {0.0838195, 0.619347, 0.0838195}, {0.0113437, 0.0838195, 0.0113437}}

In[4]:= Total[GaussianMatrix[2], 2]
Out[4]= 1.0

In[5]:= Dimensions[GaussianMatrix[2]]
Out[5]= {5, 5}

In[6]:= GaussianMatrix[{2, 1}]
Out[6]= {{0.00296902, 0.0133062, 0.0219382, 0.0133062, 0.00296902}, {0.0133062, 0.0596343, 0.0983203, 0.0596343, 0.0133062}, {0.0219382, 0.0983203, 0.162103, 0.0983203, 0.0219382}, {0.0133062, 0.0596343, 0.0983203, 0.0596343, 0.0133062}, {0.00296902, 0.0133062, 0.0219382, 0.0133062, 0.00296902}}
```

## Implementation notes

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

**Attributes:** `Protected`.

## References

- Source: [`src/imagefilter.c`](https://github.com/stblake/mathilda/blob/main/src/imagefilter.c)
- Specification: [`docs/spec/builtins/image-processing.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/image-processing.md)
- Tests: [`tests/test_image.c`](https://github.com/stblake/mathilda/blob/main/tests/test_image.c)

## Notes & additional examples

### Notes

`GaussianMatrix[r]` gives a `(2r+1) × (2r+1)` Gaussian matrix normalised to sum
`1`. `GaussianMatrix[{r, sigma}]` states the standard deviation; it defaults to
`r/2`, which puts the kernel's edge at two standard deviations.

Normalisation divides by the **realised** sum of the truncated kernel rather than
the analytic `2 pi sigma²`. The analytic constant is correct only for an infinite
kernel; using it on a truncated one leaves the sum under 1, which darkens an
image slightly on every pass — invisible once and obvious after fifty.
`GaussianFilter[image, r]` is exactly `ImageConvolve[image, GaussianMatrix[r]]`,
and `EdgeDetect`'s smoothing stage reuses this same builder.
