# GaussianMatrix

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`GaussianMatrix[r] gives a (2r+1) x (2r+1) Gaussian matrix normalised to sum 1. GaussianMatrix[{r, sigma}] states the standard deviation; it defaults to r/2, which puts the kernel's edge at two standard deviations. Normalisation divides by the realised sum rather than the analytic 2 pi sigma^2, because the analytic constant is correct only for an infinite kernel and using it on a truncated one leaves the sum under 1 -- which darkens an image slightly on every pass.`**

## Examples (2)

Every input below was run against the current Mathilda build and its output recorded.

### Basic Examples (2)

```mathematica
In[1]:= GaussianMatrix[1]
Out[1]= {{0.0113437, 0.0838195, 0.0113437}, {0.0838195, 0.619347, 0.0838195}, {0.0113437, 0.0838195, 0.0113437}}

In[2]:= Total[Flatten[GaussianMatrix[2]]]
Out[2]= 1.0
```

## Implementation notes

**Attributes:** `Protected`.

## References

- Source: [`src/imagefilter.c`](https://github.com/stblake/mathilda/blob/main/src/imagefilter.c)
- Specification: [`docs/spec/builtins/image-processing.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/image-processing.md)
- Tests: [`tests/test_image.c`](https://github.com/stblake/mathilda/blob/main/tests/test_image.c)
