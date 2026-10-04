### Worked examples

```mathematica
In[1]:= GaussianFilter[Image[{{0., 0., 0.}, {0., 1., 0.}, {0., 0., 0.}}], 1]  (* a point spreads into the kernel *)
```

```mathematica
In[1]:= Part[ImageData[GaussianFilter[Image[{{0., 0., 0.}, {0., 1., 0.}, {0., 0., 0.}}], 1]], 2, 2]  (* the centre keeps most of its weight *)
```

```mathematica
In[1]:= GaussianFilter[Image[{{0.3, 0.9}, {0.9, 0.3}}], 1]  (* a constant-sum kernel preserves overall brightness *)
```

### Notes

`GaussianFilter[image, r]` blurs with a Gaussian of radius `r`. It is exactly
`ImageConvolve[image, GaussianMatrix[r]]` — the same matrix through the same convolution, not a
second implementation — and a test asserts the identity.

The kernel is normalised by its realised sum (not the analytic `2 π σ²`), so a truncated kernel
still sums to 1 and does not darken the image on each pass. Padding replicates the border, so a
constant image comes back unchanged everywhere including the edges. The result is a `"Real"`
image, since a Gaussian of bytes is not a byte.
