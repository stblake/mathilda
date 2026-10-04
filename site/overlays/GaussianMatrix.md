### Worked examples

```mathematica
(* a radius-1 kernel is 3x3, peaked at the centre, summing to 1 *)
In[1]:= GaussianMatrix[1]
```

```mathematica
(* normalised to sum exactly 1 -- Total over both levels confirms it *)
In[1]:= Total[GaussianMatrix[2], 2]
```

```mathematica
(* radius r gives a (2r+1)x(2r+1) matrix *)
In[1]:= Dimensions[GaussianMatrix[2]]
```

```mathematica
(* the standard deviation can be stated instead of defaulting to r/2 *)
In[1]:= GaussianMatrix[{2, 1}]
```

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
