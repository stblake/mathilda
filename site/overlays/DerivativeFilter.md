### Worked examples

```mathematica
In[1]:= ramp = Image[{{0., 0.125, 0.25}, {0., 0.125, 0.25}, {0., 0.125, 0.25}}];  (* a horizontal ramp of slope 1/8 *)
```

```mathematica
In[1]:= Part[ImageData[DerivativeFilter[ramp, {0, 1}]], 2, 2]  (* the Sobel-x response reports the slope exactly: 0.125 *)
```

```mathematica
In[1]:= Part[ImageData[DerivativeFilter[Image[{{0., 0., 0.}, {0., 1., 0.}, {0., 0., 0.}}], {1, 0}]], 1, 2]  (* the Sobel-y derivative of a single bright pixel *)
```

### Notes

`DerivativeFilter[image, {n, m}]` gives the `n`-th derivative down the rows and the `m`-th
across the columns, each order 0 to 2. The kernel is a separable outer product of 1-D stencils:
order 0 is the smoothing `{1,2,1}/4`, order 1 the central difference `{-1,0,1}/2`, order 2 the
second difference `{1,-2,1}`. So `{0,1}` is Sobel-x and `{1,0}` Sobel-y.

The stencils are **normalised**, unlike the raw integer Sobel kernels that report a gradient
eight times the true slope — harmless when only the ranking of edges matters, wrong for
anything that reads the number. On `f(x) = c x` the first derivative gives exactly `c`. The
result is a `"Real"` image.
