### Worked examples

```mathematica
(* the corner response has the same dimensions as the input, one channel *)
In[1]:= ImageDimensions[CornerFilter[Image[{{0., 0, 0}, {0, 1., 0}, {0, 0, 0}}]]]
```

```mathematica
(* the result is a single-channel "Real" strength map *)
In[1]:= ImageChannels[CornerFilter[Image[{{0., 0, 0}, {0, 1., 0}, {0, 0, 0}}]]]
```

```mathematica
(* the Harris measure, selected positionally *)
In[1]:= ImageType[CornerFilter[Image[{{0., 0, 0}, {0, 1., 0}, {0, 0, 0}}], 1, "Harris"]]
```

### Notes

`CornerFilter[image]` gives the corner strength at every pixel, from the
eigenvalues of the Gaussian-weighted second-moment matrix of the gradient (the
structure tensor): both eigenvalues small is flat, one large is an edge, both
large is a corner.

`CornerFilter[image, r]` sets the window radius (default `2`), and a third
argument or a `Method ->` option selects `"MinimumEigenvalue"` (the default,
Shi–Tomasi's `lambda_min`, comparable across images) or `"Harris"`
(`det − 0.04 trace²`, cheaper and negative on edges). A straight edge scores
**zero** under both measures — every gradient in the window is parallel, so the
tensor has rank 1 and its determinant and smaller eigenvalue vanish. Colour is
reduced to luminance first. `ImageCorners` turns this response into discrete
corner positions.
