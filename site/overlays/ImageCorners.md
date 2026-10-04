### Worked examples

```mathematica
In[1]:= ImageCorners[Image[{{0., 0., 0., 0.}, {0., 1., 1., 0.}, {0., 1., 1., 0.}, {0., 0., 0., 0.}}]]  (* a bright square: positions are {row, column} *)
```

```mathematica
In[1]:= ImageCorners[Image[{{0., 0., 0., 0.}, {0., 1., 1., 0.}, {0., 1., 1., 0.}, {0., 0., 0., 0.}}], 1, 0.05, 0, 1]  (* radius, threshold, separation, max features *)
```

```mathematica
In[1]:= SeedRandom[42]; Length[ImageCorners[RandomImage[1, {32, 32}], 2, 0.2, 3]]  (* separation keeps clusters from flooding the list *)
```

### Notes

`ImageCorners[image]` gives corner positions as `{row, column}`, 1-based, so each indexes
`ImageData` directly — this is **not** Mathematica's `{x, y}` from the bottom-left.
`ImageCorners[image, r, t, d, n]` sets the window radius (default 2), the threshold as a
fraction of the largest response (0.05), the minimum separation in pixels (0), and the maximum
number of features (all, also settable with the `MaxFeatures` option).

A corner is where the structure tensor has two strong eigenvalues; along a straight edge it has
rank 1 and the response is 0 exactly. Three filters apply in order — threshold, 3×3
non-maximum suppression, then minimum separation greedy in descending response — because each
removes what the others cannot; the feature cap is applied last, and the list is sorted
strongest first.
