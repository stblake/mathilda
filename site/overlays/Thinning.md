### Worked examples

```mathematica
(* a solid 3x3 block thins to its single central pixel *)
In[1]:= ImageData[Thinning[Image[{{0, 0, 0, 0, 0}, {0, 1, 1, 1, 0}, {0, 1, 1, 1, 0}, {0, 1, 1, 1, 0}, {0, 0, 0, 0, 0}}]]]
```

```mathematica
(* the result is a "Bit" image *)
In[1]:= ImageType[Thinning[Image[{{1, 1, 1}, {1, 1, 1}, {1, 1, 1}}]]]
```

```mathematica
(* Thinning[image, n] stops after n iterations *)
In[1]:= ImageData[Thinning[Image[{{0, 0, 0, 0, 0}, {0, 1, 1, 1, 0}, {0, 1, 1, 1, 0}, {0, 1, 1, 1, 0}, {0, 0, 0, 0, 0}}], 1]]
```

### Notes

`Thinning[image]` reduces the foreground to a one-pixel-wide skeleton by
**Zhang–Suen** thinning, iterating until a pass deletes nothing;
`Thinning[image, n]` stops after `n` iterations.

The two subiterations are what preserve connectivity: deleting every
individually-removable pixel in one pass would sever a diagonal line, because two
diagonal neighbours can each be removable while removing both disconnects the
shape. Marked pixels are deleted together after each pass, never in place. A
non-binary image is thresholded at `0.5` — apply `Binarize` first for any other
rule. The result is a `"Bit"` image and is always a subset of the input. Follow
it with `Pruning` to remove the short spurs a skeleton grows at boundary
irregularities.
