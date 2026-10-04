### Worked examples

```mathematica
In[1]:= ImageData[DistanceTransform[Image[{{0, 0, 0}, {0, 1, 0}, {0, 0, 0}}]]]  (* the centre is distance 1 from the border *)
```

```mathematica
In[1]:= Part[ImageData[DistanceTransform[Image[{{0, 0, 0, 0, 0}, {0, 1, 1, 1, 0}, {0, 1, 1, 1, 0}, {0, 1, 1, 1, 0}, {0, 0, 0, 0, 0}}]]], 3, 3]  (* the blob centre is distance 2 *)
```

```mathematica
In[1]:= Part[ImageData[DistanceTransform[Image[{{0, 1, 1, 1}, {1, 1, 1, 1}, {1, 1, 1, 1}, {1, 1, 1, 1}, {1, 1, 1, 1}}]]], 5, 4]  (* three across and four down reads exactly 5 *)
```

### Notes

`DistanceTransform[image]` replaces each pixel by its **exact** Euclidean distance to the
nearest background pixel; background pixels are 0, so the value rises toward the interior of a
blob. `DistanceTransform[image, t]` takes pixels above `t` as foreground.

Exact rather than the classic two-pass chamfer approximation, which cannot represent `√2` with
integer steps and gets diagonal distances a few percent wrong — invisible on a picture, fatal
to a test. It uses Felzenszwalb and Huttenlocher's lower-envelope-of-parabolas method,
`O(n)` per row; separability is exact because squared Euclidean distance is a sum over the
axes, and the square root is taken once at the end.
