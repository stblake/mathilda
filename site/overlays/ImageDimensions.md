### Worked examples

```mathematica
In[1]:= ImageDimensions[Image[{{0., 1., 0.}}]]  (* one row of three pixels: width 3, height 1 *)
```

```mathematica
In[1]:= ImageDimensions[Image[{{0., 1., 0.}, {1., 0., 1.}}]]  (* height x width 2x3 stores, reports {3, 2} *)
```

```mathematica
In[1]:= ImageDimensions[Image3D[{{{0., 1.}, {1., 0.}}, {{1., 0.}, {0., 1.}}}]]  (* a volume reports {width, height, depth} *)
```

### Notes

`ImageDimensions` gives `{width, height}`, which is **transposed** from `ImageData`'s
`height × width` array — the convention Mathematica uses and the single most common source of
silently-wrong image code. The example above pins it with a non-square image: a 2 × 3 pixel
array (two rows of three) reports `{3, 2}`.

For an `Image3D` the result is `{width, height, depth}`, fully reversed from the
`depth × height × width` storage order.

The size is read without walking the pixels, so querying it repeatedly in a filter pipeline is
cheap — it was deliberately taken off the full-validator path that once made it `O(pixels)`.
