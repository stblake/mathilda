### Worked examples

```mathematica
In[1]:= ImageType[LocalAdaptiveBinarize[Image[{{0., 0., 1.}, {0., 1., 1.}, {1., 1., 1.}}], 1]]  (* the result is binary, typed Bit *)
```

```mathematica
In[1]:= ImageData[LocalAdaptiveBinarize[Image[{{0., 0., 1.}, {0., 1., 1.}, {1., 1., 1.}}], 1]]  (* each pixel vs its local mean *)
```

```mathematica
In[1]:= ImageData[LocalAdaptiveBinarize[Image[{{0.2, 0.9}, {0.3, 0.8}}], 1, {1, -0.2, 0.}]]  (* a negative c2 is Sauvola's rule *)
```

### Notes

`LocalAdaptiveBinarize[image, r]` binarizes by comparing each pixel to the **mean** of its own
`(2r+1) × (2r+1)` neighbourhood; `LocalAdaptiveBinarize[image, r, {c1, c2, c3}]` compares to
`c1·mean + c2·stddev + c3`. A global threshold cannot binarize unevenly lit content: if one
half of a page is darker than the other, no single number separates ink from paper in both
halves at once.

Mean alone (the default `{1, 0, 0}`) is Bradley's method; a negative `c2` is Sauvola's,
tightening the threshold where the neighbourhood is busy. Summed-area tables make the window
statistics `O(1)` per pixel regardless of `r`. The result is typed `"Bit"`, and colour is
reduced to luminance first.
