### Worked examples

```mathematica
In[1]:= BoxMatrix[1]  (* a 3x3 block of ones *)
```

```mathematica
In[1]:= BoxMatrix[0]  (* radius 0 is the 1x1 matrix {{1}} *)
```

```mathematica
In[1]:= ImageData[Dilation[Image[{{0, 0, 0}, {0, 1, 0}, {0, 0, 0}}], BoxMatrix[1]]]  (* used as a morphology element *)
```

### Notes

`BoxMatrix[r]` is a `(2r+1) × (2r+1)` matrix of `1`s. It is **not** normalised, matching
Mathematica — so `ImageConvolve[image, BoxMatrix[1]]` is nine times too bright; the normalised
version is a mean filter, which is what `MeanFilter` provides.

The result is an ordinary integer matrix, not an image. Its common uses are as a convolution
kernel for `ImageConvolve` and as the structuring element for the morphology operators
(`Dilation`, `Erosion`, `Opening`, `Closing`), where only its nonzero support matters.
