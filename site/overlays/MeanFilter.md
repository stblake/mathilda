### Worked examples

```mathematica
(* a lone bright pixel is spread over the 3x3 window: every value becomes 1/9 *)
In[1]:= ImageData[MeanFilter[Image[{{0., 0, 0}, {0, 1., 0}, {0, 0, 0}}], 1]]
```

```mathematica
(* a constant image is unchanged by averaging *)
In[1]:= ImageData[MeanFilter[Image[{{0.5, 0.5}, {0.5, 0.5}}], 1]]
```

```mathematica
(* the result is a "Real" image regardless of the input type *)
In[1]:= ImageType[MeanFilter[Image[{{0, 1}, {1, 0}}, "Bit"], 1]]
```

### Notes

`MeanFilter[image, r]` averages over a `(2r+1) × (2r+1)` neighbourhood. This *is*
a convolution with a normalised box, and it is implemented as one rather than as
a separate averaging loop — two implementations of one identity is how the
identity quietly stops holding.

Being a full rectangle, the box kernel is separable, so the cost is `kw + kh`
taps per axis rather than `kw · kh`, independent of the window area. Border
pixels use the same `"Fixed"` (replicate) padding the convolutions use, and the
result is a `"Real"` image. Contrast with `BoxMatrix[r]`, whose entries are `1`
and unnormalised — convolving with it is `(2r+1)²` times too bright.
