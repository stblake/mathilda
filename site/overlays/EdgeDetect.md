### Worked examples

```mathematica
(* the result of the Canny detector is a "Bit" image *)
In[1]:= ImageType[EdgeDetect[Image[{{0., 0, 1, 1}, {0, 0, 1, 1}}]]]
```

```mathematica
(* a step between a dark and a bright block, detected with no pre-smoothing *)
In[1]:= EdgeDetect[Image[{{0., 0, 1, 1}, {0, 0, 1, 1}, {0, 0, 1, 1}}], 0]
```

```mathematica
(* the output keeps the input's dimensions *)
In[1]:= ImageDimensions[EdgeDetect[Image[{{0., 0, 1, 1}, {0, 0, 1, 1}}]]]
```

### Notes

`EdgeDetect[image]` finds edges by the **Canny** algorithm, giving a `"Bit"`
image. `EdgeDetect[image, r]` sets the Gaussian smoothing radius (default `2`;
`0` means no smoothing), and `EdgeDetect[image, r, t]` sets the high threshold
explicitly.

Four stages: smooth (a derivative amplifies noise); gradient by the normalised
Sobel pair; non-maximum suppression along the gradient direction, which makes an
edge one pixel wide rather than a thick band; and hysteresis, keeping any pixel
above the high threshold plus any above `0.4 ×` it that is 8-connected to one.
The high threshold defaults to Otsu's method applied to the **suppressed**
magnitude, where the two classes really are edge against non-edge; run on the raw
magnitude it would be dominated by the ridge flanks.
