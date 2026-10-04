### Worked examples

```mathematica
(* enlarge a 2x2 image to 4x4 *)
In[1]:= ImageDimensions[ImageResize[Image[{{0., 1.}, {1., 0.}}], {4, 4}]]
```

```mathematica
(* shrinking to 1x1 area-averages: the mean of {0,1,1,0} is 0.5 *)
In[1]:= ImageData[ImageResize[Image[{{0., 1.}, {1., 0.}}], {1, 1}]]
```

```mathematica
(* a single width preserves the aspect ratio: 4-wide input to width 2 gives height 1 *)
In[1]:= ImageDimensions[ImageResize[Image[{{0., 1., 0., 1.}, {1., 0., 1., 0.}}], 2]]
```

### Notes

`ImageResize[image, {w, h}]` resizes to `w × h` pixels; `ImageResize[image, w]`
gives width `w` with the height following to preserve the aspect ratio.
`Resampling -> "Nearest" | "Bilinear" | "Average"` selects the method; the
default `Automatic` uses **area averaging** when either axis shrinks and bilinear
otherwise.

That default is about aliasing: point-sampling a shrinking image destroys every
frequency above half the new sampling rate, and no later interpolation can
restore it. Area averaging is a box prefilter and resample in one pass, exact for
integer reduction factors and using true fractional coverage otherwise.
Coordinates are pixel-centre aligned, avoiding the half-pixel shift the naive
`sx = i·scale` introduces at any scale other than 1:1. The result is a `"Real"`
image, and sizes must be positive integers.
