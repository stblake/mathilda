### Worked examples

```mathematica
(* a rank-3 pixel array reports its last dimension -- RGB is 3 *)
In[1]:= ImageChannels[Image[{{{1., 0, 0}, {0, 1., 0}}}]]
```

```mathematica
(* a rank-2 (grey) array is one channel *)
In[1]:= ImageChannels[Image[{{0., 1.}}]]
```

```mathematica
(* an alpha channel brings the count to 4 *)
In[1]:= ImageChannels[Image[{{{1., 0, 0, 0.5}}}]]
```

### Notes

`ImageChannels[image]` gives the number of colour channels: `1` for a grey image,
otherwise the length of each pixel's value list — `3` for RGB, `4` with an alpha
channel, `2` for grey-plus-alpha.

The channel count is the trailing (innermost, interleaved) dimension of the pixel
array, so it is read straight off the image's shape without loading any pixels.
It accepts a volumetric `Image3D` as well as a plane, applying the same rule to
the last axis of the depth × height × width (× channels) buffer.
