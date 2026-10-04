### Worked examples

```mathematica
In[1]:= ImageData[ColorConvert[Image[{{{1., 0., 0.}}}], "Grayscale"]]  (* pure red -> its Rec. 601 luminance 0.299 *)
```

```mathematica
In[1]:= ImageData[ColorConvert[Image[{{{0., 1., 0.}}}], "Grayscale"]]  (* green carries most of the perceived brightness *)
```

```mathematica
In[1]:= ImageChannels[ColorConvert[Image[{{{1., 1., 0.}, {0., 0., 1.}}}], "Gray"]]  (* the result is one channel *)
```

### Notes

`ColorConvert[image, "Grayscale"]` (or `"Gray"`) reduces an image or an `Image3D` to one
channel using the Rec. 601 weights `0.299 R + 0.587 G + 0.114 B` — the same weights every
filter here uses when it needs brightness. The weights are not a mean: green carries most of
the perceived brightness and blue almost none, so a plain average would make a saturated blue
and a saturated green look equally bright.

An already-grey image is returned unchanged, bit for bit. An image whose channels are merely
equal matches the input only to within an ulp, since those weights do not sum to exactly 1 in
binary; the weights are the standard's and are not adjusted to compensate. Only the greyscale
target is supported.
