### Worked examples

```mathematica
In[1]:= ImageDimensions[ImagePad[Image[{{0., 1.}}], 1]]  (* pad one pixel on every side: 2x1 -> 4x3 *)
```

```mathematica
In[1]:= ImageData[ImagePad[Image[{{0.5}}], 1]]  (* a 1x1 image padded with the default value 0 *)
```

```mathematica
In[1]:= ImageData[ImagePad[Image[{{0.5}}], 1, "Fixed"]]  (* the Fixed mode replicates the edge pixel instead *)
```

```mathematica
In[1]:= ImageData[ImagePad[Image[{{1., 2., 3.}}], {{1, 1}, {0, 0}}, "Reflected"]]  (* reflection mirrors without repeating the edge *)
```

### Notes

`ImagePad[image, m]` pads `m` pixels on every side; `ImagePad[image, {{left, right}, {bottom,
top}}]` names each side in Mathematica's visual order — so `top` adds rows at the *start* of
the data, since row 1 is the top of the image. Negative amounts crop but may not erase the
image.

The third argument sets the fill: a value (default 0); `"Fixed"`, which replicates the edge
pixel — the same boundary rule the filters use, so padding then filtering composes with it; or
`"Reflected"`, which mirrors *without* repeating the edge (`{1,2,3}` padded by 1 →
`{2,1,2,3,2}`, not `{1,1,2,3,3}`), because doubling the edge sample biases any later average
toward the border. `ImagePad` and `ImageCrop` are exact inverses.
