### Worked examples

```mathematica
In[1]:= ImageDimensions[ImageAssemble[{Image[{{0., 1.}}], Image[{{1., 0.}}]}]]  (* two 2x1 tiles side by side: width 4 *)
```

```mathematica
In[1]:= ImageData[ImageAssemble[{Image[{{0., 1.}}], Image[{{1., 0.}}]}]]  (* the tiles laid out left to right *)
```

```mathematica
In[1]:= ImageDimensions[ImageAssemble[{{Image[{{0.}}], Image[{{1.}}]}, {Image[{{1.}}], Image[{{0.}}]}}]]  (* a 2x2 grid *)
```

### Notes

`ImageAssemble[{{a, b}, {c, d}}]` tiles a grid; `ImageAssemble[{a, b}]` makes a single row. A
row is as tall as its tallest tile and a column as wide as its widest; each tile keeps its
natural size and any gap is left blank rather than stretched, since stretching would resample
an image the caller did not ask to resize.

A grey tile beside a colour one is promoted to colour, and the assembled sheet carries an alpha
channel if any tile did — with its gaps transparent, the property a sheet of sprites needs to
keep.
