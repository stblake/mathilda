### Worked examples

```mathematica
(* a "Byte" 255 scales out to exactly 1.0 on the unit interval *)
In[1]:= ImageData[Image[{{0, 255}, {255, 0}}, "Byte"]]
```

```mathematica
(* a second "Byte" argument asks for the stored values, unscaled *)
In[1]:= ImageData[Image[{{0, 255}}, "Byte"], "Byte"]
```

```mathematica
(* a "Real" image is already unit-scaled, so ImageData returns it as stored *)
In[1]:= ImageData[Image[{{0.25, 0.5}, {0.75, 1.}}]]
```

### Notes

`ImageData[image]` gives the pixel array as reals in `[0, 1]`, scaling out the
image's type — the array is height × width, or height × width × channels
(interleaved) for a colour image, the same shape the image was built from.

`ImageData[image, type]` gives the stored values **unscaled** instead, but only
when `type` is the image's own type: converting between types is a separate
operation with its own rounding, and this never does it silently — a mismatched
type returns unevaluated. The result shape is exactly the one `Image[...]` was
handed, which is what makes `Image[ImageData[img] ..., type]` a faithful
round-trip. It accepts a volume as well as a plane.
