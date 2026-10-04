### Worked examples

```mathematica
(* integer data in 0..255 is inferred "Byte" at construction *)
In[1]:= ImageType[Image[{{0, 255}}, "Byte"]]
```

```mathematica
(* real-valued data is "Real" *)
In[1]:= ImageType[Image[{{0., 1.}}]]
```

```mathematica
(* integer data confined to {0, 1} infers "Bit" *)
In[1]:= ImageType[Image[{{0, 1}, {1, 0}}]]
```

```mathematica
(* volumes carry a type too *)
In[1]:= ImageType[Image3D[{{{0, 1}}, {{1, 0}}}, "Bit"]]
```

### Notes

`ImageType[image]` gives the pixel type as `"Bit"`, `"Byte"`, `"Bit16"` or
`"Real"`. The type is stored on the image, set when it was constructed, not
re-derived from the pixels on each call.

The type fixes the *range* of a stored value — `"Bit"` is `{0, 1}`, `"Byte"` is
`0..255`, `"Bit16"` is `0..65535`, `"Real"` is already the unit interval — and
that is exactly what makes `ImageData`'s scaling well defined: a `"Byte"` 255
scales to `1.0`, a `"Real"` 1.0 stays `1.0`. Type inference at construction reads
only the values: all-integer in `{0, 1}` is `"Bit"`, all-integer in `0..255` is
`"Byte"`, anything else `"Real"`, so `Image[{{0, 1}}]` is a bit image while
`Image[{{0., 1.}}]` is a real one.
