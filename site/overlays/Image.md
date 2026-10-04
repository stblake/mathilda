### Worked examples

```mathematica
In[1]:= Image[{{0., 1.}, {1., 0.}}]  (* a small real-typed checkerboard, stored as a packed buffer *)
```

```mathematica
In[1]:= ImageType[Image[{{0, 1}, {1, 0}}]]  (* all values in {0,1}, so the type infers to Bit *)
```

```mathematica
In[1]:= ImageType[Image[{{0, 128, 255}}]]  (* integers to 255 infer to Byte *)
```

```mathematica
In[1]:= ImageData[Image[{{0, 300}}, "Byte"]]  (* a stated type rounds and clips: 300 -> 255, scaled out *)
```

```mathematica
In[1]:= ImageDimensions[Image[{{0., 1., 0.}}]]  (* width x height, transposed from the data's rows *)
```

### Notes

`Image[data]` normalises to the canonical `Image[data, type]`, which is also real Wolfram
syntax; `ImageQ` tests for exactly that form. The type is inferred from the values alone, so
`Image[{{0,1}}]` is a `"Bit"` image and `Image[{{0.,1.}}]` a `"Real"` one — a distinction a
caller can rely on.

The data is stored height × width (× channels), which is **transposed** relative to
`ImageDimensions`'s `{width, height}` — the single most common source of silently-wrong image
code, and Mathematica's convention. `ImageData` scales stored values back into `[0, 1]` using
the type's range, so the type is not decoration: it is what makes that scaling well defined.

Computed and constructed images share one representation — a visible packed `NDArray` — so
`===` compares two images with identical pixels as equal regardless of how each was built.
