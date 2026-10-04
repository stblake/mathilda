### Worked examples

```mathematica
(* Otsu's threshold splits dark from bright, giving a "Bit" image *)
In[1]:= Binarize[Image[{{0.1, 0.9}, {0.8, 0.2}}]]
```

```mathematica
(* the pixel data: above-threshold pixels become 1, the rest 0 *)
In[1]:= ImageData[Binarize[Image[{{0.1, 0.9}, {0.8, 0.2}}]]]
```

```mathematica
(* an explicit threshold: strictly above t becomes 1, so 0.4 at t=0.5 stays 0 *)
In[1]:= ImageData[Binarize[Image[{{0.1, 0.4, 0.9}}], 0.5]]
```

```mathematica
(* the result is typed "Bit" by construction *)
In[1]:= ImageType[Binarize[Image[{{0.2, 0.9}, {0.8, 0.1}}]]]
```

### Notes

`Binarize[image]` thresholds by Otsu's method (see `FindThreshold`), giving a
`"Bit"` image; `Binarize[image, t]` thresholds at the explicit level `t`.

A pixel **strictly above** the threshold becomes `1`, so a pixel exactly at it
becomes `0` — which matters, because "above" and "at or above" differ on exactly
the pixels a threshold was chosen to sit between. A colour image is reduced to
luminance first. The result is typed `"Bit"` rather than `"Real"` because it is
binary by construction: `ImageData` then scales nothing and the image still
reports itself as binary. With the default (Otsu) threshold a single-value image
declines, since no threshold splits one cluster.
