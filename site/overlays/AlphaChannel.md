### Worked examples

```mathematica
In[1]:= ImageData[AlphaChannel[Image[{{{1., 0., 0., 0.5}}}]]]  (* an RGBA pixel's alpha is channel 4 *)
```

```mathematica
In[1]:= AlphaChannel[Image[{{0., 1.}, {1., 0.}}]]  (* a grey image has no alpha: the result is all-opaque *)
```

```mathematica
In[1]:= ImageData[AlphaChannel[Image[{{0., 1.}, {1., 0.}}]]]  (* which reads back as 1 everywhere *)
```

### Notes

`AlphaChannel` extracts the opacity as a one-channel image. Two channels are read as
grey+alpha and four as RGB+alpha; any other channel count — a grey or an RGB image — is
treated as fully opaque, and the result is an all-`1` image rather than a decline. "How
transparent is this?" has an answer for every image.

The result is always a `"Real"` single-channel image handed back as a packed buffer, so it
feeds straight into another filter.
