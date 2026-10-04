### Worked examples

```mathematica
(* Otsu's threshold lands between the dark pair and the bright pair *)
In[1]:= FindThreshold[Image[{{0.2, 0.3, 0.7, 0.8}}]]
```

```mathematica
(* the level is on the same unit scale as ImageData, so Binarize uses it directly *)
In[1]:= Binarize[Image[{{0.2, 0.3, 0.7, 0.8}}], FindThreshold[Image[{{0.2, 0.3, 0.7, 0.8}}]]]
```

```mathematica
(* an image of one single value has no two classes, so it stays unevaluated *)
In[1]:= FindThreshold[Image[{{0.5, 0.5}, {0.5, 0.5}}]]
```

### Notes

`FindThreshold[image]` gives a threshold separating the image into two classes by
**Otsu's method**: the level maximising the between-class variance `w0 w1 (mu0 −
mu1)²`, which is algebraically the same as minimising the weighted within-class
variance but needs only one incremental pass over a 256-bin histogram.

A colour image is reduced to Rec. 601 luminance first. The returned level sits at
the upper edge of its histogram bin and is on the same unit scale as `ImageData`,
so it feeds straight into `Binarize[image, t]`. It returns unevaluated for an
image whose pixels are all identical, since no threshold splits one cluster into
two. The same routine backs `Binarize`'s default threshold and `EdgeDetect`'s
hysteresis high threshold.
