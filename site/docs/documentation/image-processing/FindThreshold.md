# FindThreshold

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`FindThreshold[image] gives a threshold separating the image into two classes, by Otsu's method: the level maximising the BETWEEN-class variance w0 w1 (mu0 - mu1)^2, which is algebraically the same as minimising the weighted within-class variance but needs only one incremental pass over a 256-bin histogram. A colour image is reduced to Rec. 601 luminance first. Returns unevaluated for an image whose pixels are all identical, since no threshold splits one cluster into two.`**

## Examples (28)

Every input below was run against the current Mathilda build and its output recorded.

### Basic Examples (9)

```mathematica
In[1]:= FindThreshold[Image[{{0., 0., 1.}, {0., 1., 1.}}]]
Out[1]= 0.00196078

In[2]:= chk = Image[Table[If[Mod[Quotient[i - 1, 2] + Quotient[j - 1, 2], 2] == 0, 0., 1.], {i, 1, 16}, {j, 1, 16}], "Real"];

In[3]:= rgb = Image[Table[{N[i/16], N[j/16], 0.5}, {i, 1, 16}, {j, 1, 16}], "Real"];

In[4]:= bit = Image[Table[Boole[Mod[i + j, 2] == 0], {i, 1, 8}, {j, 1, 8}]];

In[5]:= byte = Image[Table[Mod[i*13 + j*7, 256], {i, 1, 16}, {j, 1, 16}]];

In[6]:= FindThreshold[chk]
Out[6]= 0.00196078

In[7]:= FindThreshold[rgb]
Out[7]= 0.523529

In[8]:= FindThreshold[bit]
Out[8]= 0.00196078

In[9]:= FindThreshold[byte]
Out[9]= 0.519608
```

### Scope (12)

```mathematica
In[10]:= disk = Image[Table[N[Boole[(i - 8.5)^2 + (j - 8.5)^2 <= 25]], {i, 1, 16}, {j, 1, 16}], "Real"];

In[11]:= ramp = Image[Table[N[(j - 1)/15], {i, 1, 16}, {j, 1, 16}], "Real"];

In[12]:= zone = Image[Table[N[(1 + Cos[((i - 16)^2 + (j - 16)^2)/40.])/2], {i, 1, 32}, {j, 1, 32}], "Real"];

In[13]:= noise = Image[Table[N[Mod[i*37 + j*17, 101]]/101, {i, 1, 32}, {j, 1, 32}], "Real"];

In[14]:= rgb = Image[Table[{N[i/16], N[j/16], 0.5}, {i, 1, 16}, {j, 1, 16}], "Real"];

In[15]:= sky = Image[Table[{N[0.15 + 0.7 (16 - i)/16], N[0.35 + 0.45 (16 - i)/16], N[0.85 - 0.35 (16 - i)/16]}, {i, 1, 16}, {j, 1, 24}], "Real"];

In[16]:= FindThreshold[disk]
Out[16]= 0.00196078

In[17]:= FindThreshold[ramp]
Out[17]= 0.468627

In[18]:= FindThreshold[zone]
Out[18]= 0.488235

In[19]:= FindThreshold[noise]
Out[19]= 0.488235

In[20]:= FindThreshold[sky]
Out[20]= 0.539216

In[21]:= FindThreshold[Import[Export["/tmp/mathilda_ex.png", rgb]]]
Out[21]= 0.523529
```

### Applications (2)

```mathematica
In[22]:= chk = Image[Table[If[Mod[Quotient[i - 1, 2] + Quotient[j - 1, 2], 2] == 0, 0., 1.], {i, 1, 16}, {j, 1, 16}], "Real"];

In[23]:= Table[FindThreshold[GaussianFilter[chk, r]], {r, 1, 3}]
Out[23]= {0.194118, 0.460784, 0.382353}
```

### Neat Examples (2)

```mathematica
In[24]:= zone = Image[Table[N[(1 + Cos[((i - 16)^2 + (j - 16)^2)/40.])/2], {i, 1, 32}, {j, 1, 32}], "Real"];

In[25]:= FindThreshold[zone]
Out[25]= 0.488235
```

### Applications (3)

```mathematica
In[26]:= FindThreshold[Image[{{0.2, 0.3, 0.7, 0.8}}]]
Out[26]= 0.303922

In[27]:= Binarize[Image[{{0.2, 0.3, 0.7, 0.8}}], FindThreshold[Image[{{0.2, 0.3, 0.7, 0.8}}]]]
Out[27]= -Image-

In[28]:= FindThreshold[Image[{{0.5, 0.5}, {0.5, 0.5}}]]
Out[28]= FindThreshold[Image[NDArray[{{0.5, 0.5}, {0.5, 0.5}}], "Real"]]
```

## Implementation notes

**Algorithm.** `builtin_findthreshold` returns the single level that best
separates an image into two classes by **Otsu's method**. It reduces the image to
a Rec. 601 luminance plane (`img_grey_plane`) and calls `img_otsu`, which bins
the pixels into a 256-bin histogram (`IMG_OTSU_BINS`), with out-of-range values
clamped **into** the extreme bins rather than dropped so they still vote. It then
sweeps the threshold `t` over the bins, maintaining incrementally the class
weight `w0` and first moment `sum0`, and maximises the **between-class variance**
`sb = w0 · w1 · (mu0 − mu1)²` — algebraically equivalent to minimising the
weighted within-class variance but computable in one incremental pass. A NaN
pixel aborts (it has no bin); a strictly-greater comparison keeps the **lowest**
winning bin, so ties are deterministic. The returned threshold sits at the
**upper edge** of the winning bin, `(best_bin + 0.5)/255`, so a pixel in that bin
falls on the low side. An image whose pixels are all identical has only one
occupied cluster and no split exists, so the routine returns `false` and the head
is left unevaluated.

**Data structures.** A fixed 256-element `double` histogram on the stack; a
single flat luminance `double` plane loaded from the image. No per-threshold
allocation — the sweep is `O(bins)`.

**Complexity / limits.** `O(width · height)` to build the histogram plus
`O(bins)` for the sweep; the threshold is a `double` on the same unit scale as
`ImageData`. Returns unevaluated for a degenerate (single-value) image. The same
`img_otsu` backs `Binarize`'s default threshold and `EdgeDetect`'s high
threshold.

**Attributes:** `Protected`.

## References

**See also:** [Image](../../image-processing/Image/)

- N. Otsu, *A threshold selection method from gray-level histograms*, IEEE Trans. Systems, Man, and Cybernetics **9** (1979) 62-66.
- Source: [`src/imagefilter.c`](https://github.com/stblake/mathilda/blob/main/src/imagefilter.c)
- Specification: [`docs/spec/builtins/image-processing.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/image-processing.md)
- Tests: [`tests/test_image.c`](https://github.com/stblake/mathilda/blob/main/tests/test_image.c)

## Notes & additional examples

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
