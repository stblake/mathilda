# GaussianFilter

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`GaussianFilter[image, r] blurs image with a Gaussian of radius r. It is exactly ImageConvolve[image, GaussianMatrix[r]] -- the same matrix through the same convolution, not a second implementation -- and a test asserts the identity.`**

## Examples (50)

Every input below was run against the current Mathilda build and its output recorded.

### Basic Examples (5)

```mathematica
In[1]:= ImageData[GaussianFilter[Image[{{0., 0., 0.}, {0., 1., 0.}, {0., 0., 0.}}], 1]]
Out[1]= {{0.0113437, 0.0838195, 0.0113437}, {0.0838195, 0.619347, 0.0838195}, {0.0113437, 0.0838195, 0.0113437}}

In[2]:= chk = Image[Table[If[Mod[Quotient[i - 1, 2] + Quotient[j - 1, 2], 2] == 0, 0., 1.], {i, 1, 16}, {j, 1, 16}], "Real"];

In[3]:= GaussianFilter[chk, 2]
Out[3]= -Image-

In[4]:= ImageDimensions[GaussianFilter[chk, 2]]
Out[4]= {16, 16}

In[5]:= ImageType[GaussianFilter[chk, 1]]
Out[5]= "Real"
```

### Scope (23)

```mathematica
In[6]:= chk = Image[Table[If[Mod[Quotient[i - 1, 2] + Quotient[j - 1, 2], 2] == 0, 0., 1.], {i, 1, 16}, {j, 1, 16}], "Real"];

In[7]:= disk = Image[Table[N[Boole[(i - 8.5)^2 + (j - 8.5)^2 <= 25]], {i, 1, 16}, {j, 1, 16}], "Real"];

In[8]:= ramp = Image[Table[N[(j - 1)/15], {i, 1, 16}, {j, 1, 16}], "Real"];

In[9]:= zone = Image[Table[N[(1 + Cos[((i - 16)^2 + (j - 16)^2)/40.])/2], {i, 1, 32}, {j, 1, 32}], "Real"];

In[10]:= noise = Image[Table[N[Mod[i*37 + j*17, 101]]/101, {i, 1, 32}, {j, 1, 32}], "Real"];

In[11]:= rgb = Image[Table[{N[i/16], N[j/16], 0.5}, {i, 1, 16}, {j, 1, 16}], "Real"];

In[12]:= sky = Image[Table[{N[0.15 + 0.7 (16 - i)/16], N[0.35 + 0.45 (16 - i)/16], N[0.85 - 0.35 (16 - i)/16]}, {i, 1, 16}, {j, 1, 24}], "Real"];

In[13]:= bit = Image[Table[Boole[Mod[i + j, 2] == 0], {i, 1, 8}, {j, 1, 8}]];

In[14]:= byte = Image[Table[Mod[i*13 + j*7, 256], {i, 1, 16}, {j, 1, 16}]];

In[15]:= vol = Image3D[Table[N[Mod[z*7 + y*13 + x*3, 97]]/97, {z, 1, 8}, {y, 1, 10}, {x, 1, 12}], "Real"];

In[16]:= GaussianFilter[disk, 1]
Out[16]= -Image-

In[17]:= GaussianFilter[ramp, 2]
Out[17]= -Image-

In[18]:= GaussianFilter[zone, 2]
Out[18]= -Image-

In[19]:= GaussianFilter[noise, 3]
Out[19]= -Image-

In[20]:= GaussianFilter[rgb, 1]
Out[20]= -Image-

In[21]:= GaussianFilter[sky, 2]
Out[21]= -Image-

In[22]:= GaussianFilter[bit, 1]
Out[22]= -Image-

In[23]:= GaussianFilter[byte, 2]
Out[23]= -Image-

In[24]:= GaussianFilter[vol, 1]
Out[24]= -Image-

In[25]:= ImageChannels[GaussianFilter[rgb, 2]]
Out[25]= 3

In[26]:= ImageDimensions[GaussianFilter[vol, 1]]
Out[26]= {12, 10, 8}

In[27]:= GaussianFilter[chk, {1, 3}]
Out[27]= -Image-

In[28]:= GaussianFilter[chk, 4]
Out[28]= -Image-
```

### Applications (6)

```mathematica
In[29]:= zone = Image[Table[N[(1 + Cos[((i - 16)^2 + (j - 16)^2)/40.])/2], {i, 1, 32}, {j, 1, 32}], "Real"];

In[30]:= noise = Image[Table[N[Mod[i*37 + j*17, 101]]/101, {i, 1, 32}, {j, 1, 32}], "Real"];

In[31]:= rgb = Image[Table[{N[i/16], N[j/16], 0.5}, {i, 1, 16}, {j, 1, 16}], "Real"];

In[32]:= Binarize[GaussianFilter[noise, 2]]
Out[32]= -Image-

In[33]:= EdgeDetect[GaussianFilter[zone, 2]]
Out[33]= -Image-

In[34]:= ImageDimensions[GaussianFilter[Import[Export["/tmp/mathilda_ex.png", rgb]], 2]]
Out[34]= {16, 16}
```

### Properties & Relations (10)

```mathematica
In[35]:= chk = Image[Table[If[Mod[Quotient[i - 1, 2] + Quotient[j - 1, 2], 2] == 0, 0., 1.], {i, 1, 16}, {j, 1, 16}], "Real"];

In[36]:= ramp = Image[Table[N[(j - 1)/15], {i, 1, 16}, {j, 1, 16}], "Real"];

In[37]:= rgb = Image[Table[{N[i/16], N[j/16], 0.5}, {i, 1, 16}, {j, 1, 16}], "Real"];

In[38]:= vol = Image3D[Table[N[Mod[z*7 + y*13 + x*3, 97]]/97, {z, 1, 8}, {y, 1, 10}, {x, 1, 12}], "Real"];

In[39]:= ImageDimensions[GaussianFilter[chk, 3]] === ImageDimensions[chk]
Out[39]= True

In[40]:= ImageChannels[GaussianFilter[rgb, 2]] === ImageChannels[rgb]
Out[40]= True

In[41]:= ImageData[GaussianFilter[ramp, 0]] === ImageData[ramp]
Out[41]= True

In[42]:= Max[Flatten[ImageData[GaussianFilter[chk, 2]]]] <= 1.0
Out[42]= True

In[43]:= Min[Flatten[ImageData[GaussianFilter[chk, 2]]]] >= 0.0
Out[43]= True

In[44]:= ImageDimensions[GaussianFilter[vol, 2]] === ImageDimensions[vol]
Out[44]= True
```

### Neat Examples (3)

```mathematica
In[45]:= zone = Image[Table[N[(1 + Cos[((i - 16)^2 + (j - 16)^2)/40.])/2], {i, 1, 32}, {j, 1, 32}], "Real"];

In[46]:= GaussianFilter[zone, 4]
Out[46]= -Image-

In[47]:= GaussianFilter[zone, 1]
Out[47]= -Image-
```

### Applications (3)

A point spreads into the kernel

```mathematica
In[48]:= GaussianFilter[Image[{{0., 0., 0.}, {0., 1., 0.}, {0., 0., 0.}}], 1]
Out[48]= -Image-
```

The centre keeps most of its weight

```mathematica
In[49]:= Part[ImageData[GaussianFilter[Image[{{0., 0., 0.}, {0., 1., 0.}, {0., 0., 0.}}], 1]], 2, 2]
Out[49]= 0.619347
```

A constant-sum kernel preserves overall brightness

```mathematica
In[50]:= GaussianFilter[Image[{{0.3, 0.9}, {0.9, 0.3}}], 1]
Out[50]= -Image-
```

## Implementation notes

**Algorithm.** `builtin_gaussianfilter` blurs an image with a Gaussian of radius `r`. It is
defined as exactly `ImageConvolve[image, GaussianMatrix[r]]` and **implemented** that way —
`builtin_gaussianmatrix` builds the `(2r+1) × (2r+1)` matrix (`exp(−(dx²+dy²)/(2σ²))`, σ = r/2
by default, normalised by the *realised* sum rather than the analytic `2πσ²` so a truncated
kernel still sums to 1 and does not darken the image), and the same `convolve_dispatch` every
filter uses runs it. Two independent implementations of one identity is how the identity
quietly stops holding, so there is only one. A volume takes `gauss3_kernel` + `convolve3_run`.

The convolution is true convolution (kernel reflected), with replicate ("Fixed") padding — a
constant image convolved with a kernel summing to 1 comes back unchanged everywhere, edges
included, where zero padding would darken them. The Gaussian is symmetric, so reflection is
invisible here; it matters only on asymmetric kernels. The result is always `"Real"` — a
Gaussian of bytes is not a byte.

**Data structures.** A decoded unit buffer and a dense kernel matrix; the result through
`image_build_real` as a packed `"Real"` image (every image head returns a packed buffer —
`make check-image-packing`). `convolve_dispatch` re-derives the kernel's separability, so a
rank-1 Gaussian costs `kw + kh` taps per pixel, not `kw · kh`.

**Complexity / limits.** `O(pixels · (kw + kh))` via the separable path. The radius must be a
non-negative integer (≤ 512 planar, ≤ 32 volumetric).

**Attributes:** `Protected`.

## References

- Source: [`src/imagefilter.c`](https://github.com/stblake/mathilda/blob/main/src/imagefilter.c)
- Specification: [`docs/spec/builtins/image-processing.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/image-processing.md)
- Tests: [`tests/test_image.c`](https://github.com/stblake/mathilda/blob/main/tests/test_image.c)

## Notes & additional examples

### Notes

`GaussianFilter[image, r]` blurs with a Gaussian of radius `r`. It is exactly
`ImageConvolve[image, GaussianMatrix[r]]` — the same matrix through the same convolution, not a
second implementation — and a test asserts the identity.

The kernel is normalised by its realised sum (not the analytic `2 π σ²`), so a truncated kernel
still sums to 1 and does not darken the image on each pass. Padding replicates the border, so a
constant image comes back unchanged everywhere including the edges. The result is a `"Real"`
image, since a Gaussian of bytes is not a byte.
