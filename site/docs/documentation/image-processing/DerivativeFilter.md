# DerivativeFilter

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`DerivativeFilter[image, {n, m}] gives the n-th derivative down the rows and the m-th across the columns, each order from 0 to 2. The kernel is a separable outer product of 1-D stencils: order 0 is the smoothing {1,2,1}/4, order 1 the central difference {-1,0,1}/2, order 2 the second difference {1,-2,1}. So {0,1} is Sobel-x and {1,0} is Sobel-y. The stencils are NORMALISED, unlike the raw integer Sobel kernels, which report a gradient eight times the true slope -- harmless when only the ranking of edges matters, and wrong for anything that reads the number. On f(x) = c x the first derivative gives exactly c. The result is a "Real" image.`**

## Examples (30)

Every input below was run against the current Mathilda build and its output recorded.

### Basic Examples (6)

```mathematica
In[1]:= ImageData[DerivativeFilter[Image[{{0., 0., 1., 1.}, {0., 0., 1., 1.}}], {0, 1}]]
Out[1]= {{0.0, 0.5, 0.5, 0.0}, {0.0, 0.5, 0.5, 0.0}}

In[2]:= chk = Image[Table[If[Mod[Quotient[i - 1, 2] + Quotient[j - 1, 2], 2] == 0, 0., 1.], {i, 1, 16}, {j, 1, 16}], "Real"];

In[3]:= ramp = Image[Table[N[(j - 1)/15], {i, 1, 16}, {j, 1, 16}], "Real"];

In[4]:= DerivativeFilter[ramp, {0, 1}]
Out[4]= -Image-

In[5]:= DerivativeFilter[ramp, {1, 0}]
Out[5]= -Image-

In[6]:= ImageDimensions[DerivativeFilter[chk, {0, 1}]]
Out[6]= {16, 16}
```

### Scope (13)

```mathematica
In[7]:= chk = Image[Table[If[Mod[Quotient[i - 1, 2] + Quotient[j - 1, 2], 2] == 0, 0., 1.], {i, 1, 16}, {j, 1, 16}], "Real"];

In[8]:= disk = Image[Table[N[Boole[(i - 8.5)^2 + (j - 8.5)^2 <= 25]], {i, 1, 16}, {j, 1, 16}], "Real"];

In[9]:= zone = Image[Table[N[(1 + Cos[((i - 16)^2 + (j - 16)^2)/40.])/2], {i, 1, 32}, {j, 1, 32}], "Real"];

In[10]:= rgb = Image[Table[{N[i/16], N[j/16], 0.5}, {i, 1, 16}, {j, 1, 16}], "Real"];

In[11]:= vol = Image3D[Table[N[Mod[z*7 + y*13 + x*3, 97]]/97, {z, 1, 8}, {y, 1, 10}, {x, 1, 12}], "Real"];

In[12]:= DerivativeFilter[chk, {0, 1}]
Out[12]= -Image-

In[13]:= DerivativeFilter[chk, {1, 1}]
Out[13]= -Image-

In[14]:= DerivativeFilter[disk, {0, 2}]
Out[14]= -Image-

In[15]:= DerivativeFilter[zone, {0, 1}]
Out[15]= -Image-

In[16]:= DerivativeFilter[rgb, {0, 1}]
Out[16]= -Image-

In[17]:= DerivativeFilter[vol, {0, 0, 1}]
Out[17]= -Image-

In[18]:= ImageChannels[DerivativeFilter[rgb, {0, 1}]]
Out[18]= 3

In[19]:= ImageDimensions[DerivativeFilter[vol, {0, 1, 0}]]
Out[19]= {12, 10, 8}
```

### Applications (2)

```mathematica
In[20]:= zone = Image[Table[N[(1 + Cos[((i - 16)^2 + (j - 16)^2)/40.])/2], {i, 1, 32}, {j, 1, 32}], "Real"];

In[21]:= Binarize[DerivativeFilter[zone, {0, 1}]]
Out[21]= -Image-
```

### Properties & Relations (4)

```mathematica
In[22]:= chk = Image[Table[If[Mod[Quotient[i - 1, 2] + Quotient[j - 1, 2], 2] == 0, 0., 1.], {i, 1, 16}, {j, 1, 16}], "Real"];

In[23]:= rgb = Image[Table[{N[i/16], N[j/16], 0.5}, {i, 1, 16}, {j, 1, 16}], "Real"];

In[24]:= ImageDimensions[DerivativeFilter[chk, {0, 1}]] === ImageDimensions[chk]
Out[24]= True

In[25]:= ImageChannels[DerivativeFilter[rgb, {1, 0}]] === ImageChannels[rgb]
Out[25]= True
```

### Neat Examples (2)

```mathematica
In[26]:= zone = Image[Table[N[(1 + Cos[((i - 16)^2 + (j - 16)^2)/40.])/2], {i, 1, 32}, {j, 1, 32}], "Real"];

In[27]:= DerivativeFilter[zone, {1, 1}]
Out[27]= -Image-
```

### Applications (3)

A horizontal ramp of slope 1/8

```mathematica
In[28]:= ramp = Image[{{0., 0.125, 0.25}, {0., 0.125, 0.25}, {0., 0.125, 0.25}}];
```

The Sobel-x response reports the slope exactly: 0.125

```mathematica
In[29]:= Part[ImageData[DerivativeFilter[ramp, {0, 1}]], 2, 2]
Out[29]= 0.125
```

The Sobel-y derivative of a single bright pixel

```mathematica
In[30]:= Part[ImageData[DerivativeFilter[Image[{{0., 0., 0.}, {0., 1., 0.}, {0., 0., 0.}}], {1, 0}]], 1, 2]
Out[30]= 0.25
```

## Implementation notes

**Algorithm.** `builtin_derivativefilter` applies a separable derivative kernel.
`DerivativeFilter[image, {n, m}]` is the `n`-th derivative down the rows and the `m`-th across
the columns, each order 0–2. `deriv_kernel` builds the full 2-D kernel as the outer product of
two 1-D stencils from `deriv_stencil`:

- order 0: `{1, 2, 1}/4` — smoothing, normalised to sum 1 (so it preserves a constant, and
  being symmetric, a linear ramp exactly);
- order 1: `{+1/2, 0, −1/2}` — central difference. On `f(x) = c x` it gives exactly `c`;
- order 2: `{1, −2, 1}` — second difference. On `f(x) = c x²` it gives exactly `2c`.

So `{0, 1}` is Sobel-x and `{1, 0}` Sobel-y. Two subtleties the code pins: the order-1 stencil
is written `{+1/2, 0, −1/2}`, **pre-flipped**, because `ImageConvolve` reflects its kernel and
a reading-order stencil would compute the *negated* derivative — caught only by asserting an
exact signed value, since a gradient magnitude squares the sign away. And the stencils are
**normalised** (/4, /2), unlike the raw integer Sobel kernels that report eight times the true
slope. The full matrix is handed to the same `convolve_dispatch` every filter uses, which
re-derives the separable factorisation rather than trusting it. A volume takes `deriv3_run`
with a third stencil.

**Data structures.** One decoded unit buffer; the result through `image_build_real` as a packed
`"Real"` image (every image head returns a packed buffer — `make check-image-packing`).

**Complexity / limits.** The kernel factors to rank 1, so it costs `O(pixels · (kw + kh))`.
Orders are restricted to 0–2; the padding is replicate, as for all convolutions.

**Attributes:** `Protected`.

## References

**See also:** [ImageConvolve](../../image-processing/ImageConvolve/)

- Source: [`src/imagefilter.c`](https://github.com/stblake/mathilda/blob/main/src/imagefilter.c)
- Specification: [`docs/spec/builtins/image-processing.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/image-processing.md)
- Tests: [`tests/test_image.c`](https://github.com/stblake/mathilda/blob/main/tests/test_image.c)

## Notes & additional examples

### Notes

`DerivativeFilter[image, {n, m}]` gives the `n`-th derivative down the rows and the `m`-th
across the columns, each order 0 to 2. The kernel is a separable outer product of 1-D stencils:
order 0 is the smoothing `{1,2,1}/4`, order 1 the central difference `{-1,0,1}/2`, order 2 the
second difference `{1,-2,1}`. So `{0,1}` is Sobel-x and `{1,0}` Sobel-y.

The stencils are **normalised**, unlike the raw integer Sobel kernels that report a gradient
eight times the true slope — harmless when only the ranking of edges matters, wrong for
anything that reads the number. On `f(x) = c x` the first derivative gives exactly `c`. The
result is a `"Real"` image.
