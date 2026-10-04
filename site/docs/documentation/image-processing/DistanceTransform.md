# DistanceTransform

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`DistanceTransform[image] replaces each pixel by its EXACT Euclidean distance to the nearest background pixel; background pixels are 0, so the value rises toward the interior of a blob. DistanceTransform[image, t] takes pixels above t as foreground (default 0). Exact rather than the classic two-pass chamfer approximation, which cannot represent sqrt(2) with integer steps and so gets diagonal distances a few percent wrong -- invisible on a picture and fatal to a test. Uses Felzenszwalb and Huttenlocher's lower-envelope-of-parabolas method, O(n) per row with no sorting. Separability is EXACT here because squared Euclidean distance is a sum over the axes, so minimising it decomposes per axis; the square root is taken once at the end rather than per pass.`**

## Examples (36)

Every input below was run against the current Mathilda build and its output recorded.

### Basic Examples (6)

```mathematica
In[1]:= ImageData[DistanceTransform[Image[{{1., 1., 1., 1., 1.}, {1., 1., 1., 1., 1.}, {1., 1., 0., 1., 1.}, {1., 1., 1., 1., 1.}, {1., 1., 1., 1., 1.}}]]]
Out[1]= {{2.82843, 2.23607, 2.0, 2.23607, 2.82843}, {2.23607, 1.41421, 1.0, 1.41421, 2.23607}, {2.0, 1.0, 0.0, 1.0, 2.0}, {2.23607, 1.41421, 1.0, 1.41421, 2.23607}, {2.82843, 2.23607, 2.0, 2.23607, 2.82843}}

In[2]:= chk = Image[Table[If[Mod[Quotient[i - 1, 2] + Quotient[j - 1, 2], 2] == 0, 0., 1.], {i, 1, 16}, {j, 1, 16}], "Real"];

In[3]:= disk = Image[Table[N[Boole[(i - 8.5)^2 + (j - 8.5)^2 <= 25]], {i, 1, 16}, {j, 1, 16}], "Real"];

In[4]:= DistanceTransform[chk]
Out[4]= -Image-

In[5]:= DistanceTransform[disk]
Out[5]= -Image-

In[6]:= ImageDimensions[DistanceTransform[chk]]
Out[6]= {16, 16}
```

### Scope (20)

```mathematica
In[7]:= ramp = Image[Table[N[(j - 1)/15], {i, 1, 16}, {j, 1, 16}], "Real"];

In[8]:= zone = Image[Table[N[(1 + Cos[((i - 16)^2 + (j - 16)^2)/40.])/2], {i, 1, 32}, {j, 1, 32}], "Real"];

In[9]:= noise = Image[Table[N[Mod[i*37 + j*17, 101]]/101, {i, 1, 32}, {j, 1, 32}], "Real"];

In[10]:= rgb = Image[Table[{N[i/16], N[j/16], 0.5}, {i, 1, 16}, {j, 1, 16}], "Real"];

In[11]:= sky = Image[Table[{N[0.15 + 0.7 (16 - i)/16], N[0.35 + 0.45 (16 - i)/16], N[0.85 - 0.35 (16 - i)/16]}, {i, 1, 16}, {j, 1, 24}], "Real"];

In[12]:= bit = Image[Table[Boole[Mod[i + j, 2] == 0], {i, 1, 8}, {j, 1, 8}]];

In[13]:= byte = Image[Table[Mod[i*13 + j*7, 256], {i, 1, 16}, {j, 1, 16}]];

In[14]:= vol = Image3D[Table[N[Mod[z*7 + y*13 + x*3, 97]]/97, {z, 1, 8}, {y, 1, 10}, {x, 1, 12}], "Real"];

In[15]:= volb = Image3D[Table[N[Boole[x <= 6 && y <= 5]], {z, 1, 8}, {y, 1, 10}, {x, 1, 12}], "Real"];

In[16]:= DistanceTransform[ramp]
Out[16]= -Image-

In[17]:= DistanceTransform[zone]
Out[17]= -Image-

In[18]:= DistanceTransform[noise]
Out[18]= -Image-

In[19]:= DistanceTransform[rgb]
Out[19]= -Image-

In[20]:= DistanceTransform[sky]
Out[20]= -Image-

In[21]:= DistanceTransform[bit]
Out[21]= -Image-

In[22]:= DistanceTransform[byte]
Out[22]= -Image-

In[23]:= DistanceTransform[vol]
Out[23]= -Image-

In[24]:= DistanceTransform[volb]
Out[24]= -Image-

In[25]:= ImageChannels[DistanceTransform[rgb]]
Out[25]= 1

In[26]:= ImageDimensions[DistanceTransform[vol]]
Out[26]= {12, 10, 8}
```

### Applications (2)

```mathematica
In[27]:= disk = Image[Table[N[Boole[(i - 8.5)^2 + (j - 8.5)^2 <= 25]], {i, 1, 16}, {j, 1, 16}], "Real"];

In[28]:= Dilation[DistanceTransform[disk], 1]
Out[28]= -Image-
```

### Properties & Relations (3)

```mathematica
In[29]:= chk = Image[Table[If[Mod[Quotient[i - 1, 2] + Quotient[j - 1, 2], 2] == 0, 0., 1.], {i, 1, 16}, {j, 1, 16}], "Real"];

In[30]:= ImageDimensions[DistanceTransform[chk]] === ImageDimensions[chk]
Out[30]= True

In[31]:= Min[Flatten[ImageData[DistanceTransform[chk]]]] >= 0.0
Out[31]= True
```

### Neat Examples (2)

```mathematica
In[32]:= zone = Image[Table[N[(1 + Cos[((i - 16)^2 + (j - 16)^2)/40.])/2], {i, 1, 32}, {j, 1, 32}], "Real"];

In[33]:= DistanceTransform[zone]
Out[33]= -Image-
```

### Applications (3)

The centre is distance 1 from the border

```mathematica
In[34]:= ImageData[DistanceTransform[Image[{{0, 0, 0}, {0, 1, 0}, {0, 0, 0}}]]]
Out[34]= {{0.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 0.0}}
```

The blob centre is distance 2

```mathematica
In[35]:= Part[ImageData[DistanceTransform[Image[{{0, 0, 0, 0, 0}, {0, 1, 1, 1, 0}, {0, 1, 1, 1, 0}, {0, 1, 1, 1, 0}, {0, 0, 0, 0, 0}}]]], 3, 3]
Out[35]= 2.0
```

Three across and four down reads exactly 5

```mathematica
In[36]:= Part[ImageData[DistanceTransform[Image[{{0, 1, 1, 1}, {1, 1, 1, 1}, {1, 1, 1, 1}, {1, 1, 1, 1}, {1, 1, 1, 1}}]]], 5, 4]
Out[36]= 5.0
```

## Implementation notes

**Algorithm.** `builtin_distancetransform` replaces each pixel by its **exact** Euclidean
distance to the nearest background pixel. `DistanceTransform[image, t]` takes pixels above `t`
as foreground (default 0). The image is reduced to grey and seeded: `0` at background,
`DT_INF` at foreground, so the transform reports the distance to the nearest *seed*. It then
runs one lower-envelope pass per axis (`dt_1d`): for a line, `D(x) = min_y ((x−y)² + f(y))` is
the lower envelope of identical parabolas translated to each `y` and raised by `f(y)`. Because
all parabolas share curvature, any two intersect once, so the envelope is built in a single
sweep maintaining a stack of visible parabolas — `O(n)` per row, no sorting.

Exactness comes from choosing the envelope method over the classic two-pass **chamfer**
transform, which cannot represent `√2` with integer steps and gets diagonal distances a few
percent wrong — invisible on a picture, fatal to a test. Separability is **exact** here (unlike
the median's) because squared Euclidean distance is a *sum* over axes, so minimising it
decomposes per axis; the square root is taken once at the end, not per pass (per-axis roots
would be wrong, not merely slow). A 3-4-5 triangle reads exactly 5. A volume adds a third
`dt_1d` pass.

**Data structures.** A grey seed buffer `f` plus scratch arrays `v` (parabola indices) and `z`
(intersections) per line; the result through `image_build_real` as a packed single-channel
`"Real"` image (every image head returns a packed buffer — `make check-image-packing`).

**Complexity / limits.** `O(pixels)` — linear in the number of pixels, independent of how far
the nearest seed is. An all-foreground image has no seed, so every pixel reads `DT_INF`.

**Attributes:** `Protected`.

## References

**See also:** [Image3D](../../image-processing/Image3D/)

- P. F. Felzenszwalb and D. P. Huttenlocher, *Distance Transforms of Sampled Functions*, Theory of Computing **8** (2012) 415-428.
- Source: [`src/imagefilter.c`](https://github.com/stblake/mathilda/blob/main/src/imagefilter.c)
- Specification: [`docs/spec/builtins/image-processing.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/image-processing.md)
- Tests: [`tests/test_image.c`](https://github.com/stblake/mathilda/blob/main/tests/test_image.c)

## Notes & additional examples

### Notes

`DistanceTransform[image]` replaces each pixel by its **exact** Euclidean distance to the
nearest background pixel; background pixels are 0, so the value rises toward the interior of a
blob. `DistanceTransform[image, t]` takes pixels above `t` as foreground.

Exact rather than the classic two-pass chamfer approximation, which cannot represent `√2` with
integer steps and gets diagonal distances a few percent wrong — invisible on a picture, fatal
to a test. It uses Felzenszwalb and Huttenlocher's lower-envelope-of-parabolas method,
`O(n)` per row; separability is exact because squared Euclidean distance is a sum over the
axes, and the square root is taken once at the end.
