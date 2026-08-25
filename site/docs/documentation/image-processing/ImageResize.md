# ImageResize

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`ImageResize[image, {w, h}] resizes to w x h pixels; ImageResize[image, w] gives width w with the height following to preserve the aspect ratio. Resampling -> "Nearest" | "Bilinear" | "Average" selects the method; the default Automatic uses AREA AVERAGING when either axis shrinks and bilinear otherwise. That default is about aliasing: point-sampling a shrinking image destroys every frequency above half the new sampling rate -- a fine checkerboard reduced by nearest-neighbour comes back a flat field -- and no interpolation afterwards can restore what point-sampling discarded. Area averaging is a box prefilter and a resample in one pass, exact for integer reduction factors, using true fractional coverage so a 3 -> 2 reduction is as correct as 4 -> 2. Enlarging has no frequencies to remove, so bilinear is used there; area averaging on an enlargement would degenerate to nearest. Coordinates are centre-aligned, avoiding the half-pixel shift that sx = i * scale introduces at any scale other than 1:1. The result is a "Real" image; sizes must be positive integers.`**

## Examples (35)

Every input below was run against the current Mathilda build and its output recorded.

### Basic Examples (5)

```mathematica
In[1]:= chk = Image[Table[If[Mod[Quotient[i - 1, 2] + Quotient[j - 1, 2], 2] == 0, 0., 1.], {i, 1, 16}, {j, 1, 16}], "Real"];

In[2]:= disk = Image[Table[N[Boole[(i - 8.5)^2 + (j - 8.5)^2 <= 25]], {i, 1, 16}, {j, 1, 16}], "Real"];

In[3]:= ImageResize[chk, {8, 8}]
Out[3]= -Image-

In[4]:= ImageDimensions[ImageResize[chk, {8, 8}]]
Out[4]= {8, 8}

In[5]:= ImageResize[disk, {12, 12}]
Out[5]= -Image-
```

### Scope (19)

```mathematica
In[6]:= chk = Image[{{0.,1.,0.,1.},{1.,0.,1.,0.},{0.,1.,0.,1.},{1.,0.,1.,0.}}];
```

Area averaging

```mathematica
In[7]:= ImageData[ImageResize[chk, {2, 2}]]
Out[7]= {{0.5, 0.5}, {0.5, 0.5}}
```

```mathematica
In[8]:= ImageData[ImageResize[chk, {2, 2}, Resampling -> "Nearest"]]
Out[8]= {{0.0, 0.0}, {0.0, 0.0}}

In[9]:= ramp = Image[Table[N[(j - 1)/15], {i, 1, 16}, {j, 1, 16}], "Real"];

In[10]:= zone = Image[Table[N[(1 + Cos[((i - 16)^2 + (j - 16)^2)/40.])/2], {i, 1, 32}, {j, 1, 32}], "Real"];

In[11]:= noise = Image[Table[N[Mod[i*37 + j*17, 101]]/101, {i, 1, 32}, {j, 1, 32}], "Real"];

In[12]:= rgb = Image[Table[{N[i/16], N[j/16], 0.5}, {i, 1, 16}, {j, 1, 16}], "Real"];

In[13]:= sky = Image[Table[{N[0.15 + 0.7 (16 - i)/16], N[0.35 + 0.45 (16 - i)/16], N[0.85 - 0.35 (16 - i)/16]}, {i, 1, 16}, {j, 1, 24}], "Real"];

In[14]:= bit = Image[Table[Boole[Mod[i + j, 2] == 0], {i, 1, 8}, {j, 1, 8}]];

In[15]:= byte = Image[Table[Mod[i*13 + j*7, 256], {i, 1, 16}, {j, 1, 16}]];

In[16]:= ImageResize[rgb, {8, 8}]
Out[16]= -Image-

In[17]:= ImageResize[sky, {12, 8}]
Out[17]= -Image-

In[18]:= ImageResize[bit, {4, 4}]
Out[18]= -Image-

In[19]:= ImageResize[byte, {8, 8}]
Out[19]= -Image-

In[20]:= ImageResize[zone, {16, 16}]
Out[20]= -Image-

In[21]:= ImageResize[noise, {16, 24}]
Out[21]= -Image-

In[22]:= ImageResize[ramp, {8, 16}]
Out[22]= -Image-

In[23]:= ImageChannels[ImageResize[rgb, {8, 8}]]
Out[23]= 3

In[24]:= ImageDimensions[ImageResize[zone, {20, 10}]]
Out[24]= {20, 10}
```

### Applications (4)

```mathematica
In[25]:= disk = Image[Table[N[Boole[(i - 8.5)^2 + (j - 8.5)^2 <= 25]], {i, 1, 16}, {j, 1, 16}], "Real"];

In[26]:= zone = Image[Table[N[(1 + Cos[((i - 16)^2 + (j - 16)^2)/40.])/2], {i, 1, 32}, {j, 1, 32}], "Real"];

In[27]:= Binarize[ImageResize[zone, {16, 16}]]
Out[27]= -Image-

In[28]:= EdgeDetect[ImageResize[disk, {12, 12}]]
Out[28]= -Image-
```

### Properties & Relations (5)

```mathematica
In[29]:= chk = Image[Table[If[Mod[Quotient[i - 1, 2] + Quotient[j - 1, 2], 2] == 0, 0., 1.], {i, 1, 16}, {j, 1, 16}], "Real"];

In[30]:= rgb = Image[Table[{N[i/16], N[j/16], 0.5}, {i, 1, 16}, {j, 1, 16}], "Real"];

In[31]:= ImageDimensions[ImageResize[chk, {8, 8}]] === {8, 8}
Out[31]= True

In[32]:= ImageChannels[ImageResize[rgb, {4, 4}]] === ImageChannels[rgb]
Out[32]= True

In[33]:= ImageDimensions[ImageResize[chk, {16, 16}]] === ImageDimensions[chk]
Out[33]= True
```

### Neat Examples (2)

```mathematica
In[34]:= zone = Image[Table[N[(1 + Cos[((i - 16)^2 + (j - 16)^2)/40.])/2], {i, 1, 32}, {j, 1, 32}], "Real"];

In[35]:= ImageResize[zone, {24, 8}]
Out[35]= -Image-
```

## Implementation notes

**Attributes:** `Protected`.

## References

**See also:** [NDArray](../../linear-algebra/NDArray/), [ImageData](../../image-processing/ImageData/), [List](../../other-advanced/List/), [ImageConvolve](../../image-processing/ImageConvolve/), [GaussianFilter](../../image-processing/GaussianFilter/), [Table](../../lists-and-iteration/Table/)

- Source: [`src/imagegeom.c`](https://github.com/stblake/mathilda/blob/main/src/imagegeom.c)
- Specification: [`docs/spec/builtins/image-processing.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/image-processing.md)
- Tests: [`tests/test_image.c`](https://github.com/stblake/mathilda/blob/main/tests/test_image.c)
