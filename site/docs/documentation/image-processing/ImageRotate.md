# ImageRotate

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`ImageRotate[image] rotates a quarter turn counterclockwise; ImageRotate[image, angle] rotates by angle in radians (use n Degree for degrees). A multiple of a right angle takes an EXACT index-permutation path -- every pixel lands on another pixel's position, nothing is interpolated, and four quarter turns are exactly the identity. An odd number of quarter turns swaps the dimensions. Any other angle interpolates bilinearly, sampling the source per destination pixel (inverse mapping, so every output is filled exactly once; forward mapping leaves holes wherever the rotation stretches). Area rotated in from outside reads as 0 rather than the replicated edge, because that area was never photographed and smearing the border across it would invent content.`**

## Examples (39)

Every input below was run against the current Mathilda build and its output recorded.

### Basic Examples (8)

```mathematica
In[1]:= ImageData[ImageRotate[Image[{{1., 2.}, {3., 4.}}]]]
Out[1]= {{3.0, 1.0}, {4.0, 2.0}}

In[2]:= Module[{img = Image[{{1., 2.}, {3., 4.}}]}, Nest[ImageRotate, img, 4] === img]
Out[2]= True

In[3]:= chk = Image[Table[If[Mod[Quotient[i - 1, 2] + Quotient[j - 1, 2], 2] == 0, 0., 1.], {i, 1, 16}, {j, 1, 16}], "Real"];

In[4]:= disk = Image[Table[N[Boole[(i - 8.5)^2 + (j - 8.5)^2 <= 25]], {i, 1, 16}, {j, 1, 16}], "Real"];

In[5]:= ImageRotate[chk]
Out[5]= -Image-

In[6]:= ImageRotate[disk]
Out[6]= -Image-

In[7]:= ImageDimensions[ImageRotate[chk]]
Out[7]= {16, 16}

In[8]:= ImageRotate[chk, 0.4]
Out[8]= -Image-
```

### Scope (19)

```mathematica
In[9]:= disk = Image[Table[N[Boole[(i - 8.5)^2 + (j - 8.5)^2 <= 25]], {i, 1, 16}, {j, 1, 16}], "Real"];

In[10]:= ramp = Image[Table[N[(j - 1)/15], {i, 1, 16}, {j, 1, 16}], "Real"];

In[11]:= zone = Image[Table[N[(1 + Cos[((i - 16)^2 + (j - 16)^2)/40.])/2], {i, 1, 32}, {j, 1, 32}], "Real"];

In[12]:= noise = Image[Table[N[Mod[i*37 + j*17, 101]]/101, {i, 1, 32}, {j, 1, 32}], "Real"];

In[13]:= rgb = Image[Table[{N[i/16], N[j/16], 0.5}, {i, 1, 16}, {j, 1, 16}], "Real"];

In[14]:= sky = Image[Table[{N[0.15 + 0.7 (16 - i)/16], N[0.35 + 0.45 (16 - i)/16], N[0.85 - 0.35 (16 - i)/16]}, {i, 1, 16}, {j, 1, 24}], "Real"];

In[15]:= bit = Image[Table[Boole[Mod[i + j, 2] == 0], {i, 1, 8}, {j, 1, 8}]];

In[16]:= byte = Image[Table[Mod[i*13 + j*7, 256], {i, 1, 16}, {j, 1, 16}]];

In[17]:= ImageRotate[rgb]
Out[17]= -Image-

In[18]:= ImageRotate[sky]
Out[18]= -Image-

In[19]:= ImageRotate[bit]
Out[19]= -Image-

In[20]:= ImageRotate[byte]
Out[20]= -Image-

In[21]:= ImageRotate[zone]
Out[21]= -Image-

In[22]:= ImageRotate[ramp]
Out[22]= -Image-

In[23]:= ImageRotate[noise, 0.8]
Out[23]= -Image-

In[24]:= ImageRotate[disk, 1.2]
Out[24]= -Image-

In[25]:= ImageRotate[zone, 0.3]
Out[25]= -Image-

In[26]:= ImageChannels[ImageRotate[rgb]]
Out[26]= 3

In[27]:= ImageDimensions[ImageRotate[sky]]
Out[27]= {16, 24}
```

### Applications (4)

```mathematica
In[28]:= chk = Image[Table[If[Mod[Quotient[i - 1, 2] + Quotient[j - 1, 2], 2] == 0, 0., 1.], {i, 1, 16}, {j, 1, 16}], "Real"];

In[29]:= zone = Image[Table[N[(1 + Cos[((i - 16)^2 + (j - 16)^2)/40.])/2], {i, 1, 32}, {j, 1, 32}], "Real"];

In[30]:= EdgeDetect[ImageRotate[chk]]
Out[30]= -Image-

In[31]:= Binarize[ImageRotate[zone, 0.5]]
Out[31]= -Image-
```

### Properties & Relations (6)

```mathematica
In[32]:= chk = Image[Table[If[Mod[Quotient[i - 1, 2] + Quotient[j - 1, 2], 2] == 0, 0., 1.], {i, 1, 16}, {j, 1, 16}], "Real"];

In[33]:= disk = Image[Table[N[Boole[(i - 8.5)^2 + (j - 8.5)^2 <= 25]], {i, 1, 16}, {j, 1, 16}], "Real"];

In[34]:= rgb = Image[Table[{N[i/16], N[j/16], 0.5}, {i, 1, 16}, {j, 1, 16}], "Real"];

In[35]:= ImageData[ImageRotate[ImageRotate[ImageRotate[ImageRotate[chk]]]]] === ImageData[chk]
Out[35]= True

In[36]:= ImageChannels[ImageRotate[rgb]] === ImageChannels[rgb]
Out[36]= True

In[37]:= ImageData[ImageRotate[disk, 0.]] === ImageData[disk]
Out[37]= True
```

### Neat Examples (2)

```mathematica
In[38]:= zone = Image[Table[N[(1 + Cos[((i - 16)^2 + (j - 16)^2)/40.])/2], {i, 1, 32}, {j, 1, 32}], "Real"];

In[39]:= ImageRotate[zone, 0.7]
Out[39]= -Image-
```

## Implementation notes

**Attributes:** `Protected`.

## References

- Source: [`src/imagegeom.c`](https://github.com/stblake/mathilda/blob/main/src/imagegeom.c)
- Specification: [`docs/spec/builtins/image-processing.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/image-processing.md)
- Tests: [`tests/test_image.c`](https://github.com/stblake/mathilda/blob/main/tests/test_image.c)
