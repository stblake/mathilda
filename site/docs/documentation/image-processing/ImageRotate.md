# ImageRotate

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`ImageRotate[image] rotates a quarter turn counterclockwise; ImageRotate[image, angle] rotates by angle in radians counterclockwise (use n Degree for degrees; a negative angle turns clockwise); ImageRotate[image, side] turns the top of the image to face side (Left, Right, Top or Bottom) and ImageRotate[image, side1 -> side2] turns side1 onto side2. A multiple of a right angle takes an EXACT index-permutation path -- every pixel lands on another pixel's position, nothing is interpolated, and four quarter turns are exactly the identity. An odd number of quarter turns swaps the dimensions. Any other angle interpolates bilinearly, sampling the source per destination pixel (inverse mapping, so every output is filled exactly once; forward mapping leaves holes wherever the rotation stretches). Area rotated in from outside reads as 0 rather than the replicated edge, because that area was never photographed and smearing the border across it would invent content. Unlike Mathematica, a free angle keeps the input's dimensions (Mathematica's size Full) rather than enlarging to enclose the rotated image.`**

## Examples (42)

Every input below was run against the current Mathilda build and its output recorded.

### Basic Examples (8)

```mathematica
In[1]:= ImageData[ImageRotate[Image[{{1., 2.}, {3., 4.}}]]]
Out[1]= {{2.0, 4.0}, {1.0, 3.0}}

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

### Scope (22)

```mathematica
In[9]:= ImageData[ImageRotate[Image[{{1., 2., 3.}, {4., 5., 6.}}]]]
Out[9]= {{3.0, 6.0}, {2.0, 5.0}, {1.0, 4.0}}

In[10]:= ImageData[ImageRotate[Image[{{1., 2., 3.}, {4., 5., 6.}}], -Pi/2]]
Out[10]= {{4.0, 1.0}, {5.0, 2.0}, {6.0, 3.0}}

In[11]:= ImageData[ImageRotate[Image[{{1., 2., 3.}, {4., 5., 6.}}], Right]]
Out[11]= {{4.0, 1.0}, {5.0, 2.0}, {6.0, 3.0}}

In[12]:= disk = Image[Table[N[Boole[(i - 8.5)^2 + (j - 8.5)^2 <= 25]], {i, 1, 16}, {j, 1, 16}], "Real"];

In[13]:= ramp = Image[Table[N[(j - 1)/15], {i, 1, 16}, {j, 1, 16}], "Real"];

In[14]:= zone = Image[Table[N[(1 + Cos[((i - 16)^2 + (j - 16)^2)/40.])/2], {i, 1, 32}, {j, 1, 32}], "Real"];

In[15]:= noise = Image[Table[N[Mod[i*37 + j*17, 101]]/101, {i, 1, 32}, {j, 1, 32}], "Real"];

In[16]:= rgb = Image[Table[{N[i/16], N[j/16], 0.5}, {i, 1, 16}, {j, 1, 16}], "Real"];

In[17]:= sky = Image[Table[{N[0.15 + 0.7 (16 - i)/16], N[0.35 + 0.45 (16 - i)/16], N[0.85 - 0.35 (16 - i)/16]}, {i, 1, 16}, {j, 1, 24}], "Real"];

In[18]:= bit = Image[Table[Boole[Mod[i + j, 2] == 0], {i, 1, 8}, {j, 1, 8}]];

In[19]:= byte = Image[Table[Mod[i*13 + j*7, 256], {i, 1, 16}, {j, 1, 16}]];

In[20]:= ImageRotate[rgb]
Out[20]= -Image-

In[21]:= ImageRotate[sky]
Out[21]= -Image-

In[22]:= ImageRotate[bit]
Out[22]= -Image-

In[23]:= ImageRotate[byte]
Out[23]= -Image-

In[24]:= ImageRotate[zone]
Out[24]= -Image-

In[25]:= ImageRotate[ramp]
Out[25]= -Image-

In[26]:= ImageRotate[noise, 0.8]
Out[26]= -Image-

In[27]:= ImageRotate[disk, 1.2]
Out[27]= -Image-

In[28]:= ImageRotate[zone, 0.3]
Out[28]= -Image-

In[29]:= ImageChannels[ImageRotate[rgb]]
Out[29]= 3

In[30]:= ImageDimensions[ImageRotate[sky]]
Out[30]= {16, 24}
```

### Applications (4)

```mathematica
In[31]:= chk = Image[Table[If[Mod[Quotient[i - 1, 2] + Quotient[j - 1, 2], 2] == 0, 0., 1.], {i, 1, 16}, {j, 1, 16}], "Real"];

In[32]:= zone = Image[Table[N[(1 + Cos[((i - 16)^2 + (j - 16)^2)/40.])/2], {i, 1, 32}, {j, 1, 32}], "Real"];

In[33]:= EdgeDetect[ImageRotate[chk]]
Out[33]= -Image-

In[34]:= Binarize[ImageRotate[zone, 0.5]]
Out[34]= -Image-
```

### Properties & Relations (6)

```mathematica
In[35]:= chk = Image[Table[If[Mod[Quotient[i - 1, 2] + Quotient[j - 1, 2], 2] == 0, 0., 1.], {i, 1, 16}, {j, 1, 16}], "Real"];

In[36]:= disk = Image[Table[N[Boole[(i - 8.5)^2 + (j - 8.5)^2 <= 25]], {i, 1, 16}, {j, 1, 16}], "Real"];

In[37]:= rgb = Image[Table[{N[i/16], N[j/16], 0.5}, {i, 1, 16}, {j, 1, 16}], "Real"];

In[38]:= ImageData[ImageRotate[ImageRotate[ImageRotate[ImageRotate[chk]]]]] === ImageData[chk]
Out[38]= True

In[39]:= ImageChannels[ImageRotate[rgb]] === ImageChannels[rgb]
Out[39]= True

In[40]:= ImageData[ImageRotate[disk, 0.]] === ImageData[disk]
Out[40]= True
```

### Neat Examples (2)

```mathematica
In[41]:= zone = Image[Table[N[(1 + Cos[((i - 16)^2 + (j - 16)^2)/40.])/2], {i, 1, 32}, {j, 1, 32}], "Real"];

In[42]:= ImageRotate[zone, 0.7]
Out[42]= -Image-
```

## Implementation notes

**Attributes:** `Protected`.

## References

- Source: [`src/imagegeom.c`](https://github.com/stblake/mathilda/blob/main/src/imagegeom.c)
- Specification: [`docs/spec/builtins/image-processing.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/image-processing.md)
- Tests: [`tests/test_image.c`](https://github.com/stblake/mathilda/blob/main/tests/test_image.c)
