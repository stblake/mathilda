# Image3D

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Image3D[data] is a volumetric image, normalising to Image3D[data, type]. The data is a depth x height x width array of voxels, or depth x height x width x channels for colour, so it is indexed data[[z, y, x]] with slices outermost. ImageDimensions reports {width, height, depth} -- FULLY REVERSED from that order, which is Mathematica's convention. Type inference and the accessors match Image: ImageQ is False for a volume (use Image3DQ), while ImageDimensions, ImageChannels, ImageType and ImageData all accept either rank.`**

## Examples (67)

Every input below was run against the current Mathilda build and its output recorded.

### Basic Examples (14)

```mathematica
In[1]:= SeedRandom[5];
```

Uniform noise fills the volume to its boundary, so the block shows what it holds

```mathematica
In[2]:= Image3D[RandomReal[1, {24, 24, 24}], "Real"]
Out[2]= -Image-
```

Three channels, one per axis: the faces of the block are three colour gradients

```mathematica
In[3]:= Image3D[Table[{N[x/16], N[y/16], N[z/16]}, {z, 1, 16}, {y, 1, 16}, {x, 1, 16}], "Real"]
Out[3]= -Image-
```

```mathematica
In[4]:= ball = Image3D[Table[N[Boole[(x - 16)^2 + (y - 16)^2 + (z - 16)^2 <= 100]], {z, 1, 32}, {y, 1, 32}, {x, 1, 32}], "Real"];

In[5]:= shell = Image3D[Table[N[Mod[Floor[Sqrt[(x - 16.)^2 + (y - 16.)^2 + (z - 16.)^2]/3], 2]], {z, 1, 32}, {y, 1, 32}, {x, 1, 32}], "Real"];

In[6]:= cube = Image3D[Table[If[Mod[Quotient[x - 1, 4] + Quotient[y - 1, 4] + Quotient[z - 1, 4], 2] == 0, 0., 1.], {z, 1, 32}, {y, 1, 32}, {x, 1, 32}], "Real"];

In[7]:= ball
Out[7]= -Image-

In[8]:= ImageDimensions[ball]
Out[8]= {32, 32, 32}

In[9]:= Dimensions[ImageData[ball]]
Out[9]= {32, 32, 32}

In[10]:= ImageChannels[ball]
Out[10]= 1

In[11]:= ImageType[ball]
Out[11]= "Real"

In[12]:= Image3DQ[ball]
Out[12]= True

In[13]:= shell
Out[13]= -Image-

In[14]:= cube
Out[14]= -Image-
```

### Scope (28)

```mathematica
In[15]:= colvol = Image3D[Table[{N[x/16], N[y/16], N[z/16]}, {z, 1, 16}, {y, 1, 16}, {x, 1, 16}], "Real"];

In[16]:= bitvol = Image3D[Table[Boole[Mod[x + y + z, 2] == 0], {z, 1, 8}, {y, 1, 8}, {x, 1, 8}]];

In[17]:= rampz = Image3D[Table[N[(z - 1)/31], {z, 1, 32}, {y, 1, 24}, {x, 1, 24}], "Real"];

In[18]:= ball = Image3D[Table[N[Boole[(x - 16)^2 + (y - 16)^2 + (z - 16)^2 <= 100]], {z, 1, 32}, {y, 1, 32}, {x, 1, 32}], "Real"];

In[19]:= shell = Image3D[Table[N[Mod[Floor[Sqrt[(x - 16.)^2 + (y - 16.)^2 + (z - 16.)^2]/3], 2]], {z, 1, 32}, {y, 1, 32}, {x, 1, 32}], "Real"];

In[20]:= cube = Image3D[Table[If[Mod[Quotient[x - 1, 4] + Quotient[y - 1, 4] + Quotient[z - 1, 4], 2] == 0, 0., 1.], {z, 1, 32}, {y, 1, 32}, {x, 1, 32}], "Real"];

In[21]:= colvol
Out[21]= -Image-

In[22]:= ImageChannels[colvol]
Out[22]= 3

In[23]:= bitvol
Out[23]= -Image-

In[24]:= ImageType[bitvol]
Out[24]= "Bit"

In[25]:= rampz
Out[25]= -Image-

In[26]:= ImageDimensions[rampz]
Out[26]= {24, 24, 32}

In[27]:= GaussianFilter[ball, 2]
Out[27]= -Image-

In[28]:= MeanFilter[ball, 1]
Out[28]= -Image-

In[29]:= MedianFilter[ball, 1]
Out[29]= -Image-

In[30]:= Dilation[ball, 2]
Out[30]= -Image-

In[31]:= Erosion[ball, 2]
Out[31]= -Image-

In[32]:= Opening[ball, 2]
Out[32]= -Image-

In[33]:= Closing[ball, 2]
Out[33]= -Image-

In[34]:= Binarize[shell]
Out[34]= -Image-

In[35]:= LocalAdaptiveBinarize[shell, 3]
Out[35]= -Image-

In[36]:= DistanceTransform[ball]
Out[36]= -Image-

In[37]:= ImagePad[cube, 2]
Out[37]= -Image-

In[38]:= ImageReflect[cube, Front]
Out[38]= -Image-

In[39]:= ColorConvert[colvol, "Grayscale"]
Out[39]= -Image-

In[40]:= DerivativeFilter[cube, {0, 0, 1}]
Out[40]= -Image-

In[41]:= CornerFilter[ball, 2]
Out[41]= -Image-

In[42]:= ImageCorners[ball]
Out[42]= {{9, 16, 9}, {16, 9, 9}, {16, 23, 9}, {23, 16, 9}, {9, 16, 23}, {9, 23, 16}, {16, 9, 23}, {16, 23, 23}, {23, 16, 23}, {23, 23, 16}, {9, 9, 16}, {23, 9, 16}, {7, 12, 16}, {25, 12, 16}, {12, 16, 7}, {12, 16, 25}, {16, 12, 7}, {16, 12, 25}, {16, 25, 20}, {20, 16, 7}, {20, 16, 25}, {7, 20, 16}, {12, 7, 16}, {12, 25, 16}, {16, 7, 12}, {16, 20, 7}, {16, 20, 25}, {16, 25, 12}, {20, 7, 16}, {20, 25, 16}, {25, 20, 16}, {7, 16, 12}, {7, 16, 20}, {16, 7, 20}, {25, 16, 12}, {25, 16, 20}, {7, 19, 13}, {13, 13, 7}, {13, 13, 25}, {13, 25, 13}, {19, 13, 7}, {19, 13, 25}, {19, 25, 13}, {25, 19, 13}, {7, 13, 13}, {25, 13, 13}, {7, 13, 19}, {7, 19, 19}, {25, 13, 19}, {25, 19, 19}, {13, 7, 19}, {13, 19, 7}, {13, 19, 25}, {13, 25, 19}, {19, 7, 19}, {19, 19, 7}, {19, 19, 25}, {19, 25, 19}, {13, 7, 13}, {19, 7, 13}, {6, 16, 16}, {16, 6, 16}, {16, 16, 6}, {16, 16, 26}, {16, 26, 16}, {26, 16, 16}, {10, 10, 10}, {10, 22, 10}, {22, 10, 10}, {22, 22, 10}, {10, 10, 22}, {10, 22, 22}, {22, 10, 22}, {22, 22, 22}}
```

### Applications (8)

```mathematica
In[43]:= ball = Image3D[Table[N[Boole[(x - 16)^2 + (y - 16)^2 + (z - 16)^2 <= 100]], {z, 1, 32}, {y, 1, 32}, {x, 1, 32}], "Real"];

In[44]:= shell = Image3D[Table[N[Mod[Floor[Sqrt[(x - 16.)^2 + (y - 16.)^2 + (z - 16.)^2]/3], 2]], {z, 1, 32}, {y, 1, 32}, {x, 1, 32}], "Real"];

In[45]:= Image[Part[ImageData[ball], 16]]
Out[45]= -Image-

In[46]:= Image[Part[ImageData[GaussianFilter[ball, 2]], 16]]
Out[46]= -Image-

In[47]:= Image[Join[Part[ImageData[shell], 8], Part[ImageData[shell], 16], Part[ImageData[shell], 24], 2]]
Out[47]= -Image-

In[48]:= Total[Flatten[ImageData[ball]]]
Out[48]= 4169.0

In[49]:= Total[Flatten[ImageData[ball] - ImageData[Erosion[ball, 1]]]]
Out[49]= 1640.0

In[50]:= Max[Flatten[ImageData[DistanceTransform[ball]]]]
Out[50]= 10.0499
```

### Properties & Relations (12)

```mathematica
In[51]:= ball = Image3D[Table[N[Boole[(x - 16)^2 + (y - 16)^2 + (z - 16)^2 <= 100]], {z, 1, 32}, {y, 1, 32}, {x, 1, 32}], "Real"];

In[52]:= cube = Image3D[Table[If[Mod[Quotient[x - 1, 4] + Quotient[y - 1, 4] + Quotient[z - 1, 4], 2] == 0, 0., 1.], {z, 1, 32}, {y, 1, 32}, {x, 1, 32}], "Real"];

In[53]:= colvol = Image3D[Table[{N[x/16], N[y/16], N[z/16]}, {z, 1, 16}, {y, 1, 16}, {x, 1, 16}], "Real"];

In[54]:= ImageDimensions[ball] === {32, 32, 32}
Out[54]= True

In[55]:= Reverse[ImageDimensions[ball]] === Dimensions[ImageData[ball]]
Out[55]= True

In[56]:= ImageData[ImageReflect[ImageReflect[cube, Front], Front]] === ImageData[cube]
Out[56]= True

In[57]:= ImageData[Dilation[ball, 0]] === ImageData[ball]
Out[57]= True

In[58]:= ImageData[Opening[Opening[ball, 2], 2]] === ImageData[Opening[ball, 2]]
Out[58]= True

In[59]:= ImageDimensions[GaussianFilter[ball, 2]] === ImageDimensions[ball]
Out[59]= True

In[60]:= ImageChannels[ColorConvert[colvol, "Grayscale"]] === 1
Out[60]= True

In[61]:= Image3DQ[ball] && Not[ImageQ[ball]]
Out[61]= True

In[62]:= Max[Flatten[ImageData[CornerFilter[Image3D[Table[If[x <= 16, 0., 1.], {z, 1, 16}, {y, 1, 16}, {x, 1, 16}], "Real"], 2]]]]
Out[62]= 0.0
```

### Neat Examples (5)

```mathematica
In[63]:= shell = Image3D[Table[N[Mod[Floor[Sqrt[(x - 16.)^2 + (y - 16.)^2 + (z - 16.)^2]/3], 2]], {z, 1, 32}, {y, 1, 32}, {x, 1, 32}], "Real"];

In[64]:= cube = Image3D[Table[If[Mod[Quotient[x - 1, 4] + Quotient[y - 1, 4] + Quotient[z - 1, 4], 2] == 0, 0., 1.], {z, 1, 32}, {y, 1, 32}, {x, 1, 32}], "Real"];

In[65]:= shell
Out[65]= -Image-

In[66]:= Image[Part[ImageData[shell], 16]]
Out[66]= -Image-

In[67]:= Image[Part[ImageData[Dilation[cube, 1]], 16]]
Out[67]= -Image-
```

## Implementation notes

**Attributes:** `Protected`.

## References

**See also:** [ImageDimensions](../../image-processing/ImageDimensions/), [ImageQ](../../image-processing/ImageQ/), [Image3DQ](../../image-processing/Image3DQ/), [Image](../../image-processing/Image/), [ImageData](../../image-processing/ImageData/), [ImageChannels](../../image-processing/ImageChannels/), [ImageType](../../image-processing/ImageType/)

- Source: [`src/image.c`](https://github.com/stblake/mathilda/blob/main/src/image.c)
- Specification: [`docs/spec/builtins/image-processing.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/image-processing.md)
- Tests: [`tests/test_image.c`](https://github.com/stblake/mathilda/blob/main/tests/test_image.c)
