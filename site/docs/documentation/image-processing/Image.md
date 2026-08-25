# Image

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Image[data] is a raster image, normalising to the canonical Image[data, type]. The data is a rectangular height x width array of pixel values, or height x width x channels for a colour image, so it is indexed data[[y, x]] with rows running down the image -- note that ImageDimensions reports {width, height}, transposed relative to this. The type is inferred from the values: all-integer data in {0, 1} is "Bit", all-integer in 0..255 is "Byte", anything else is "Real". Image[data, type] states the type instead, and declines if the data does not fit it. Ragged data declines rather than being padded.`**

## Examples (38)

Every input below was run against the current Mathilda build and its output recorded.

### Basic Examples (6)

```mathematica
In[1]:= Image[{{0., 1.}, {1., 0.}}]
Out[1]= -Image-

In[2]:= Image[{{0., 0.25, 0.5, 0.75, 1.}}]
Out[2]= -Image-

In[3]:= ImageDimensions[Image[{{0., 1., 0.5}, {1., 0., 0.25}}]]
Out[3]= {3, 2}

In[4]:= ImageChannels[Image[{{0., 1.}, {1., 0.}}]]
Out[4]= 1

In[5]:= ImageType[Image[{{0., 1.}, {1., 0.}}]]
Out[5]= "Real"

In[6]:= ImageData[Image[{{0., 1.}, {1., 0.}}]]
Out[6]= {{0.0, 1.0}, {1.0, 0.0}}
```

### Scope (32)

The type is inferred from the values

```mathematica
In[7]:= ImageType[Image[{{0, 1}, {1, 0}}]]
Out[7]= "Bit"
```

```mathematica
In[8]:= ImageType[Image[{{0, 128}, {255, 7}}]]
Out[8]= "Byte"

In[9]:= ImageType[Image[{{0., 0.5}}]]
Out[9]= "Real"
```

A stated type must fit the data, or the call declines

```mathematica
In[10]:= Head[Image[{{0, 300}}, "Byte"]]
Out[10]= Image
```

```mathematica
In[11]:= Head[Image[{{0, 2}}, "Bit"]]
Out[11]= Image

In[12]:= ImageType[Image[{{0, 1}, {1, 0}}, "Byte"]]
Out[12]= "Byte"
```

Ragged data declines rather than being padded

```mathematica
In[13]:= Head[Image[{{1., 2.}, {3.}}]]
Out[13]= Image
```

```mathematica
In[14]:= Head[Image[{}]]
Out[14]= Image
```

Data is indexed [[y, x]] while ImageDimensions reports {width, height}

```mathematica
In[15]:= Dimensions[ImageData[Image[{{1., 2., 3.}, {4., 5., 6.}}]]]
Out[15]= {2, 3}
```

```mathematica
In[16]:= ImageDimensions[Image[{{1., 2., 3.}, {4., 5., 6.}}]]
Out[16]= {3, 2}

In[17]:= ImageData[Image[{{1., 2., 3.}, {4., 5., 6.}}]][[1, 3]]
Out[17]= 3.0
```

Three channels make a colour image

```mathematica
In[18]:= ImageChannels[Image[{{{1., 0., 0.}, {0., 1., 0.}}}]]
Out[18]= 3
```

```mathematica
In[19]:= Image[{{{1., 0., 0.}, {0., 1., 0.}, {0., 0., 1.}}}]
Out[19]= -Image-

In[20]:= ImageDimensions[Image[{{{1., 0., 0.}, {0., 1., 0.}}}]]
Out[20]= {2, 1}
```

Two or four channels carry alpha

```mathematica
In[21]:= ImageChannels[Image[{{{1., 0.5}, {0., 1.}}}]]
Out[21]= 2
```

```mathematica
In[22]:= ImageChannels[Image[{{{1., 0., 0., 0.5}}}]]
Out[22]= 4
```

ImageData reports STORED values with an explicit type

```mathematica
In[23]:= ImageData[Image[{{0, 255}}, "Byte"], "Byte"]
Out[23]= {{0, 255}}
```

```mathematica
In[24]:= ImageData[Image[{{0, 255}}, "Byte"]]
Out[24]= {{0.0, 1.0}}

In[25]:= ImageData[Image[{{1, 0}, {0, 1}}, "Bit"], "Bit"]
Out[25]= {{1, 0}, {0, 1}}
```

An already-canonical image is left alone, so evaluation reaches a fixed point

```mathematica
In[26]:= Image[Image[{{0., 1.}}]] === Image[{{0., 1.}}]
Out[26]= False
```

Pixels survive a round trip through ImageData exactly

```mathematica
In[27]:= Module[{d = {{0.1, 0.2}, {0.3, 0.4}}}, ImageData[Image[d]] === d]
Out[27]= True
```

```mathematica
In[28]:= Image[Table[N[i j]/9, {i, 3}, {j, 3}]]
Out[28]= -Image-

In[29]:= Image[Table[N[Mod[i + j, 2]], {i, 8}, {j, 8}]]
Out[29]= -Image-

In[30]:= Image[Table[N[(i - 1)/7], {i, 8}, {j, 8}]]
Out[30]= -Image-
```

A colour ramp

```mathematica
In[31]:= Image[Table[{N[(j - 1)/7], N[(i - 1)/7], 0.5}, {i, 8}, {j, 8}]]
Out[31]= -Image-
```

```mathematica
In[32]:= Image[Table[If[Mod[Quotient[i - 1, 2] + Quotient[j - 1, 2], 2] == 0, 0., 1.], {i, 8}, {j, 8}]]
Out[32]= -Image-

In[33]:= Image[Table[N[Boole[(i - 4.5)^2 + (j - 4.5)^2 <= 9]], {i, 8}, {j, 8}]]
Out[33]= -Image-

In[34]:= ImageQ[Image[{{0., 1.}}]]
Out[34]= True

In[35]:= ImageQ[{{0., 1.}}]
Out[35]= False

In[36]:= Image3DQ[Image[{{0., 1.}}]]
Out[36]= False
```

The storage is a packed buffer, not a tree of Expr nodes

```mathematica
In[37]:= Head[Part[Image[Table[N[i j]/64, {i, 8}, {j, 8}]], 1]]
Out[37]= List
```

```mathematica
In[38]:= Part[Image[{{0., 1.}}], 2]
Out[38]= "Real"
```

## Implementation notes

**Attributes:** `Protected`.

## References

**See also:** [ImageQ](../../image-processing/ImageQ/), [ImageDimensions](../../image-processing/ImageDimensions/)

- Source: [`src/image.c`](https://github.com/stblake/mathilda/blob/main/src/image.c)
- Specification: [`docs/spec/builtins/image-processing.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/image-processing.md)
- Tests: [`tests/test_image.c`](https://github.com/stblake/mathilda/blob/main/tests/test_image.c)
