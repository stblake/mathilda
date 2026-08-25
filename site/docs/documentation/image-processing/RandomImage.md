# RandomImage

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`RandomImage[] gives a 150x150 grey image of uniform noise on [0, 1]. RandomImage[max] scales the range to [0, max]; RandomImage[max, {w, h}] sets the size, and a single n means {n, n}. ColorSpace -> "RGB" gives three independent channels. Samples are drawn from the same stream as RandomReal, so SeedRandom makes the result reproducible.`**

## Examples (20)

Every input below was run against the current Mathilda build and its output recorded.

### Basic Examples (5)

```mathematica
In[1]:= SeedRandom[42];

In[2]:= RandomImage[1, {32, 32}]
Out[2]= -Image-

In[3]:= ImageDimensions[RandomImage[]]
Out[3]= {150, 150}

In[4]:= ImageType[RandomImage[1, {8, 8}]]
Out[4]= "Real"

In[5]:= ImageChannels[RandomImage[1, {8, 8}]]
Out[5]= 1
```

### Scope (6)

```mathematica
In[6]:= SeedRandom[7];

In[7]:= ImageDimensions[RandomImage[1, {64, 16}]]
Out[7]= {64, 16}

In[8]:= ImageDimensions[RandomImage[1, 24]]
Out[8]= {24, 24}

In[9]:= RandomImage[1, {24, 24}, ColorSpace -> "RGB"]
Out[9]= -Image-

In[10]:= ImageChannels[RandomImage[1, {8, 8}, ColorSpace -> "RGB"]]
Out[10]= 3

In[11]:= Max[Flatten[ImageData[RandomImage[255, {16, 16}]]]] > 200
Out[11]= True
```

### Applications (5)

```mathematica
In[12]:= SeedRandom[3];
```

Noise is what shows a smoothing filter doing anything at all

```mathematica
In[13]:= GaussianFilter[RandomImage[1, {48, 48}], 2]
Out[13]= -Image-
```

A median filter removes salt-and-pepper noise a mean filter would only spread

```mathematica
In[14]:= MedianFilter[RandomImage[1, {48, 48}], 2]
Out[14]= -Image-
```

```mathematica
In[15]:= Binarize[RandomImage[1, {32, 32}]]
Out[15]= -Image-
```

And it is the honest input for a timing comparison, having no structure to exploit

```mathematica
In[16]:= ImageDimensions[Dilation[RandomImage[1, {64, 64}], 3]]
Out[16]= {64, 64}
```

### Properties & Relations (4)

The same seed gives the same image

```mathematica
In[17]:= Module[{a, b}, SeedRandom[7]; a = ImageData[RandomImage[1, {4, 4}]]; SeedRandom[7]; b = ImageData[RandomImage[1, {4, 4}]]; a === b]
Out[17]= True
```

And a different seed does not

```mathematica
In[18]:= Module[{a, b}, SeedRandom[7]; a = ImageData[RandomImage[1, {4, 4}]]; SeedRandom[8]; b = ImageData[RandomImage[1, {4, 4}]]; a =!= b]
Out[18]= True
```

The result is packed, like every other image-returning head

```mathematica
In[19]:= Head[Part[RandomImage[1, {16, 16}], 1]]
Out[19]= NDArray
```

An unsupported colour space declines

```mathematica
In[20]:= Head[RandomImage[1, {4, 4}, ColorSpace -> "CMYK"]]
Out[20]= RandomImage
```

## Algorithm

imageio.c -- Import and Export for raster image files.

Until this landed, every image in the system had to be typed out as an array of numbers, which makes the whole subsystem a demonstration rather than a tool: a filter is judged on photographs, and a synthetic checkerboard cannot show what a bilateral filter does that a Gaussian does not.

WHY A VENDORED DECODER. JPEG decoding is a baseline-Huffman-plus-IDCT project of its own and PNG needs an inflate, so the choice is between vendoring or making libpng and libjpeg hard build requirements. Two dependency-free public-domain headers cost less than either, and -- unlike a system library -- they cannot be missing at a user's site, which for an `Import` is the whole point. The headers are included HERE AND NOWHERE ELSE so that this is the only object file carrying third-party code.

WHAT A SAMPLE MEANS. A decoded 8-bit sample is scaled by 1/255 into the unit interval, because that is what the rest of the subsystem means by a brightness (see `image_load`) and the type a filter answers with is always "Real". So `Import` produces a "Real" image, not a "Byte" one: an image whose stored range depended on the file's bit depth would make every downstream kernel's scale depend on it too.

## Implementation notes

- `Protected`.
- Samples are drawn from the **same stream as `RandomReal`**, so `SeedRandom` makes a
  random image reproducible. A private generator would have made this the one random
  builtin that ignores the seed.
- The range is scaled, not clamped: a `"Real"` image may hold values above 1, and
  clamping belongs in `Export`, where 8 bits actually run out.
- Noise is the input a filter is most often judged on — a smoothing radius means nothing
  on a checkerboard and everything on a noise field.
- An unsupported colour space declines rather than silently returning grey.

**Attributes:** `Protected`.

## References

**See also:** [RandomReal](../../random-number-generation/RandomReal/), [SeedRandom](../../random-number-generation/SeedRandom/), [Export](../../file-io/Export/)

- Source: [`src/imageio.c`](https://github.com/stblake/mathilda/blob/main/src/imageio.c)
- Specification: [`docs/spec/builtins/image-processing.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/image-processing.md)
- Tests: [`tests/test_image.c`](https://github.com/stblake/mathilda/blob/main/tests/test_image.c)
