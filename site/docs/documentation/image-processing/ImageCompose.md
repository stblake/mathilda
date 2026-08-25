# ImageCompose

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`ImageCompose[base, over] alpha-composites over onto base, centred, keeping base's size and clipping whatever falls outside. ImageCompose[base, over, {x, y}] centres the overlay at {x, y} in image coordinates -- x from the left, y from the BOTTOM. ImageCompose[base, {over, a}] scales the overlay's opacity by a. A grey image composed with a colour one produces colour: grey means the same value in every channel, so it is replicated rather than zero-padded.`**

## Examples (15)

Every input below was run against the current Mathilda build and its output recorded.

### Basic Examples (7)

```mathematica
In[1]:= a = Image[Table[N[(i + j)/32], {i, 1, 16}, {j, 1, 16}], "Real"];

In[2]:= red = Image[Table[{1., 0., 0.}, {i, 1, 6}, {j, 1, 6}], "Real"];

In[3]:= ImageCompose[a, red]
Out[3]= -Image-

In[4]:= ImageDimensions[ImageCompose[a, red]]
Out[4]= {16, 16}

In[5]:= ImageChannels[ImageCompose[a, red]]
Out[5]= 3

In[6]:= ImageCompose[a, red, {4, 4}]
Out[6]= -Image-

In[7]:= ImageCompose[a, {red, 0.4}]
Out[7]= -Image-
```

### Applications (3)

```mathematica
In[8]:= a = Image[Table[N[(i + j)/32], {i, 1, 32}, {j, 1, 32}], "Real"];
```

An edge map laid over the image it came from

```mathematica
In[9]:= ImageCompose[a, {EdgeDetect[a], 0.6}]
Out[9]= -Image-
```

A blurred copy blended halfway: the classic soft-focus composite

```mathematica
In[10]:= ImageCompose[a, {GaussianFilter[a, 3], 0.5}]
Out[10]= -Image-
```

### Properties & Relations (5)

```mathematica
In[11]:= a = Image[Table[N[(i + j)/32], {i, 1, 16}, {j, 1, 16}], "Real"];

In[12]:= red = Image[Table[{1., 0., 0.}, {i, 1, 6}, {j, 1, 6}], "Real"];
```

The size is the base's, whichever way round the two are given

```mathematica
In[13]:= {ImageDimensions[ImageCompose[a, red]], ImageDimensions[ImageCompose[red, a]]}
Out[13]= {{16, 16}, {6, 6}}
```

Outside the overlay, the grey base is replicated across all three channels

```mathematica
In[14]:= Module[{d = ImageData[ImageCompose[a, red]]}, d[[1, 1, 1]] === d[[1, 1, 2]] && d[[1, 1, 2]] === d[[1, 1, 3]]]
Out[14]= True
```

At zero opacity the overlay contributes nothing, even where it covers

```mathematica
In[15]:= Module[{d = ImageData[ImageCompose[a, {red, 0.}]]}, d[[8, 8, 1]] === d[[8, 8, 2]]]
Out[15]= True
```

## Implementation notes

- `Protected`.
- The result keeps the **base's** size and clips whatever falls outside: composition is "draw on
  this", not "make something bigger".
- A grey image composed with a colour one produces colour. Grey means the same value in every
  channel, so it is **replicated**, never zero-padded — padding would turn a grey pixel red.
- The result carries alpha only if the base did: compositing onto an opaque image gives an opaque
  image.

**Attributes:** `Protected`.

## References

- Source: [`src/imagecompose.c`](https://github.com/stblake/mathilda/blob/main/src/imagecompose.c)
- Specification: [`docs/spec/builtins/image-processing.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/image-processing.md)
- Tests: [`tests/test_image.c`](https://github.com/stblake/mathilda/blob/main/tests/test_image.c)
