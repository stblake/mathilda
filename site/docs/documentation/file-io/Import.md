# Import

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Import["file"] reads a raster image file (PNG, JPEG, BMP, GIF, TGA, PSD, HDR, PNM) and returns an Image. Import["file", "Image"] is the same. Samples are scaled by 1/255 into the unit interval, so the result is a "Real" image whatever the file's bit depth, and the file's channel count is preserved -- grey stays 1 channel, RGBA keeps its alpha. Gives $Failed for a missing or malformed file.`**

## Examples (21)

Every input below was run against the current Mathilda build and its output recorded.

### Basic Examples (5)

```mathematica
In[1]:= img = Image[Table[{N[i/24], N[j/32], N[Mod[i + j, 8]]/8}, {i, 1, 24}, {j, 1, 32}], "Real"];

In[2]:= Export["/tmp/mathilda_doc.png", img]
Out[2]= "/tmp/mathilda_doc.png"

In[3]:= Import["/tmp/mathilda_doc.png"]
Out[3]= -Image-

In[4]:= ImageDimensions[Import["/tmp/mathilda_doc.png"]]
Out[4]= {32, 24}

In[5]:= ImageType[Import["/tmp/mathilda_doc.png"]]
Out[5]= "Real"
```

### Scope (7)

```mathematica
In[6]:= img = Image[Table[{N[i/24], N[j/32], N[Mod[i + j, 8]]/8}, {i, 1, 24}, {j, 1, 32}], "Real"];

In[7]:= ImageChannels[Import[Export["/tmp/mathilda_doc_g.png", Image[Table[N[i/16], {i, 1, 16}, {j, 1, 16}], "Real"]]]]
Out[7]= 1

In[8]:= ImageChannels[Import[Export["/tmp/mathilda_doc_a.png", Image[Table[{0.2, 0.4, 0.6, 0.8}, {i, 1, 8}, {j, 1, 8}], "Real"]]]]
Out[8]= 4

In[9]:= ImageDimensions[Import[Export["/tmp/mathilda_doc.jpg", img]]]
Out[9]= {32, 24}

In[10]:= ImageDimensions[Import[Export["/tmp/mathilda_doc.bmp", img]]]
Out[10]= {32, 24}

In[11]:= ImageDimensions[Import[Export["/tmp/mathilda_doc.tga", img]]]
Out[11]= {32, 24}

In[12]:= Import["/tmp/mathilda_doc_missing.png"]
Out[12]= $Failed
```

### Applications (4)

```mathematica
In[13]:= img = Image[Table[{N[i/24], N[j/32], N[Mod[i + j, 8]]/8}, {i, 1, 24}, {j, 1, 32}], "Real"];
```

Filters compose with an imported image exactly as with a constructed one

```mathematica
In[14]:= ImageDimensions[GaussianFilter[Import["/tmp/mathilda_doc.png"], 2]]
Out[14]= {32, 24}
```

A pipeline written end to end: read, edge-detect, write

```mathematica
In[15]:= Export["/tmp/mathilda_doc_edges.png", EdgeDetect[Import["/tmp/mathilda_doc.png"]]]
Out[15]= "/tmp/mathilda_doc_edges.png"
```

```mathematica
In[16]:= ImageDimensions[Import["/tmp/mathilda_doc_edges.png"]]
Out[16]= {32, 24}
```

### Properties & Relations (5)

```mathematica
In[17]:= img = Image[Table[{N[i/24], N[j/32], N[Mod[i + j, 8]]/8}, {i, 1, 24}, {j, 1, 32}], "Real"];
```

PNG is lossless, so a round trip is exact to within half a quantisation level

```mathematica
In[18]:= Max[Abs[Flatten[ImageData[Import[Export["/tmp/mathilda_doc.png", img]]] - ImageData[img]]]] <= 1/510. + 1.*^-12
Out[18]= True
```

An imported image is packed, like every filter's result

```mathematica
In[19]:= Head[Part[Import["/tmp/mathilda_doc.png"], 1]]
Out[19]= NDArray
```

JPEG is lossy: the same round trip is bounded, not exact

```mathematica
In[20]:= 0 < Mean[Flatten[Abs[ImageData[Import[Export["/tmp/mathilda_doc.jpg", img]]] - ImageData[img]]]] < 0.1
Out[20]= True
```

A format nothing here claims stays unevaluated, which is not the same failure as a missing file

```mathematica
In[21]:= Head[Import["/tmp/mathilda_doc.xyz"]]
Out[21]= Import
```

## Algorithm

imageio.c -- Import and Export for raster image files.

Until this landed, every image in the system had to be typed out as an array of numbers, which makes the whole subsystem a demonstration rather than a tool: a filter is judged on photographs, and a synthetic checkerboard cannot show what a bilateral filter does that a Gaussian does not.

WHY A VENDORED DECODER. JPEG decoding is a baseline-Huffman-plus-IDCT project of its own and PNG needs an inflate, so the choice is between vendoring or making libpng and libjpeg hard build requirements. Two dependency-free public-domain headers cost less than either, and -- unlike a system library -- they cannot be missing at a user's site, which for an `Import` is the whole point. The headers are included HERE AND NOWHERE ELSE so that this is the only object file carrying third-party code.

WHAT A SAMPLE MEANS. A decoded 8-bit sample is scaled by 1/255 into the unit interval, because that is what the rest of the subsystem means by a brightness (see `image_load`) and the type a filter answers with is always "Real". So `Import` produces a "Real" image, not a "Byte" one: an image whose stored range depended on the file's bit depth would make every downstream kernel's scale depend on it too.

## Implementation notes

- `Protected`.
- Samples are scaled by `1/255` into the unit interval, so the result is a `"Real"`
  image whatever the file's bit depth. The stored range of an image fixes what every
  downstream kernel's arithmetic means (see `ImageData`), and a range that depended on
  the file would make a Gaussian's scale depend on it too.
- The file's channel count is preserved: a grey file stays 1-channel and an RGBA file
  keeps its alpha. Forcing 3 channels would invent two for the first and silently
  discard transparency from the second.
- The result is packed and canonical — the same representation a filter produces, so an
  imported photograph needs no special-casing downstream.
- `$Failed` for a missing or malformed file. A path whose format is not handled at all
  stays unevaluated instead, which keeps `Import` from appearing to implement every
  format in existence.
- Decoding is by the vendored `stb_image` (public domain), so no system image library is
  a build requirement.

**Attributes:** `Protected`.

## References

**See also:** [Image](../../image-processing/Image/), [ImageData](../../image-processing/ImageData/)

- Source: [`src/imageio.c`](https://github.com/stblake/mathilda/blob/main/src/imageio.c)
- Specification: [`docs/spec/builtins/file-io.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/file-io.md)
- Tests: [`tests/test_image.c`](https://github.com/stblake/mathilda/blob/main/tests/test_image.c)
