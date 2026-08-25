# Export

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Export["file", obj] writes obj to a file, choosing the format from the file extension; Export["file", obj, "FMT"] states it explicitly. An Image writes to a raster file (PNG, JPEG, BMP, TGA); its samples outside the unit interval are clamped, since 8-bit output has no room for them. A Graphics object (the result of Plot, ListPlot, Graphics, ...) writes to PDF, PNG, or JPEG: PDF is a resolution-independent vector file produced without any external library and works headless, while PNG and JPEG render through the graphics backend and so need graphics support compiled in and a display. A Graphics3D object (Plot3D, ParametricPlot3D, ...) exports to PNG or JPEG the same way; it has no vector-PDF form. Returns the file name, so Import[Export[f, img]] round-trips.`**

## Examples (15)

Every input below was run against the current Mathilda build and its output recorded.

### Basic Examples (4)

```mathematica
In[1]:= img = Image[Table[{N[i/24], N[j/32], N[Mod[i + j, 8]]/8}, {i, 1, 24}, {j, 1, 32}], "Real"];

In[2]:= Export["/tmp/mathilda_doc_e.png", img]
Out[2]= "/tmp/mathilda_doc_e.png"

In[3]:= FileExistsQ["/tmp/mathilda_doc_e.png"]
Out[3]= True

In[4]:= ImageDimensions[Import["/tmp/mathilda_doc_e.png"]]
Out[4]= {32, 24}
```

### Scope (5)

```mathematica
In[5]:= img = Image[Table[{N[i/24], N[j/32], N[Mod[i + j, 8]]/8}, {i, 1, 24}, {j, 1, 32}], "Real"];

In[6]:= Export["/tmp/mathilda_doc_e.jpg", img]
Out[6]= "/tmp/mathilda_doc_e.jpg"

In[7]:= Export["/tmp/mathilda_doc_e.bmp", img]
Out[7]= "/tmp/mathilda_doc_e.bmp"

In[8]:= Export["/tmp/mathilda_doc_e.tga", img]
Out[8]= "/tmp/mathilda_doc_e.tga"
```

A grey image writes a 1-channel file

```mathematica
In[9]:= ImageChannels[Import[Export["/tmp/mathilda_doc_eg.png", Image[Table[N[i/16], {i, 1, 16}, {j, 1, 16}], "Real"]]]]
Out[9]= 1
```

### Options (3)

```mathematica
In[10]:= img = Image[Table[{N[i/24], N[j/32], N[Mod[i + j, 8]]/8}, {i, 1, 24}, {j, 1, 32}], "Real"];
```

The format may be stated rather than inferred, which is the only way to write a file with no extension

```mathematica
In[11]:= Export["/tmp/mathilda_doc_noext", img, "PNG"]
Out[11]= "/tmp/mathilda_doc_noext"
```

```mathematica
In[12]:= ImageDimensions[Import["/tmp/mathilda_doc_noext", "Image"]]
Out[12]= {32, 24}
```

### Properties & Relations (3)

```mathematica
In[13]:= img = Image[Table[{N[i/24], N[j/32], N[Mod[i + j, 8]]/8}, {i, 1, 24}, {j, 1, 32}], "Real"];
```

Out-of-range samples clamp to the ends rather than wrapping

```mathematica
In[14]:= ImageData[Import[Export["/tmp/mathilda_doc_clamp.png", Image[{{2.0, -1.0}, {1.0, 0.0}}, "Real"]]]]
Out[14]= {{1.0, 0.0}, {1.0, 0.0}}
```

A volume is declined rather than silently reduced to a slice

```mathematica
In[15]:= Head[Export["/tmp/mathilda_doc_vol.png", Image3D[Table[0.5, {z, 1, 2}, {y, 1, 2}, {x, 1, 2}], "Real"]]]
Out[15]= Export
```

## Options & behaviour

### Graphics export

- **PDF** is a resolution-independent **vector** file written by a small built-in PDF
  emitter — no external library, no display, so it works headless and in every build. It
  walks the graphics primitives directly (`Line`, `Point`, `Polygon`, `Disk`/`Circle`,
  `Rectangle`, `Arrow`, `Text`) with the `RGBColor`/`GrayLevel`/`Hue`/`CMYKColor`,
  `Opacity`, `Thickness` and `PointSize` directives, and draws a framed set of axes with
  "nice" ticks and numeric labels. Text uses the PDF base-14 Helvetica, so no font is
  embedded. This is the recommended format for print and for the book.
- **PNG** and **JPEG** render through the graphics backend into an offscreen buffer, so the
  file is pixel-identical to the on-screen plot (the same axes, ticks, labels and text).
  They therefore need graphics support compiled in (`USE_GRAPHICS`) **and** a usable GUI
  session; with none (a headless box, `ssh`, cron) they return `$Failed` gracefully rather
  than crashing, while PDF still works. Resolution follows the `ImageSize` option: a width
  (default 800) with the height derived from the plot's `AspectRatio`, sized exactly as the
  on-screen window is, so an aspect-driven plot (`ArrayPlot`, `DensityPlot`, `ContourPlot`,
  …) fills its frame edge-to-edge instead of letterboxing inside a fixed canvas
  `ImageSize -> {w, h}` pins both dimensions (then `AspectRatio` shapes the data inside that
  box). The pixels are encoded by the vendored `stb_image_write`, so JPEG output does not
  depend on which formats the Raylib build happens to support.
- A `Graphics3D` object (`Plot3D`, `ParametricPlot3D`, `ComplexPlot3D`, ...) exports to
  **PNG or JPEG** through the 3D renderer, with the same graphics-support/display
  requirement; it has no vector-PDF form (PDF of a 3D scene returns `$Failed`).
- Samples outside the unit interval are **clamped**, not wrapped. An unsharp mask
  legitimately overshoots and 8-bit output has nowhere to put the overshoot; wrapping
  would turn a bright highlight black, which reads as a bug in the filter rather than in
  the writer. `NaN` clamps to 0.
- JPEG is written at quality 90 — a documented constant rather than a silent one. Use PNG
  when the bytes must survive.
- An `Image3D` is declined (the expression stays unevaluated): a volume has no single
  raster, and quietly writing its middle slice would misreport what was exported.
- Writing is by the vendored `stb_image_write` (public domain).

## Algorithm

imageio.c -- Import and Export for raster image files.

Until this landed, every image in the system had to be typed out as an array of numbers, which makes the whole subsystem a demonstration rather than a tool: a filter is judged on photographs, and a synthetic checkerboard cannot show what a bilateral filter does that a Gaussian does not.

WHY A VENDORED DECODER. JPEG decoding is a baseline-Huffman-plus-IDCT project of its own and PNG needs an inflate, so the choice is between vendoring or making libpng and libjpeg hard build requirements. Two dependency-free public-domain headers cost less than either, and -- unlike a system library -- they cannot be missing at a user's site, which for an `Import` is the whole point. The headers are included HERE AND NOWHERE ELSE so that this is the only object file carrying third-party code.

WHAT A SAMPLE MEANS. A decoded 8-bit sample is scaled by 1/255 into the unit interval, because that is what the rest of the subsystem means by a brightness (see `image_load`) and the type a filter answers with is always "Real". So `Import` produces a "Real" image, not a "Byte" one: an image whose stored range depended on the file's bit depth would make every downstream kernel's scale depend on it too.

## Implementation notes

- `Protected`.
- Returns the file name, so `Import[Export[f, img]]` is a round trip that can be written
  as a single expression.

**Attributes:** `Protected`.

## References

**See also:** [Image](../../image-processing/Image/), [Plot](../../graphics/Plot/), [ListPlot](../../graphics/ListPlot/), [CMYKColor](../../graphics/CMYKColor/), [ImageSize](../../other-advanced/ImageSize/), [AspectRatio](../../other-advanced/AspectRatio/), [ArrayPlot](../../graphics/ArrayPlot/), [DensityPlot](../../graphics/DensityPlot/)

- Source: [`src/imageio.c`](https://github.com/stblake/mathilda/blob/main/src/imageio.c)
- Specification: [`docs/spec/builtins/file-io.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/file-io.md)
- Tests: [`tests/test_graphics.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graphics.c)
- Tests: [`tests/test_image.c`](https://github.com/stblake/mathilda/blob/main/tests/test_image.c)
