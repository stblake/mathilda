# Export

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Export["file", obj] writes obj to a file, choosing the format from the file extension; Export["file", obj, "FMT"] states it explicitly. An Image writes to a raster file (PNG, JPEG, BMP, TGA); its samples outside the unit interval are clamped, since 8-bit output has no room for them. A Graphics object (the result of Plot, ListPlot, Graphics, ...) writes to PDF, PNG, or JPEG: PDF is a resolution-independent vector file produced without any external library and works headless, while PNG and JPEG render through the graphics backend and so need graphics support compiled in and a display. A Graphics3D object (Plot3D, ParametricPlot3D, ...) exports to PNG or JPEG the same way; it has no vector-PDF form. Returns the file name, so Import[Export[f, img]] round-trips.`**

## Examples (19)

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

### Applications (4)

An 8x12 grey image

```mathematica
In[16]:= img = Image[Table[N[(i + j)/16], {i, 8}, {j, 12}], "Real"];
```

Returns the file name

```mathematica
In[17]:= Export["/tmp/mathilda_export.png", img]
Out[17]= "/tmp/mathilda_export.png"
```

```mathematica
In[18]:= FileExistsQ["/tmp/mathilda_export.png"]
Out[18]= True
```

{width, height} = {12, 8}

```mathematica
In[19]:= ImageDimensions[Import["/tmp/mathilda_export.png"]]
Out[19]= {12, 8}
```

## Options & behaviour

### Graphics export

- **PDF** is a resolution-independent **vector** file written by a small built-in PDF
  emitter — no external library, no display, so it works headless and in every build. It
  walks the graphics primitives directly (`Line`, `Point`, `Polygon`, `Disk`/`Circle`,
  `Rectangle`, `Arrow`, `Text`) with the `RGBColor`/`GrayLevel`/`Hue`/`CMYKColor`,
  `Opacity`, `Thickness`/`AbsoluteThickness`, `PointSize`/`AbsolutePointSize`,
  `Dashing`/`AbsoluteDashing` (`Dashed`, `Dotted`, `DotDashed`) and `Directive`
  directives, each `{...}` list scoping the directives set inside it, and draws a framed
  set of axes with "nice" ticks and numeric labels. It honours `PlotStyle` (the style a
  single-curve `Plot`/`ListPlot`/`ParametricPlot` is drawn in), `Prolog`/`Epilog`,
  `PlotLabel`, `AxesLabel` and `PlotLegends` (a swatch-and-label legend to the right of the
  plot). Coordinates may be exact (`1/2`, `Pi/4`, `Sqrt[2]`): they are converted the same
  way the on-screen renderer converts them. Text uses the PDF base-14 Helvetica, so no font
  is embedded. This is the recommended format for print and for the book.
- In the PDF, `Text[s, pos, {ox, oy}]` aligns as Mathematica does (`{-1, 0}` puts the left
  end of `s` at `pos`, `{0, 0}` centres it, using the Helvetica advance widths), and
  `Text[Style[s, n | FontSize -> n | colour, ...], ...]` sets that string's size and
  colour. `Arrowheads[s]` fixes the arrowhead length at `s` times the plot width, and an
  arrow's shaft stops inside its head rather than poking past the point.
  `AspectRatio -> Automatic` maps x and y with one scale (the page height follows the
  data unless `ImageSize -> {w, h}` fixes both, in which case the picture is centred).
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
  depend on which formats the Raylib build happens to support. Directive sizes mean the
  same thing in both formats: `Thickness`/`PointSize`/`Dashing` are fractions of the plot
  width and the `Absolute*` forms are printer's points (one pixel in a raster).
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

**Algorithm.** `builtin_export` writes an object to a file and **returns the file name** (so
`Import[Export[f, img]]` is a one-expression round trip). The format comes from an explicit
third-argument string, else from the path's extension.

A `Graphics`/`Graphics3D` second argument is handled first. **PDF** of 2D graphics goes through
`graphics_export_pdf`, a dependency-free built-in vector emitter — no external library and no
display, so it works headless; a 3D scene has no vector projection and is declined. **PNG/JPEG**
render the scene to an RGBA buffer through the Raylib backend (`graphics_render_rgba` /
`graphics3d_render_rgba`, sized by `ImageSize`/`AspectRatio`) and encode it with `stb_image_write`
— so they need `USE_GRAPHICS` **and** a usable GUI session, and return `$Failed` gracefully
otherwise while PDF still works.

Otherwise the second argument is loaded as an image by `image_load` (which declines an
`Image3D`, returning `NULL`, rather than silently writing a middle slice). Each `double` sample
is converted to a byte with **clamping, not wrapping**: `v <= 0` (and `NaN`, via `!(v > 0)`) → `0`,
`v >= 1` → `255`, else `v·255 + 0.5`. The bytes are written by extension/format with
`stbi_write_png` / `stbi_write_jpg` (quality 90) / `stbi_write_bmp` / `stbi_write_tga`; an
unclaimed format returns `NULL` (unevaluated) and a failed write returns `$Failed`.

**Data structures.** The vendored public-domain **`stb_image_write`** encoder
(`src/external/stb/stb_image_write.h`), a transient `unsigned char` byte buffer, and — for
graphics — the Raylib RGBA render buffer or the PDF emitter's primitive walk. Encoding the
bytes with `stb` rather than Raylib's own writers means JPEG output does not depend on which
formats a given Raylib build supports.

**Complexity / limits.** `O(w·h·channels)` for an image; a plot's cost is dominated by
rendering. Clamping is deliberate: 8-bit output has nowhere to put an unsharp mask's legitimate
overshoot, and wrapping would turn a bright highlight black. `ATTR_PROTECTED`.

- `Protected`.
- Returns the file name, so `Import[Export[f, img]]` is a round trip that can be written
  as a single expression.

**Attributes:** `Protected`.

## References

**See also:** [Image](../../image-processing/Image/), [Plot](../../graphics/Plot/), [ListPlot](../../graphics/ListPlot/), [CMYKColor](../../graphics/CMYKColor/), [ParametricPlot](../../graphics/ParametricPlot/), [ImageSize](../../other-advanced/ImageSize/), [AspectRatio](../../other-advanced/AspectRatio/), [ArrayPlot](../../graphics/ArrayPlot/)

- Source: [`src/imageio.c`](https://github.com/stblake/mathilda/blob/main/src/imageio.c)
- Specification: [`docs/spec/builtins/file-io.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/file-io.md)
- Tests: [`tests/test_graphics.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graphics.c)
- Tests: [`tests/test_graphplot.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graphplot.c)
- Tests: [`tests/test_image.c`](https://github.com/stblake/mathilda/blob/main/tests/test_image.c)

## Notes & additional examples

### Notes

`Export["file", obj]` writes `obj` to a file, choosing the format from the
extension; `Export["file", obj, "FMT"]` states it. An `Image` writes to a raster
file (PNG, JPEG, BMP, TGA), and a `Graphics`/`Graphics3D` object writes to a graphic
(PDF, PNG, or JPEG). It **returns the file name**, which is what makes
`Import[Export[f, img]]` a single-expression round trip.

PDF of a 2D `Graphics` is a resolution-independent vector file from a built-in
emitter — no external library and no display, so it works headless and is the
recommended print format. PNG/JPEG render through the graphics backend, so they
need `USE_GRAPHICS` and a GUI session and otherwise return `$Failed` gracefully;
PDF still works. On raster export, samples outside the unit interval are **clamped**
(not wrapped), and `NaN` clamps to `0`. An `Image3D` is declined rather than
silently reduced to a slice. Writing is by the vendored `stb_image_write`.

There is no in-memory file target, so `Export` necessarily touches the filesystem;
these examples write under `/tmp`.
