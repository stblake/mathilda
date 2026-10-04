### Worked examples

```mathematica
In[1]:= img = Image[Table[N[(i + j)/16], {i, 8}, {j, 12}], "Real"];  (* an 8x12 grey image *)
In[2]:= Export["/tmp/mathilda_export.png", img]  (* returns the file name *)
In[3]:= FileExistsQ["/tmp/mathilda_export.png"]
In[4]:= ImageDimensions[Import["/tmp/mathilda_export.png"]]  (* {width, height} = {12, 8} *)
```

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
