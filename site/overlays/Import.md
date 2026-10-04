### Worked examples

```mathematica
In[1]:= img = Image[Table[{N[i/8], N[j/8], 0.5}, {i, 8}, {j, 8}], "Real"];  (* an 8x8 RGB image *)
In[2]:= Export["/tmp/mathilda_import.png", img];  (* write it so there is something to read *)
In[3]:= ImageDimensions[Import["/tmp/mathilda_import.png"]]  (* decode it back *)
In[4]:= ImageType[Import["/tmp/mathilda_import.png"]]  (* always a "Real" image *)
```

### Notes

`Import["file"]` decodes a raster image (PNG, JPEG, BMP, GIF, TGA, PSD, HDR, PNM) by
content and returns an `Image`; `Import["file", "Image"]` states it explicitly.
Samples are scaled by `1/255` into the unit interval, so the result is a `"Real"`
image whatever the file's bit depth, and the file's channel count is **preserved**
— a grey file stays 1-channel and an RGBA file keeps its alpha.

The result is packed and canonical, the same representation a filter produces, so
`Import[Export[f, img]]` round-trips and the image needs no special-casing
downstream. A missing or malformed file gives `$Failed`; a path whose format is not
handled at all stays unevaluated (so `Import` does not appear to implement every
format). Decoding is by the vendored `stb_image`, so no system image library is a
build requirement.

These examples write to and read from `$TemporaryDirectory`-style paths under
`/tmp`; there is no in-memory image source, so `Import`/`Export` necessarily touch
the filesystem.
