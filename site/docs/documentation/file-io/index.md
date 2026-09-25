# File I/O

23 built-in function(s) in this category.

- [`Close`](Close.md) — Close[stream]  _(Stable)_
- [`Export`](Export.md) — Export["file", obj] writes obj to a file, choosing the format from the file extension; Export["file", obj, "FMT"] states it explicitly. An Image writes to a raster file (PNG, JPEG, BMP, TGA); its samples outside the unit interval are clamped, since 8-bit output has no room for them. A Graphics object (the result of Plot, ListPlot, Graphics, ...) writes to PDF, PNG, or JPEG: PDF is a resolution-independent vector file produced without any external library and works headless, while PNG and JPEG render through the graphics backend and so need graphics support compiled in and a display. A Graphics3D object (Plot3D, ParametricPlot3D, ...) exports to PNG or JPEG the same way; it has no vector-PDF form. Returns the file name, so Import[Export[f, img]] round-trips.  _(Stable)_
- [`FileBaseName`](FileBaseName.md) — FileBaseName["file"]  _(Stable)_
- [`FileExistsQ`](FileExistsQ.md) — FileExistsQ["name"]  _(Stable)_
- [`FileExtension`](FileExtension.md) — FileExtension["file"]  _(Stable)_
- [`FileNameJoin`](FileNameJoin.md) — FileNameJoin[{"name1", "name2", ...}]  _(Stable)_
- [`FileNameSplit`](FileNameSplit.md) — FileNameSplit["name"]  _(Stable)_
- [`FilePrint`](FilePrint.md) — FilePrint["file"]  _(Stable)_
- [`Get`](Get.md) — Get["filename"]  _(Stable)_
- [`Import`](Import.md) — Import["file"] reads a raster image file (PNG, JPEG, BMP, GIF, TGA, PSD, HDR, PNM) and returns an Image. Import["file", "Image"] is the same. Samples are scaled by 1/255 into the unit interval, so the result is a "Real" image whatever the file's bit depth, and the file's channel count is preserved -- grey stays 1 channel, RGBA keeps its alpha. Gives $Failed for a missing or malformed file.  _(Stable)_
- [`LoadModule`](LoadModule.md) — LoadModule["relpath"]  _(Experimental)_
- [`OpenAppend`](OpenAppend.md) — OpenAppend["file"]  _(Stable)_
- [`OpenRead`](OpenRead.md) — OpenRead["file"]  _(Stable)_
- [`OpenWrite`](OpenWrite.md) — OpenWrite["file"]  _(Stable)_
- [`Put`](Put.md) — Put[expr, "filename"] or expr >> "filename"  _(Stable)_
- [`PutAppend`](PutAppend.md) — PutAppend[expr, "filename"] or expr >>> "filename"  _(Stable)_
- [`Read`](Read.md) — Read[stream]  _(Stable)_
- [`ReadList`](ReadList.md) — ReadList["file"]  _(Stable)_
- [`SetStreamPosition`](SetStreamPosition.md) — SetStreamPosition[stream, n]  _(Stable)_
- [`StreamPosition`](StreamPosition.md) — StreamPosition[stream]  _(Stable)_
- [`Streams`](Streams.md) — Streams[]  _(Stable)_
- [`Write`](Write.md) — Write[stream, expr1, expr2, ...]  _(Stable)_
- [`WriteString`](WriteString.md) — WriteString[stream, str1, str2, ...]  _(Stable)_
