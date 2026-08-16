# File I/O

12 built-in function(s) in this category.

- [`Export`](Export.md) — Export["file", image] writes an Image to a raster file, choosing the format from the file extension (PNG, JPEG, BMP, TGA). Export["file", image, "PNG"] states the format explicitly. Samples outside the unit interval are clamped, since 8-bit output has no room for them. Returns the file name, so Import[Export[f, img]] is a round trip.  _(Stable)_
- [`FileBaseName`](FileBaseName.md) — FileBaseName["file"]  _(Stable)_
- [`FileExistsQ`](FileExistsQ.md) — FileExistsQ["name"]  _(Stable)_
- [`FileExtension`](FileExtension.md) — FileExtension["file"]  _(Stable)_
- [`FileNameJoin`](FileNameJoin.md) — FileNameJoin[{"name1", "name2", ...}]  _(Stable)_
- [`FileNameSplit`](FileNameSplit.md) — FileNameSplit["name"]  _(Stable)_
- [`FilePrint`](FilePrint.md) — FilePrint["file"]  _(Stable)_
- [`Get`](Get.md) — Get["filename"]  _(Stable)_
- [`Import`](Import.md) — Import["file"] reads a raster image file (PNG, JPEG, BMP, GIF, TGA, PSD, HDR, PNM) and returns an Image. Import["file", "Image"] is the same. Samples are scaled by 1/255 into the unit interval, so the result is a "Real" image whatever the file's bit depth, and the file's channel count is preserved -- grey stays 1 channel, RGBA keeps its alpha. Gives $Failed for a missing or malformed file.  _(Stable)_
- [`LoadModule`](LoadModule.md) — LoadModule["relpath"]  _(Experimental)_
- [`Put`](Put.md) — Put[expr, "filename"] or expr >> "filename"  _(Stable)_
- [`PutAppend`](PutAppend.md) — PutAppend[expr, "filename"] or expr >>> "filename"  _(Stable)_
