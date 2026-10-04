---
source: src/imageio.c
---
**Algorithm.** `builtin_import` decodes a raster file into an `Image`. It takes a path and an
optional second argument that must name an image format (`"Image"`, `"PNG"`, `"JPEG"`/`"JPG"`,
`"BMP"`, `"GIF"`, `"TGA"`, `"PNM"`); with no second argument the path must be one `is_decodable`
recognises by extension, otherwise the call is left unevaluated (`NULL`) — so `Import` does not
pretend to implement formats it has no decoder for, and the door stays open for non-image
imports elsewhere. Decoding is `stbi_load(path, &w, &h, &n, 0)`: the `0` requested-channel count
means "**however many the file has**", so a grey file stays 1-channel and an RGBA file keeps its
alpha rather than being forced to RGB.

The `unsigned char` samples are scaled by `1/255` into a `double` buffer, and
`image_build_real(buf, w, h, n)` produces a **canonical, packed `"Real"` image** — the same
representation every filter yields, so an imported photograph is not a second-class citizen
downstream. A missing or malformed file returns `$Failed` (`failed()`); a zero/invalid
dimension from the decoder also returns `$Failed`.

**Data structures.** The vendored public-domain **`stb_image`** decoder
(`src/external/stb/stb_image.h`), a transient `double` sample buffer, and the packed `Image`
(an `NDArray`-backed real tensor) that `image_build_real` returns. No system image library is a
build requirement; `stbi_load` frees its own pixel buffer via `stbi_image_free`.

**Complexity / limits.** `O(w·h·channels)` to decode and rescale. The fixed `1/255` scaling
means the stored range is always `[0,1]` whatever the file's bit depth, which is what keeps
every downstream kernel's arithmetic meaning stable (see `ImageData`). `ATTR_PROTECTED`.
