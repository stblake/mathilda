---
source: src/imageio.c
---
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
