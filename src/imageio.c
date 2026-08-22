/* imageio.c -- Import and Export for raster image files.
 *
 * Until this landed, every image in the system had to be typed out as an array of numbers, which
 * makes the whole subsystem a demonstration rather than a tool: a filter is judged on photographs,
 * and a synthetic checkerboard cannot show what a bilateral filter does that a Gaussian does not.
 *
 * WHY A VENDORED DECODER. JPEG decoding is a baseline-Huffman-plus-IDCT project of its own and PNG
 * needs an inflate, so the choice is between vendoring or making libpng and libjpeg hard build
 * requirements. Two dependency-free public-domain headers cost less than either, and -- unlike a
 * system library -- they cannot be missing at a user's site, which for an `Import` is the whole
 * point. The headers are included HERE AND NOWHERE ELSE so that this is the only object file
 * carrying third-party code.
 *
 * WHAT A SAMPLE MEANS. A decoded 8-bit sample is scaled by 1/255 into the unit interval, because
 * that is what the rest of the subsystem means by a brightness (see `image_load`) and the type a
 * filter answers with is always "Real". So `Import` produces a "Real" image, not a "Byte" one:
 * an image whose stored range depended on the file's bit depth would make every downstream
 * kernel's scale depend on it too.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "expr.h"
#include "symtab.h"
#include "attr.h"
#include "image.h"
#include "imageio.h"
#include "random.h"
#include "options.h"
#include "sym_names.h"
#include "graphics/graphics_export.h"

/* The vendored decoders. Their diagnostics are suppressed rather than fixed: the files are
 * upstream text, byte-identical to the release, so they can be re-fetched without a merge -- see
 * src/external/stb/README.md. Suppression is scoped to this include and popped immediately. */
#if defined(__GNUC__) || defined(__clang__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-parameter"
#pragma GCC diagnostic ignored "-Wsign-compare"
#pragma GCC diagnostic ignored "-Wmissing-field-initializers"
#pragma GCC diagnostic ignored "-Wcast-qual"
/* With internal linkage (below) any stb entry point imageio.c never calls is an "unused static
 * function"; harmless for a vendored header we deliberately include whole. */
#pragma GCC diagnostic ignored "-Wunused-function"
/* Because we define STBI_FAILURE_USERMSG (below), stb's stbi__err(x,y) expands to stbi__err(y) and
 * discards x. One PNG error site builds an x it then never reads (invalid_chunk), which GCC flags
 * as set-but-not-used. Upstream, unfixable without editing the vendored header. */
#pragma GCC diagnostic ignored "-Wunused-but-set-variable"
#endif

/* STB_IMAGE_STATIC / STB_IMAGE_WRITE_STATIC give every stb function internal linkage, so its code
 * lives in imageio.o alone and exports no `stbi_*` symbols. raylib's libraylib.a bakes its OWN copy
 * of stb_image into rtextures.c.o with external linkage; without this scoping the two copies collide
 * at link time when USE_GRAPHICS is on -- "multiple definition of stbi_load ..." (issue #65). Making
 * ours file-local removes the clash at the root, with no need for the `-z muldefs` linker hack, which
 * would merely paper over an ODR violation by silently keeping whichever copy the linker saw first.
 * This is sound precisely because -- as the file header notes -- imageio.c is the ONLY translation
 * unit that both defines and calls these functions. */
#define STB_IMAGE_STATIC
#define STB_IMAGE_IMPLEMENTATION
#define STBI_FAILURE_USERMSG      /* readable reasons, e.g. "unknown image type" */
#include "external/stb/stb_image.h"

#define STB_IMAGE_WRITE_STATIC
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "external/stb/stb_image_write.h"

#if defined(__GNUC__) || defined(__clang__)
#pragma GCC diagnostic pop
#endif

/* JPEG is lossy, so the quality is a documented constant rather than a silent one. 90 is the usual
 * "visually indistinguishable" setting; a caller who cares should Export to PNG. */
#define JPEG_QUALITY 90

/* ------------------------------------------------------------------ helpers */

/* Case-insensitive extension test. `path` may have no extension at all. */
static bool has_ext(const char* path, const char* ext)
{
    const char* dot = strrchr(path, '.');
    size_t i;
    if (!dot) return false;
    dot++;
    for (i = 0; ext[i] && dot[i]; i++) {
        int a = dot[i], b = ext[i];
        if (a >= 'A' && a <= 'Z') a += 'a' - 'A';
        if (b >= 'A' && b <= 'Z') b += 'a' - 'A';
        if (a != b) return false;
    }
    return ext[i] == '\0' && dot[i] == '\0';
}

/* The formats stb can decode. Listed so that a file we cannot read is distinguishable from a file
 * we can read and failed to: the first leaves Import unevaluated (a future format handler may
 * claim it), the second is $Failed. */
static bool is_decodable(const char* path)
{
    static const char* exts[] = {"png", "jpg", "jpeg", "bmp", "gif", "tga", "psd",
                                 "hdr", "pic", "ppm", "pgm", "pnm", NULL};
    size_t i;
    for (i = 0; exts[i]; i++)
        if (has_ext(path, exts[i])) return true;
    return false;
}

static Expr* failed(void) { return expr_new_symbol("$Failed"); }

#ifdef USE_GRAPHICS
/* Pixel dimensions for a rasterised Graphics export: honour ImageSize
 * (a number w -> {w, 3w/4}, or an explicit {w, h}), else a readable default.
 * Only needed for the PNG/JPEG (Raylib) path, hence the USE_GRAPHICS guard. */
static void graphics_raster_size(const Expr* g, int* w, int* h)
{
    *w = 720; *h = 540;                 /* default 4:3 */
    if (!g || g->type != EXPR_FUNCTION) return;
    for (size_t i = 1; i < g->data.function.arg_count; i++) {
        const Expr* a = g->data.function.args[i];
        if (!a || a->type != EXPR_FUNCTION || !a->data.function.head
            || a->data.function.head->type != EXPR_SYMBOL
            || a->data.function.head->data.symbol.name != SYM_Rule
            || a->data.function.arg_count != 2) continue;
        const Expr* lhs = a->data.function.args[0];
        const Expr* rhs = a->data.function.args[1];
        if (!lhs || lhs->type != EXPR_SYMBOL || lhs->data.symbol.name != SYM_ImageSize)
            continue;
        if (rhs->type == EXPR_INTEGER) {
            *w = (int)rhs->data.integer; *h = (int)(rhs->data.integer * 3 / 4);
        } else if (rhs->type == EXPR_REAL) {
            *w = (int)rhs->data.real; *h = (int)(rhs->data.real * 0.75);
        } else if (rhs->type == EXPR_FUNCTION && rhs->data.function.head
                   && rhs->data.function.head->type == EXPR_SYMBOL
                   && rhs->data.function.head->data.symbol.name == SYM_List
                   && rhs->data.function.arg_count == 2) {
            const Expr* ew = rhs->data.function.args[0];
            const Expr* eh = rhs->data.function.args[1];
            if (ew->type == EXPR_INTEGER)   *w = (int)ew->data.integer;
            else if (ew->type == EXPR_REAL) *w = (int)ew->data.real;
            if (eh->type == EXPR_INTEGER)   *h = (int)eh->data.integer;
            else if (eh->type == EXPR_REAL) *h = (int)eh->data.real;
        }
    }
    if (*w < 16)   *w = 16;
    if (*h < 16)   *h = 16;
    if (*w > 8000) *w = 8000;
    if (*h > 8000) *h = 8000;
}
#endif /* USE_GRAPHICS */

/* A borrowed string argument, or NULL when the argument is not a string. */
static const char* str_arg(const Expr* e, size_t i)
{
    const Expr* a;
    if (e->type != EXPR_FUNCTION || i >= e->data.function.arg_count) return NULL;
    a = e->data.function.args[i];
    return (a && a->type == EXPR_STRING) ? a->data.string : NULL;
}

/* ------------------------------------------------------------------- Import */

/* Import[file] / Import[file, "Image"] -- decode a raster file into an Image.
 *
 * Returns $Failed when the file is missing or malformed, and NULL (unevaluated) when the path is
 * not a format this handles at all, which keeps the door open for Import of non-image formats
 * without this function having to claim them.
 */
static Expr* builtin_import(Expr* res)
{
    const char* path;
    const char* fmt;
    unsigned char* px;
    double* buf;
    Expr* out;
    int w = 0, h = 0, n = 0;
    size_t count, i;

    if (res->type != EXPR_FUNCTION) return NULL;
    if (res->data.function.arg_count < 1 || res->data.function.arg_count > 2) return NULL;

    path = str_arg(res, 0);
    if (!path) return NULL;

    /* A second argument is an element or format request. Only the ones that mean "the image
     * itself" are handled; "Data" or a channel request would be a different return type and is
     * better absent than wrong. */
    if (res->data.function.arg_count == 2) {
        fmt = str_arg(res, 1);
        if (!fmt) return NULL;
        if (strcmp(fmt, "Image") != 0 && strcmp(fmt, "PNG") != 0 && strcmp(fmt, "JPEG") != 0 &&
            strcmp(fmt, "JPG") != 0 && strcmp(fmt, "BMP") != 0 && strcmp(fmt, "GIF") != 0 &&
            strcmp(fmt, "TGA") != 0 && strcmp(fmt, "PNM") != 0)
            return NULL;
    } else if (!is_decodable(path)) {
        return NULL;
    }

    /* 0 channels means "however many the file has", so a grey file stays grey and an RGBA file
     * keeps its alpha -- forcing RGB would invent two channels for the first and discard
     * transparency for the second. */
    px = stbi_load(path, &w, &h, &n, 0);
    if (!px) return failed();
    if (w <= 0 || h <= 0 || n <= 0) { stbi_image_free(px); return failed(); }

    count = (size_t)w * (size_t)h * (size_t)n;
    buf = (double*)malloc(count * sizeof(double));
    if (!buf) { stbi_image_free(px); return failed(); }
    for (i = 0; i < count; i++) buf[i] = (double)px[i] / 255.0;
    stbi_image_free(px);

    /* image_build_real hands back a canonical, PACKED image -- the same representation every
     * filter produces, so an imported photograph is not a second-class citizen. */
    out = image_build_real(buf, (size_t)w, (size_t)h, (size_t)n);
    free(buf);
    return out ? out : failed();
}

/* ------------------------------------------------------------------- Export */

/* Export[file, image] / Export[file, image, "PNG"] -- write an Image to a raster file.
 *
 * Returns the file name, which is what makes `Import[Export[f, img]]` a round trip that can be
 * written as one expression, and is Mathematica's return value too.
 */
static Expr* builtin_export(Expr* res)
{
    const char* path;
    const char* fmt = NULL;
    const Expr* img;
    double* pix = NULL;
    unsigned char* bytes;
    size_t w = 0, h = 0, ch = 0, count, i;
    int ok = 0;
    Expr* out;

    if (res->type != EXPR_FUNCTION) return NULL;
    if (res->data.function.arg_count < 2 || res->data.function.arg_count > 3) return NULL;

    path = str_arg(res, 0);
    if (!path) return NULL;
    if (res->data.function.arg_count == 3) {
        fmt = str_arg(res, 2);
        if (!fmt) return NULL;
    }

    /* Graphics[...] (Plot/ListPlot/Graphics output) -> PDF/PNG/JPEG. Handled
     * before image_load, which only understands Image[]. PDF is a
     * dependency-free vector writer; PNG/JPEG reuse the Raylib renderer. */
    {
        const Expr* gobj = res->data.function.args[1];
        if (gobj && gobj->type == EXPR_FUNCTION && gobj->data.function.head
            && gobj->data.function.head->type == EXPR_SYMBOL
            && (gobj->data.function.head->data.symbol.name == SYM_Graphics
                || gobj->data.function.head->data.symbol.name == SYM_Graphics3D)) {
            int is3d   = (gobj->data.function.head->data.symbol.name == SYM_Graphics3D);
            int is_pdf = (fmt && strcmp(fmt, "PDF") == 0) || (!fmt && has_ext(path, "pdf"));
            int is_png = (fmt && strcmp(fmt, "PNG") == 0) || (!fmt && has_ext(path, "png"));
            int is_jpg = (fmt && (strcmp(fmt, "JPEG") == 0 || strcmp(fmt, "JPG") == 0))
                         || (!fmt && (has_ext(path, "jpg") || has_ext(path, "jpeg")));
            if (!is_pdf && !is_png && !is_jpg) return NULL;  /* unclaimed format */

            int gok = 0;
            if (is_pdf) {
                /* PDF is a 2D vector format; 3D has no vector projection here. */
                gok = is3d ? 0 : graphics_export_pdf(gobj, path);
            } else {
#ifdef USE_GRAPHICS
                /* Raylib renders the plot to RGBA pixels (the 3D renderer for a
                 * Graphics3D); stb encodes them, so JPEG works regardless of
                 * what this Raylib build supports. */
                int gw, gh, ow = 0, oh = 0;
                graphics_raster_size(gobj, &gw, &gh);
                unsigned char* rgba = is3d
                    ? graphics3d_render_rgba(gobj, gw, gh, &ow, &oh)
                    : graphics_render_rgba(gobj, gw, gh, &ow, &oh);
                if (rgba) {
                    if (is_jpg) gok = stbi_write_jpg(path, ow, oh, 4, rgba, JPEG_QUALITY);
                    else        gok = stbi_write_png(path, ow, oh, 4, rgba, ow * 4);
                    free(rgba);
                }
#else
                gok = 0;        /* PNG/JPEG of graphics needs the Raylib renderer */
#endif
            }
            if (!gok) return failed();
            out = expr_new_string(path);
            return out ? out : failed();
        }
    }

    img = res->data.function.args[1];
    if (!image_load(img, &w, &h, &ch, &pix)) {
        /* A volume has no single raster to write, and silently writing its middle slice would be
         * a lie about what was exported. */
        return NULL;
    }

    count = w * h * ch;
    bytes = (unsigned char*)malloc(count ? count : 1);
    if (!bytes) { free(pix); return failed(); }
    for (i = 0; i < count; i++) {
        /* Clamp rather than wrap. A filter can legitimately overshoot the unit interval (an
         * unsharp mask does), and 8-bit output has nowhere to put that; wrapping would turn a
         * bright highlight black, which looks like a bug in the filter. The !(v > 0) form also
         * catches NaN. */
        double v = pix[i];
        if (!(v > 0.0)) bytes[i] = 0;
        else if (v >= 1.0) bytes[i] = 255;
        else bytes[i] = (unsigned char)(v * 255.0 + 0.5);
    }
    free(pix);

    if ((fmt && (strcmp(fmt, "PNG") == 0)) || (!fmt && has_ext(path, "png")))
        ok = stbi_write_png(path, (int)w, (int)h, (int)ch, bytes, (int)(w * ch));
    else if ((fmt && (strcmp(fmt, "JPEG") == 0 || strcmp(fmt, "JPG") == 0)) ||
             (!fmt && (has_ext(path, "jpg") || has_ext(path, "jpeg"))))
        ok = stbi_write_jpg(path, (int)w, (int)h, (int)ch, bytes, JPEG_QUALITY);
    else if ((fmt && strcmp(fmt, "BMP") == 0) || (!fmt && has_ext(path, "bmp")))
        ok = stbi_write_bmp(path, (int)w, (int)h, (int)ch, bytes);
    else if ((fmt && strcmp(fmt, "TGA") == 0) || (!fmt && has_ext(path, "tga")))
        ok = stbi_write_tga(path, (int)w, (int)h, (int)ch, bytes);
    else {
        free(bytes);
        return NULL;              /* an unclaimed format, not a failed write */
    }
    free(bytes);
    if (!ok) return failed();

    out = expr_new_string(path);
    return out ? out : failed();
}


/* -------------------------------------------------------------- RandomImage */

/* The default size, matching Mathematica's. Large enough that a filter's effect is visible and
 * small enough to print, resize and Export without thinking about it. */
#define RANDIMG_DEFAULT 150

/* A dimension pair {w, h}, or a single n meaning {n, n}. Returns false when the argument is
 * neither -- a symbolic size has to stay unevaluated rather than become a guess. */
static bool read_size(const Expr* e, size_t* w, size_t* h)
{
    if (!e) return false;
    if (e->type == EXPR_INTEGER && e->data.integer > 0) {
        *w = *h = (size_t)e->data.integer;
        return true;
    }
    if (e->type == EXPR_FUNCTION && e->data.function.arg_count == 2 &&
        e->data.function.head && e->data.function.head->type == EXPR_SYMBOL &&
        strcmp(e->data.function.head->data.symbol.name, "List") == 0) {
        const Expr* a = e->data.function.args[0];
        const Expr* b = e->data.function.args[1];
        if (a && b && a->type == EXPR_INTEGER && b->type == EXPR_INTEGER &&
            a->data.integer > 0 && b->data.integer > 0) {
            *w = (size_t)a->data.integer;
            *h = (size_t)b->data.integer;
            return true;
        }
    }
    return false;
}

/* RandomImage[] / RandomImage[max] / RandomImage[max, {w, h}] / ColorSpace -> "RGB"
 *
 * Noise is the input a filter is most often judged on -- a smoothing radius means nothing on a
 * checkerboard and everything on a noise field -- and until this existed the only way to get one
 * was to write a deterministic Mod expression and call it random.
 *
 * Samples come from `random_uniform_01`, the SAME stream RandomReal draws on, so SeedRandom makes
 * a random image reproducible. A private generator would have been simpler and would have quietly
 * made this the one random builtin that ignores the seed.
 */
static Expr* builtin_random_image(Expr* res)
{
    const Expr* colorspace = NULL;
    bool cs_given = false;
    OptEntry opts[] = { { "ColorSpace", &colorspace, &cs_given } };
    size_t argc = 0;
    size_t w = RANDIMG_DEFAULT, h = RANDIMG_DEFAULT, ch = 1, n, i;
    double max = 1.0;
    double* buf;
    Expr* out;

    if (res->type != EXPR_FUNCTION) return NULL;
    if (!options_extract(res, "RandomImage", opts, 1, &argc)) return NULL;
    if (argc > 2) return NULL;

    if (argc >= 1) {
        const Expr* a = res->data.function.args[0];
        if (!a) return NULL;
        if (a->type == EXPR_INTEGER)   max = (double)a->data.integer;
        else if (a->type == EXPR_REAL) max = a->data.real;
        else return NULL;
        if (!(max > 0.0)) return NULL;
    }
    if (argc == 2 && !read_size(res->data.function.args[1], &w, &h)) return NULL;

    if (cs_given && colorspace && colorspace->type == EXPR_STRING) {
        const char* cs = colorspace->data.string;
        if (strcmp(cs, "RGB") == 0)            ch = 3;
        else if (strcmp(cs, "Grayscale") == 0) ch = 1;
        else return NULL;   /* an unsupported space is not silently grey */
    }

    n = w * h * ch;
    buf = (double*)malloc(n * sizeof(double));
    if (!buf) return NULL;
    /* Scaled by `max` and NOT clamped: the caller asked for that range, and a "Real" image is
     * free to hold values above 1 -- Export is where clamping belongs. */
    for (i = 0; i < n; i++) buf[i] = random_uniform_01() * max;

    out = image_build_real(buf, w, h, ch);
    free(buf);
    return out;
}

void imageio_init(void)
{
    symtab_add_builtin("Import", builtin_import);
    symtab_get_def("Import")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("Import",
        "Import[\"file\"] reads a raster image file (PNG, JPEG, BMP, GIF, TGA, PSD, HDR, PNM) and "
        "returns an Image. Import[\"file\", \"Image\"] is the same. Samples are scaled by 1/255 "
        "into the unit interval, so the result is a \"Real\" image whatever the file's bit depth, "
        "and the file's channel count is preserved -- grey stays 1 channel, RGBA keeps its alpha. "
        "Gives $Failed for a missing or malformed file.");

    symtab_add_builtin("Export", builtin_export);
    symtab_get_def("Export")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("Export",
        "Export[\"file\", obj] writes obj to a file, choosing the format from the file extension; "
        "Export[\"file\", obj, \"FMT\"] states it explicitly. An Image writes to a raster file "
        "(PNG, JPEG, BMP, TGA); its samples outside the unit interval are clamped, since 8-bit "
        "output has no room for them. A Graphics object (the result of Plot, ListPlot, Graphics, "
        "...) writes to PDF, PNG, or JPEG: PDF is a resolution-independent vector file produced "
        "without any external library and works headless, while PNG and JPEG render through the "
        "graphics backend and so need graphics support compiled in and a display. A Graphics3D "
        "object (Plot3D, ParametricPlot3D, ...) exports to PNG or JPEG the same way; it has no "
        "vector-PDF form. Returns the file name, so Import[Export[f, img]] round-trips.");
    symtab_add_builtin("RandomImage", builtin_random_image);
    symtab_get_def("RandomImage")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("RandomImage",
        "RandomImage[] gives a 150x150 grey image of uniform noise on [0, 1]. RandomImage[max] "
        "scales the range to [0, max]; RandomImage[max, {w, h}] sets the size, and a single n "
        "means {n, n}. ColorSpace -> \"RGB\" gives three independent channels. Samples are drawn "
        "from the same stream as RandomReal, so SeedRandom makes the result reproducible.");
}
