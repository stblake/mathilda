/* image.c -- Image, ImageQ, ImageDimensions, ImageChannels, ImageType, ImageData.
 *
 * The image subsystem starts here, and starts small on purpose: a validated representation with
 * accessors, and nothing that filters. ImageConvolve, GaussianFilter, Binarize and the rest all
 * need to ask an image its size, its channel count and its pixel type before they can do
 * anything, so those questions get answered first and answered once.
 *
 * THE CANONICAL FORM IS Image[data, type], WHICH IS ALSO REAL WOLFRAM SYNTAX. Image[data]
 * normalises to it by inferring the type. That normalisation is what makes validity decidable:
 * a builtin that returns NULL leaves its expression alone, so a valid Image[data] and a
 * nonsensical Image["hello"] would be indistinguishable -- both unevaluated. Normalising means
 * a canonical two-argument form is exactly what passed validation, and ImageQ tests for it.
 *
 * ROW-MAJOR, AND THE TRANSPOSITION IS THE WHOLE TRAP. ImageData is a HEIGHT x WIDTH array --
 * data[y][x], rows down the image -- while ImageDimensions reports {WIDTH, HEIGHT}. The two are
 * transposed relative to each other, which is Wolfram's convention and the single most common
 * source of silently-wrong image code. A test pins a non-square image precisely because a
 * square one cannot tell the two apart.
 *
 * WHY THE TYPE MATTERS RATHER THAN BEING DECORATION. It fixes the RANGE of a stored value, and
 * ImageData returns reals in the unit interval by scaling out that range: a "Byte" 255 becomes
 * exactly 1.0 and a "Bit" 1 becomes exactly 1.0. Storing the original integers and scaling on
 * read (rather than scaling on construction) keeps Image[...] round-trippable -- what went in
 * is what FullForm shows -- and makes the scaling an assertion about one function instead of a
 * property spread across two.
 */
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "expr.h"
#include "symtab.h"
#include "attr.h"
#include "sym_names.h"
#include "linalg/numarray.h"
#include "image.h"

/* ---- shape and type inspection of a nested List ---- */

static bool img_is_list(const Expr* e) {
    return e && e->type == EXPR_FUNCTION && e->data.function.head
        && e->data.function.head->type == EXPR_SYMBOL
        && e->data.function.head->data.symbol.name == SYM_List;
}

/* Walk a rank-2 or rank-3 numeric array, checking it is RECTANGULAR and collecting its shape.
 *
 * Raggedness is rejected rather than padded or truncated. A ragged array is not an image, and
 * every downstream filter indexes it as if it were rectangular -- so accepting one would turn a
 * clear refusal here into an out-of-bounds read somewhere far away. */
static bool img_shape(const Expr* d, size_t* h, size_t* w, size_t* c,
                      bool* all_int, double* lo, double* hi) {
    if (!img_is_list(d) || d->data.function.arg_count == 0) return false;
    size_t rows = d->data.function.arg_count;
    size_t cols = 0, chan = 0;
    bool rank3 = false;
    *all_int = true;
    *lo = 0.0; *hi = 0.0;
    bool first = true;

    for (size_t y = 0; y < rows; y++) {
        const Expr* row = d->data.function.args[y];
        if (!img_is_list(row) || row->data.function.arg_count == 0) return false;
        if (y == 0) cols = row->data.function.arg_count;
        else if (row->data.function.arg_count != cols) return false;   /* ragged */

        for (size_t x = 0; x < cols; x++) {
            const Expr* px = row->data.function.args[x];
            if (img_is_list(px)) {
                /* A list-valued pixel: this is a multichannel image. Every pixel must then be
                 * a list of the SAME length, or the channel count is not a property of the
                 * image at all. */
                if (y == 0 && x == 0) { rank3 = true; chan = px->data.function.arg_count; }
                if (!rank3 || px->data.function.arg_count != chan || chan == 0) return false;
                for (size_t k = 0; k < chan; k++) {
                    double re = 0.0, im = 0.0;
                    const Expr* v = px->data.function.args[k];
                    if (!na_read_scalar(v, &re, &im) || im != 0.0) return false;
                    if (v->type != EXPR_INTEGER) *all_int = false;
                    if (first) { *lo = *hi = re; first = false; }
                    else { if (re < *lo) *lo = re; if (re > *hi) *hi = re; }
                }
            } else {
                if (rank3) return false;         /* mixed rank-2 and rank-3 rows */
                double re = 0.0, im = 0.0;
                if (!na_read_scalar(px, &re, &im) || im != 0.0) return false;
                if (px->type != EXPR_INTEGER) *all_int = false;
                if (first) { *lo = *hi = re; first = false; }
                else { if (re < *lo) *lo = re; if (re > *hi) *hi = re; }
            }
        }
    }
    *h = rows; *w = cols; *c = rank3 ? chan : 1;
    return true;
}

static const char* img_type_name(ImgType t) {
    switch (t) {
        case IMG_BIT:  return "Bit";
        case IMG_BYTE: return "Byte";
        default:       return "Real";
    }
}

static bool img_type_from_name(const char* s, ImgType* out) {
    if (strcmp(s, "Bit") == 0)  { *out = IMG_BIT;  return true; }
    if (strcmp(s, "Byte") == 0) { *out = IMG_BYTE; return true; }
    if (strcmp(s, "Real") == 0) { *out = IMG_REAL; return true; }
    return false;
}

/* Shape of a canonical image WITHOUT walking every pixel.
 *
 * Added because a benchmark caught the obvious implementation being far too slow: routing
 * image_info through the full validator made ImageDimensions cost 0.59 ms on a 512x512 image,
 * because asking how wide an image is touched all 262144 pixels. A filter pipeline queries
 * dimensions constantly, so that compounds.
 *
 * The canonical two-argument form can only be produced by Image[], which validates every pixel,
 * so re-validating them here is redundant for any image the system built. What is still checked
 * is RECTANGULARITY -- O(height) row-length comparisons, ~500x cheaper than O(pixels) -- because
 * a hand-typed Image[{{1,2},{3}}, "Bit"] must not pass as an image when every filter downstream
 * will index it as rectangular. Per-pixel numeric checking is the part that stays in
 * construction, where it is paid once. */
static bool img_shape_fast(const Expr* d, size_t* h, size_t* w, size_t* c) {
    if (!img_is_list(d) || d->data.function.arg_count == 0) return false;
    size_t rows = d->data.function.arg_count;
    const Expr* r0 = d->data.function.args[0];
    if (!img_is_list(r0) || r0->data.function.arg_count == 0) return false;
    size_t cols = r0->data.function.arg_count;
    for (size_t y = 1; y < rows; y++) {
        const Expr* row = d->data.function.args[y];
        if (!img_is_list(row) || row->data.function.arg_count != cols) return false;
    }
    const Expr* p0 = r0->data.function.args[0];
    size_t chan = 1;
    if (img_is_list(p0)) {
        chan = p0->data.function.arg_count;
        if (chan == 0) return false;
    }
    *h = rows; *w = cols; *c = chan;
    return true;
}

bool image_info(const Expr* e, size_t* width, size_t* height,
                size_t* channels, ImgType* type) {
    if (!e || e->type != EXPR_FUNCTION || e->data.function.arg_count != 2) return false;
    const Expr* hd = e->data.function.head;
    if (!hd || hd->type != EXPR_SYMBOL || strcmp(hd->data.symbol.name, "Image") != 0)
        return false;
    const Expr* ty = e->data.function.args[1];
    if (!ty || ty->type != EXPR_STRING) return false;
    ImgType t;
    if (!img_type_from_name(ty->data.string, &t)) return false;

    size_t h = 0, w = 0, c = 0;
    if (!img_shape_fast(e->data.function.args[0], &h, &w, &c)) return false;
    if (width) *width = w;
    if (height) *height = h;
    if (channels) *channels = c;
    if (type) *type = t;
    return true;
}

/* Image[data] -- validate and normalise to Image[data, type].
 * Image[data, "type"] -- validate against the stated type.
 *
 * TYPE INFERENCE, and where it deliberately stops. All-integer data whose values are only 0 and
 * 1 is "Bit"; all-integer within 0..255 is "Byte"; anything else is "Real". Inference looks at
 * the values and nothing else, so it is predictable from the data alone -- and it is why
 * Image[{{0, 1}, {1, 0}}] is a bit image while Image[{{0., 1.}, {1., 0.}}] is a real one, which
 * is a distinction a caller can rely on rather than having to guess. */
static Expr* builtin_image(Expr* res) {
    size_t argc = res->data.function.arg_count;
    if (argc != 1 && argc != 2) return NULL;

    /* Already canonical: leave it alone, or the evaluator would never reach a fixed point. */
    if (argc == 2 && image_info(res, NULL, NULL, NULL, NULL)) return NULL;

    Expr* data = res->data.function.args[0];
    size_t h = 0, w = 0, c = 0; bool all_int = false; double lo = 0.0, hi = 0.0;
    if (!img_shape(data, &h, &w, &c, &all_int, &lo, &hi)) return NULL;

    ImgType t;
    if (argc == 2) {
        Expr* ty = res->data.function.args[1];
        if (!ty || ty->type != EXPR_STRING || !img_type_from_name(ty->data.string, &t))
            return NULL;
        /* A stated type must be consistent with the data. Silently reinterpreting 300 as a byte
         * would corrupt every later scaling, so it declines instead. */
        if (t == IMG_BIT && !(all_int && lo >= 0.0 && hi <= 1.0)) return NULL;
        if (t == IMG_BYTE && !(all_int && lo >= 0.0 && hi <= 255.0)) return NULL;
    } else {
        if (all_int && lo >= 0.0 && hi <= 1.0)        t = IMG_BIT;
        else if (all_int && lo >= 0.0 && hi <= 255.0) t = IMG_BYTE;
        else                                          t = IMG_REAL;
    }

    Expr* two[2];
    two[0] = expr_copy(data);
    two[1] = expr_new_string(img_type_name(t));
    if (!two[0] || !two[1]) { expr_free(two[0]); expr_free(two[1]); return NULL; }
    return expr_new_function(expr_new_symbol("Image"), two, 2);
}

static Expr* builtin_imageq(Expr* res) {
    if (res->data.function.arg_count != 1) return NULL;
    bool q = image_info(res->data.function.args[0], NULL, NULL, NULL, NULL);
    return expr_new_symbol(q ? SYM_True : SYM_False);
}

/* ImageDimensions -- {WIDTH, HEIGHT}, transposed relative to ImageData's height x width. */
static Expr* builtin_imagedimensions(Expr* res) {
    if (res->data.function.arg_count != 1) return NULL;
    size_t w = 0, h = 0;
    if (!image_info(res->data.function.args[0], &w, &h, NULL, NULL)) return NULL;
    Expr* two[2];
    two[0] = expr_new_integer((int64_t)w);
    two[1] = expr_new_integer((int64_t)h);
    if (!two[0] || !two[1]) { expr_free(two[0]); expr_free(two[1]); return NULL; }
    return expr_new_function(expr_new_symbol(SYM_List), two, 2);
}

static Expr* builtin_imagechannels(Expr* res) {
    if (res->data.function.arg_count != 1) return NULL;
    size_t c = 0;
    if (!image_info(res->data.function.args[0], NULL, NULL, &c, NULL)) return NULL;
    return expr_new_integer((int64_t)c);
}

static Expr* builtin_imagetype(Expr* res) {
    if (res->data.function.arg_count != 1) return NULL;
    ImgType t;
    if (!image_info(res->data.function.args[0], NULL, NULL, NULL, &t)) return NULL;
    return expr_new_string(img_type_name(t));
}

/* Scale one stored value into the unit interval according to the image's type. */
static double img_to_unit(double v, ImgType t) {
    if (t == IMG_BYTE) return v / 255.0;
    return v;                       /* Bit is already 0 or 1; Real is already unit-scaled */
}

/* Recursively rebuild the pixel array, scaling leaves. The structure is copied rather than
 * flattened because ImageData must give back the SAME shape it was handed -- height x width for
 * grey, height x width x channels for colour, interleaved. */
static Expr* img_scale_tree(const Expr* e, ImgType t, bool raw) {
    if (img_is_list(e)) {
        size_t n = e->data.function.arg_count;
        Expr** kids = malloc(sizeof(Expr*) * n);
        if (!kids) return NULL;
        bool ok = true;
        for (size_t i = 0; i < n; i++) {
            kids[i] = img_scale_tree(e->data.function.args[i], t, raw);
            if (!kids[i]) ok = false;
        }
        if (!ok) {
            for (size_t i = 0; i < n; i++) expr_free(kids[i]);
            free(kids);
            return NULL;
        }
        Expr* out = expr_new_function(expr_new_symbol(SYM_List), kids, n);
        free(kids);
        return out;
    }
    double re = 0.0, im = 0.0;
    if (!na_read_scalar(e, &re, &im) || im != 0.0) return NULL;
    /* expr_copy takes a mutable pointer; the cast is safe because it only reads. */
    if (raw) return expr_copy((Expr*)e);
    return expr_new_real(img_to_unit(re, t));
}

/* ImageData[img] -- pixels as reals in [0, 1].
 * ImageData[img, "Byte"] / [img, "Bit"] / [img, "Real"] -- the stored values, unscaled.
 *
 * The default SCALES, which is Wolfram's behaviour and the reason the type is stored at all: a
 * "Byte" 255 comes back as exactly 1.0 and a 0 as exactly 0.0. Requesting a type returns the
 * raw stored values instead, which is the round-trip a caller needs to get back what it put in.
 */
static Expr* builtin_imagedata(Expr* res) {
    size_t argc = res->data.function.arg_count;
    if (argc != 1 && argc != 2) return NULL;
    Expr* img = res->data.function.args[0];
    ImgType t;
    if (!image_info(img, NULL, NULL, NULL, &t)) return NULL;

    bool raw = false;
    if (argc == 2) {
        Expr* ty = res->data.function.args[1];
        ImgType want;
        if (!ty || ty->type != EXPR_STRING || !img_type_from_name(ty->data.string, &want))
            return NULL;
        /* Only the image's OWN type may be requested. Converting between types is a separate
         * operation with its own rounding decisions, and quietly doing it here would hide them. */
        if (want != t) return NULL;
        raw = true;
    }
    return img_scale_tree(img->data.function.args[0], t, raw);
}

void image_init(void) {
    symtab_add_builtin("Image", builtin_image);
    symtab_get_def("Image")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("Image",
        "Image[data] is a raster image, normalising to the canonical Image[data, type]. "
        "The data is a rectangular height x width array of pixel values, or "
        "height x width x channels for a colour image, so it is indexed data[[y, x]] with "
        "rows running down the image -- note that ImageDimensions reports {width, height}, "
        "transposed relative to this. The type is inferred from the values: all-integer data "
        "in {0, 1} is \"Bit\", all-integer in 0..255 is \"Byte\", anything else is \"Real\". "
        "Image[data, type] states the type instead, and declines if the data does not fit it. "
        "Ragged data declines rather than being padded.");

    symtab_add_builtin("ImageQ", builtin_imageq);
    symtab_get_def("ImageQ")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("ImageQ",
        "ImageQ[expr] gives True if expr is a valid image in canonical form, and False "
        "otherwise. Malformed input to Image stays unevaluated, so ImageQ is how validity is "
        "tested.");

    symtab_add_builtin("ImageDimensions", builtin_imagedimensions);
    symtab_get_def("ImageDimensions")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("ImageDimensions",
        "ImageDimensions[image] gives {width, height}. This is TRANSPOSED relative to "
        "ImageData, which returns a height x width array -- the same convention Mathematica "
        "uses.");

    symtab_add_builtin("ImageChannels", builtin_imagechannels);
    symtab_get_def("ImageChannels")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("ImageChannels",
        "ImageChannels[image] gives the number of colour channels: 1 for a grey image, "
        "otherwise the length of each pixel's value list (3 for RGB, 4 with an alpha "
        "channel).");

    symtab_add_builtin("ImageType", builtin_imagetype);
    symtab_get_def("ImageType")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("ImageType",
        "ImageType[image] gives the pixel type as \"Bit\", \"Byte\" or \"Real\". The type "
        "fixes the range of a stored value, which is what makes ImageData's scaling to the "
        "unit interval well defined.");

    symtab_add_builtin("ImageData", builtin_imagedata);
    symtab_get_def("ImageData")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("ImageData",
        "ImageData[image] gives the pixel array as reals in [0, 1], scaling out the image's "
        "type -- a \"Byte\" 255 comes back as exactly 1.0. The array is height x width, or "
        "height x width x channels for a colour image, interleaved. "
        "ImageData[image, type] gives the stored values unscaled instead, where type must be "
        "the image's own type; converting between types is a separate operation with its own "
        "rounding, not something this does silently.");
}
