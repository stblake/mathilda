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
#include "ndarray.h"
#include "pack.h"
#include "image.h"

/* ---- shape and type inspection of a nested List ---- */

static bool img_is_list(const Expr* e) {
    return e && e->type == EXPR_FUNCTION && e->data.function.head
        && e->data.function.head->type == EXPR_SYMBOL
        && e->data.function.head->data.symbol.name == SYM_List;
}

/* Defined below; img_shape needs it for the packed case. */
static bool img_shape_fast(const Expr* d, size_t* h, size_t* w, size_t* c);

/* Walk a rank-2 or rank-3 numeric array, checking it is RECTANGULAR and collecting its shape.
 *
 * Raggedness is rejected rather than padded or truncated. A ragged array is not an image, and
 * every downstream filter indexes it as if it were rectangular -- so accepting one would turn a
 * clear refusal here into an out-of-bounds read somewhere far away. */
static bool img_shape(const Expr* d, size_t* h, size_t* w, size_t* c,
                      bool* all_int, double* lo, double* hi) {
    /* A packed buffer is machine reals or machine integers by construction, so there is nothing
     * per-element to validate: every element is already numeric and non-complex. Read the range
     * straight from the buffer, which is what makes constructing an Image from a filter's output
     * O(pixels) of arithmetic rather than O(pixels) of Expr inspection. */
    if (is_ndarray(d)) {
        const NDArrayData* a = &d->data.ndarray;
        if (!img_shape_fast(d, h, w, c)) return false;
        if (a->dtype == NDT_COMPLEX64 || a->dtype == NDT_COMPLEX32) return false;
        size_t n = (*h) * (*w) * (*c);
        *all_int = (a->dtype == NDT_INT64 || a->dtype == NDT_BOOL);
        *lo = *hi = 0.0;
        for (size_t i = 0; i < n; i++) {
            double re = 0.0, imv = 0.0;
            ndt_get(a->data, i, a->dtype, &re, &imv);
            if (imv != 0.0) return false;
            if (i == 0) { *lo = *hi = re; }
            else { if (re < *lo) *lo = re; if (re > *hi) *hi = re; }
        }
        return true;
    }
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
    /* A PACKED image answers in O(1) from its dims, with no walking at all -- the point of
     * storing one. rank 2 is grey, rank 3 is height x width x channels. */
    if (is_ndarray(d)) {
        const NDArrayData* a = &d->data.ndarray;
        if (a->rank == 2) { *h = (size_t)a->dims[0]; *w = (size_t)a->dims[1]; *c = 1; return true; }
        if (a->rank == 3) {
            *h = (size_t)a->dims[0]; *w = (size_t)a->dims[1]; *c = (size_t)a->dims[2];
            return *c > 0;
        }
        return false;
    }
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
    /* If the caller handed us a PACKED LIST -- which Table and Range produce for any sizeable
     * array -- keep the buffer by storing it on the visible NDArray surface. Left as a packed
     * List it would be materialised into Expr nodes the moment this node came to rest, by the
     * same post-gate described in image_build_real, so an image built from Table[...] would
     * arrive already un-packed and every later filter would pay the walk. */
    if (two[0] && is_packed_list(two[0]))
        two[0]->data.ndarray.present_as = NDA_HEAD_NDARRAY;
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
    size_t w = 0, h = 0, dp = 0;
    /* A 3-D image reports {width, height, depth} -- three elements, and REVERSED from the
     * depth x height x width storage order. Six orderings are possible with three axes and a cubic
     * volume validates none of them, so every test uses distinct extents. */
    if (image3d_info(res->data.function.args[0], &w, &h, &dp, NULL, NULL)) {
        Expr* three[3];
        three[0] = expr_new_integer((int64_t)w);
        three[1] = expr_new_integer((int64_t)h);
        three[2] = expr_new_integer((int64_t)dp);
        if (!three[0] || !three[1] || !three[2]) {
            expr_free(three[0]); expr_free(three[1]); expr_free(three[2]); return NULL;
        }
        return expr_new_function(expr_new_symbol(SYM_List), three, 3);
    }
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
    if (image3d_info(res->data.function.args[0], NULL, NULL, NULL, &c, NULL))
        return expr_new_integer((int64_t)c);
    if (!image_info(res->data.function.args[0], NULL, NULL, &c, NULL)) return NULL;
    return expr_new_integer((int64_t)c);
}

static Expr* builtin_imagetype(Expr* res) {
    if (res->data.function.arg_count != 1) return NULL;
    ImgType t;
    if (image3d_info(res->data.function.args[0], NULL, NULL, NULL, NULL, &t))
        return expr_new_string(img_type_name(t));
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
/* Rebuild nested Lists from a flat row-major buffer of ANY rank, scaling leaves.
 *
 * The 2-D and 3-D cases are unrolled below for directness, but a colour VOLUME is rank 4 and
 * unrolling a fourth level would be the point at which the pattern should have been a recursion
 * from the start. One axis per call, with the stride being the product of the remaining dims. */
static Expr* nd_nest(const void* data, NDType dt, const int64_t* dims, int rank,
                     size_t offset, ImgType t, bool raw) {
    size_t n = (size_t)dims[0];
    Expr** kids = malloc(sizeof(Expr*) * n);
    if (!kids) return NULL;
    bool ok = true;
    for (size_t i = 0; i < n; i++) kids[i] = NULL;

    if (rank == 1) {
        for (size_t i = 0; i < n && ok; i++) {
            if (raw) {
                kids[i] = ndarray_buffer_element_to_expr(data, offset + i, dt);
            } else {
                double re = 0.0, imv = 0.0;
                ndt_get(data, offset + i, dt, &re, &imv);
                kids[i] = expr_new_real(img_to_unit(re, t));
            }
            if (!kids[i]) ok = false;
        }
    } else {
        size_t stride = 1;
        for (int r = 1; r < rank; r++) stride *= (size_t)dims[r];
        for (size_t i = 0; i < n && ok; i++) {
            kids[i] = nd_nest(data, dt, dims + 1, rank - 1, offset + i * stride, t, raw);
            if (!kids[i]) ok = false;
        }
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

static Expr* img_scale_tree(const Expr* e, ImgType t, bool raw) {
    /* A buffer-backed image: rebuild the nested List from the buffer, scaling as we go.
     *
     * ImageData must always answer with a List whatever the storage is -- a caller writing
     * Part[ImageData[img], y, x] cannot be asked to care whether the image happened to exceed the
     * packing threshold. This branch is what makes the storage an implementation detail rather
     * than a visible fork in the API, and its absence is what the first version got wrong: every
     * test used a small image, stayed on the nested path, and so never touched it. */
    if (is_ndarray(e)) {
        /* Rank 4 is a colour VOLUME. Rebuilding it as nested Lists needs a general recursion over
         * dims rather than the two-and-a-bit levels the 2-D case hard-codes, so ranks above 3
         * recurse one axis at a time. */
        const NDArrayData* a = &e->data.ndarray;
        if (a->rank > 3) {
            return nd_nest(a->data, a->dtype, a->dims, a->rank, 0, t, raw);
        }
        size_t h = (size_t)a->dims[0];
        size_t w = (a->rank >= 2) ? (size_t)a->dims[1] : 1;
        size_t c = (a->rank >= 3) ? (size_t)a->dims[2] : 1;
        Expr** rows = malloc(sizeof(Expr*) * h);
        if (!rows) return NULL;
        bool ok = true;
        for (size_t y = 0; y < h; y++) rows[y] = NULL;
        for (size_t y = 0; y < h && ok; y++) {
            Expr** cols = malloc(sizeof(Expr*) * w);
            if (!cols) { ok = false; break; }
            for (size_t x = 0; x < w; x++) cols[x] = NULL;
            for (size_t x = 0; x < w && ok; x++) {
                if (c == 1) {
                    double re = 0.0, imv = 0.0;
                    ndt_get(a->data, y * w + x, a->dtype, &re, &imv);
                    cols[x] = raw ? ndarray_buffer_element_to_expr(a->data, y * w + x, a->dtype)
                                  : expr_new_real(img_to_unit(re, t));
                } else {
                    Expr** ch = malloc(sizeof(Expr*) * c);
                    if (!ch) { ok = false; break; }
                    bool okc = true;
                    for (size_t k = 0; k < c; k++) {
                        size_t idx = (y * w + x) * c + k;
                        double re = 0.0, imv = 0.0;
                        ndt_get(a->data, idx, a->dtype, &re, &imv);
                        ch[k] = raw ? ndarray_buffer_element_to_expr(a->data, idx, a->dtype)
                                    : expr_new_real(img_to_unit(re, t));
                        if (!ch[k]) okc = false;
                    }
                    if (okc) cols[x] = expr_new_function(expr_new_symbol(SYM_List), ch, c);
                    else for (size_t k = 0; k < c; k++) expr_free(ch[k]);
                    free(ch);
                }
                if (!cols[x]) ok = false;
            }
            if (ok) rows[y] = expr_new_function(expr_new_symbol(SYM_List), cols, w);
            else for (size_t x = 0; x < w; x++) expr_free(cols[x]);
            free(cols);
            if (!rows[y]) ok = false;
        }
        if (!ok) {
            for (size_t y = 0; y < h; y++) expr_free(rows[y]);
            free(rows);
            return NULL;
        }
        Expr* out = expr_new_function(expr_new_symbol(SYM_List), rows, h);
        free(rows);
        return out;
    }
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
    /* Works for both ranks: the scaling walk is shape-agnostic, so a volume needs only the type. */
    if (!image3d_info(img, NULL, NULL, NULL, NULL, &t)
        && !image_info(img, NULL, NULL, NULL, &t)) return NULL;

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

bool image_load(const Expr* img, size_t* width, size_t* height, size_t* channels,
                double** buf) {
    size_t w = 0, h = 0, c = 0; ImgType t;
    if (!image_info(img, &w, &h, &c, &t)) return false;
    double* out = malloc(sizeof(double) * w * h * c);
    if (!out) return false;

    const Expr* d = img->data.function.args[0];

    /* THE fast path this whole change exists for: a packed image is already a flat row-major
     * buffer in exactly the order needed, so loading it is a scale over contiguous memory rather
     * than a walk over height*width*channels Expr nodes. Measurements put the walk at ~4.6 ms per
     * 262144 pixels, which dominated every filter at small kernel radii. */
    if (is_ndarray(d)) {
        const NDArrayData* a = &d->data.ndarray;
        size_t n = h * w * c;
        for (size_t i = 0; i < n; i++) {
            double re = 0.0, imv = 0.0;
            ndt_get(a->data, i, a->dtype, &re, &imv);
            if (imv != 0.0) { free(out); return false; }
            out[i] = img_to_unit(re, t);
        }
        if (width) *width = w;
        if (height) *height = h;
        if (channels) *channels = c;
        *buf = out;
        return true;
    }

    for (size_t y = 0; y < h; y++) {
        const Expr* row = d->data.function.args[y];
        for (size_t x = 0; x < w; x++) {
            const Expr* px = row->data.function.args[x];
            for (size_t k = 0; k < c; k++) {
                const Expr* v = (c == 1) ? px : px->data.function.args[k];
                double re = 0.0, im = 0.0;
                if (!na_read_scalar(v, &re, &im) || im != 0.0) { free(out); return false; }
                out[(y * w + x) * c + k] = img_to_unit(re, t);
            }
        }
    }
    if (width) *width = w;
    if (height) *height = h;
    if (channels) *channels = c;
    *buf = out;
    return true;
}

Expr* image_build_real(const double* buf, size_t width, size_t height, size_t channels) {
    if (!buf || width == 0 || height == 0 || channels == 0) return NULL;

    /* Build a PACKED array when packing is available: one memcpy-shaped loop into a machine
     * buffer instead of height*width*channels expr_new_real calls. ndbuild_open returns NULL when
     * packing is disabled or the array is under the size threshold, and the nested-Expr path below
     * then runs exactly as before -- so small images keep their old representation and every
     * existing test keeps its exact answer. */
    {
        int64_t dims[3];
        int rank;
        if (channels == 1) { rank = 2; dims[0] = (int64_t)height; dims[1] = (int64_t)width; }
        else { rank = 3; dims[0] = (int64_t)height; dims[1] = (int64_t)width;
               dims[2] = (int64_t)channels; }
        void* raw = NULL;
        Expr* nd = ndbuild_open(rank, dims, NDT_FLOAT64, &raw);
        if (nd && raw) {
            memcpy(raw, buf, sizeof(double) * width * height * channels);
            /* Present as a VISIBLE NDArray, not a packed List, and this is the crux.
             *
             * eval.c's POST-GATE materialises any packed LIST still sitting inside a node that
             * has come to rest -- unconditionally, with no aware check, because for an ordinary
             * head a resting buffer means "a fast path declined it" and an inert Mod[buffer]
             * would then behave differently from Mod[plain list]. That reasoning is right, and it
             * is fatal here: Image[...] is a CONTAINER, so coming to rest holding its pixels is
             * exactly what it is for. A packed-List image was materialised on every single
             * evaluation, which is why marking the image heads AWARE changed nothing.
             *
             * The gate never touches a visible NDArray -- the asymmetry SPEC.md documents -- so
             * that is the surface an image's storage has to use. */
            nd->data.ndarray.present_as = NDA_HEAD_NDARRAY;
            Expr* two[2];
            two[0] = nd;
            two[1] = expr_new_string("Real");
            if (!two[1]) { expr_free(nd); return NULL; }
            return expr_new_function(expr_new_symbol("Image"), two, 2);
        }
        if (nd) expr_free(nd);
    }

    Expr** rows = malloc(sizeof(Expr*) * height);
    if (!rows) return NULL;
    bool ok = true;
    for (size_t y = 0; y < height; y++) rows[y] = NULL;
    for (size_t y = 0; y < height && ok; y++) {
        Expr** cols = malloc(sizeof(Expr*) * width);
        if (!cols) { ok = false; break; }
        for (size_t x = 0; x < width; x++) cols[x] = NULL;
        for (size_t x = 0; x < width && ok; x++) {
            if (channels == 1) {
                cols[x] = expr_new_real(buf[(y * width + x)]);
            } else {
                Expr** ch = malloc(sizeof(Expr*) * channels);
                if (!ch) { ok = false; break; }
                bool okc = true;
                for (size_t k = 0; k < channels; k++) {
                    ch[k] = expr_new_real(buf[(y * width + x) * channels + k]);
                    if (!ch[k]) okc = false;
                }
                if (okc) cols[x] = expr_new_function(expr_new_symbol(SYM_List), ch, channels);
                else for (size_t k = 0; k < channels; k++) expr_free(ch[k]);
                free(ch);
            }
            if (!cols[x]) ok = false;
        }
        if (ok) rows[y] = expr_new_function(expr_new_symbol(SYM_List), cols, width);
        else for (size_t x = 0; x < width; x++) expr_free(cols[x]);
        free(cols);
        if (!rows[y]) ok = false;
    }
    if (!ok) {
        for (size_t y = 0; y < height; y++) expr_free(rows[y]);
        free(rows);
        return NULL;
    }
    Expr* data = expr_new_function(expr_new_symbol(SYM_List), rows, height);
    free(rows);
    if (!data) return NULL;
    Expr* two[2];
    two[0] = data;
    two[1] = expr_new_string("Real");
    if (!two[1]) { expr_free(data); return NULL; }
    return expr_new_function(expr_new_symbol("Image"), two, 2);
}

/* ---- Image3D: volumetric images ------------------------------------------
 *
 * A volume is depth x height x width (slices outermost, indexed data[[z, y, x]]) for grey, or
 * depth x height x width x channels for colour. ImageDimensions reports {width, height, depth} --
 * FULLY REVERSED from the storage order. That is Mathematica's convention and it is the 3-D version
 * of the trap the 2-D accessors already carry, only worse: with three axes there are six possible
 * orderings and a cubic test volume validates none of them. Every test here therefore uses
 * DISTINCT depth, height and width.
 *
 * The validation and storage follow the 2-D path exactly: a buffer answers its shape in O(1) from
 * its dims, a nested List is walked once at construction, and computed volumes are stored on the
 * visible NDArray surface so the post-gate cannot flatten them.
 */

/* Shape of a canonical 3-D image without walking the voxels. */
static bool img3_shape_fast(const Expr* d, size_t* dp, size_t* h, size_t* w, size_t* c) {
    if (is_ndarray(d)) {
        const NDArrayData* a = &d->data.ndarray;
        if (a->rank == 3) {
            *dp = (size_t)a->dims[0]; *h = (size_t)a->dims[1]; *w = (size_t)a->dims[2]; *c = 1;
            return true;
        }
        if (a->rank == 4) {
            *dp = (size_t)a->dims[0]; *h = (size_t)a->dims[1]; *w = (size_t)a->dims[2];
            *c = (size_t)a->dims[3];
            return *c > 0;
        }
        return false;
    }
    /* Nested: check every slice is rectangular and of the same shape. O(depth * height), not
     * O(voxels) -- the same trade the 2-D accessor makes, and for the same reason. */
    if (!img_is_list(d) || d->data.function.arg_count == 0) return false;
    size_t nd = d->data.function.arg_count;
    size_t nh = 0, nw = 0, nc = 1;
    for (size_t z = 0; z < nd; z++) {
        const Expr* sl = d->data.function.args[z];
        if (!img_is_list(sl) || sl->data.function.arg_count == 0) return false;
        if (z == 0) nh = sl->data.function.arg_count;
        else if (sl->data.function.arg_count != nh) return false;
        for (size_t y = 0; y < nh; y++) {
            const Expr* row = sl->data.function.args[y];
            if (!img_is_list(row) || row->data.function.arg_count == 0) return false;
            if (z == 0 && y == 0) {
                nw = row->data.function.arg_count;
                const Expr* p0 = row->data.function.args[0];
                if (img_is_list(p0)) {
                    nc = p0->data.function.arg_count;
                    if (nc == 0) return false;
                }
            } else if (row->data.function.arg_count != nw) return false;
        }
    }
    *dp = nd; *h = nh; *w = nw; *c = nc;
    return true;
}

/* Full validation: every voxel numeric and real, plus the range and integer-ness for type
 * inference. */
static bool img3_shape(const Expr* d, size_t* dp, size_t* h, size_t* w, size_t* c,
                       bool* all_int, double* lo, double* hi) {
    if (!img3_shape_fast(d, dp, h, w, c)) return false;
    *all_int = true; *lo = 0.0; *hi = 0.0;
    bool first = true;

    if (is_ndarray(d)) {
        const NDArrayData* a = &d->data.ndarray;
        if (a->dtype == NDT_COMPLEX64 || a->dtype == NDT_COMPLEX32) return false;
        size_t n = (*dp) * (*h) * (*w) * (*c);
        *all_int = (a->dtype == NDT_INT64 || a->dtype == NDT_BOOL);
        for (size_t i = 0; i < n; i++) {
            double re = 0.0, imv = 0.0;
            ndt_get(a->data, i, a->dtype, &re, &imv);
            if (imv != 0.0) return false;
            if (first) { *lo = *hi = re; first = false; }
            else { if (re < *lo) *lo = re; if (re > *hi) *hi = re; }
        }
        return true;
    }

    for (size_t z = 0; z < *dp; z++) {
        const Expr* sl = d->data.function.args[z];
        for (size_t y = 0; y < *h; y++) {
            const Expr* row = sl->data.function.args[y];
            for (size_t x = 0; x < *w; x++) {
                const Expr* px = row->data.function.args[x];
                if (*c > 1) {
                    if (!img_is_list(px) || px->data.function.arg_count != *c) return false;
                    for (size_t k = 0; k < *c; k++) {
                        const Expr* v = px->data.function.args[k];
                        double re = 0.0, imv = 0.0;
                        if (!na_read_scalar(v, &re, &imv) || imv != 0.0) return false;
                        if (v->type != EXPR_INTEGER) *all_int = false;
                        if (first) { *lo = *hi = re; first = false; }
                        else { if (re < *lo) *lo = re; if (re > *hi) *hi = re; }
                    }
                } else {
                    if (img_is_list(px)) return false;   /* mixed rank */
                    double re = 0.0, imv = 0.0;
                    if (!na_read_scalar(px, &re, &imv) || imv != 0.0) return false;
                    if (px->type != EXPR_INTEGER) *all_int = false;
                    if (first) { *lo = *hi = re; first = false; }
                    else { if (re < *lo) *lo = re; if (re > *hi) *hi = re; }
                }
            }
        }
    }
    return true;
}

bool image3d_info(const Expr* e, size_t* width, size_t* height, size_t* depth,
                  size_t* channels, ImgType* type) {
    if (!e || e->type != EXPR_FUNCTION || e->data.function.arg_count != 2) return false;
    const Expr* hd = e->data.function.head;
    if (!hd || hd->type != EXPR_SYMBOL || strcmp(hd->data.symbol.name, "Image3D") != 0)
        return false;
    const Expr* ty = e->data.function.args[1];
    if (!ty || ty->type != EXPR_STRING) return false;
    ImgType t;
    if (!img_type_from_name(ty->data.string, &t)) return false;
    size_t dp = 0, h = 0, w = 0, c = 0;
    if (!img3_shape_fast(e->data.function.args[0], &dp, &h, &w, &c)) return false;
    if (width) *width = w;
    if (height) *height = h;
    if (depth) *depth = dp;
    if (channels) *channels = c;
    if (type) *type = t;
    return true;
}

bool image3d_load(const Expr* img, size_t* width, size_t* height, size_t* depth,
                  size_t* channels, double** buf) {
    size_t w = 0, h = 0, dp = 0, c = 0; ImgType t;
    if (!image3d_info(img, &w, &h, &dp, &c, &t)) return false;
    double* out = malloc(sizeof(double) * dp * h * w * c);
    if (!out) return false;
    const Expr* d = img->data.function.args[0];

    if (is_ndarray(d)) {
        const NDArrayData* a = &d->data.ndarray;
        size_t n = dp * h * w * c;
        for (size_t i = 0; i < n; i++) {
            double re = 0.0, imv = 0.0;
            ndt_get(a->data, i, a->dtype, &re, &imv);
            if (imv != 0.0) { free(out); return false; }
            out[i] = img_to_unit(re, t);
        }
    } else {
        size_t i = 0;
        for (size_t z = 0; z < dp; z++) {
            const Expr* sl = d->data.function.args[z];
            for (size_t y = 0; y < h; y++) {
                const Expr* row = sl->data.function.args[y];
                for (size_t x = 0; x < w; x++) {
                    const Expr* px = row->data.function.args[x];
                    for (size_t k = 0; k < c; k++) {
                        const Expr* v = (c == 1) ? px : px->data.function.args[k];
                        double re = 0.0, imv = 0.0;
                        if (!na_read_scalar(v, &re, &imv) || imv != 0.0) { free(out); return false; }
                        out[i++] = img_to_unit(re, t);
                    }
                }
            }
        }
    }
    if (width) *width = w;
    if (height) *height = h;
    if (depth) *depth = dp;
    if (channels) *channels = c;
    *buf = out;
    return true;
}

/* Nested Lists from a flat double buffer of any rank -- the fallback when packing declines.
 *
 * Distinct from nd_nest, which reads an NDArray buffer through a dtype and applies unit scaling.
 * This one takes plain doubles already in unit scale, which is what a filter produces. */
static Expr* dbuf_nest(const double* buf, const int64_t* dims, int rank, size_t offset) {
    size_t n = (size_t)dims[0];
    Expr** kids = malloc(sizeof(Expr*) * n);
    if (!kids) return NULL;
    bool ok = true;
    for (size_t i = 0; i < n; i++) kids[i] = NULL;
    if (rank == 1) {
        for (size_t i = 0; i < n && ok; i++) {
            kids[i] = expr_new_real(buf[offset + i]);
            if (!kids[i]) ok = false;
        }
    } else {
        size_t stride = 1;
        for (int r = 1; r < rank; r++) stride *= (size_t)dims[r];
        for (size_t i = 0; i < n && ok; i++) {
            kids[i] = dbuf_nest(buf, dims + 1, rank - 1, offset + i * stride);
            if (!kids[i]) ok = false;
        }
    }
    if (!ok) {
        for (size_t i = 0; i < n; i++) expr_free(kids[i]);
        free(kids);
        return NULL;
    }
    Expr* o = expr_new_function(expr_new_symbol(SYM_List), kids, n);
    free(kids);
    return o;
}

Expr* image3d_build_real(const double* buf, size_t width, size_t height, size_t depth,
                         size_t channels) {
    if (!buf || width == 0 || height == 0 || depth == 0 || channels == 0) return NULL;
    int64_t dims[4];
    int rank;
    if (channels == 1) { rank = 3; dims[0] = (int64_t)depth; dims[1] = (int64_t)height;
                         dims[2] = (int64_t)width; }
    else { rank = 4; dims[0] = (int64_t)depth; dims[1] = (int64_t)height;
           dims[2] = (int64_t)width; dims[3] = (int64_t)channels; }
    void* raw = NULL;
    Expr* nd = ndbuild_open(rank, dims, NDT_FLOAT64, &raw);
    if (nd && raw) {
        memcpy(raw, buf, sizeof(double) * depth * height * width * channels);
        /* Visible surface, for the reason image_build_real documents: a packed List inside a
         * resting container is materialised by eval.c's post-gate on every evaluation. */
        nd->data.ndarray.present_as = NDA_HEAD_NDARRAY;
        Expr* two[2];
        two[0] = nd;
        two[1] = expr_new_string("Real");
        if (!two[1]) { expr_free(nd); return NULL; }
        return expr_new_function(expr_new_symbol("Image3D"), two, 2);
    }
    if (nd) expr_free(nd);

    /* NESTED FALLBACK, and its absence was a live bug.
     *
     * ndbuild_open declines an array under the packing threshold, so a small volume -- a 3-voxel
     * z-line, an 8-voxel cube -- got NULL here and the whole convolution DECLINED rather than
     * returning a nested result. The 2-D builder always had this fallback; the 3-D one was written
     * without it and only large volumes were tried, so it looked correct. Building the nested form
     * for the small case is what makes the two ranks behave alike. */
    Expr* out = dbuf_nest(buf, dims, rank, 0);
    if (!out) return NULL;
    Expr* two[2];
    two[0] = out;
    two[1] = expr_new_string("Real");
    if (!two[1]) { expr_free(out); return NULL; }
    return expr_new_function(expr_new_symbol("Image3D"), two, 2);
}

/* Image3D[data] / Image3D[data, type] -- same normalisation and inference as Image. */
static Expr* builtin_image3d(Expr* res) {
    size_t argc = res->data.function.arg_count;
    if (argc != 1 && argc != 2) return NULL;
    if (argc == 2 && image3d_info(res, NULL, NULL, NULL, NULL, NULL)) return NULL;

    Expr* data = res->data.function.args[0];
    size_t dp = 0, h = 0, w = 0, c = 0; bool all_int = false; double lo = 0.0, hi = 0.0;
    if (!img3_shape(data, &dp, &h, &w, &c, &all_int, &lo, &hi)) return NULL;

    ImgType t;
    if (argc == 2) {
        Expr* ty = res->data.function.args[1];
        if (!ty || ty->type != EXPR_STRING || !img_type_from_name(ty->data.string, &t))
            return NULL;
        if (t == IMG_BIT && !(all_int && lo >= 0.0 && hi <= 1.0)) return NULL;
        if (t == IMG_BYTE && !(all_int && lo >= 0.0 && hi <= 255.0)) return NULL;
    } else {
        if (all_int && lo >= 0.0 && hi <= 1.0)        t = IMG_BIT;
        else if (all_int && lo >= 0.0 && hi <= 255.0) t = IMG_BYTE;
        else                                          t = IMG_REAL;
    }
    Expr* two[2];
    two[0] = expr_copy(data);
    if (two[0] && is_packed_list(two[0]))
        two[0]->data.ndarray.present_as = NDA_HEAD_NDARRAY;
    two[1] = expr_new_string(img_type_name(t));
    if (!two[0] || !two[1]) { expr_free(two[0]); expr_free(two[1]); return NULL; }
    return expr_new_function(expr_new_symbol("Image3D"), two, 2);
}

static Expr* builtin_image3dq(Expr* res) {
    if (res->data.function.arg_count != 1) return NULL;
    bool q = image3d_info(res->data.function.args[0], NULL, NULL, NULL, NULL, NULL);
    return expr_new_symbol(q ? SYM_True : SYM_False);
}

void image_init(void) {
    symtab_add_builtin("Image3D", builtin_image3d);
    symtab_get_def("Image3D")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("Image3D",
        "Image3D[data] is a volumetric image, normalising to Image3D[data, type]. The data is a "
        "depth x height x width array of voxels, or depth x height x width x channels for colour, "
        "so it is indexed data[[z, y, x]] with slices outermost. ImageDimensions reports "
        "{width, height, depth} -- FULLY REVERSED from that order, which is Mathematica's "
        "convention. Type inference and the accessors match Image: ImageQ is False for a volume "
        "(use Image3DQ), while ImageDimensions, ImageChannels, ImageType and ImageData all accept "
        "either rank.");

    symtab_add_builtin("Image3DQ", builtin_image3dq);
    symtab_get_def("Image3DQ")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("Image3DQ",
        "Image3DQ[expr] gives True if expr is a valid volumetric image in canonical form. "
        "Malformed input to Image3D stays unevaluated, so this is how validity is tested.");

    symtab_set_packed_aware("Image3D");
    symtab_set_packed_aware("Image3DQ");
    /* The image heads are listed in src/pack.c's AWARE table, which is where
     * tools/check_packed_aware.py looks. Setting the flag here as well because pack_init() runs
     * before image_init() in core_init, and the table's effect must survive whatever order the
     * modules happen to initialise in. */
    symtab_set_packed_aware("Image");
    symtab_set_packed_aware("ImageData");
    symtab_set_packed_aware("ImageDimensions");
    symtab_set_packed_aware("ImageChannels");
    symtab_set_packed_aware("ImageType");
    symtab_set_packed_aware("ImageQ");
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
