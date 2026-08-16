/* imagecompose.c -- the alpha channel, and putting images together.
 *
 * Five heads that share one problem: two images almost never agree on channel count, and every
 * operation here has to decide what to do about it. The decision is made ONCE, in `promote`, and
 * every head calls it: a 1-channel image standing next to a 3-channel one is replicated into
 * three, never padded with zeros, because grey means "the same in every channel" and zero-padding
 * would turn a grey pixel red.
 *
 * COORDINATES. `ImageCompose[a, b, {x, y}]` places b's CENTRE at {x, y} in Mathematica's image
 * coordinates: x from the left, y from the BOTTOM, both in pixels. Our storage is row-major from
 * the top, so the y axis is flipped on the way in. That flip is the whole reason this is written
 * down -- a composition that silently mirrored its placement would look plausible on a symmetric
 * test image and be wrong on every real one.
 */

#include <stdlib.h>
#include <string.h>

#include "expr.h"
#include "symtab.h"
#include "attr.h"
#include "image.h"
#include "imagecompose.h"

/* ------------------------------------------------------------------ helpers */

/* One decoded image: unit-interval samples, row-major, channels innermost. */
typedef struct {
    size_t w, h, c;
    double* buf;
} Img;

static void img_free(Img* im) { free(im->buf); im->buf = NULL; }

static bool img_read(const Expr* e, Img* out)
{
    out->buf = NULL;
    return image_load(e, &out->w, &out->h, &out->c, &out->buf);
}

/* Sample of channel `k` at (x, y), replicating a grey channel across colour channels and
 * reporting alpha as 1 for an image that has none.
 *
 * `nc` is the channel count the CALLER wants to read in: for an image with `c` channels, colour
 * channel k >= c reads channel 0 (grey replication) and the alpha slot reads 1 (fully opaque). */
static double sample(const Img* im, size_t x, size_t y, size_t k, size_t colour_channels)
{
    size_t base = (y * im->w + x) * im->c;
    if (k < colour_channels) {
        /* A grey source read as colour: the same value in each channel. */
        size_t src = (im->c == 1 || im->c == 2) ? 0 : (k < im->c ? k : im->c - 1);
        return im->buf[base + src];
    }
    /* The alpha slot. 2 channels is grey+alpha, 4 is RGB+alpha; anything else has no alpha. */
    if (im->c == 2) return im->buf[base + 1];
    if (im->c == 4) return im->buf[base + 3];
    return 1.0;
}

/* How many COLOUR channels an image has, alpha excluded. */
static size_t colour_count(size_t c) { return (c == 2) ? 1 : (c == 4 ? 3 : c); }
static bool has_alpha(size_t c) { return c == 2 || c == 4; }

/* The colour channel count two images should be combined in: the larger of the two, so a grey
 * image composed with a colour one produces colour rather than discarding hue. */
static size_t promote(size_t ca, size_t cb)
{
    size_t a = colour_count(ca), b = colour_count(cb);
    return a > b ? a : b;
}

static const Expr* arg(const Expr* e, size_t i)
{
    if (e->type != EXPR_FUNCTION || i >= e->data.function.arg_count) return NULL;
    return e->data.function.args[i];
}

/* A real number argument, in [0, 1] or not -- the caller decides the range. */
static bool as_double(const Expr* e, double* out)
{
    if (!e) return false;
    if (e->type == EXPR_INTEGER) { *out = (double)e->data.integer; return true; }
    if (e->type == EXPR_REAL)    { *out = e->data.real; return true; }
    return false;
}

/* An integer pair {a, b}. */
static bool as_pair(const Expr* e, double* a, double* b)
{
    if (!e || e->type != EXPR_FUNCTION || e->data.function.arg_count != 2) return false;
    if (!e->data.function.head || e->data.function.head->type != EXPR_SYMBOL ||
        strcmp(e->data.function.head->data.symbol.name, "List") != 0) return false;
    return as_double(e->data.function.args[0], a) && as_double(e->data.function.args[1], b);
}

/* ------------------------------------------------------------ AlphaChannel */

/* AlphaChannel[image] -- the opacity, as a grey image.
 *
 * An image with no alpha channel answers with an all-opaque one rather than declining: "how
 * transparent is this?" has an answer for every image, and it is "not at all". */
static Expr* builtin_alphachannel(Expr* res)
{
    Img im;
    double* out;
    size_t n, i;
    Expr* r;

    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 1) return NULL;
    if (!img_read(arg(res, 0), &im)) return NULL;

    n = im.w * im.h;
    out = (double*)malloc(n * sizeof(double));
    if (!out) { img_free(&im); return NULL; }
    for (i = 0; i < n; i++) {
        if (im.c == 2)      out[i] = im.buf[i * 2 + 1];
        else if (im.c == 4) out[i] = im.buf[i * 4 + 3];
        else                out[i] = 1.0;
    }
    img_free(&im);
    r = image_build_real(out, im.w, im.h, 1);
    free(out);
    return r;
}

/* ---------------------------------------------------------- SetAlphaChannel */

/* SetAlphaChannel[image] / [image, a] -- attach or replace the opacity.
 *
 * `a` may be a number (one opacity everywhere) or an image (per-pixel, read as grey and required
 * to match in size -- a mismatched mask is a mistake, not something to resample silently).
 */
static Expr* builtin_setalphachannel(Expr* res)
{
    Img im, mask;
    bool have_mask = false;
    double a = 1.0;
    double* out;
    size_t cc, oc, x, y, k;
    Expr* r;

    if (res->type != EXPR_FUNCTION) return NULL;
    if (res->data.function.arg_count < 1 || res->data.function.arg_count > 2) return NULL;
    if (!img_read(arg(res, 0), &im)) return NULL;

    if (res->data.function.arg_count == 2) {
        const Expr* second = arg(res, 1);
        if (as_double(second, &a)) {
            if (a < 0.0 || a > 1.0) { img_free(&im); return NULL; }
        } else if (img_read(second, &mask)) {
            if (mask.w != im.w || mask.h != im.h) { img_free(&im); img_free(&mask); return NULL; }
            have_mask = true;
        } else {
            img_free(&im);
            return NULL;
        }
    }

    cc = colour_count(im.c);
    oc = cc + 1;                       /* the result always carries alpha */
    out = (double*)malloc(im.w * im.h * oc * sizeof(double));
    if (!out) { img_free(&im); if (have_mask) img_free(&mask); return NULL; }

    for (y = 0; y < im.h; y++) {
        for (x = 0; x < im.w; x++) {
            size_t d = (y * im.w + x) * oc;
            for (k = 0; k < cc; k++) out[d + k] = sample(&im, x, y, k, cc);
            if (have_mask) {
                /* Grey of the mask: its own channels averaged, so a colour mask is not silently
                 * read as its red channel. */
                size_t mb = (y * mask.w + x) * mask.c, j;
                double s = 0.0;
                size_t mcc = colour_count(mask.c);
                for (j = 0; j < mcc; j++) s += mask.buf[mb + j];
                out[d + cc] = s / (double)mcc;
            } else {
                out[d + cc] = a;
            }
        }
    }
    r = image_build_real(out, im.w, im.h, oc);
    free(out);
    img_free(&im);
    if (have_mask) img_free(&mask);
    return r;
}

/* ------------------------------------------------------- RemoveAlphaChannel */

/* RemoveAlphaChannel[image] / [image, background]
 *
 * Without a background the alpha is simply dropped. With one, the image is COMPOSITED over it,
 * which is the difference between "forget the transparency" and "resolve it": a half-transparent
 * white pixel over black is grey, and dropping alpha would leave it white.
 */
static Expr* builtin_removealphachannel(Expr* res)
{
    Img im;
    double bg = 0.0;
    bool composite = false;
    double* out;
    size_t cc, x, y, k;
    Expr* r;

    if (res->type != EXPR_FUNCTION) return NULL;
    if (res->data.function.arg_count < 1 || res->data.function.arg_count > 2) return NULL;
    if (!img_read(arg(res, 0), &im)) return NULL;

    if (res->data.function.arg_count == 2) {
        if (!as_double(arg(res, 1), &bg)) { img_free(&im); return NULL; }
        composite = true;
    }

    cc = colour_count(im.c);
    out = (double*)malloc(im.w * im.h * cc * sizeof(double));
    if (!out) { img_free(&im); return NULL; }
    for (y = 0; y < im.h; y++) {
        for (x = 0; x < im.w; x++) {
            size_t d = (y * im.w + x) * cc;
            double al = has_alpha(im.c) ? sample(&im, x, y, cc, cc) : 1.0;
            for (k = 0; k < cc; k++) {
                double v = sample(&im, x, y, k, cc);
                out[d + k] = composite ? v * al + bg * (1.0 - al) : v;
            }
        }
    }
    r = image_build_real(out, im.w, im.h, cc);
    free(out);
    img_free(&im);
    return r;
}

/* ------------------------------------------------------------- ImageCompose */

/* ImageCompose[base, over] / [base, over, {x, y}] / [base, {over, opacity}]
 *
 * `over` is alpha-composited onto `base`, centred by default and at {x, y} otherwise. The result
 * keeps the base's size -- composition is "draw on this", not "make something bigger" -- and
 * parts of `over` that fall outside are clipped.
 *
 * {x, y} is Mathematica's image coordinate: x from the left, y from the BOTTOM.
 */
static Expr* builtin_imagecompose(Expr* res)
{
    Img base, over;
    double opacity = 1.0;
    bool have_pos = false;
    double px = 0.0, py = 0.0;
    const Expr* second;
    double* out;
    size_t cc, oc, x, y, k;
    long ox, oy;                        /* top-left of `over` in base pixels, may be negative */
    Expr* r;

    if (res->type != EXPR_FUNCTION) return NULL;
    if (res->data.function.arg_count < 2 || res->data.function.arg_count > 3) return NULL;
    if (!img_read(arg(res, 0), &base)) return NULL;

    /* {over, opacity} -- a constant transparency for the whole overlay, which is how a caller
     * fades one image into another without building a mask. */
    second = arg(res, 1);
    if (second && second->type == EXPR_FUNCTION && second->data.function.arg_count == 2 &&
        second->data.function.head && second->data.function.head->type == EXPR_SYMBOL &&
        strcmp(second->data.function.head->data.symbol.name, "List") == 0 &&
        as_double(second->data.function.args[1], &opacity)) {
        if (!img_read(second->data.function.args[0], &over)) { img_free(&base); return NULL; }
        if (opacity < 0.0 || opacity > 1.0) { img_free(&base); img_free(&over); return NULL; }
    } else if (!img_read(second, &over)) {
        img_free(&base);
        return NULL;
    }

    if (res->data.function.arg_count == 3) {
        if (!as_pair(arg(res, 2), &px, &py)) { img_free(&base); img_free(&over); return NULL; }
        have_pos = true;
    }

    if (have_pos) {
        /* Centre at {x, y}, with y measured up from the bottom -- hence h - y. */
        ox = (long)(px - (double)over.w / 2.0 + 0.5);
        oy = (long)((double)base.h - py - (double)over.h / 2.0 + 0.5);
    } else {
        ox = (long)(((double)base.w - (double)over.w) / 2.0 + 0.5);
        oy = (long)(((double)base.h - (double)over.h) / 2.0 + 0.5);
    }

    cc = promote(base.c, over.c);
    /* The result carries alpha only if the BASE did: compositing onto an opaque image produces an
     * opaque image, and inventing an alpha channel would change the type of the thing drawn on. */
    oc = has_alpha(base.c) ? cc + 1 : cc;

    out = (double*)malloc(base.w * base.h * oc * sizeof(double));
    if (!out) { img_free(&base); img_free(&over); return NULL; }

    for (y = 0; y < base.h; y++) {
        for (x = 0; x < base.w; x++) {
            size_t d = (y * base.w + x) * oc;
            long sx = (long)x - ox, sy = (long)y - oy;
            bool inside = (sx >= 0 && sy >= 0 && (size_t)sx < over.w && (size_t)sy < over.h);
            double al = 0.0;
            if (inside)
                al = sample(&over, (size_t)sx, (size_t)sy, cc, cc) * opacity;
            for (k = 0; k < cc; k++) {
                double b = sample(&base, x, y, k, cc);
                double o = inside ? sample(&over, (size_t)sx, (size_t)sy, k, cc) : 0.0;
                out[d + k] = o * al + b * (1.0 - al);
            }
            if (oc > cc) {
                /* Alpha of a composite: the base's own opacity, raised toward opaque wherever the
                 * overlay covered it. */
                double ba = sample(&base, x, y, cc, cc);
                out[d + cc] = ba + al * (1.0 - ba);
            }
        }
    }
    r = image_build_real(out, base.w, base.h, oc);
    free(out);
    img_free(&base);
    img_free(&over);
    return r;
}

/* ------------------------------------------------------------ ImageAssemble */

/* ImageAssemble[{{a, b}, {c, d}}] / ImageAssemble[{a, b}]
 *
 * Tiles images into one. A flat list is a single ROW, which is the reading order a list already
 * implies. Rows may differ in height and columns in width; each tile is placed at its natural
 * size and any gap is left black rather than stretched, because stretching would silently resample
 * an image the caller did not ask to resize.
 */
static Expr* builtin_imageassemble(Expr* res)
{
    const Expr* grid;
    Img* tiles = NULL;
    size_t rows = 0, cols = 0, nt = 0, i, j, x, y, k;
    size_t *row_h = NULL, *col_w = NULL;
    size_t total_w = 0, total_h = 0, cc = 1, oc;
    double* out = NULL;
    Expr* r = NULL;
    bool flat, any_alpha = false;

    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 1) return NULL;
    grid = arg(res, 0);
    if (!grid || grid->type != EXPR_FUNCTION || !grid->data.function.head ||
        grid->data.function.head->type != EXPR_SYMBOL ||
        strcmp(grid->data.function.head->data.symbol.name, "List") != 0 ||
        grid->data.function.arg_count == 0)
        return NULL;

    /* Flat or nested? An element that is itself a List of images means nested. */
    {
        const Expr* first = grid->data.function.args[0];
        Img probe;
        flat = img_read(first, &probe);
        if (flat) img_free(&probe);
    }

    if (flat) {
        rows = 1;
        cols = grid->data.function.arg_count;
    } else {
        rows = grid->data.function.arg_count;
        for (i = 0; i < rows; i++) {
            const Expr* row = grid->data.function.args[i];
            if (!row || row->type != EXPR_FUNCTION || !row->data.function.head ||
                row->data.function.head->type != EXPR_SYMBOL ||
                strcmp(row->data.function.head->data.symbol.name, "List") != 0)
                return NULL;
            if (row->data.function.arg_count > cols) cols = row->data.function.arg_count;
        }
    }
    if (!rows || !cols) return NULL;

    nt = rows * cols;
    tiles = (Img*)calloc(nt, sizeof(Img));
    row_h = (size_t*)calloc(rows, sizeof(size_t));
    col_w = (size_t*)calloc(cols, sizeof(size_t));
    if (!tiles || !row_h || !col_w) goto done;

    for (i = 0; i < rows; i++) {
        const Expr* row = flat ? grid : grid->data.function.args[i];
        size_t n = row->data.function.arg_count;
        for (j = 0; j < n && j < cols; j++) {
            Img* t = &tiles[i * cols + j];
            if (!img_read(row->data.function.args[j], t)) goto done;
            if (t->h > row_h[i]) row_h[i] = t->h;
            if (t->w > col_w[j]) col_w[j] = t->w;
            if (colour_count(t->c) > cc) cc = colour_count(t->c);
            if (has_alpha(t->c)) any_alpha = true;
        }
    }
    for (i = 0; i < rows; i++) total_h += row_h[i];
    for (j = 0; j < cols; j++) total_w += col_w[j];
    if (!total_w || !total_h) goto done;

    /* Alpha survives if ANY tile had it: dropping it would make an assembled sheet of sprites
     * opaque, which is the one property such a sheet needs to keep. */
    oc = any_alpha ? cc + 1 : cc;
    out = (double*)calloc(total_w * total_h * oc, sizeof(double));
    if (!out) goto done;
    if (any_alpha) {
        /* Gaps are transparent rather than opaque black when the sheet has alpha. */
        for (i = 0; i < total_w * total_h; i++) out[i * oc + cc] = 0.0;
    }

    {
        size_t y0 = 0;
        for (i = 0; i < rows; i++) {
            size_t x0 = 0;
            for (j = 0; j < cols; j++) {
                Img* t = &tiles[i * cols + j];
                if (t->buf) {
                    for (y = 0; y < t->h; y++) {
                        for (x = 0; x < t->w; x++) {
                            size_t d = ((y0 + y) * total_w + (x0 + x)) * oc;
                            for (k = 0; k < cc; k++) out[d + k] = sample(t, x, y, k, cc);
                            if (oc > cc) out[d + cc] = sample(t, x, y, cc, cc);
                        }
                    }
                }
                x0 += col_w[j];
            }
            y0 += row_h[i];
        }
    }
    r = image_build_real(out, total_w, total_h, oc);

done:
    if (tiles) {
        for (i = 0; i < nt; i++) img_free(&tiles[i]);
        free(tiles);
    }
    free(row_h);
    free(col_w);
    free(out);
    return r;
}

void imagecompose_init(void)
{
    symtab_add_builtin("AlphaChannel", builtin_alphachannel);
    symtab_get_def("AlphaChannel")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("AlphaChannel",
        "AlphaChannel[image] gives the image's opacity as a one-channel image. An image with no "
        "alpha channel answers with an all-opaque one rather than declining: \"how transparent is "
        "this?\" has an answer for every image, and it is \"not at all\". Two channels are read as "
        "grey+alpha and four as RGB+alpha.");

    symtab_add_builtin("SetAlphaChannel", builtin_setalphachannel);
    symtab_get_def("SetAlphaChannel")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("SetAlphaChannel",
        "SetAlphaChannel[image] attaches a fully opaque alpha channel. SetAlphaChannel[image, a] "
        "sets one opacity everywhere when a is a number in [0, 1], or per pixel when a is an image "
        "of the same dimensions (read as grey, so a colour mask is not taken as its red channel "
        "alone). A mask of the wrong size is declined rather than resampled.");

    symtab_add_builtin("RemoveAlphaChannel", builtin_removealphachannel);
    symtab_get_def("RemoveAlphaChannel")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("RemoveAlphaChannel",
        "RemoveAlphaChannel[image] drops the alpha channel. RemoveAlphaChannel[image, b] instead "
        "COMPOSITES over a background of brightness b, which is the difference between forgetting "
        "the transparency and resolving it: a half-transparent white pixel over black is grey, "
        "where dropping alpha would leave it white.");

    symtab_add_builtin("ImageCompose", builtin_imagecompose);
    symtab_get_def("ImageCompose")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("ImageCompose",
        "ImageCompose[base, over] alpha-composites over onto base, centred, keeping base's size "
        "and clipping whatever falls outside. ImageCompose[base, over, {x, y}] centres the overlay "
        "at {x, y} in image coordinates -- x from the left, y from the BOTTOM. "
        "ImageCompose[base, {over, a}] scales the overlay's opacity by a. A grey image composed "
        "with a colour one produces colour: grey means the same value in every channel, so it is "
        "replicated rather than zero-padded.");

    symtab_add_builtin("ImageAssemble", builtin_imageassemble);
    symtab_get_def("ImageAssemble")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("ImageAssemble",
        "ImageAssemble[{{a, b}, {c, d}}] tiles a grid of images into one; ImageAssemble[{a, b}] "
        "makes a single row. Each tile keeps its natural size -- a row is as tall as its tallest "
        "tile and a column as wide as its widest, and any gap is left blank rather than stretched, "
        "since stretching would resample an image the caller did not ask to resize. Alpha survives "
        "if any tile had it.");
}
