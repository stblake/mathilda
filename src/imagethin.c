/* imagethin.c -- Thinning and Pruning: reducing a shape to its skeleton, and tidying it.
 *
 * Both are ITERATIVE and both are defined on a BINARY image, which is what separates them from
 * everything in imagefilter.c: a filter reads a window and writes one value, where these two delete
 * pixels in passes and stop when a pass changes nothing. A neighbourhood kernel cannot express that
 * — the result of one pass is the input to the next.
 *
 * WHY ZHANG-SUEN. It is the standard two-subiteration thinning, and its two subiterations exist for
 * a reason worth stating: deleting every deletable pixel in ONE pass severs a diagonal line, because
 * two diagonal neighbours can each be individually removable while removing both disconnects the
 * shape. Alternating the two conditions removes from opposite sides on alternating passes, which is
 * what preserves connectivity. A single-pass "delete if removable" thinning looks correct on a thick
 * blob and quietly breaks every diagonal stroke.
 *
 * THRESHOLD. A non-binary image is thresholded at 0.5 rather than at "nonzero". Nonzero is the right
 * rule for MorphologicalComponents, where the caller has usually binarised already, but for these
 * two it would make almost every grey image entirely foreground and the skeleton would be a frame
 * around the border. Callers wanting another rule should apply Binarize first, which is a decision
 * they can see.
 */

#include <stdlib.h>
#include <string.h>

#include "expr.h"
#include "symtab.h"
#include "attr.h"
#include "image.h"
#include "imagethin.h"

/* ------------------------------------------------------------------ helpers */

/* Foreground mask from an image: channels averaged, thresholded at 0.5. Caller frees. */
static bool mask_from_image(const Expr* e, size_t* w, size_t* h, unsigned char** out)
{
    double* buf = NULL;
    size_t ww = 0, hh = 0, cc = 0, i, k;
    unsigned char* m;

    if (!image_load(e, &ww, &hh, &cc, &buf)) return false;
    size_t npix = ww * hh;
    m = (unsigned char*)malloc(npix ? npix : 1);
    if (!m) { free(buf); return false; }
    for (i = 0; i < npix; i++) {
        double s = 0.0;
        for (k = 0; k < cc; k++) s += buf[i * cc + k];
        m[i] = (s / (double)cc) >= 0.5 ? 1 : 0;
    }
    free(buf);
    *w = ww; *h = hh; *out = m;
    return true;
}

/* The eight neighbours in Zhang-Suen's order: P2 is north, then clockwise. Outside the image counts
 * as background, which makes a shape touching the border thin as though it ended there. */
static void neighbours(const unsigned char* m, size_t w, size_t h,
                       size_t x, size_t y, unsigned char p[8])
{
    static const int DX[8] = { 0, 1, 1, 1, 0, -1, -1, -1 };
    static const int DY[8] = { -1, -1, 0, 1, 1, 1, 0, -1 };
    int i;
    for (i = 0; i < 8; i++) {
        long nx = (long)x + DX[i], ny = (long)y + DY[i];
        p[i] = (nx < 0 || ny < 0 || (size_t)nx >= w || (size_t)ny >= h)
             ? 0 : m[(size_t)ny * w + (size_t)nx];
    }
}

static int transitions(const unsigned char p[8])
{
    /* 0 -> 1 transitions around the ring, P2..P9,P2. */
    int i, a = 0;
    for (i = 0; i < 8; i++)
        if (!p[i] && p[(i + 1) % 8]) a++;
    return a;
}

static int neighbour_count(const unsigned char p[8])
{
    int i, b = 0;
    for (i = 0; i < 8; i++) b += p[i];
    return b;
}

/* One Zhang-Suen subiteration. `second` selects the second condition pair. Returns the number of
 * pixels deleted, which is what tells the caller whether the skeleton has settled. */
static size_t thin_pass(unsigned char* m, size_t w, size_t h, bool second)
{
    size_t npix = w * h;
    unsigned char* del = (unsigned char*)calloc(npix ? npix : 1, 1);
    size_t x, y, n = 0;
    if (!del) return 0;

    for (y = 0; y < h; y++) {
        for (x = 0; x < w; x++) {
            unsigned char p[8];
            int b, a;
            if (!m[y * w + x]) continue;
            neighbours(m, w, h, x, y, p);
            b = neighbour_count(p);
            if (b < 2 || b > 6) continue;          /* an end point or an interior pixel */
            a = transitions(p);
            if (a != 1) continue;                  /* deleting it would break connectivity */
            /* p[0]=N(P2) p[2]=E(P4) p[4]=S(P6) p[6]=W(P8) */
            if (!second) {
                if (p[0] && p[2] && p[4]) continue;
                if (p[2] && p[4] && p[6]) continue;
            } else {
                if (p[0] && p[2] && p[6]) continue;
                if (p[0] && p[4] && p[6]) continue;
            }
            del[y * w + x] = 1;
        }
    }
    /* Deleted TOGETHER, after the whole pass. Deleting in place would let a pixel's removal change
     * the verdict on its neighbour mid-pass, which is exactly the connectivity break the two
     * subiterations exist to avoid. */
    for (y = 0; y < npix; y++)
        if (del[y]) { m[y] = 0; n++; }
    free(del);
    return n;
}

/* An end point: foreground with exactly one foreground neighbour. An ISOLATED pixel has none and is
 * therefore not an end point -- pruning shortens branches, it does not erase specks, and a rule that
 * deleted isolated pixels would quietly remove every one-pixel component. */
static size_t prune_pass(unsigned char* m, size_t w, size_t h)
{
    size_t npix = w * h;
    unsigned char* del = (unsigned char*)calloc(npix ? npix : 1, 1);
    size_t x, y, i, n = 0;
    if (!del) return 0;
    for (y = 0; y < h; y++) {
        for (x = 0; x < w; x++) {
            unsigned char p[8];
            if (!m[y * w + x]) continue;
            neighbours(m, w, h, x, y, p);
            if (neighbour_count(p) == 1) del[y * w + x] = 1;
        }
    }
    for (i = 0; i < npix; i++)
        if (del[i]) { m[i] = 0; n++; }
    free(del);
    return n;
}

/* An optional non-negative integer second argument. */
static bool read_count(const Expr* res, size_t argc, long* n, bool* given)
{
    *given = false;
    if (argc < 2) return true;
    if (argc > 2) return false;
    {
        const Expr* a = res->data.function.args[1];
        if (!a || a->type != EXPR_INTEGER || a->data.integer < 0) return false;
        *n = (long)a->data.integer;
        *given = true;
    }
    return true;
}

/* ------------------------------------------------------------------ builtins */

/* Thinning[image] / Thinning[image, n] */
static Expr* builtin_thinning(Expr* res)
{
    unsigned char* m = NULL;
    size_t w = 0, h = 0;
    long limit = -1;
    bool given = false;
    Expr* out;

    if (res->type != EXPR_FUNCTION) return NULL;
    if (!read_count(res, res->data.function.arg_count, &limit, &given)) return NULL;
    if (res->data.function.arg_count < 1) return NULL;
    if (!mask_from_image(res->data.function.args[0], &w, &h, &m)) return NULL;

    {
        long iter = 0;
        for (;;) {
            size_t changed;
            if (given && iter >= limit) break;
            changed = thin_pass(m, w, h, false);
            changed += thin_pass(m, w, h, true);
            iter++;
            if (!changed) break;      /* settled: a further pass would delete nothing */
        }
    }
    out = image_build_bit(m, w, h);
    free(m);
    return out;
}

/* Pruning[image] / Pruning[image, n] */
static Expr* builtin_pruning(Expr* res)
{
    unsigned char* m = NULL;
    size_t w = 0, h = 0;
    long n = 1;
    bool given = false;
    Expr* out;
    long i;

    if (res->type != EXPR_FUNCTION) return NULL;
    if (!read_count(res, res->data.function.arg_count, &n, &given)) return NULL;
    if (res->data.function.arg_count < 1) return NULL;
    if (!mask_from_image(res->data.function.args[0], &w, &h, &m)) return NULL;

    for (i = 0; i < n; i++)
        if (!prune_pass(m, w, h)) break;    /* nothing left to prune */

    out = image_build_bit(m, w, h);
    free(m);
    return out;
}

void imagethin_init(void)
{
    symtab_add_builtin("Thinning", builtin_thinning);
    symtab_get_def("Thinning")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("Thinning",
        "Thinning[image] reduces the foreground to a one-pixel-wide skeleton by Zhang-Suen "
        "thinning, iterating until a pass deletes nothing. Thinning[image, n] stops after n "
        "iterations. The two subiterations are what preserve connectivity: deleting every "
        "individually-removable pixel in one pass severs a diagonal line, since two diagonal "
        "neighbours can each be removable while removing both disconnects the shape. A non-binary "
        "image is thresholded at 0.5 -- apply Binarize first for any other rule. The result is a "
        "\"Bit\" image, and it is always a subset of the input.");

    symtab_add_builtin("Pruning", builtin_pruning);
    symtab_get_def("Pruning")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("Pruning",
        "Pruning[image] removes one pixel from every free end of the foreground; "
        "Pruning[image, n] repeats that n times, which shortens each branch by up to n and deletes "
        "any branch shorter than that. Used after Thinning to remove the short spurs a skeleton "
        "grows at boundary irregularities. An end point has exactly one foreground neighbour, so an "
        "ISOLATED pixel is not one and survives: pruning shortens branches rather than erasing "
        "specks. Pruning[image, 0] is the image unchanged. The result is a \"Bit\" image.");
}
