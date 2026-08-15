/* imagefilter.c -- ImageConvolve, GaussianMatrix, BoxMatrix, GaussianFilter.
 *
 * CONVOLUTION, NOT CORRELATION, AND THE DIFFERENCE IS NOT PEDANTRY. Convolution reflects the
 * kernel before summing; correlation does not. For a symmetric kernel -- a Gaussian, a box --
 * they agree exactly, which is why the distinction is so easy to get wrong and so hard to
 * notice: every smoothing filter anyone tries first will look right either way. It shows up only
 * on an asymmetric kernel, where the two answers are mirror images. Mathematica draws the same
 * line, with ImageConvolve reflecting and ImageCorrelate not, so a test here pins the reflection
 * with an asymmetric kernel over a delta image, where convolution gives {1, 2, 3} and
 * correlation gives {3, 2, 1}.
 *
 * PADDING IS "Fixed": out-of-range reads clamp to the nearest edge pixel, replicating the border.
 * This is Mathematica's default for ImageConvolve, and it is the right one for smoothing: zero
 * padding would darken every edge, which looks exactly like a real vignetting bug. Clamping also
 * makes an exact test available -- a constant image convolved with a kernel summing to 1 comes
 * back as the SAME constant, everywhere including the border, which zero padding would break.
 *
 * The result is always a "Real" image. A Gaussian of bytes is not a byte, and rounding back into
 * the input type would discard precision the caller never asked to lose.
 *
 * ON SPEED, HONESTLY. This is the direct O(w * h * kw * kh) form. Two things would make it
 * dramatically faster and both are deliberately not here yet: a SEPARABLE kernel (a Gaussian is
 * an outer product of two 1-D Gaussians) turns kw * kh into kw + kh, and Accelerate's vImage --
 * already linked -- has hand-tuned SIMD convolution. The direct form lands first because it is
 * the reference the fast paths have to agree with, and there is nothing to check a fast path
 * against until it exists.
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

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ---- kernels as plain matrices ---- */

static bool ker_is_list(const Expr* e) {
    return e && e->type == EXPR_FUNCTION && e->data.function.head
        && e->data.function.head->type == EXPR_SYMBOL
        && e->data.function.head->data.symbol.name == SYM_List;
}

/* Read a rank-2 numeric matrix into a flat buffer. Rejects ragged input, for the same reason
 * the image reader does: every index below assumes rectangularity. */
static bool ker_load(const Expr* e, size_t* kh, size_t* kw, double** buf) {
    if (!ker_is_list(e) || e->data.function.arg_count == 0) return false;
    size_t h = e->data.function.arg_count;
    const Expr* r0 = e->data.function.args[0];
    if (!ker_is_list(r0) || r0->data.function.arg_count == 0) return false;
    size_t w = r0->data.function.arg_count;
    double* k = malloc(sizeof(double) * h * w);
    if (!k) return false;
    for (size_t i = 0; i < h; i++) {
        const Expr* row = e->data.function.args[i];
        if (!ker_is_list(row) || row->data.function.arg_count != w) { free(k); return false; }
        for (size_t j = 0; j < w; j++) {
            double re = 0.0, im = 0.0;
            if (!na_read_scalar(row->data.function.args[j], &re, &im) || im != 0.0) {
                free(k); return false;
            }
            k[i * w + j] = re;
        }
    }
    *kh = h; *kw = w; *buf = k;
    return true;
}

/* Clamp an index to [0, n) -- the "Fixed" padding rule. Written on int64_t because the
 * expression y - i + ci goes negative near the top edge, and doing that arithmetic in size_t
 * would wrap to an enormous positive index instead. */
static size_t clampi(int64_t v, size_t n) {
    if (v < 0) return 0;
    if ((uint64_t)v >= (uint64_t)n) return n - 1;
    return (size_t)v;
}

/* The convolution itself. `dst` and `src` are h*w*c unit-scaled buffers. */
static void convolve_planes(const double* src, double* dst,
                            size_t w, size_t h, size_t c,
                            const double* k, size_t kw, size_t kh) {
    /* Centre of the kernel. Integer division puts it at the true middle for an odd size and
     * just past the middle for an even one, which is the usual convention. */
    int64_t ci = (int64_t)(kh / 2), cj = (int64_t)(kw / 2);

    for (size_t y = 0; y < h; y++) {
        for (size_t x = 0; x < w; x++) {
            for (size_t ch = 0; ch < c; ch++) {
                double acc = 0.0;
                for (size_t i = 0; i < kh; i++) {
                    /* MINUS i, plus the centre: this is the reflection that makes it a
                     * convolution. Correlation would be `+ i - ci`. */
                    size_t sy = clampi((int64_t)y - (int64_t)i + ci, h);
                    for (size_t j = 0; j < kw; j++) {
                        size_t sx = clampi((int64_t)x - (int64_t)j + cj, w);
                        acc += k[i * kw + j] * src[(sy * w + sx) * c + ch];
                    }
                }
                dst[(y * w + x) * c + ch] = acc;
            }
        }
    }
}

/* ImageConvolve[image, kernel] */
static Expr* builtin_imageconvolve(Expr* res) {
    if (res->data.function.arg_count != 2) return NULL;
    size_t w = 0, h = 0, c = 0; double* src = NULL;
    if (!image_load(res->data.function.args[0], &w, &h, &c, &src)) return NULL;
    size_t kw = 0, kh = 0; double* k = NULL;
    if (!ker_load(res->data.function.args[1], &kh, &kw, &k)) { free(src); return NULL; }

    double* dst = malloc(sizeof(double) * w * h * c);
    Expr* out = NULL;
    if (dst) {
        convolve_planes(src, dst, w, h, c, k, kw, kh);
        out = image_build_real(dst, w, h, c);
    }
    free(src); free(k); free(dst);
    return out;
}

/* GaussianMatrix[r] -- a (2r+1) x (2r+1) Gaussian, normalised to sum 1.
 * GaussianMatrix[{r, sigma}] -- with the standard deviation stated.
 *
 * SIGMA DEFAULTS TO r/2, which is Mathematica's convention and worth stating because it is not
 * the only reasonable one: it puts the kernel's edge at two standard deviations, where the
 * Gaussian has fallen to about 13% of its peak, so truncating there loses little.
 *
 * Normalised to sum EXACTLY 1 by dividing by the realised sum rather than by the analytic
 * 2 pi sigma^2. The analytic constant is only correct for an infinite kernel; using it on a
 * truncated one leaves the sum slightly under 1, which darkens an image a little on every pass
 * -- invisible once and obvious after fifty. */
static Expr* builtin_gaussianmatrix(Expr* res) {
    if (res->data.function.arg_count != 1) return NULL;
    Expr* a = res->data.function.args[0];
    double r = 0.0, sigma = -1.0, im = 0.0;

    if (ker_is_list(a)) {
        if (a->data.function.arg_count != 2) return NULL;
        if (!na_read_scalar(a->data.function.args[0], &r, &im) || im != 0.0) return NULL;
        if (!na_read_scalar(a->data.function.args[1], &sigma, &im) || im != 0.0) return NULL;
        if (!(sigma > 0.0)) return NULL;
    } else {
        if (!na_read_scalar(a, &r, &im) || im != 0.0) return NULL;
    }
    /* A radius must be a non-negative integer: a fractional radius has no matrix size. */
    if (!(r >= 0.0) || r != floor(r) || r > 512.0) return NULL;
    size_t rr = (size_t)r;
    if (sigma < 0.0) sigma = (rr == 0) ? 1.0 : (double)rr / 2.0;

    size_t n = 2 * rr + 1;
    double* k = malloc(sizeof(double) * n * n);
    if (!k) return NULL;
    double sum = 0.0;
    for (size_t i = 0; i < n; i++) {
        double dy = (double)i - (double)rr;
        for (size_t j = 0; j < n; j++) {
            double dx = (double)j - (double)rr;
            double v = exp(-(dx * dx + dy * dy) / (2.0 * sigma * sigma));
            k[i * n + j] = v;
            sum += v;
        }
    }
    if (!(sum > 0.0)) { free(k); return NULL; }
    for (size_t i = 0; i < n * n; i++) k[i] /= sum;

    Expr** rows = malloc(sizeof(Expr*) * n);
    Expr* out = NULL;
    if (rows) {
        bool ok = true;
        for (size_t i = 0; i < n; i++) rows[i] = NULL;
        for (size_t i = 0; i < n && ok; i++) {
            Expr** cols = malloc(sizeof(Expr*) * n);
            if (!cols) { ok = false; break; }
            bool okc = true;
            for (size_t j = 0; j < n; j++) {
                cols[j] = expr_new_real(k[i * n + j]);
                if (!cols[j]) okc = false;
            }
            if (okc) rows[i] = expr_new_function(expr_new_symbol(SYM_List), cols, n);
            else for (size_t j = 0; j < n; j++) expr_free(cols[j]);
            free(cols);
            if (!rows[i]) ok = false;
        }
        if (ok) out = expr_new_function(expr_new_symbol(SYM_List), rows, n);
        else for (size_t i = 0; i < n; i++) expr_free(rows[i]);
        free(rows);
    }
    free(k);
    return out;
}

/* BoxMatrix[r] -- a (2r+1) x (2r+1) matrix of 1s.
 *
 * NOT normalised, which is Mathematica's definition and a trap worth naming: convolving with it
 * multiplies brightness by the element count, so ImageConvolve[img, BoxMatrix[1]] is nine times
 * too bright. The normalised version is a mean filter. Kept faithful rather than helpfully
 * rescaled, because a caller reaching for BoxMatrix in an arithmetic expression needs the ones. */
static Expr* builtin_boxmatrix(Expr* res) {
    if (res->data.function.arg_count != 1) return NULL;
    double r = 0.0, im = 0.0;
    if (!na_read_scalar(res->data.function.args[0], &r, &im) || im != 0.0) return NULL;
    if (!(r >= 0.0) || r != floor(r) || r > 512.0) return NULL;
    size_t n = 2 * (size_t)r + 1;

    Expr** rows = malloc(sizeof(Expr*) * n);
    if (!rows) return NULL;
    bool ok = true;
    for (size_t i = 0; i < n; i++) rows[i] = NULL;
    for (size_t i = 0; i < n && ok; i++) {
        Expr** cols = malloc(sizeof(Expr*) * n);
        if (!cols) { ok = false; break; }
        bool okc = true;
        for (size_t j = 0; j < n; j++) {
            cols[j] = expr_new_integer(1);
            if (!cols[j]) okc = false;
        }
        if (okc) rows[i] = expr_new_function(expr_new_symbol(SYM_List), cols, n);
        else for (size_t j = 0; j < n; j++) expr_free(cols[j]);
        free(cols);
        if (!rows[i]) ok = false;
    }
    Expr* out = NULL;
    if (ok) out = expr_new_function(expr_new_symbol(SYM_List), rows, n);
    else for (size_t i = 0; i < n; i++) expr_free(rows[i]);
    free(rows);
    return out;
}

/* GaussianFilter[image, r] -- exactly ImageConvolve[image, GaussianMatrix[r]].
 *
 * Implemented by building the same matrix and calling the same convolution rather than by a
 * separate code path, because Mathematica DOCUMENTS the two as equal and a test asserts it. Two
 * independent implementations of the same identity is how the identity quietly stops holding. */
static Expr* builtin_gaussianfilter(Expr* res) {
    if (res->data.function.arg_count != 2) return NULL;
    Expr* rad[1];
    rad[0] = expr_copy(res->data.function.args[1]);
    if (!rad[0]) return NULL;
    Expr* call = expr_new_function(expr_new_symbol("GaussianMatrix"), rad, 1);
    if (!call) { expr_free(rad[0]); return NULL; }
    Expr* ker = builtin_gaussianmatrix(call);
    expr_free(call);
    if (!ker) return NULL;

    size_t w = 0, h = 0, c = 0; double* src = NULL;
    if (!image_load(res->data.function.args[0], &w, &h, &c, &src)) { expr_free(ker); return NULL; }
    size_t kw = 0, kh = 0; double* k = NULL;
    if (!ker_load(ker, &kh, &kw, &k)) { expr_free(ker); free(src); return NULL; }
    expr_free(ker);

    double* dst = malloc(sizeof(double) * w * h * c);
    Expr* out = NULL;
    if (dst) {
        convolve_planes(src, dst, w, h, c, k, kw, kh);
        out = image_build_real(dst, w, h, c);
    }
    free(src); free(k); free(dst);
    return out;
}

void imagefilter_init(void) {
    symtab_add_builtin("ImageConvolve", builtin_imageconvolve);
    symtab_get_def("ImageConvolve")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("ImageConvolve",
        "ImageConvolve[image, kernel] convolves image with the rank-2 numeric kernel. This is "
        "true convolution: the kernel is REFLECTED before summing, so it differs from "
        "correlation on an asymmetric kernel (the two agree exactly on a symmetric one such as "
        "a Gaussian or a box). Out-of-range reads clamp to the nearest edge pixel, replicating "
        "the border, so a constant image convolved with a kernel summing to 1 comes back "
        "unchanged everywhere including the edges -- zero padding would darken them. The result "
        "is always a \"Real\" image of the same dimensions, since a filtered byte is not "
        "generally a byte. Each colour channel is convolved independently.");

    symtab_add_builtin("GaussianMatrix", builtin_gaussianmatrix);
    symtab_get_def("GaussianMatrix")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("GaussianMatrix",
        "GaussianMatrix[r] gives a (2r+1) x (2r+1) Gaussian matrix normalised to sum 1. "
        "GaussianMatrix[{r, sigma}] states the standard deviation; it defaults to r/2, which "
        "puts the kernel's edge at two standard deviations. Normalisation divides by the "
        "realised sum rather than the analytic 2 pi sigma^2, because the analytic constant is "
        "correct only for an infinite kernel and using it on a truncated one leaves the sum "
        "under 1 -- which darkens an image slightly on every pass.");

    symtab_add_builtin("BoxMatrix", builtin_boxmatrix);
    symtab_get_def("BoxMatrix")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("BoxMatrix",
        "BoxMatrix[r] gives a (2r+1) x (2r+1) matrix of 1s. It is NOT normalised, matching "
        "Mathematica, so ImageConvolve[image, BoxMatrix[1]] is nine times too bright; the "
        "normalised version is a mean filter. Kept faithful rather than helpfully rescaled, "
        "since a caller using BoxMatrix in arithmetic needs the ones.");

    symtab_add_builtin("GaussianFilter", builtin_gaussianfilter);
    symtab_get_def("GaussianFilter")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("GaussianFilter",
        "GaussianFilter[image, r] blurs image with a Gaussian of radius r. It is exactly "
        "ImageConvolve[image, GaussianMatrix[r]] -- the same matrix through the same "
        "convolution, not a second implementation -- and a test asserts the identity.");
}
