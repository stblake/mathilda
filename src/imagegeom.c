/* imagegeom.c -- ImageResize, and the resampling that makes it correct.
 *
 * ALIASING IS THE WHOLE PROBLEM WITH DOWNSAMPLING, and it is invisible until it is catastrophic.
 * Picking every other pixel out of a smooth photograph looks fine, so nearest-neighbour
 * downsampling passes every casual test. Hand it a fine checkerboard and it returns a flat field:
 * the pattern is sampled at exactly the frequency that makes it vanish. Nyquist says every
 * frequency above half the new sampling rate has to be REMOVED BEFORE resampling, not after, and
 * no interpolation applied afterwards can put back what point-sampling threw away.
 *
 * So the default for shrinking is AREA AVERAGING: each destination pixel is the mean of the exact
 * source region it covers, which is a box prefilter and a resample in one pass. It is not the best
 * possible antialiasing filter -- a windowed sinc or a Gaussian prefilter has better stopband
 * behaviour -- but it is the one that is exactly right for integer reduction factors, cheap, and
 * impossible to get subtly wrong.
 *
 * That difference is what the tests pin. A 4x4 alternating 0/1 pattern reduced to 2x2 gives a
 * constant 0.5 under area averaging and a constant 0 or 1 under nearest -- exact values, and no way
 * for a broken implementation to land between them.
 *
 * FRACTIONAL COVERAGE, not integer blocks. The area path weights each source pixel by how much of
 * it the destination pixel actually overlaps, so a 3 -> 2 reduction is handled as correctly as
 * 4 -> 2. Restricting it to integer factors would have been simpler and would have quietly fallen
 * back to something worse on the sizes people actually ask for.
 *
 * PIXEL CENTRES, not corners. A destination pixel i covers source coordinates
 * [i * s, (i + 1) * s) and its centre sits at (i + 0.5) * s, so the bilinear map is
 * sx = (i + 0.5) * s - 0.5. The naive sx = i * s is the classic half-pixel shift: it looks right
 * at 1:1, drifts the image half a pixel at any other scale, and is asymmetric -- the left edge
 * gains a border the right edge does not.
 */
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "expr.h"
#include "symtab.h"
#include "attr.h"
#include "sym_names.h"
#include "linalg/numarray.h"
#include "eval.h"
#include "image.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

typedef enum { RS_AUTO, RS_NEAREST, RS_BILINEAR, RS_AVERAGE } Resample;

static size_t clamp_idx(int64_t v, size_t n) {
    if (v < 0) return 0;
    if ((uint64_t)v >= (uint64_t)n) return n - 1;
    return (size_t)v;
}

static void rs_nearest(const double* src, double* dst, size_t sw, size_t sh,
                       size_t dw, size_t dh, size_t c) {
    double xs = (double)sw / (double)dw, ys = (double)sh / (double)dh;
    for (size_t y = 0; y < dh; y++) {
        size_t sy = clamp_idx((int64_t)floor(((double)y + 0.5) * ys), sh);
        for (size_t x = 0; x < dw; x++) {
            size_t sx = clamp_idx((int64_t)floor(((double)x + 0.5) * xs), sw);
            for (size_t k = 0; k < c; k++)
                dst[(y * dw + x) * c + k] = src[(sy * sw + sx) * c + k];
        }
    }
}

static void rs_bilinear(const double* src, double* dst, size_t sw, size_t sh,
                        size_t dw, size_t dh, size_t c) {
    double xs = (double)sw / (double)dw, ys = (double)sh / (double)dh;
    for (size_t y = 0; y < dh; y++) {
        /* Centre-aligned, so the image does not drift half a pixel at scales other than 1:1. */
        double fy = ((double)y + 0.5) * ys - 0.5;
        int64_t y0 = (int64_t)floor(fy);
        double wy = fy - (double)y0;
        size_t ya = clamp_idx(y0, sh), yb = clamp_idx(y0 + 1, sh);
        for (size_t x = 0; x < dw; x++) {
            double fx = ((double)x + 0.5) * xs - 0.5;
            int64_t x0 = (int64_t)floor(fx);
            double wx = fx - (double)x0;
            size_t xa = clamp_idx(x0, sw), xb = clamp_idx(x0 + 1, sw);
            for (size_t k = 0; k < c; k++) {
                double p00 = src[(ya * sw + xa) * c + k];
                double p01 = src[(ya * sw + xb) * c + k];
                double p10 = src[(yb * sw + xa) * c + k];
                double p11 = src[(yb * sw + xb) * c + k];
                double top = p00 + (p01 - p00) * wx;
                double bot = p10 + (p11 - p10) * wx;
                dst[(y * dw + x) * c + k] = top + (bot - top) * wy;
            }
        }
    }
}

/* Area averaging with exact fractional coverage.
 *
 * Destination pixel i spans source [i * s, (i + 1) * s). Every source pixel overlapping that span
 * contributes in proportion to the overlap LENGTH, so the weights sum to the span width and the
 * normalised result is a true mean over the covered region. For an integer factor the fractional
 * parts vanish and this reduces to an exact block mean -- which is what makes the 4 -> 2 test an
 * equality rather than a tolerance. */
static void rs_average(const double* src, double* dst, size_t sw, size_t sh,
                       size_t dw, size_t dh, size_t c) {
    double xs = (double)sw / (double)dw, ys = (double)sh / (double)dh;
    for (size_t y = 0; y < dh; y++) {
        double y_lo = (double)y * ys, y_hi = y_lo + ys;
        size_t iy0 = (size_t)floor(y_lo);
        size_t iy1 = (size_t)ceil(y_hi);
        if (iy1 > sh) iy1 = sh;
        if (iy0 >= sh) iy0 = sh - 1;
        for (size_t x = 0; x < dw; x++) {
            double x_lo = (double)x * xs, x_hi = x_lo + xs;
            size_t ix0 = (size_t)floor(x_lo);
            size_t ix1 = (size_t)ceil(x_hi);
            if (ix1 > sw) ix1 = sw;
            if (ix0 >= sw) ix0 = sw - 1;

            for (size_t k = 0; k < c; k++) {
                double acc = 0.0, wsum = 0.0;
                for (size_t sy = iy0; sy < iy1; sy++) {
                    double top = (double)sy > y_lo ? (double)sy : y_lo;
                    double bot = (double)(sy + 1) < y_hi ? (double)(sy + 1) : y_hi;
                    double why = bot - top;
                    if (!(why > 0.0)) continue;
                    for (size_t sx = ix0; sx < ix1; sx++) {
                        double lft = (double)sx > x_lo ? (double)sx : x_lo;
                        double rgt = (double)(sx + 1) < x_hi ? (double)(sx + 1) : x_hi;
                        double whx = rgt - lft;
                        if (!(whx > 0.0)) continue;
                        double wgt = whx * why;
                        acc += wgt * src[(sy * sw + sx) * c + k];
                        wsum += wgt;
                    }
                }
                dst[(y * dw + x) * c + k] = (wsum > 0.0) ? acc / wsum : 0.0;
            }
        }
    }
}

/* Parse Resampling -> "..." from the trailing options. */
static bool rs_parse_option(Expr* o, Resample* out) {
    if (!o || o->type != EXPR_FUNCTION || o->data.function.arg_count != 2) return false;
    Expr* h = o->data.function.head;
    if (!h || h->type != EXPR_SYMBOL) return false;
    if (h->data.symbol.name != SYM_Rule && h->data.symbol.name != SYM_RuleDelayed) return false;
    Expr* lhs = o->data.function.args[0];
    if (!lhs || lhs->type != EXPR_SYMBOL
        || strcmp(lhs->data.symbol.name, "Resampling") != 0) return false;
    Expr* rhs = o->data.function.args[1];
    if (rhs && rhs->type == EXPR_SYMBOL && strcmp(rhs->data.symbol.name, "Automatic") == 0) {
        *out = RS_AUTO; return true;
    }
    if (!rhs || rhs->type != EXPR_STRING) return false;
    const char* s = rhs->data.string;
    if (strcmp(s, "Nearest") == 0)  { *out = RS_NEAREST;  return true; }
    if (strcmp(s, "Bilinear") == 0) { *out = RS_BILINEAR; return true; }
    if (strcmp(s, "Average") == 0 || strcmp(s, "Mean") == 0) { *out = RS_AVERAGE; return true; }
    return false;
}

/* ---- volumetric resampling -------------------------------------------------
 *
 * The 2-D reasoning carries over unchanged, and so does the trap: aliasing is what makes area
 * averaging the right default when shrinking, and a 3-D checkerboard is what proves it. A
 * 2x2x2-periodic pattern halved on every axis comes back a constant 0.5 under area averaging and a
 * FLAT FIELD under nearest -- the pattern sampled at exactly the frequency that annihilates it.
 *
 * Fractional coverage in three axes: a destination voxel's weight for a source voxel is the product
 * of the three per-axis overlaps, so an integer reduction factor gives an exact block mean and a
 * non-integer one is handled correctly rather than approximated.
 *
 * ONE THING GETS HARDER IN 3-D, and it is not the arithmetic. The size specification is
 * {width, height, depth} while the storage is depth x height x width -- fully reversed -- so a resize
 * has to reverse the spec before it indexes anything. A CUBIC test volume validates none of that, and
 * neither does a cubic TARGET, so the tests resize a non-cubic volume to a non-cubic size.
 */
static void rs3_nearest(const double* src, double* dst,
                        size_t sw, size_t sh, size_t sd,
                        size_t dw, size_t dh, size_t dd, size_t c) {
    double xs = (double)sw / (double)dw, ys = (double)sh / (double)dh,
           zs = (double)sd / (double)dd;
    for (size_t z = 0; z < dd; z++) {
        size_t sz = clamp_idx((int64_t)floor(((double)z + 0.5) * zs), sd);
        for (size_t y = 0; y < dh; y++) {
            size_t sy = clamp_idx((int64_t)floor(((double)y + 0.5) * ys), sh);
            for (size_t x = 0; x < dw; x++) {
                size_t sx = clamp_idx((int64_t)floor(((double)x + 0.5) * xs), sw);
                for (size_t k = 0; k < c; k++)
                    dst[((z * dh + y) * dw + x) * c + k] =
                        src[((sz * sh + sy) * sw + sx) * c + k];
            }
        }
    }
}

static void rs3_trilinear(const double* src, double* dst,
                          size_t sw, size_t sh, size_t sd,
                          size_t dw, size_t dh, size_t dd, size_t c) {
    double xs = (double)sw / (double)dw, ys = (double)sh / (double)dh,
           zs = (double)sd / (double)dd;
    for (size_t z = 0; z < dd; z++) {
        double fz = ((double)z + 0.5) * zs - 0.5;
        int64_t z0 = (int64_t)floor(fz);
        double wz = fz - (double)z0;
        size_t za = clamp_idx(z0, sd), zb = clamp_idx(z0 + 1, sd);
        for (size_t y = 0; y < dh; y++) {
            double fy = ((double)y + 0.5) * ys - 0.5;
            int64_t y0 = (int64_t)floor(fy);
            double wy = fy - (double)y0;
            size_t ya = clamp_idx(y0, sh), yb = clamp_idx(y0 + 1, sh);
            for (size_t x = 0; x < dw; x++) {
                double fx = ((double)x + 0.5) * xs - 0.5;
                int64_t x0 = (int64_t)floor(fx);
                double wx = fx - (double)x0;
                size_t xa = clamp_idx(x0, sw), xb = clamp_idx(x0 + 1, sw);
                for (size_t k = 0; k < c; k++) {
                    /* Eight corners, collapsed axis by axis: x, then y, then z. */
                    double c000 = src[((za * sh + ya) * sw + xa) * c + k];
                    double c001 = src[((za * sh + ya) * sw + xb) * c + k];
                    double c010 = src[((za * sh + yb) * sw + xa) * c + k];
                    double c011 = src[((za * sh + yb) * sw + xb) * c + k];
                    double c100 = src[((zb * sh + ya) * sw + xa) * c + k];
                    double c101 = src[((zb * sh + ya) * sw + xb) * c + k];
                    double c110 = src[((zb * sh + yb) * sw + xa) * c + k];
                    double c111 = src[((zb * sh + yb) * sw + xb) * c + k];
                    double a0 = c000 + (c001 - c000) * wx;
                    double a1 = c010 + (c011 - c010) * wx;
                    double b0 = c100 + (c101 - c100) * wx;
                    double b1 = c110 + (c111 - c110) * wx;
                    double p = a0 + (a1 - a0) * wy;
                    double q = b0 + (b1 - b0) * wy;
                    dst[((z * dh + y) * dw + x) * c + k] = p + (q - p) * wz;
                }
            }
        }
    }
}

static void rs3_average(const double* src, double* dst,
                        size_t sw, size_t sh, size_t sd,
                        size_t dw, size_t dh, size_t dd, size_t c) {
    double xs = (double)sw / (double)dw, ys = (double)sh / (double)dh,
           zs = (double)sd / (double)dd;
    for (size_t z = 0; z < dd; z++) {
        double z_lo = (double)z * zs, z_hi = z_lo + zs;
        size_t iz0 = (size_t)floor(z_lo), iz1 = (size_t)ceil(z_hi);
        if (iz1 > sd) iz1 = sd;
        if (iz0 >= sd) iz0 = sd - 1;
        for (size_t y = 0; y < dh; y++) {
            double y_lo = (double)y * ys, y_hi = y_lo + ys;
            size_t iy0 = (size_t)floor(y_lo), iy1 = (size_t)ceil(y_hi);
            if (iy1 > sh) iy1 = sh;
            if (iy0 >= sh) iy0 = sh - 1;
            for (size_t x = 0; x < dw; x++) {
                double x_lo = (double)x * xs, x_hi = x_lo + xs;
                size_t ix0 = (size_t)floor(x_lo), ix1 = (size_t)ceil(x_hi);
                if (ix1 > sw) ix1 = sw;
                if (ix0 >= sw) ix0 = sw - 1;
                for (size_t k = 0; k < c; k++) {
                    double acc = 0.0, wsum = 0.0;
                    for (size_t sz = iz0; sz < iz1; sz++) {
                        double zt = (double)sz > z_lo ? (double)sz : z_lo;
                        double zb2 = (double)(sz + 1) < z_hi ? (double)(sz + 1) : z_hi;
                        double whz = zb2 - zt;
                        if (!(whz > 0.0)) continue;
                        for (size_t sy = iy0; sy < iy1; sy++) {
                            double yt = (double)sy > y_lo ? (double)sy : y_lo;
                            double yb2 = (double)(sy + 1) < y_hi ? (double)(sy + 1) : y_hi;
                            double why = yb2 - yt;
                            if (!(why > 0.0)) continue;
                            for (size_t sx = ix0; sx < ix1; sx++) {
                                double xl = (double)sx > x_lo ? (double)sx : x_lo;
                                double xr = (double)(sx + 1) < x_hi ? (double)(sx + 1) : x_hi;
                                double whx = xr - xl;
                                if (!(whx > 0.0)) continue;
                                double wgt = whx * why * whz;
                                acc += wgt * src[((sz * sh + sy) * sw + sx) * c + k];
                                wsum += wgt;
                            }
                        }
                    }
                    dst[((z * dh + y) * dw + x) * c + k] = (wsum > 0.0) ? acc / wsum : 0.0;
                }
            }
        }
    }
}

/* Resize a volume. `spec` is {width, height, depth} or a single width. */
static Expr* resize3_run(Expr* vol, Expr* spec, Resample how) {
    size_t sw = 0, sh = 0, sd = 0, c = 0; double* src = NULL;
    if (!image3d_load(vol, &sw, &sh, &sd, &c, &src)) return NULL;

    double a = 0.0, b = 0.0, d3 = 0.0, im = 0.0;
    bool full = false;
    if (spec && spec->type == EXPR_FUNCTION && spec->data.function.head
        && spec->data.function.head->type == EXPR_SYMBOL
        && spec->data.function.head->data.symbol.name == SYM_List) {
        if (spec->data.function.arg_count != 3) { free(src); return NULL; }
        if (!na_read_scalar(spec->data.function.args[0], &a, &im) || im != 0.0
         || !na_read_scalar(spec->data.function.args[1], &b, &im) || im != 0.0
         || !na_read_scalar(spec->data.function.args[2], &d3, &im) || im != 0.0) {
            free(src); return NULL;
        }
        full = true;
    } else {
        if (!na_read_scalar(spec, &a, &im) || im != 0.0) { free(src); return NULL; }
    }
    if (!(a >= 1.0) || a != floor(a) || a > 4096.0) { free(src); return NULL; }
    if (full && (!(b >= 1.0) || b != floor(b) || b > 4096.0
              || !(d3 >= 1.0) || d3 != floor(d3) || d3 > 4096.0)) { free(src); return NULL; }

    /* The spec is {width, height, depth}; the buffers are depth x height x width. Reversing here,
     * once, rather than at each index is what keeps the three loops readable -- and a non-cubic
     * TARGET is what a test needs to catch getting it wrong. */
    size_t dw = (size_t)a;
    size_t dh = full ? (size_t)b : (size_t)(((double)sh * (double)dw / (double)sw) + 0.5);
    size_t dd = full ? (size_t)d3 : (size_t)(((double)sd * (double)dw / (double)sw) + 0.5);
    if (dh == 0) dh = 1;
    if (dd == 0) dd = 1;

    Resample use = how;
    if (use == RS_AUTO) use = (dw < sw || dh < sh || dd < sd) ? RS_AVERAGE : RS_BILINEAR;

    double* dst = malloc(sizeof(double) * dw * dh * dd * c);
    Expr* out = NULL;
    if (dst) {
        if (use == RS_NEAREST)       rs3_nearest(src, dst, sw, sh, sd, dw, dh, dd, c);
        else if (use == RS_BILINEAR) rs3_trilinear(src, dst, sw, sh, sd, dw, dh, dd, c);
        else                         rs3_average(src, dst, sw, sh, sd, dw, dh, dd, c);
        out = image3d_build_real(dst, dw, dh, dd, c);
    }
    free(src); free(dst);
    return out;
}

/* ImageResize[image, {w, h}] / ImageResize[image, w] / with Resampling -> method.
 *
 * A single number is a WIDTH, and the height follows to preserve the aspect ratio -- Mathematica's
 * convention. Rounding the derived height rather than truncating keeps a 3:2 image from becoming
 * 3:1.999 on the way down. */
static Expr* builtin_imageresize(Expr* res) {
    size_t argc = res->data.function.arg_count;
    if (argc < 2) return NULL;

    Resample how = RS_AUTO;
    for (size_t i = 2; i < argc; i++)
        if (!rs_parse_option(res->data.function.args[i], &how)) return NULL;

    /* A VOLUME takes the rank-3 path. Dispatching on the image rather than on the spec, so a
     * two-element spec handed to a volume is a mistake to decline rather than a plane to guess at. */
    if (image3d_info(res->data.function.args[0], NULL, NULL, NULL, NULL, NULL))
        return resize3_run(res->data.function.args[0], res->data.function.args[1], how);

    size_t sw = 0, sh = 0, c = 0; double* src = NULL;
    if (!image_load(res->data.function.args[0], &sw, &sh, &c, &src)) return NULL;

    /* Target size. */
    Expr* spec = res->data.function.args[1];
    double dwd = 0.0, dhd = 0.0, im = 0.0;
    bool have_h = false;
    if (spec && spec->type == EXPR_FUNCTION && spec->data.function.head
        && spec->data.function.head->type == EXPR_SYMBOL
        && spec->data.function.head->data.symbol.name == SYM_List) {
        if (spec->data.function.arg_count != 2) { free(src); return NULL; }
        if (!na_read_scalar(spec->data.function.args[0], &dwd, &im) || im != 0.0) {
            free(src); return NULL;
        }
        if (!na_read_scalar(spec->data.function.args[1], &dhd, &im) || im != 0.0) {
            free(src); return NULL;
        }
        have_h = true;
    } else {
        if (!na_read_scalar(spec, &dwd, &im) || im != 0.0) { free(src); return NULL; }
    }
    /* Sizes must be positive integers: a fractional pixel count has no meaning, and rounding one
     * silently would make ImageResize[img, 10.5] quietly mean something the caller did not say. */
    if (!(dwd >= 1.0) || dwd != floor(dwd) || dwd > 1e6) { free(src); return NULL; }
    if (have_h && (!(dhd >= 1.0) || dhd != floor(dhd) || dhd > 1e6)) { free(src); return NULL; }
    size_t dw = (size_t)dwd;
    size_t dh = have_h ? (size_t)dhd
                       : (size_t)(( (double)sh * (double)dw / (double)sw ) + 0.5);
    if (dh == 0) dh = 1;

    /* Automatic: AREA AVERAGING when either axis shrinks, bilinear otherwise.
     *
     * Shrinking is where aliasing lives, and area averaging is the choice that cannot alias.
     * Enlarging has no frequencies to remove, so bilinear -- smooth and cheap -- is right there;
     * area averaging on an enlargement would degenerate to nearest, since each destination pixel
     * would fall inside a single source pixel. */
    Resample use = how;
    if (use == RS_AUTO) use = (dw < sw || dh < sh) ? RS_AVERAGE : RS_BILINEAR;

    double* dst = malloc(sizeof(double) * dw * dh * c);
    Expr* out = NULL;
    if (dst) {
        if (use == RS_NEAREST)       rs_nearest(src, dst, sw, sh, dw, dh, c);
        else if (use == RS_BILINEAR) rs_bilinear(src, dst, sw, sh, dw, dh, c);
        else                         rs_average(src, dst, sw, sh, dw, dh, c);
        out = image_build_real(dst, dw, dh, c);
    }
    free(src); free(dst);
    return out;
}

/* ---- rotation and reflection ----------------------------------------------
 *
 * RIGHT ANGLES ARE A PURE INDEX PERMUTATION, and keeping them on that path rather than sending them
 * through the general resampler is the whole design. A quarter turn moves every pixel to another pixel's
 * exact position -- no interpolation, no rounding, nothing lost -- so four of them are EXACTLY the
 * identity, and a test asserts that with `===` rather than a tolerance. Route 90 degrees through a
 * bilinear resampler and the identity becomes approximate for no reason at all: the sample points land
 * on pixel centres, so the interpolation weights are 1 and 0, but the arithmetic still runs and the
 * half-pixel convention still has to be exactly right for it to come out clean.
 *
 * A quarter turn also SWAPS THE DIMENSIONS, which is the other thing a test has to pin, and it needs a
 * non-square image to say anything at all.
 *
 * ARBITRARY ANGLES cannot be exact -- a rotated pixel grid does not land on a pixel grid -- so they
 * interpolate, and the honest test is different in kind: rotating by theta and then by -theta returns
 * the image approximately, not exactly, and the interior is what recovers while the corners have rotated
 * out of frame and back as background. Testing the interior and saying so is better than a loose
 * whole-image tolerance that hides where the error actually is.
 *
 * Out-of-frame samples read as 0 (black) rather than replicating the edge. Replication is right for a
 * FILTER, where the border is a boundary condition on an operation happening inside the image; it is
 * wrong for a ROTATION, where the area outside genuinely was not photographed and smearing the edge
 * across it invents content.
 */
static void rot90_run(const double* src, double* dst, size_t w, size_t h, size_t c, int quarters) {
    /* quarters is 1, 2 or 3. Destination dimensions are swapped for the odd turns; the caller has
     * already allocated accordingly. */
    for (size_t y = 0; y < h; y++)
        for (size_t x = 0; x < w; x++)
            for (size_t k = 0; k < c; k++) {
                size_t dx, dy, dw;
                if (quarters == 1)      { dx = h - 1 - y; dy = x;         dw = h; }
                else if (quarters == 2) { dx = w - 1 - x; dy = h - 1 - y; dw = w; }
                else                    { dx = y;         dy = w - 1 - x; dw = h; }
                dst[(dy * dw + dx) * c + k] = src[(y * w + x) * c + k];
            }
}

/* Bilinear rotation about the image centre, sampling the SOURCE for each destination pixel.
 *
 * Inverse mapping, not forward: iterating over destination pixels and asking where each came from
 * fills every output exactly once. Forward mapping -- iterating over the source and writing where each
 * pixel lands -- leaves holes wherever the rotation stretches, which is the classic artefact and is why
 * inverse mapping is universal. */
static void rot_free_run(const double* src, double* dst, size_t w, size_t h, size_t c, double rad) {
    double ca = cos(rad), sa = sin(rad);
    double cx = 0.5 * (double)(w - 1), cy = 0.5 * (double)(h - 1);
    for (size_t y = 0; y < h; y++)
        for (size_t x = 0; x < w; x++) {
            double dx = (double)x - cx, dy = (double)y - cy;
            /* Rotate the destination offset BACKWARDS to find the source point. */
            double sxf =  ca * dx + sa * dy + cx;
            double syf = -sa * dx + ca * dy + cy;
            int64_t x0 = (int64_t)floor(sxf), y0 = (int64_t)floor(syf);
            double fx = sxf - (double)x0, fy = syf - (double)y0;
            for (size_t k = 0; k < c; k++) {
                /* Outside the frame is 0, not the clamped edge: that area was never photographed. */
                double p00 = 0.0, p01 = 0.0, p10 = 0.0, p11 = 0.0;
                if (x0 >= 0 && (size_t)x0 < w && y0 >= 0 && (size_t)y0 < h)
                    p00 = src[((size_t)y0 * w + (size_t)x0) * c + k];
                if (x0 + 1 >= 0 && (size_t)(x0 + 1) < w && y0 >= 0 && (size_t)y0 < h)
                    p01 = src[((size_t)y0 * w + (size_t)(x0 + 1)) * c + k];
                if (x0 >= 0 && (size_t)x0 < w && y0 + 1 >= 0 && (size_t)(y0 + 1) < h)
                    p10 = src[((size_t)(y0 + 1) * w + (size_t)x0) * c + k];
                if (x0 + 1 >= 0 && (size_t)(x0 + 1) < w && y0 + 1 >= 0 && (size_t)(y0 + 1) < h)
                    p11 = src[((size_t)(y0 + 1) * w + (size_t)(x0 + 1)) * c + k];
                double top = p00 + (p01 - p00) * fx;
                double bot = p10 + (p11 - p10) * fx;
                dst[(y * w + x) * c + k] = top + (bot - top) * fy;
            }
        }
}

/* ImageRotate[image] (a quarter turn) / [image, angle] / [image, n Degree] */
static Expr* builtin_imagerotate(Expr* res) {
    size_t argc = res->data.function.arg_count;
    if (argc != 1 && argc != 2) return NULL;

    double rad = M_PI / 2.0;               /* Mathematica's default is a quarter turn */
    if (argc == 2) {
        /* NUMERICALISE the angle before reading it. `Pi`, `Pi/2` and `90 Degree` are exact symbolic
         * values, not machine reals, so na_read_scalar refuses them -- ImageRotate[img, Pi] declined
         * outright until this was added, which is a poor answer to the most natural way of writing a
         * half turn. Wrapping in N and evaluating handles every exact form at once rather than
         * special-casing Pi. */
        Expr* one[1];
        one[0] = expr_copy(res->data.function.args[1]);
        if (!one[0]) return NULL;
        Expr* nexpr = expr_new_function(expr_new_symbol("N"), one, 1);
        if (!nexpr) { expr_free(one[0]); return NULL; }
        Expr* num = evaluate(nexpr);
        double v = 0.0, im = 0.0;
        bool ok = num && na_read_scalar(num, &v, &im) && im == 0.0;
        expr_free(num);
        if (!ok) return NULL;
        rad = v;
    }
    size_t w = 0, h = 0, c = 0; double* src = NULL;
    if (!image_load(res->data.function.args[0], &w, &h, &c, &src)) return NULL;

    /* Is the angle a multiple of a right angle, to within the rounding of Pi/2 arithmetic? If so take
     * the exact permutation path. The tolerance is on the ANGLE, not on the pixels: a caller writing
     * Pi/2 gets the exact path, and one writing 1.5707963 gets the resampler, which is the honest
     * reading of what they asked for. */
    double q = rad / (M_PI / 2.0);
    double qr = floor(q + 0.5);
    bool right = fabs(q - qr) < 1e-9;

    Expr* out = NULL;
    if (right) {
        int quarters = (int)(((int64_t)qr % 4 + 4) % 4);
        if (quarters == 0) {
            out = image_build_real(src, w, h, c);      /* full turn: the image itself */
        } else {
            size_t dw = (quarters == 2) ? w : h;
            size_t dh = (quarters == 2) ? h : w;
            double* dst = malloc(sizeof(double) * dw * dh * c);
            if (dst) {
                rot90_run(src, dst, w, h, c, quarters);
                out = image_build_real(dst, dw, dh, c);
            }
            free(dst);
        }
    } else {
        double* dst = malloc(sizeof(double) * w * h * c);
        if (dst) {
            rot_free_run(src, dst, w, h, c, rad);
            out = image_build_real(dst, w, h, c);
        }
        free(dst);
    }
    free(src);
    return out;
}

/* ImageReflect[image] / [image, side] -- a pure index permutation, hence exact and self-inverse. */
static Expr* builtin_imagereflect(Expr* res) {
    size_t argc = res->data.function.arg_count;
    if (argc != 1 && argc != 2) return NULL;

    /* Mathematica's default reflects top-to-bottom. The named sides are given as symbols. */
    int mode = 0;                          /* 0 = vertical flip, 1 = horizontal, 2 = both */
    if (argc == 2) {
        Expr* sd = res->data.function.args[1];
        if (!sd || sd->type != EXPR_SYMBOL) return NULL;
        const char* n = sd->data.symbol.name;
        if      (strcmp(n, "Top") == 0 || strcmp(n, "Bottom") == 0) mode = 0;
        else if (strcmp(n, "Left") == 0 || strcmp(n, "Right") == 0) mode = 1;
        else return NULL;
    }
    size_t w = 0, h = 0, c = 0; double* src = NULL;
    if (!image_load(res->data.function.args[0], &w, &h, &c, &src)) return NULL;
    double* dst = malloc(sizeof(double) * w * h * c);
    Expr* out = NULL;
    if (dst) {
        for (size_t y = 0; y < h; y++)
            for (size_t x = 0; x < w; x++) {
                size_t sy = (mode == 1) ? y : h - 1 - y;
                size_t sx = (mode == 0) ? x : w - 1 - x;
                for (size_t k = 0; k < c; k++)
                    dst[(y * w + x) * c + k] = src[(sy * w + sx) * c + k];
            }
        out = image_build_real(dst, w, h, c);
    }
    free(src); free(dst);
    return out;
}

/* ---- padding and cropping -------------------------------------------------
 *
 * PAD THEN CROP IS THE IDENTITY, exactly, and that is the property the pair exists to satisfy. Both are
 * index arithmetic with no interpolation, so cropping back to the original size after padding returns
 * the original pixels bit for bit -- asserted with `===`. It is also the test that catches an
 * off-by-one on either side independently: a pad that adds one row too many at the top and one too few
 * at the bottom still has the right total size, and only a round trip notices.
 *
 * THE VERTICAL CONVENTION IS THE TRAP, as everywhere else in this subsystem. Mathematica's
 * ImagePad[image, {{left, right}, {bottom, top}}] names the pair in VISUAL order -- bottom before top --
 * while the data's first row is the TOP of the image. So `bottom` padding adds rows at the END of the
 * array and `top` at the beginning, which is the reverse of how the spec reads. A test pads
 * asymmetrically and checks which end grew, since a symmetric pad cannot tell.
 *
 * PADDING MODES. A constant fills with that value. "Fixed" replicates the edge pixel, which is the
 * boundary rule the filters use, so padding by r and filtering is equivalent to filtering with Fixed
 * padding -- worth having for exactly that composition. "Reflected" mirrors WITHOUT repeating the edge
 * pixel: {1,2,3} padded by 1 gives {2,1,2,3,2}, not {1,1,2,3,3}. That distinction matters because the
 * repeating variant doubles the edge sample and so biases any subsequent average toward the border.
 *
 * ImageCrop with no size TRIMS A UNIFORM BORDER, which is a different operation from cropping to a
 * size: it asks how much of the frame carries no information. The border colour is taken from a corner
 * rather than assumed to be black, because a scanned page's margin is white and assuming black would
 * trim nothing.
 */
typedef enum { PAD_VALUE, PAD_FIXED, PAD_REFLECT } PadMode;

/* Map a padded coordinate back into the source, or report that it lies outside. */
static bool pad_src_index(int64_t v, size_t n, PadMode mode, size_t* out) {
    if (v >= 0 && (size_t)v < n) { *out = (size_t)v; return true; }
    if (mode == PAD_VALUE) return false;
    if (mode == PAD_FIXED) { *out = (v < 0) ? 0 : n - 1; return true; }
    /* Reflect without repeating the edge: index -1 maps to 1, index n maps to n-2. A period of
     * 2n - 2 covers arbitrarily deep padding rather than only one width. */
    if (n == 1) { *out = 0; return true; }
    int64_t period = 2 * (int64_t)n - 2;
    int64_t m = v % period;
    if (m < 0) m += period;
    if (m >= (int64_t)n) m = period - m;
    *out = (size_t)m;
    return true;
}

/* ImagePad[image, m] / [image, {{l,r},{b,t}}] / with a value or "Fixed"/"Reflected" */
static Expr* builtin_imagepad(Expr* res) {
    size_t argc = res->data.function.arg_count;
    if (argc != 2 && argc != 3) return NULL;

    PadMode mode = PAD_VALUE;
    double fill = 0.0;
    if (argc == 3) {
        Expr* pv = res->data.function.args[2];
        double im = 0.0;
        if (pv && pv->type == EXPR_STRING) {
            if (strcmp(pv->data.string, "Fixed") == 0) mode = PAD_FIXED;
            else if (strcmp(pv->data.string, "Reflected") == 0) mode = PAD_REFLECT;
            else return NULL;
        } else if (na_read_scalar(pv, &fill, &im) && im == 0.0) {
            mode = PAD_VALUE;
        } else return NULL;
    }

    /* The amounts. A single number pads all four sides; {{l,r},{b,t}} is Mathematica's form. */
    int64_t pl = 0, pr = 0, pb = 0, pt = 0;
    Expr* sp = res->data.function.args[1];
    double a = 0.0, im = 0.0;
    if (na_read_scalar(sp, &a, &im) && im == 0.0) {
        if (a != floor(a)) return NULL;
        pl = pr = pb = pt = (int64_t)a;
    } else if (sp && sp->type == EXPR_FUNCTION && sp->data.function.head
               && sp->data.function.head->type == EXPR_SYMBOL
               && sp->data.function.head->data.symbol.name == SYM_List
               && sp->data.function.arg_count == 2) {
        Expr* hx = sp->data.function.args[0];
        Expr* hy = sp->data.function.args[1];
        double l = 0, r = 0, b = 0, t = 0;
        if (!hx || hx->type != EXPR_FUNCTION || hx->data.function.arg_count != 2) return NULL;
        if (!hy || hy->type != EXPR_FUNCTION || hy->data.function.arg_count != 2) return NULL;
        if (!na_read_scalar(hx->data.function.args[0], &l, &im) || im != 0.0) return NULL;
        if (!na_read_scalar(hx->data.function.args[1], &r, &im) || im != 0.0) return NULL;
        if (!na_read_scalar(hy->data.function.args[0], &b, &im) || im != 0.0) return NULL;
        if (!na_read_scalar(hy->data.function.args[1], &t, &im) || im != 0.0) return NULL;
        if (l != floor(l) || r != floor(r) || b != floor(b) || t != floor(t)) return NULL;
        pl = (int64_t)l; pr = (int64_t)r; pb = (int64_t)b; pt = (int64_t)t;
    } else return NULL;

    size_t w = 0, h = 0, c = 0; double* src = NULL;
    if (!image_load(res->data.function.args[0], &w, &h, &c, &src)) return NULL;

    int64_t nw = (int64_t)w + pl + pr;
    int64_t nh = (int64_t)h + pb + pt;
    if (nw < 1 || nh < 1) { free(src); return NULL; }   /* negative padding may not erase the image */

    double* dst = malloc(sizeof(double) * (size_t)nw * (size_t)nh * c);
    Expr* out = NULL;
    if (dst) {
        /* THE INTERIOR IS A ROW-BLOCK COPY, not a per-pixel map. The first version called the
         * coordinate mapper once per pixel per channel, which put a 512x512 constant pad at 0.57 ms
         * against numpy's 0.055 -- 10x off, on an operation that only moves memory. Every mode agrees
         * on the interior (it is the source, unshifted), so the modes differ only over the border,
         * which is a thin frame. Copying each row's in-range span with memcpy and computing only the
         * frame per-pixel is the whole difference.
         *
         * Written as a span per destination row, which also handles NEGATIVE padding without a second
         * code path: a negative offset simply makes the in-range span start later and end sooner. */
        size_t rowbytes_all = sizeof(double) * (size_t)nw * c;
        if (mode == PAD_VALUE) {
            if (fill == 0.0) {
                memset(dst, 0, rowbytes_all * (size_t)nh);      /* the common case, and the fastest */
            } else {
                size_t n = (size_t)nw * (size_t)nh * c;
                for (size_t i = 0; i < n; i++) dst[i] = fill;
            }
        }
        for (int64_t y = 0; y < nh; y++) {
            int64_t sy = y - pt;
            double* drow = dst + (size_t)y * (size_t)nw * c;
            if (sy >= 0 && sy < (int64_t)h) {
                /* The destination columns that map inside the source: [x0, x1). */
                int64_t x0 = pl > 0 ? pl : 0;
                int64_t x1 = pl + (int64_t)w;
                if (x1 > nw) x1 = nw;
                if (x1 > x0)
                    memcpy(drow + (size_t)x0 * c,
                           src + ((size_t)sy * w + (size_t)(x0 - pl)) * c,
                           sizeof(double) * (size_t)(x1 - x0) * c);
                if (mode != PAD_VALUE) {                        /* only the left/right frame remains */
                    for (int64_t x = 0; x < nw; x++) {
                        if (x >= x0 && x < x1) continue;
                        size_t ux;
                        pad_src_index(x - pl, w, mode, &ux);
                        for (size_t k = 0; k < c; k++)
                            drow[(size_t)x * c + k] = src[((size_t)sy * w + ux) * c + k];
                    }
                }
            } else if (mode != PAD_VALUE) {
                size_t uy;
                pad_src_index(sy, h, mode, &uy);                /* a whole row from outside the frame */
                for (int64_t x = 0; x < nw; x++) {
                    size_t ux;
                    pad_src_index(x - pl, w, mode, &ux);
                    for (size_t k = 0; k < c; k++)
                        drow[(size_t)x * c + k] = src[(uy * w + ux) * c + k];
                }
            }
        }
        out = image_build_real(dst, (size_t)nw, (size_t)nh, c);
    }
    free(src); free(dst);
    return out;
}

/* Shrink the frame while its outer ring matches the corner colour, one edge at a time.
 *
 * A plain function rather than the macro this first was: the macro needed a GNU statement expression
 * to return a value, which is exactly the extension this tree forbids -- it compiles on clang and is
 * not C99. Each edge is tested independently, so a border uniform on three sides and not the fourth
 * trims the three. */
static bool span_uniform(const double* p, size_t w, size_t c, const double* ref,
                         size_t a0, size_t a1, size_t fixed, bool row) {
    for (size_t i = a0; i < a1; i++) {
        size_t idx = row ? (fixed * w + i) : (i * w + fixed);
        for (size_t k = 0; k < c; k++)
            if (p[idx * c + k] != ref[k]) return false;
    }
    return true;
}

static void border_trim(const double* p, size_t w, size_t h, size_t c,
                        size_t* l, size_t* r, size_t* t, size_t* b) {
    const double* ref = p;                       /* the top-left pixel is the border colour */
    *l = 0; *r = w; *t = 0; *b = h;
    while (*t < *b && span_uniform(p, w, c, ref, *l, *r, *t, true)) (*t)++;
    while (*b > *t && span_uniform(p, w, c, ref, *l, *r, *b - 1, true)) (*b)--;
    while (*l < *r && span_uniform(p, w, c, ref, *t, *b, *l, false)) (*l)++;
    while (*r > *l && span_uniform(p, w, c, ref, *t, *b, *r - 1, false)) (*r)--;
}

/* ImageCrop[image, {w, h}] -- centred crop. ImageCrop[image] -- trim a uniform border. */
static Expr* builtin_imagecrop(Expr* res) {
    size_t argc = res->data.function.arg_count;
    if (argc != 1 && argc != 2) return NULL;
    size_t w = 0, h = 0, c = 0; double* src = NULL;
    if (!image_load(res->data.function.args[0], &w, &h, &c, &src)) return NULL;

    size_t x0 = 0, y0 = 0, cw = w, ch = h;
    if (argc == 2) {
        Expr* sp = res->data.function.args[1];
        double a = 0.0, b = 0.0, im = 0.0;
        if (!sp || sp->type != EXPR_FUNCTION || sp->data.function.arg_count != 2
            || sp->data.function.head->type != EXPR_SYMBOL
            || sp->data.function.head->data.symbol.name != SYM_List) { free(src); return NULL; }
        if (!na_read_scalar(sp->data.function.args[0], &a, &im) || im != 0.0
         || !na_read_scalar(sp->data.function.args[1], &b, &im) || im != 0.0) {
            free(src); return NULL;
        }
        if (!(a >= 1.0) || a != floor(a) || !(b >= 1.0) || b != floor(b)) { free(src); return NULL; }
        cw = (size_t)a; ch = (size_t)b;
        if (cw > w || ch > h) { free(src); return NULL; }   /* a crop cannot enlarge */
        /* Centred, with any odd remainder going to the RIGHT and BOTTOM. That is the same
         * floor-division convention the kernel centres use, so a pad of m followed by a crop back to
         * the original size lands exactly where it started. */
        x0 = (w - cw) / 2;
        y0 = (h - ch) / 2;
    } else {
        /* Trim a uniform border. The colour comes from the top-left corner rather than being assumed
         * black: a scanned page's margin is white, and assuming black would trim nothing. */
        size_t l = 0, r = w, t = 0, bo = h;
        border_trim(src, w, h, c, &l, &r, &t, &bo);
        /* An entirely uniform image has no content to keep; returning nothing is not an image, so it
         * comes back unchanged rather than as a zero-sized one. */
        if (r <= l || bo <= t) { x0 = 0; y0 = 0; cw = w; ch = h; }
        else { x0 = l; y0 = t; cw = r - l; ch = bo - t; }
    }

    double* dst = malloc(sizeof(double) * cw * ch * c);
    Expr* out = NULL;
    if (dst) {
        for (size_t y = 0; y < ch; y++)
            for (size_t x = 0; x < cw; x++)
                for (size_t k = 0; k < c; k++)
                    dst[(y * cw + x) * c + k] = src[((y0 + y) * w + (x0 + x)) * c + k];
        out = image_build_real(dst, cw, ch, c);
    }
    free(src); free(dst);
    return out;
}

void imagegeom_init(void) {
    symtab_add_builtin("ImagePad", builtin_imagepad);
    symtab_get_def("ImagePad")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("ImagePad",
        "ImagePad[image, m] pads m pixels on every side; ImagePad[image, {{left, right}, "
        "{bottom, top}}] pads each side separately, in Mathematica's VISUAL order -- so `top` adds rows "
        "at the start of the data, since row 1 is the top of the image. Negative amounts crop, but may "
        "not erase the image. ImagePad[image, m, v] fills with the value v (default 0); "
        "ImagePad[image, m, \"Fixed\"] replicates the edge pixel, the same boundary rule the filters "
        "use, so padding then filtering composes with it; ImagePad[image, m, \"Reflected\"] mirrors "
        "WITHOUT repeating the edge -- {1,2,3} padded by 1 gives {2,1,2,3,2}, not {1,1,2,3,3}, because "
        "doubling the edge sample biases any later average toward the border. Reflection uses a period "
        "of 2n-2, so padding deeper than the image still works.");

    symtab_add_builtin("ImageCrop", builtin_imagecrop);
    symtab_get_def("ImageCrop")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("ImageCrop",
        "ImageCrop[image, {w, h}] crops to w x h about the centre, any odd remainder going to the right "
        "and bottom -- the same floor-division convention the kernel centres use, which is what makes "
        "ImageCrop[ImagePad[image, m], ImageDimensions[image]] exactly the original image. A crop may "
        "not enlarge. ImageCrop[image] instead TRIMS A UNIFORM BORDER, asking how much of the frame "
        "carries no information; the border colour is read from a corner rather than assumed black, "
        "since a scanned page's margin is white. An entirely uniform image comes back unchanged, there "
        "being no content to keep and a zero-sized image not being one.");
    symtab_add_builtin("ImageRotate", builtin_imagerotate);
    symtab_get_def("ImageRotate")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("ImageRotate",
        "ImageRotate[image] rotates a quarter turn counterclockwise; ImageRotate[image, angle] rotates "
        "by angle in radians (use n Degree for degrees). A multiple of a right angle takes an EXACT "
        "index-permutation path -- every pixel lands on another pixel's position, nothing is "
        "interpolated, and four quarter turns are exactly the identity. An odd number of quarter turns "
        "swaps the dimensions. Any other angle interpolates bilinearly, sampling the source per "
        "destination pixel (inverse mapping, so every output is filled exactly once; forward mapping "
        "leaves holes wherever the rotation stretches). Area rotated in from outside reads as 0 rather "
        "than the replicated edge, because that area was never photographed and smearing the border "
        "across it would invent content.");

    symtab_add_builtin("ImageReflect", builtin_imagereflect);
    symtab_get_def("ImageReflect")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("ImageReflect",
        "ImageReflect[image] reflects top to bottom; ImageReflect[image, Left] or Right reflects left "
        "to right, and Top or Bottom reflects vertically. A pure index permutation, so it is exact and "
        "self-inverse: reflecting twice is exactly the identity.");
    symtab_add_builtin("ImageResize", builtin_imageresize);
    symtab_get_def("ImageResize")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("ImageResize",
        "ImageResize[image, {w, h}] resizes to w x h pixels; ImageResize[image, w] gives width w "
        "with the height following to preserve the aspect ratio. "
        "Resampling -> \"Nearest\" | \"Bilinear\" | \"Average\" selects the method; the default "
        "Automatic uses AREA AVERAGING when either axis shrinks and bilinear otherwise. That "
        "default is about aliasing: point-sampling a shrinking image destroys every frequency "
        "above half the new sampling rate -- a fine checkerboard reduced by nearest-neighbour "
        "comes back a flat field -- and no interpolation afterwards can restore what "
        "point-sampling discarded. Area averaging is a box prefilter and a resample in one pass, "
        "exact for integer reduction factors, using true fractional coverage so a 3 -> 2 "
        "reduction is as correct as 4 -> 2. Enlarging has no frequencies to remove, so bilinear "
        "is used there; area averaging on an enlargement would degenerate to nearest. Coordinates "
        "are centre-aligned, avoiding the half-pixel shift that sx = i * scale introduces at any "
        "scale other than 1:1. The result is a \"Real\" image; sizes must be positive integers.");
}
