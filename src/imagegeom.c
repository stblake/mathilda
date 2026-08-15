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
#include "image.h"

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

void imagegeom_init(void) {
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
