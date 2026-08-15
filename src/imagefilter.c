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

/* Is the kernel SEPARABLE -- that is, rank 1, an outer product u (x) v?
 *
 * A separable kernel turns kw * kh multiply-adds per pixel into kw + kh, which at radius 4 is 18
 * instead of 81. Gaussians are separable by construction, since
 * exp(-(dx^2 + dy^2)/2s^2) = exp(-dx^2/2s^2) exp(-dy^2/2s^2), and so is any box.
 *
 * THE TOLERANCE IS THE WHOLE RISK, and it is deliberately tight. Treating a non-separable kernel
 * as separable does not make it slightly wrong -- it computes a completely different filter. So the
 * check is a relative one against the kernel's own magnitude at 1e-12: a Gaussian factorises to
 * within about 1e-16 relative, so real separable kernels pass with room to spare, while anything
 * genuinely rank 2 or higher (an identity matrix, a rotation-sensitive edge kernel) fails long
 * before it could be mistaken for rank 1. A test convolves with {{1,0},{0,1}} specifically to pin
 * that the direct path still runs for it.
 *
 * The factorisation itself is the textbook one: pick the largest-magnitude entry as pivot, take its
 * column as u and its row (scaled) as v, then verify every entry. Pivoting on the largest entry
 * rather than on K[0][0] matters, because a kernel with a zero in the corner -- which a
 * derivative kernel has -- would otherwise divide by it. */
static bool ker_separable(const double* k, size_t kh, size_t kw, double* u, double* v) {
    size_t pi = 0, pj = 0;
    double best = 0.0;
    for (size_t i = 0; i < kh; i++)
        for (size_t j = 0; j < kw; j++) {
            double a = fabs(k[i * kw + j]);
            if (a > best) { best = a; pi = i; pj = j; }
        }
    if (!(best > 0.0)) return false;              /* an all-zero kernel: nothing to factor */

    double piv = k[pi * kw + pj];
    for (size_t i = 0; i < kh; i++) u[i] = k[i * kw + pj];
    for (size_t j = 0; j < kw; j++) v[j] = k[pi * kw + j] / piv;

    double tol = 1e-12 * best;
    for (size_t i = 0; i < kh; i++)
        for (size_t j = 0; j < kw; j++)
            if (fabs(k[i * kw + j] - u[i] * v[j]) > tol) return false;
    return true;
}

/* Two-pass separable convolution: horizontal with v, then vertical with u.
 *
 * EXACTLY equal to the direct form, and not merely close, because "Fixed" padding is itself
 * separable -- the direct read src[clamp(y)][clamp(x)] clamps the two axes independently, which is
 * precisely what doing one axis and then the other does. Only the summation ORDER differs, so the
 * two agree to floating-point rounding; and because ImageConvolve routes BOTH itself and
 * GaussianFilter through this same path when the kernel is separable, the documented
 * GaussianFilter == ImageConvolve[.., GaussianMatrix[..]] identity stays bit-exact rather than
 * becoming approximate. */
static void convolve_separable(const double* src, double* dst, double* tmp,
                               size_t w, size_t h, size_t c,
                               const double* u, size_t kh, const double* v, size_t kw) {
    int64_t ci = (int64_t)(kh / 2), cj = (int64_t)(kw / 2);

    for (size_t y = 0; y < h; y++)
        for (size_t x = 0; x < w; x++)
            for (size_t ch = 0; ch < c; ch++) {
                double acc = 0.0;
                for (size_t j = 0; j < kw; j++) {
                    size_t sx = clampi((int64_t)x - (int64_t)j + cj, w);
                    acc += v[j] * src[(y * w + sx) * c + ch];
                }
                tmp[(y * w + x) * c + ch] = acc;
            }

    for (size_t y = 0; y < h; y++)
        for (size_t x = 0; x < w; x++)
            for (size_t ch = 0; ch < c; ch++) {
                double acc = 0.0;
                for (size_t i = 0; i < kh; i++) {
                    size_t sy = clampi((int64_t)y - (int64_t)i + ci, h);
                    acc += u[i] * tmp[(sy * w + x) * c + ch];
                }
                dst[(y * w + x) * c + ch] = acc;
            }
}

/* Convolve, taking the separable path when the kernel allows it. */
static bool convolve_dispatch(const double* src, double* dst, size_t w, size_t h, size_t c,
                              const double* k, size_t kw, size_t kh) {
    double* u = malloc(sizeof(double) * kh);
    double* v = malloc(sizeof(double) * kw);
    bool sep = false;
    if (u && v) sep = ker_separable(k, kh, kw, u, v);
    if (sep) {
        double* tmp = malloc(sizeof(double) * w * h * c);
        if (tmp) {
            convolve_separable(src, dst, tmp, w, h, c, u, kh, v, kw);
            free(tmp); free(u); free(v);
            return true;
        }
        sep = false;                              /* no scratch: fall back rather than fail */
    }
    free(u); free(v);
    convolve_planes(src, dst, w, h, c, k, kw, kh);
    return true;
}

/* ---- volumetric convolution -----------------------------------------------
 *
 * SEPARABILITY MATTERS MORE IN THREE DIMENSIONS THAN IN TWO, and by a widening margin: a rank-1
 * kernel costs kw + kh + kd taps instead of kw * kh * kd, so radius 1 goes from 27 to 9 and radius 4
 * from 729 to 27. In two dimensions skipping it costs a factor of the radius; here it costs the
 * SQUARE of it, which is the difference between a volume filter being usable and not.
 *
 * The rank-3 factorisation follows the rank-2 one: pivot on the largest-magnitude entry so a kernel
 * with a zero at a corner still factors, take the three axis-lines through that pivot as the
 * candidate factors, then VERIFY every entry against their product at a tight relative tolerance.
 * Verifying is the whole point -- treating a non-separable kernel as separable computes a different
 * filter, not a slightly wrong one.
 */
static bool ker3_load(const Expr* e, size_t* kd, size_t* kh, size_t* kw, double** buf) {
    if (!ker_is_list(e) || e->data.function.arg_count == 0) return false;
    size_t d = e->data.function.arg_count;
    const Expr* s0 = e->data.function.args[0];
    if (!ker_is_list(s0) || s0->data.function.arg_count == 0) return false;
    size_t hh = s0->data.function.arg_count;
    const Expr* r0 = s0->data.function.args[0];
    if (!ker_is_list(r0) || r0->data.function.arg_count == 0) return false;
    size_t ww = r0->data.function.arg_count;

    double* k = malloc(sizeof(double) * d * hh * ww);
    if (!k) return false;
    for (size_t z = 0; z < d; z++) {
        const Expr* sl = e->data.function.args[z];
        if (!ker_is_list(sl) || sl->data.function.arg_count != hh) { free(k); return false; }
        for (size_t y = 0; y < hh; y++) {
            const Expr* row = sl->data.function.args[y];
            if (!ker_is_list(row) || row->data.function.arg_count != ww) { free(k); return false; }
            for (size_t x = 0; x < ww; x++) {
                double re = 0.0, im = 0.0;
                if (!na_read_scalar(row->data.function.args[x], &re, &im) || im != 0.0) {
                    free(k); return false;
                }
                k[(z * hh + y) * ww + x] = re;
            }
        }
    }
    *kd = d; *kh = hh; *kw = ww; *buf = k;
    return true;
}

static bool ker3_separable(const double* k, size_t kd, size_t kh, size_t kw,
                           double* u, double* v, double* t) {
    size_t pz = 0, py = 0, px = 0;
    double best = 0.0;
    for (size_t z = 0; z < kd; z++)
        for (size_t y = 0; y < kh; y++)
            for (size_t x = 0; x < kw; x++) {
                double a = fabs(k[(z * kh + y) * kw + x]);
                if (a > best) { best = a; pz = z; py = y; px = x; }
            }
    if (!(best > 0.0)) return false;
    double piv = k[(pz * kh + py) * kw + px];

    /* The three axis-lines through the pivot, each scaled so its pivot entry is 1. The pivot value
     * is then folded into u alone, so the product u v t reproduces K rather than K/piv^2. */
    for (size_t z = 0; z < kd; z++) u[z] = k[(z * kh + py) * kw + px] / piv;
    for (size_t y = 0; y < kh; y++) v[y] = k[(pz * kh + y) * kw + px] / piv;
    for (size_t x = 0; x < kw; x++) t[x] = k[(pz * kh + py) * kw + x] / piv;

    double tol = 1e-12 * best;
    for (size_t z = 0; z < kd; z++)
        for (size_t y = 0; y < kh; y++)
            for (size_t x = 0; x < kw; x++)
                if (fabs(k[(z * kh + y) * kw + x] - piv * u[z] * v[y] * t[x]) > tol) return false;
    for (size_t z = 0; z < kd; z++) u[z] *= piv;
    return true;
}

/* Three passes: x, then y, then z. Clamping is per-axis and therefore separable, exactly as in the
 * 2-D case, so this equals the direct form up to summation order. */
static void convolve3_separable(const double* src, double* dst, double* a, double* b,
                                size_t w, size_t h, size_t d, size_t c,
                                const double* u, size_t kd, const double* v, size_t kh,
                                const double* t, size_t kw) {
    int64_t cz = (int64_t)(kd / 2), cy = (int64_t)(kh / 2), cx = (int64_t)(kw / 2);
    size_t plane = w * h * c;

    for (size_t z = 0; z < d; z++)
      for (size_t y = 0; y < h; y++)
        for (size_t x = 0; x < w; x++)
          for (size_t ch = 0; ch < c; ch++) {
            double acc = 0.0;
            for (size_t j = 0; j < kw; j++) {
                size_t sx = clampi((int64_t)x - (int64_t)j + cx, w);
                acc += t[j] * src[z * plane + (y * w + sx) * c + ch];
            }
            a[z * plane + (y * w + x) * c + ch] = acc;
          }

    for (size_t z = 0; z < d; z++)
      for (size_t y = 0; y < h; y++)
        for (size_t x = 0; x < w; x++)
          for (size_t ch = 0; ch < c; ch++) {
            double acc = 0.0;
            for (size_t i = 0; i < kh; i++) {
                size_t sy = clampi((int64_t)y - (int64_t)i + cy, h);
                acc += v[i] * a[z * plane + (sy * w + x) * c + ch];
            }
            b[z * plane + (y * w + x) * c + ch] = acc;
          }

    for (size_t z = 0; z < d; z++)
      for (size_t y = 0; y < h; y++)
        for (size_t x = 0; x < w; x++)
          for (size_t ch = 0; ch < c; ch++) {
            double acc = 0.0;
            for (size_t m = 0; m < kd; m++) {
                size_t sz = clampi((int64_t)z - (int64_t)m + cz, d);
                acc += u[m] * b[sz * plane + (y * w + x) * c + ch];
            }
            dst[z * plane + (y * w + x) * c + ch] = acc;
          }
}

static void convolve3_direct(const double* src, double* dst,
                             size_t w, size_t h, size_t d, size_t c,
                             const double* k, size_t kd, size_t kh, size_t kw) {
    int64_t cz = (int64_t)(kd / 2), cy = (int64_t)(kh / 2), cx = (int64_t)(kw / 2);
    size_t plane = w * h * c;
    for (size_t z = 0; z < d; z++)
      for (size_t y = 0; y < h; y++)
        for (size_t x = 0; x < w; x++)
          for (size_t ch = 0; ch < c; ch++) {
            double acc = 0.0;
            for (size_t m = 0; m < kd; m++) {
                size_t sz = clampi((int64_t)z - (int64_t)m + cz, d);
                for (size_t i = 0; i < kh; i++) {
                    size_t sy = clampi((int64_t)y - (int64_t)i + cy, h);
                    for (size_t j = 0; j < kw; j++) {
                        size_t sx = clampi((int64_t)x - (int64_t)j + cx, w);
                        acc += k[(m * kh + i) * kw + j] * src[sz * plane + (sy * w + sx) * c + ch];
                    }
                }
            }
            dst[z * plane + (y * w + x) * c + ch] = acc;
          }
}

/* Convolve a volume, separably when the kernel allows. Returns the built Image3D or NULL. */
static Expr* convolve3_run(Expr* vol, const double* k, size_t kd, size_t kh, size_t kw) {
    size_t w = 0, h = 0, d = 0, c = 0; double* src = NULL;
    if (!image3d_load(vol, &w, &h, &d, &c, &src)) return NULL;
    size_t n = w * h * d * c;
    double* dst = malloc(sizeof(double) * n);
    double* u = malloc(sizeof(double) * kd);
    double* v = malloc(sizeof(double) * kh);
    double* t = malloc(sizeof(double) * kw);
    Expr* out = NULL;
    if (dst && u && v && t) {
        if (ker3_separable(k, kd, kh, kw, u, v, t)) {
            double* a = malloc(sizeof(double) * n);
            double* b = malloc(sizeof(double) * n);
            if (a && b) convolve3_separable(src, dst, a, b, w, h, d, c, u, kd, v, kh, t, kw);
            else convolve3_direct(src, dst, w, h, d, c, k, kd, kh, kw);
            free(a); free(b);
        } else {
            convolve3_direct(src, dst, w, h, d, c, k, kd, kh, kw);
        }
        out = image3d_build_real(dst, w, h, d, c);
    }
    free(src); free(dst); free(u); free(v); free(t);
    return out;
}

/* A separable 3-D Gaussian of radius r, normalised to sum 1. Built as a full rank-3 kernel and
 * handed to the same dispatcher, so its separability is re-derived and verified rather than
 * assumed -- the same discipline the 2-D derivative kernels follow. */
static double* gauss3_kernel(size_t r, size_t* kd, size_t* kh, size_t* kw) {
    size_t n = 2 * r + 1;
    double sigma = (r == 0) ? 1.0 : (double)r / 2.0;
    double* k = malloc(sizeof(double) * n * n * n);
    if (!k) return NULL;
    double sum = 0.0;
    for (size_t z = 0; z < n; z++) {
        double dz = (double)z - (double)r;
        for (size_t y = 0; y < n; y++) {
            double dy = (double)y - (double)r;
            for (size_t x = 0; x < n; x++) {
                double dx = (double)x - (double)r;
                double val = exp(-(dx * dx + dy * dy + dz * dz) / (2.0 * sigma * sigma));
                k[(z * n + y) * n + x] = val;
                sum += val;
            }
        }
    }
    if (!(sum > 0.0)) { free(k); return NULL; }
    for (size_t i = 0; i < n * n * n; i++) k[i] /= sum;
    *kd = *kh = *kw = n;
    return k;
}

/* ImageConvolve[image, kernel] */
static Expr* builtin_imageconvolve(Expr* res) {
    if (res->data.function.arg_count != 2) return NULL;
    /* A VOLUME takes the rank-3 path. Dispatching on the image rather than on the kernel is
     * deliberate: a rank-3 kernel handed to a plane is a mistake worth declining, not something to
     * reinterpret as a stack of 2-D kernels. */
    if (image3d_info(res->data.function.args[0], NULL, NULL, NULL, NULL, NULL)) {
        size_t kd = 0, kh3 = 0, kw3 = 0; double* k3 = NULL;
        if (!ker3_load(res->data.function.args[1], &kd, &kh3, &kw3, &k3)) return NULL;
        Expr* o = convolve3_run(res->data.function.args[0], k3, kd, kh3, kw3);
        free(k3);
        return o;
    }
    size_t w = 0, h = 0, c = 0; double* src = NULL;
    if (!image_load(res->data.function.args[0], &w, &h, &c, &src)) return NULL;
    size_t kw = 0, kh = 0; double* k = NULL;
    if (!ker_load(res->data.function.args[1], &kh, &kw, &k)) { free(src); return NULL; }

    double* dst = malloc(sizeof(double) * w * h * c);
    Expr* out = NULL;
    if (dst) {
        convolve_dispatch(src, dst, w, h, c, k, kw, kh);
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
    if (image3d_info(res->data.function.args[0], NULL, NULL, NULL, NULL, NULL)) {
        double rr = 0.0, im = 0.0;
        if (!na_read_scalar(res->data.function.args[1], &rr, &im) || im != 0.0) return NULL;
        if (!(rr >= 0.0) || rr != floor(rr) || rr > 32.0) return NULL;
        size_t kd = 0, kh3 = 0, kw3 = 0;
        double* k3 = gauss3_kernel((size_t)rr, &kd, &kh3, &kw3);
        if (!k3) return NULL;
        Expr* o = convolve3_run(res->data.function.args[0], k3, kd, kh3, kw3);
        free(k3);
        return o;
    }
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
        convolve_dispatch(src, dst, w, h, c, k, kw, kh);
        out = image_build_real(dst, w, h, c);
    }
    free(src); free(k); free(dst);
    return out;
}

/* ---- greyscale, thresholding, segmentation ---------------------------------
 *
 * REC. 601 LUMINANCE, 0.299 R + 0.587 G + 0.114 B, which is what Mathematica's conversion to
 * greyscale uses. The weights are not a mean because the eye is not equally sensitive across the
 * spectrum -- green carries most of the perceived brightness and blue almost none -- so a plain
 * average would make a saturated blue and a saturated green look equally bright, which they are
 * not. This matters for thresholding specifically: an unweighted average puts pure red and pure
 * blue on the same side of any threshold, when perceptually they are far apart.
 */
static void img_to_grey(const double* src, double* dst, size_t w, size_t h, size_t c) {
    static const double R = 0.299, G = 0.587, B = 0.114;
    for (size_t i = 0; i < w * h; i++) {
        if (c == 1) { dst[i] = src[i]; continue; }
        if (c >= 3) {
            dst[i] = R * src[i * c + 0] + G * src[i * c + 1] + B * src[i * c + 2];
        } else {
            /* Two channels is grey-plus-alpha; the alpha is not brightness. */
            dst[i] = src[i * c + 0];
        }
    }
}

#define IMG_OTSU_BINS 256

/* Otsu's threshold: the level maximising BETWEEN-CLASS variance.
 *
 * The classic 1979 result, and the reason it is the default everywhere: maximising the variance
 * *between* the two classes is algebraically the same as minimising the weighted variance
 * *within* them, so one cheap pass over a histogram optimises the thing you actually want -- how
 * well separated the two groups are -- without ever computing a within-class variance.
 *
 *   sigma_b^2(t) = w0(t) w1(t) (mu0(t) - mu1(t))^2
 *
 * Both class weights and both means update INCREMENTALLY as t advances by one bin, so the whole
 * search is O(bins) after the histogram, not O(bins^2). Written with running sums for exactly
 * that reason: recomputing the means per candidate is the obvious implementation and is
 * quadratic.
 *
 * Values are binned over [0, 1] because that is the range image_load guarantees. Anything
 * outside -- a Real image may legitimately hold it, since Image stores faithfully rather than
 * clamping -- is clamped INTO the histogram rather than dropped, so an out-of-range pixel still
 * votes for the extreme it belongs to instead of vanishing from the statistics.
 *
 * Returns false for a degenerate image where every pixel is identical: there is no threshold
 * that separates one cluster into two, and inventing one would be a fiction. */
static bool img_otsu(const double* grey, size_t n, double* thresh) {
    if (!grey || n == 0) return false;
    double hist[IMG_OTSU_BINS];
    for (size_t i = 0; i < IMG_OTSU_BINS; i++) hist[i] = 0.0;

    for (size_t i = 0; i < n; i++) {
        double v = grey[i];
        if (!(v == v)) return false;                    /* a NaN pixel has no bin */
        if (v < 0.0) v = 0.0;
        if (v > 1.0) v = 1.0;
        size_t b = (size_t)(v * (double)(IMG_OTSU_BINS - 1) + 0.5);
        if (b >= IMG_OTSU_BINS) b = IMG_OTSU_BINS - 1;
        hist[b] += 1.0;
    }

    /* One occupied bin means one cluster: no split exists. */
    size_t occupied = 0;
    for (size_t i = 0; i < IMG_OTSU_BINS; i++) if (hist[i] > 0.0) occupied++;
    if (occupied < 2) return false;

    double total = (double)n, sum_all = 0.0;
    for (size_t i = 0; i < IMG_OTSU_BINS; i++) sum_all += (double)i * hist[i];

    double w0 = 0.0, sum0 = 0.0, best = -1.0;
    size_t best_bin = 0;
    for (size_t t = 0; t + 1 < IMG_OTSU_BINS; t++) {
        w0 += hist[t];
        sum0 += (double)t * hist[t];
        double w1 = total - w0;
        if (!(w0 > 0.0) || !(w1 > 0.0)) continue;       /* an empty class has no mean */
        double mu0 = sum0 / w0;
        double mu1 = (sum_all - sum0) / w1;
        double d = mu0 - mu1;
        double sb = w0 * w1 * d * d;
        /* Strictly greater, so a tie keeps the LOWEST bin -- deterministic, and the same
         * lowest-index tie-break the rest of the tree uses. */
        if (sb > best) { best = sb; best_bin = t; }
    }
    if (!(best > 0.0)) return false;

    /* The threshold sits at the UPPER EDGE of the winning bin, so a pixel in that bin falls on
     * the low side. Bin b covers values centred on b/255, so its upper edge is halfway to the
     * next centre. */
    *thresh = ((double)best_bin + 0.5) / (double)(IMG_OTSU_BINS - 1);
    if (*thresh > 1.0) *thresh = 1.0;
    return true;
}

/* Load an image and reduce it to a greyscale plane. */
static bool img_grey_plane(Expr* img, size_t* w, size_t* h, double** grey) {
    size_t ww = 0, hh = 0, cc = 0; double* src = NULL;
    if (!image_load(img, &ww, &hh, &cc, &src)) return false;
    double* g = malloc(sizeof(double) * ww * hh);
    if (!g) { free(src); return false; }
    img_to_grey(src, g, ww, hh, cc);
    free(src);
    *w = ww; *h = hh; *grey = g;
    return true;
}

/* FindThreshold[image] -- Otsu's threshold as a real number. */
static Expr* builtin_findthreshold(Expr* res) {
    if (res->data.function.arg_count != 1) return NULL;
    size_t w = 0, h = 0; double* g = NULL;
    if (!img_grey_plane(res->data.function.args[0], &w, &h, &g)) return NULL;
    double t = 0.0;
    bool ok = img_otsu(g, w * h, &t);
    free(g);
    return ok ? expr_new_real(t) : NULL;
}

/* Binarize[image] -- threshold by Otsu.
 * Binarize[image, t] -- threshold at t.
 *
 * A pixel STRICTLY ABOVE the threshold becomes 1, so a pixel exactly at it becomes 0. That
 * matches Mathematica and it is the boundary a test pins, because "above" and "at or above"
 * differ on exactly the pixels a threshold was chosen to sit between.
 *
 * Colour is reduced to luminance first, which is also what Mathematica does. */
static Expr* builtin_binarize(Expr* res) {
    size_t argc = res->data.function.arg_count;
    if (argc != 1 && argc != 2) return NULL;
    size_t w = 0, h = 0; double* g = NULL;
    if (!img_grey_plane(res->data.function.args[0], &w, &h, &g)) return NULL;

    double t = 0.0;
    if (argc == 2) {
        double im = 0.0;
        if (!na_read_scalar(res->data.function.args[1], &t, &im) || im != 0.0) {
            free(g); return NULL;
        }
    } else if (!img_otsu(g, w * h, &t)) {
        free(g); return NULL;
    }

    /* Built as a Bit image directly rather than through image_build_real: the result is 0/1 by
     * construction, and typing it "Real" would lose that -- a later ImageData would then scale
     * nothing but the caller could no longer tell it was binary. */
    Expr** rows = malloc(sizeof(Expr*) * h);
    Expr* out = NULL;
    if (rows) {
        bool ok = true;
        for (size_t y = 0; y < h; y++) rows[y] = NULL;
        for (size_t y = 0; y < h && ok; y++) {
            Expr** cols = malloc(sizeof(Expr*) * w);
            if (!cols) { ok = false; break; }
            bool okc = true;
            for (size_t x = 0; x < w; x++) {
                cols[x] = expr_new_integer(g[y * w + x] > t ? 1 : 0);
                if (!cols[x]) okc = false;
            }
            if (okc) rows[y] = expr_new_function(expr_new_symbol(SYM_List), cols, w);
            else for (size_t x = 0; x < w; x++) expr_free(cols[x]);
            free(cols);
            if (!rows[y]) ok = false;
        }
        if (ok) {
            Expr* data = expr_new_function(expr_new_symbol(SYM_List), rows, h);
            if (data) {
                Expr* two[2];
                two[0] = data;
                two[1] = expr_new_string("Bit");
                if (two[1]) out = expr_new_function(expr_new_symbol("Image"), two, 2);
                else expr_free(data);
            }
        } else {
            for (size_t y = 0; y < h; y++) expr_free(rows[y]);
        }
        free(rows);
    }
    free(g);
    return out;
}

/* ColorConvert[image, "Grayscale"] -- Rec. 601 luminance.
 *
 * Only "Grayscale" is accepted. The other colour spaces Mathematica supports (LAB, HSB, XYZ, ...)
 * each carry their own white point and transfer-function decisions, and accepting the name while
 * doing something approximate would be worse than declining it. */
static Expr* builtin_colorconvert(Expr* res) {
    if (res->data.function.arg_count != 2) return NULL;
    Expr* sp = res->data.function.args[1];
    if (!sp || sp->type != EXPR_STRING) return NULL;
    if (strcmp(sp->data.string, "Grayscale") != 0 && strcmp(sp->data.string, "Gray") != 0)
        return NULL;
    size_t w = 0, h = 0; double* g = NULL;
    if (!img_grey_plane(res->data.function.args[0], &w, &h, &g)) return NULL;
    Expr* out = image_build_real(g, w, h, 1);
    free(g);
    return out;
}

/* ---- derivative and gradient filters --------------------------------------
 *
 * THE STENCILS, and why these normalisations. A derivative filter is a separable outer product of
 * two 1-D stencils, one per axis:
 *
 *   order 0   {1, 2, 1}/4     smooth. Normalised to sum 1, so it preserves a constant -- and,
 *                             being symmetric, it preserves a LINEAR ramp exactly too.
 *   order 1   {-1, 0, 1}/2    central difference. On f(x) = c x this gives
 *                             (c(x+1) - c(x-1))/2 = c EXACTLY, which is what makes "the gradient
 *                             of a ramp is its slope" an equality rather than an approximation.
 *   order 2   {1, -2, 1}      second difference. On f(x) = c x^2 it gives 2c exactly.
 *
 * So DerivativeFilter[img, {0, 1}] is the classic Sobel-x, {1, 0} Sobel-y. Dividing by 4 and 2
 * rather than leaving the raw integer stencils matters: unnormalised Sobel reports a gradient
 * eight times the true slope, which is harmless for edge DETECTION -- where only the ranking
 * matters -- and wrong for anything that reads the number, including a physical gradient or a
 * threshold carried over from another tool.
 *
 * Both stencils are separable and one of them has a zero in the middle, so this is exactly the
 * case the separable path's largest-magnitude pivot exists for: pivoting on K[[1,1]] would divide
 * by that zero. The kernels are built as full 2-D matrices and handed to the same
 * convolve_dispatch every other filter uses, which re-derives the factorisation rather than
 * assuming it -- the factorisation is then verified by the same tolerance as any other kernel,
 * instead of being trusted because the author knew it was separable.
 */
static bool deriv_stencil(int order, double* st, size_t* n) {
    switch (order) {
        case 0: st[0] = 0.25; st[1] = 0.5;  st[2] = 0.25; *n = 3; return true;
        /* NOTE THE ORDER: {+1/2, 0, -1/2}, not {-1/2, 0, +1/2}.
         *
         * ImageConvolve REFLECTS its kernel, so a stencil written in the natural reading order
         * computes the NEGATED derivative: convolution with {-1/2, 0, 1/2} gives
         * (s[x-1] - s[x+1])/2. A derivative filter is really a CORRELATION, which is why every CV
         * library's Sobel is defined with correlate semantics -- so the stencil is pre-flipped here
         * and convolution then yields the true derivative, positive where brightness increases to
         * the right.
         *
         * Caught by asserting an exact value rather than an absolute one: the gradient MAGNITUDE
         * squares the sign away and looked perfectly correct, and so would any edge-detection
         * result, because only the ranking matters there. Only "the derivative of a ramp of slope
         * 1/8 is exactly +1/8" could see it. The order-2 stencil is symmetric and so unaffected. */
        case 1: st[0] = 0.5;  st[1] = 0.0;  st[2] = -0.5; *n = 3; return true;
        case 2: st[0] = 1.0;  st[1] = -2.0; st[2] = 1.0;  *n = 3; return true;
        default: return false;
    }
}

/* Build the (ny x nx) kernel for a {y-order, x-order} derivative. Caller frees. */
static double* deriv_kernel(int oy, int ox, size_t* kh, size_t* kw) {
    double sy[3], sx[3]; size_t ny = 0, nx = 0;
    if (!deriv_stencil(oy, sy, &ny) || !deriv_stencil(ox, sx, &nx)) return NULL;
    double* k = malloc(sizeof(double) * ny * nx);
    if (!k) return NULL;
    for (size_t i = 0; i < ny; i++)
        for (size_t j = 0; j < nx; j++)
            k[i * nx + j] = sy[i] * sx[j];
    *kh = ny; *kw = nx;
    return k;
}

/* Read a {n, m} derivative-order specification. */
static bool deriv_orders(Expr* spec, int* oy, int* ox) {
    if (!ker_is_list(spec) || spec->data.function.arg_count != 2) return false;
    double a = 0.0, b = 0.0, im = 0.0;
    if (!na_read_scalar(spec->data.function.args[0], &a, &im) || im != 0.0) return false;
    if (!na_read_scalar(spec->data.function.args[1], &b, &im) || im != 0.0) return false;
    if (a != floor(a) || b != floor(b) || a < 0.0 || b < 0.0 || a > 2.0 || b > 2.0) return false;
    *oy = (int)a; *ox = (int)b;
    return true;
}

/* DerivativeFilter[image, {n, m}] -- n-th derivative down the rows, m-th across the columns. */
static Expr* builtin_derivativefilter(Expr* res) {
    if (res->data.function.arg_count != 2) return NULL;
    int oy = 0, ox = 0;
    if (!deriv_orders(res->data.function.args[1], &oy, &ox)) return NULL;
    /* {0, 0} is a plain smoothing, which is a legitimate request but not a derivative; it is
     * allowed because the stencil table defines it and refusing would be arbitrary. */
    size_t kh = 0, kw = 0;
    double* k = deriv_kernel(oy, ox, &kh, &kw);
    if (!k) return NULL;

    size_t w = 0, h = 0, c = 0; double* src = NULL;
    if (!image_load(res->data.function.args[0], &w, &h, &c, &src)) { free(k); return NULL; }
    double* dst = malloc(sizeof(double) * w * h * c);
    Expr* out = NULL;
    if (dst) {
        convolve_dispatch(src, dst, w, h, c, k, kw, kh);
        out = image_build_real(dst, w, h, c);
    }
    free(src); free(dst); free(k);
    return out;
}

/* GradientFilter[image] -- the gradient MAGNITUDE, Sqrt[dx^2 + dy^2].
 *
 * ROTATION INVARIANCE is the property that makes this the right combination rather than, say,
 * |dx| + |dy|: the magnitude of a vector does not depend on which way the axes point, so an edge
 * at 45 degrees reports the same strength as one at 0. The absolute-sum alternative is cheaper and
 * reports an edge at 45 degrees as sqrt(2) times stronger, which biases every downstream threshold
 * by orientation.
 *
 * Colour is reduced to luminance FIRST, then differentiated once -- not differentiated per channel
 * and combined. Per-channel gradients would have to be combined by some rule (max? sum? norm?) and
 * every choice is arbitrary; taking the gradient of brightness is the one interpretation that
 * needs no such choice. */
static Expr* builtin_gradientfilter(Expr* res) {
    if (res->data.function.arg_count != 1) return NULL;
    size_t w = 0, h = 0; double* g = NULL;
    if (!img_grey_plane(res->data.function.args[0], &w, &h, &g)) return NULL;

    size_t khx = 0, kwx = 0, khy = 0, kwy = 0;
    double* kx = deriv_kernel(0, 1, &khx, &kwx);
    double* ky = deriv_kernel(1, 0, &khy, &kwy);
    double* dx = malloc(sizeof(double) * w * h);
    double* dy = malloc(sizeof(double) * w * h);
    Expr* out = NULL;
    if (kx && ky && dx && dy) {
        convolve_dispatch(g, dx, w, h, 1, kx, kwx, khx);
        convolve_dispatch(g, dy, w, h, 1, ky, kwy, khy);
        for (size_t i = 0; i < w * h; i++)
            dx[i] = sqrt(dx[i] * dx[i] + dy[i] * dy[i]);
        out = image_build_real(dx, w, h, 1);
    }
    free(g); free(kx); free(ky); free(dx); free(dy);
    return out;
}

/* ---- EdgeDetect: Canny --------------------------------------------------------
 *
 * Four stages, and each exists to fix a specific failure of the previous one:
 *
 *   1. SMOOTH. A derivative amplifies noise -- differencing doubles the noise amplitude while a
 *      real edge keeps its step -- so gradients are taken of a blurred image, never a raw one.
 *   2. GRADIENT. Magnitude and direction from the normalised Sobel pair.
 *   3. NON-MAXIMUM SUPPRESSION. A gradient magnitude ridge is several pixels wide; thresholding it
 *      alone gives a thick band. NMS keeps only pixels that are maximal ALONG THE GRADIENT
 *      DIRECTION, which is what makes a Canny edge one pixel wide.
 *   4. HYSTERESIS. One threshold either breaks long edges wherever they weaken, or admits noise
 *      everywhere. Two thresholds plus connectivity keeps a weak pixel only if it is reachable
 *      from a strong one, so a genuine edge survives its own faint stretches while isolated weak
 *      responses do not.
 *
 * THINNING TO EXACTLY ONE PIXEL DEPENDS ON AN ASYMMETRIC COMPARISON, and this is the subtle part.
 * A clean step edge does not give a single-pixel gradient peak: with the central difference, a step
 * at column k responds 0.5 at BOTH k-1 and k -- a two-wide plateau of exactly equal values. Testing
 * `mag >= both neighbours` keeps both and yields a two-pixel edge. So the test is
 * `mag > backward && mag >= forward`: on a tie the lower-index side wins, deterministically, and
 * the ridge thins to one pixel. That is the difference between an edge detector and a thick mask,
 * and a test pins the width at exactly 1.
 *
 * NMS IS INSENSITIVE TO THE GRADIENT SIGN, which is worth stating after the last change: it compares
 * both neighbours along the direction, and direction and direction+180 degrees select the same pair.
 * So the sign bug fixed in DerivativeFilter would not have been caught here either -- another reason
 * that exact signed value had to be asserted where it was.
 *
 * The direction is quantised into four bins rather than interpolated along the exact angle.
 * Interpolation is more accurate on curved edges; quantisation is exactly reproducible and is what
 * the original algorithm specifies, and reproducibility is what lets a test pin a pixel pattern.
 */
static void canny_nms(const double* mag, const double* dx, const double* dy,
                      double* out, size_t w, size_t h) {
    for (size_t y = 0; y < h; y++) {
        for (size_t x = 0; x < w; x++) {
            size_t i = y * w + x;
            double gx = dx[i], gy = dy[i];
            double m = mag[i];
            if (!(m > 0.0)) { out[i] = 0.0; continue; }

            /* Quantise the direction to one of four neighbour pairs. atan2 is avoided: comparing
             * |gy| against |gx| times tan(22.5) and tan(67.5) picks the same bin without a
             * transcendental call per pixel. */
            double agx = fabs(gx), agy = fabs(gy);
            int64_t ox, oy;
            if (agy <= 0.41421356237309503 * agx)      { ox = 1;  oy = 0;  }   /* tan 22.5 */
            else if (agy >= 2.414213562373095 * agx)   { ox = 0;  oy = 1;  }   /* tan 67.5 */
            else if ((gx > 0.0) == (gy > 0.0))         { ox = 1;  oy = 1;  }
            else                                       { ox = 1;  oy = -1; }

            /* Off the edge counts as zero: a border pixel has no neighbour to lose to, so it is
             * kept if it beats the one neighbour it has. Clamping instead would compare a pixel
             * against itself and keep every border pixel unconditionally. */
            double fwd = 0.0, bwd = 0.0;
            int64_t fx = (int64_t)x + ox, fy = (int64_t)y + oy;
            int64_t bx = (int64_t)x - ox, by = (int64_t)y - oy;
            if (fx >= 0 && (size_t)fx < w && fy >= 0 && (size_t)fy < h)
                fwd = mag[(size_t)fy * w + (size_t)fx];
            if (bx >= 0 && (size_t)bx < w && by >= 0 && (size_t)by < h)
                bwd = mag[(size_t)by * w + (size_t)bx];

            /* ASYMMETRIC on purpose -- see the note above. Strictly greater one way, at-least the
             * other, so an even plateau resolves to one pixel rather than two. */
            out[i] = (m > bwd && m >= fwd) ? m : 0.0;
        }
    }
}

/* Hysteresis: keep every pixel above `hi`, and every pixel above `lo` reachable from one.
 *
 * Iterative with an explicit stack rather than recursive. A long edge across a large image is a
 * connected component thousands of pixels deep, and recursion there is a stack overflow on a real
 * photograph -- not a hypothetical one. */
static void canny_hysteresis(const double* nms, double* out, size_t w, size_t h,
                             double lo, double hi, size_t* stack) {
    size_t n = w * h, top = 0;
    for (size_t i = 0; i < n; i++) {
        if (nms[i] >= hi) { out[i] = 1.0; stack[top++] = i; }
        else out[i] = 0.0;
    }
    while (top > 0) {
        size_t i = stack[--top];
        size_t y = i / w, x = i % w;
        for (int64_t dyy = -1; dyy <= 1; dyy++) {
            for (int64_t dxx = -1; dxx <= 1; dxx++) {
                if (dxx == 0 && dyy == 0) continue;
                int64_t nx = (int64_t)x + dxx, ny = (int64_t)y + dyy;
                if (nx < 0 || (size_t)nx >= w || ny < 0 || (size_t)ny >= h) continue;
                size_t j = (size_t)ny * w + (size_t)nx;
                /* 8-connectivity, because a diagonal edge is connected in any sensible reading and
                 * 4-connectivity would break every 45-degree line into dots. */
                if (out[j] == 0.0 && nms[j] >= lo) { out[j] = 1.0; stack[top++] = j; }
            }
        }
    }
}

/* EdgeDetect[image] / EdgeDetect[image, r] / EdgeDetect[image, r, t] */
static Expr* builtin_edgedetect(Expr* res) {
    size_t argc = res->data.function.arg_count;
    if (argc < 1 || argc > 3) return NULL;

    double rad = 2.0, thi = -1.0, im = 0.0;
    if (argc >= 2) {
        if (!na_read_scalar(res->data.function.args[1], &rad, &im) || im != 0.0) return NULL;
        if (!(rad >= 0.0) || rad != floor(rad) || rad > 64.0) return NULL;
    }
    if (argc == 3) {
        if (!na_read_scalar(res->data.function.args[2], &thi, &im) || im != 0.0) return NULL;
        if (!(thi >= 0.0)) return NULL;
    }

    size_t w = 0, h = 0; double* g = NULL;
    if (!img_grey_plane(res->data.function.args[0], &w, &h, &g)) return NULL;
    size_t n = w * h;

    double* sm  = malloc(sizeof(double) * n);
    double* dx  = malloc(sizeof(double) * n);
    double* dy  = malloc(sizeof(double) * n);
    double* mag = malloc(sizeof(double) * n);
    double* nms = malloc(sizeof(double) * n);
    size_t* stk = malloc(sizeof(size_t) * n);
    size_t khx = 0, kwx = 0, khy = 0, kwy = 0;
    double* kx = deriv_kernel(0, 1, &khx, &kwx);
    double* ky = deriv_kernel(1, 0, &khy, &kwy);
    Expr* out = NULL;

    if (sm && dx && dy && mag && nms && stk && kx && ky) {
        /* 1. Smooth. Radius 0 means no smoothing, which is a legitimate request on already-clean
         * synthetic input and is what makes an exact single-step test possible. */
        if (rad > 0.0) {
            size_t gkh = 0, gkw = 0;
            double* gk = NULL;
            {
                Expr* r1[1]; r1[0] = expr_new_integer((int64_t)rad);
                Expr* call = r1[0] ? expr_new_function(expr_new_symbol("GaussianMatrix"), r1, 1)
                                   : NULL;
                Expr* kexpr = call ? builtin_gaussianmatrix(call) : NULL;
                if (call) expr_free(call);
                if (kexpr) { ker_load(kexpr, &gkh, &gkw, &gk); expr_free(kexpr); }
            }
            if (gk) { convolve_dispatch(g, sm, w, h, 1, gk, gkw, gkh); free(gk); }
            else memcpy(sm, g, sizeof(double) * n);
        } else {
            memcpy(sm, g, sizeof(double) * n);
        }

        /* 2. Gradient. */
        convolve_dispatch(sm, dx, w, h, 1, kx, kwx, khx);
        convolve_dispatch(sm, dy, w, h, 1, ky, kwy, khy);
        for (size_t i = 0; i < n; i++) mag[i] = sqrt(dx[i] * dx[i] + dy[i] * dy[i]);

        /* 3. Thin. */
        canny_nms(mag, dx, dy, nms, w, h);

        /* 4. Threshold. The high threshold defaults to Otsu ON THE SUPPRESSED MAGNITUDE, not on the
         * raw magnitude: after thinning, the histogram really is "edge pixels against the rest",
         * which is the two-class problem Otsu solves. Run on the raw magnitude it would be
         * dominated by the wide ridge flanks. Low is 0.4 * high, the conventional ratio, stated
         * because it is a choice rather than a derivation. */
        double hi = thi;
        if (hi < 0.0 && !img_otsu(nms, n, &hi)) {
            /* No two classes -- a blank or perfectly uniform gradient field. No edges is the honest
             * answer, not an arbitrary threshold. */
            hi = 2.0;
        }
        double lo = 0.4 * hi;
        double* bits = malloc(sizeof(double) * n);
        if (bits) {
            canny_hysteresis(nms, bits, w, h, lo, hi, stk);
            /* Built as a Bit image: the result is 0/1 by construction. */
            Expr** rows = malloc(sizeof(Expr*) * h);
            if (rows) {
                bool ok = true;
                for (size_t y = 0; y < h; y++) rows[y] = NULL;
                for (size_t y = 0; y < h && ok; y++) {
                    Expr** cols = malloc(sizeof(Expr*) * w);
                    if (!cols) { ok = false; break; }
                    bool okc = true;
                    for (size_t x = 0; x < w; x++) {
                        cols[x] = expr_new_integer(bits[y * w + x] > 0.5 ? 1 : 0);
                        if (!cols[x]) okc = false;
                    }
                    if (okc) rows[y] = expr_new_function(expr_new_symbol(SYM_List), cols, w);
                    else for (size_t x = 0; x < w; x++) expr_free(cols[x]);
                    free(cols);
                    if (!rows[y]) ok = false;
                }
                if (ok) {
                    Expr* data = expr_new_function(expr_new_symbol(SYM_List), rows, h);
                    if (data) {
                        Expr* two[2];
                        two[0] = data;
                        two[1] = expr_new_string("Bit");
                        if (two[1]) out = expr_new_function(expr_new_symbol("Image"), two, 2);
                        else expr_free(data);
                    }
                } else for (size_t y = 0; y < h; y++) expr_free(rows[y]);
                free(rows);
            }
            free(bits);
        }
    }
    free(g); free(sm); free(dx); free(dy); free(mag); free(nms); free(stk); free(kx); free(ky);
    return out;
}

/* ---- morphology ------------------------------------------------------------
 *
 * FLAT morphology: only the SUPPORT of the structuring element matters, not its values. Dilation is
 * the maximum over the neighbourhood the element marks, erosion the minimum. A greyscale element
 * (adding its values before taking the max) is a different operator and is not what Dilation[r]
 * means anywhere; using the support keeps `Dilation[img, BoxMatrix[1]]` and `Dilation[img, 1]` the
 * same operation, which is what a caller expects.
 *
 * WHY THESE FOUR ARE WORTH TESTING TOGETHER: morphology has algebraic laws, and they are exact.
 *
 *   duality      Erosion[f, k] == 1 - Dilation[1 - f, k]      for a symmetric k
 *   ordering     Erosion <= Opening <= f <= Closing <= Dilation   pointwise, everywhere
 *   idempotence  Opening[Opening[f]] == Opening[f]             and likewise Closing
 *
 * Idempotence is the DEFINING property of an opening: it is why "open then open again" is not a
 * sharpening loop. Each of the three would fail for a different bug -- duality for a swapped
 * min/max, ordering for a wrong element centre, idempotence for a mis-composed pair -- so the three
 * together pin more than any accuracy figure could.
 *
 * PADDING IS REPLICATE, the same rule the convolutions use, and it is what makes the laws hold AT
 * THE BORDER. Zero padding would let a dilation at the edge see a black neighbour that is not
 * there, breaking `Dilation >= f` on the boundary; and it is not self-dual, so duality would fail
 * there too. Replicate is self-dual, so both survive.
 *
 * A FULL RECTANGLE IS SEPARABLE FOR MAX AND MIN, exactly as it is for a sum: the maximum over a
 * rectangle is the maximum over rows of the maxima over columns. That gives kw + kh comparisons
 * instead of kw * kh, the same win the convolution gets, and it applies to the common case since
 * BoxMatrix and an integer radius both give full rectangles.
 */
typedef enum { MORPH_DILATE, MORPH_ERODE } MorphOp;

/* Support of a structuring element: 1 where nonzero. Returns false on a malformed element. */
static bool morph_support(const Expr* e, size_t* kh, size_t* kw, unsigned char** sup,
                          bool* full) {
    size_t h = 0, w = 0; double* k = NULL;
    if (!ker_load(e, &h, &w, &k)) return false;
    unsigned char* m = malloc(h * w);
    if (!m) { free(k); return false; }
    bool all = true, any = false;
    for (size_t i = 0; i < h * w; i++) {
        m[i] = (k[i] != 0.0) ? 1 : 0;
        if (!m[i]) all = false; else any = true;
    }
    free(k);
    if (!any) { free(m); return false; }   /* an empty element has no neighbourhood */
    *kh = h; *kw = w; *sup = m; *full = all;
    return true;
}

static void morph_direct(const double* src, double* dst, size_t w, size_t h, size_t c,
                         const unsigned char* sup, size_t kh, size_t kw, MorphOp op) {
    int64_t ci = (int64_t)(kh / 2), cj = (int64_t)(kw / 2);
    for (size_t y = 0; y < h; y++)
      for (size_t x = 0; x < w; x++)
        for (size_t ch = 0; ch < c; ch++) {
            double acc = (op == MORPH_DILATE) ? -INFINITY : INFINITY;
            for (size_t i = 0; i < kh; i++) {
                size_t sy = clampi((int64_t)y + (int64_t)i - ci, h);
                for (size_t j = 0; j < kw; j++) {
                    if (!sup[i * kw + j]) continue;
                    size_t sx = clampi((int64_t)x + (int64_t)j - cj, w);
                    double v = src[(sy * w + sx) * c + ch];
                    if (op == MORPH_DILATE) { if (v > acc) acc = v; }
                    else                    { if (v < acc) acc = v; }
                }
            }
            dst[(y * w + x) * c + ch] = acc;
        }
}

/* Separable max/min over a full rectangle: rows then columns. */
static void morph_separable(const double* src, double* dst, double* tmp,
                            size_t w, size_t h, size_t c,
                            size_t kh, size_t kw, MorphOp op) {
    int64_t ci = (int64_t)(kh / 2), cj = (int64_t)(kw / 2);
    for (size_t y = 0; y < h; y++)
      for (size_t x = 0; x < w; x++)
        for (size_t ch = 0; ch < c; ch++) {
            double acc = (op == MORPH_DILATE) ? -INFINITY : INFINITY;
            for (size_t j = 0; j < kw; j++) {
                size_t sx = clampi((int64_t)x + (int64_t)j - cj, w);
                double v = src[(y * w + sx) * c + ch];
                if (op == MORPH_DILATE) { if (v > acc) acc = v; }
                else                    { if (v < acc) acc = v; }
            }
            tmp[(y * w + x) * c + ch] = acc;
        }
    for (size_t y = 0; y < h; y++)
      for (size_t x = 0; x < w; x++)
        for (size_t ch = 0; ch < c; ch++) {
            double acc = (op == MORPH_DILATE) ? -INFINITY : INFINITY;
            for (size_t i = 0; i < kh; i++) {
                size_t sy = clampi((int64_t)y + (int64_t)i - ci, h);
                double v = tmp[(sy * w + x) * c + ch];
                if (op == MORPH_DILATE) { if (v > acc) acc = v; }
                else                    { if (v < acc) acc = v; }
            }
            dst[(y * w + x) * c + ch] = acc;
        }
}

/* Run one morphological pass over a loaded buffer. */
static void morph_run(const double* src, double* dst, size_t w, size_t h, size_t c,
                      const unsigned char* sup, size_t kh, size_t kw, bool full, MorphOp op) {
    if (full) {
        double* tmp = malloc(sizeof(double) * w * h * c);
        if (tmp) { morph_separable(src, dst, tmp, w, h, c, kh, kw, op); free(tmp); return; }
    }
    morph_direct(src, dst, w, h, c, sup, kh, kw, op);
}

/* Read the second argument: an integer radius (a full (2r+1) square) or an explicit element. */
static bool morph_element(Expr* spec, size_t* kh, size_t* kw, unsigned char** sup, bool* full) {
    double r = 0.0, im = 0.0;
    if (na_read_scalar(spec, &r, &im) && im == 0.0) {
        if (!(r >= 0.0) || r != floor(r) || r > 256.0) return false;
        size_t n = 2 * (size_t)r + 1;
        unsigned char* m = malloc(n * n);
        if (!m) return false;
        memset(m, 1, n * n);
        *kh = *kw = n; *sup = m; *full = true;
        return true;
    }
    return morph_support(spec, kh, kw, sup, full);
}

/* One-pass operators (Dilation, Erosion) and two-pass ones (Opening, Closing). */
static Expr* morph_builtin(Expr* res, MorphOp first, bool two_pass) {
    if (res->data.function.arg_count != 2) return NULL;
    size_t kh = 0, kw = 0; unsigned char* sup = NULL; bool full = false;
    if (!morph_element(res->data.function.args[1], &kh, &kw, &sup, &full)) return NULL;

    size_t w = 0, h = 0, c = 0; double* src = NULL;
    if (!image_load(res->data.function.args[0], &w, &h, &c, &src)) { free(sup); return NULL; }
    size_t n = w * h * c;
    double* a = malloc(sizeof(double) * n);
    Expr* out = NULL;
    if (a) {
        morph_run(src, a, w, h, c, sup, kh, kw, full, first);
        if (!two_pass) {
            out = image_build_real(a, w, h, c);
        } else {
            /* Opening is erode-then-dilate, closing dilate-then-erode: the SAME element both
             * times, which is what makes the pair idempotent. Using a different element for the
             * second pass would still smooth, and would no longer be an opening. */
            double* b = malloc(sizeof(double) * n);
            if (b) {
                MorphOp second = (first == MORPH_ERODE) ? MORPH_DILATE : MORPH_ERODE;
                morph_run(a, b, w, h, c, sup, kh, kw, full, second);
                out = image_build_real(b, w, h, c);
                free(b);
            }
        }
    }
    free(src); free(a); free(sup);
    return out;
}

static Expr* builtin_dilation(Expr* res) { return morph_builtin(res, MORPH_DILATE, false); }
static Expr* builtin_erosion(Expr* res)  { return morph_builtin(res, MORPH_ERODE,  false); }
static Expr* builtin_opening(Expr* res)  { return morph_builtin(res, MORPH_ERODE,  true);  }
static Expr* builtin_closing(Expr* res)  { return morph_builtin(res, MORPH_DILATE, true);  }

void imagefilter_init(void) {
    symtab_add_builtin("Dilation", builtin_dilation);
    symtab_get_def("Dilation")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("Dilation",
        "Dilation[image, r] gives the maximum over a (2r+1) x (2r+1) square neighbourhood; "
        "Dilation[image, elem] uses the SUPPORT of the matrix elem -- its nonzero positions -- as "
        "the neighbourhood. This is flat morphology: the element's values do not enter the "
        "maximum, which is what keeps Dilation[img, BoxMatrix[1]] and Dilation[img, 1] the same "
        "operation. Padding replicates the border, the same rule the convolutions use, which is "
        "what makes Dilation >= image hold at the edges too. A full rectangle is separable for the "
        "maximum exactly as for a sum, so it costs kw + kh comparisons rather than kw * kh.");

    symtab_add_builtin("Erosion", builtin_erosion);
    symtab_get_def("Erosion")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("Erosion",
        "Erosion[image, r] gives the minimum over a (2r+1) x (2r+1) square neighbourhood; "
        "Erosion[image, elem] uses the support of elem. Dual to Dilation: for a symmetric element, "
        "Erosion[f, k] equals 1 - Dilation[1 - f, k] exactly, which holds at the border only "
        "because the replicate padding is itself self-dual.");

    symtab_add_builtin("Opening", builtin_opening);
    symtab_get_def("Opening")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("Opening",
        "Opening[image, r] erodes then dilates with the same element, removing bright features "
        "smaller than it while leaving larger ones close to their original size. IDEMPOTENT: "
        "Opening[Opening[f]] equals Opening[f], which is the defining property and the reason "
        "opening twice is not a sharpening loop.");

    symtab_add_builtin("Closing", builtin_closing);
    symtab_get_def("Closing")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("Closing",
        "Closing[image, r] dilates then erodes with the same element, filling dark features "
        "smaller than it. Idempotent, like Opening, and the two bracket the image: "
        "Erosion <= Opening <= image <= Closing <= Dilation pointwise everywhere.");
    symtab_add_builtin("EdgeDetect", builtin_edgedetect);
    symtab_get_def("EdgeDetect")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("EdgeDetect",
        "EdgeDetect[image] finds edges by the Canny algorithm, giving a \"Bit\" image. "
        "EdgeDetect[image, r] sets the Gaussian smoothing radius (default 2; 0 means no "
        "smoothing). EdgeDetect[image, r, t] sets the high threshold explicitly. Four stages: "
        "smooth, because a derivative amplifies noise; gradient by the normalised Sobel pair; "
        "non-maximum suppression along the gradient direction, which is what makes an edge ONE "
        "pixel wide rather than a thick band; and hysteresis, keeping any pixel above the high "
        "threshold plus any above 0.4 of it that is 8-connected to one, so a real edge survives "
        "its faint stretches while isolated weak responses do not. The high threshold defaults to "
        "Otsu's method applied to the SUPPRESSED magnitude, where the two classes really are edge "
        "against non-edge; on the raw magnitude it would be dominated by the ridge flanks.");
    symtab_add_builtin("DerivativeFilter", builtin_derivativefilter);
    symtab_get_def("DerivativeFilter")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("DerivativeFilter",
        "DerivativeFilter[image, {n, m}] gives the n-th derivative down the rows and the m-th "
        "across the columns, each order from 0 to 2. The kernel is a separable outer product of "
        "1-D stencils: order 0 is the smoothing {1,2,1}/4, order 1 the central difference "
        "{-1,0,1}/2, order 2 the second difference {1,-2,1}. So {0,1} is Sobel-x and {1,0} is "
        "Sobel-y. The stencils are NORMALISED, unlike the raw integer Sobel kernels, which report "
        "a gradient eight times the true slope -- harmless when only the ranking of edges matters, "
        "and wrong for anything that reads the number. On f(x) = c x the first derivative gives "
        "exactly c. The result is a \"Real\" image.");

    symtab_add_builtin("GradientFilter", builtin_gradientfilter);
    symtab_get_def("GradientFilter")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("GradientFilter",
        "GradientFilter[image] gives the gradient magnitude Sqrt[dx^2 + dy^2], using the "
        "normalised Sobel derivatives of DerivativeFilter. The magnitude rather than |dx| + |dy| "
        "because it is ROTATION INVARIANT: an edge at 45 degrees reports the same strength as one "
        "at 0, where the absolute sum would report it sqrt(2) times stronger and so bias every "
        "downstream threshold by orientation. A colour image is reduced to luminance first and "
        "differentiated once, rather than differentiated per channel and combined by some "
        "arbitrary rule.");
    symtab_add_builtin("FindThreshold", builtin_findthreshold);
    symtab_get_def("FindThreshold")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("FindThreshold",
        "FindThreshold[image] gives a threshold separating the image into two classes, by "
        "Otsu's method: the level maximising the BETWEEN-class variance "
        "w0 w1 (mu0 - mu1)^2, which is algebraically the same as minimising the weighted "
        "within-class variance but needs only one incremental pass over a 256-bin histogram. "
        "A colour image is reduced to Rec. 601 luminance first. Returns unevaluated for an "
        "image whose pixels are all identical, since no threshold splits one cluster into two.");

    symtab_add_builtin("Binarize", builtin_binarize);
    symtab_get_def("Binarize")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("Binarize",
        "Binarize[image] thresholds image by Otsu's method (see FindThreshold), giving a "
        "\"Bit\" image. Binarize[image, t] thresholds at t. A pixel STRICTLY ABOVE the "
        "threshold becomes 1, so a pixel exactly at it becomes 0 -- which matters, because "
        "\"above\" and \"at or above\" differ on exactly the pixels a threshold was chosen "
        "to sit between. A colour image is reduced to luminance first.");

    symtab_add_builtin("ColorConvert", builtin_colorconvert);
    symtab_get_def("ColorConvert")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("ColorConvert",
        "ColorConvert[image, \"Grayscale\"] converts to greyscale by Rec. 601 luminance, "
        "0.299 R + 0.587 G + 0.114 B -- weighted rather than averaged because the eye is far "
        "more sensitive to green than to blue, so an unweighted mean would put a saturated blue "
        "and a saturated green at the same brightness. Only \"Grayscale\" is accepted; other "
        "colour spaces carry their own white-point and transfer-function decisions, and "
        "accepting the name while doing something approximate would be worse than declining.");
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
