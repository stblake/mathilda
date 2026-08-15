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
#ifdef USE_FFTW
#include <fftw3.h>
#endif

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
/* ---- convolution through the transform -------------------------------------
 *
 * A direct convolution costs w*h*kw*kh multiply-adds, which at 512x512 with a 32x32 kernel is 275
 * million and measured 197 ms -- against SciPy's 143 for the same dense loop. Neither is doing
 * anything clever, and the convolution theorem is the thing to do instead: O(n log n) regardless of
 * kernel size.
 *
 * THE BORDER IS THE WHOLE DIFFICULTY. A transform gives CIRCULAR convolution, while every filter here
 * replicates the edge -- so the wrap has to be made impossible rather than accepted. Padding the
 * image to PH = h + kh - 1 by replication does exactly that: the outputs that matter then sit at
 * indices kh-1 .. kh-2+h of the linear result, and each of them reads only padded rows that exist, so
 * no output touches a wrapped one. Rounding the transform up to a 5-smooth size adds zeros beyond the
 * pad, which cannot reach those indices either.
 *
 * The index algebra, written out because an off-by-one here is a shifted image rather than an error.
 * The direct loop computes dst[y] = sum_i src[clamp(y - i + ci)] k[i] with ci = kh/2. Define
 *
 *     P[a] = src[clamp(a - (kh-1) + ci)]
 *
 * then sum_i P[(y + kh - 1) - i] k[i] = sum_i src[clamp(y - i + ci)] k[i] = dst[y], and the left side
 * is the linear convolution of P with k at index y + kh - 1. So no kernel flip is needed: this is a
 * convolution, and the output is a shifted window of the transform's result.
 *
 * The kernel's transform is computed ONCE and reused across channels; only the image transform and the
 * inverse are per channel.
 */
#ifdef USE_FFTW
/* The next size >= n that is a product of 2, 3 and 5, where FFTW's algorithms are the good ones. */
static size_t fft_smooth(size_t n) {
    if (n < 2) return 2;
    for (;;) {
        size_t m = n;
        while (m % 2 == 0) m /= 2;
        while (m % 3 == 0) m /= 3;
        while (m % 5 == 0) m /= 5;
        if (m == 1) return n;
        n++;
    }
}

static bool convolve_fft(const double* src, double* dst, size_t w, size_t h, size_t c,
                         const double* k, size_t kw, size_t kh) {
    int64_t ci = (int64_t)(kh / 2), cj = (int64_t)(kw / 2);
    size_t PH = h + kh - 1, PW = w + kw - 1;
    size_t NH = fft_smooth(PH), NW = fft_smooth(PW);
    size_t NC = NW / 2 + 1;                       /* r2c packs the last axis */
    double scale = 1.0 / (double)(NH * NW);

    double* rbuf = fftw_malloc(sizeof(double) * NH * NW);
    double* kbuf = fftw_malloc(sizeof(double) * NH * NW);
    fftw_complex* F = fftw_malloc(sizeof(fftw_complex) * NH * NC);
    fftw_complex* G = fftw_malloc(sizeof(fftw_complex) * NH * NC);
    if (!rbuf || !kbuf || !F || !G) {
        fftw_free(rbuf); fftw_free(kbuf); fftw_free(F); fftw_free(G);
        return false;                             /* caller falls back to the direct loop */
    }

    /* Plans are built per call on the actual pointers, matching fourier.c: FFTW_ESTIMATE does not
     * touch the arrays while planning, and a plan built on the buffer it will run on cannot meet a
     * different alignment later. */
    int nn[2]; nn[0] = (int)NH; nn[1] = (int)NW;
    fftw_plan pf = fftw_plan_dft_r2c(2, nn, rbuf, F, FFTW_ESTIMATE);
    fftw_plan pk = fftw_plan_dft_r2c(2, nn, kbuf, G, FFTW_ESTIMATE);
    fftw_plan pb = fftw_plan_dft_c2r(2, nn, F, rbuf, FFTW_ESTIMATE);
    if (!pf || !pk || !pb) {
        if (pf) fftw_destroy_plan(pf);
        if (pk) fftw_destroy_plan(pk);
        if (pb) fftw_destroy_plan(pb);
        fftw_free(rbuf); fftw_free(kbuf); fftw_free(F); fftw_free(G);
        return false;
    }

    /* The kernel at the origin, once for every channel. */
    memset(kbuf, 0, sizeof(double) * NH * NW);
    for (size_t i = 0; i < kh; i++)
        for (size_t j = 0; j < kw; j++) kbuf[i * NW + j] = k[i * kw + j];
    fftw_execute(pk);

    for (size_t ch = 0; ch < c; ch++) {
        memset(rbuf, 0, sizeof(double) * NH * NW);
        for (size_t y = 0; y < PH; y++) {
            size_t sy = clampi((int64_t)y - (int64_t)(kh - 1) + ci, h);
            for (size_t x = 0; x < PW; x++) {
                size_t sx = clampi((int64_t)x - (int64_t)(kw - 1) + cj, w);
                rbuf[y * NW + x] = src[(sy * w + sx) * c + ch];
            }
        }
        fftw_execute(pf);
        for (size_t i = 0; i < NH * NC; i++) {
            double ar = F[i][0], ai = F[i][1], br = G[i][0], bi = G[i][1];
            F[i][0] = ar * br - ai * bi;
            F[i][1] = ar * bi + ai * br;
        }
        fftw_execute(pb);
        for (size_t y = 0; y < h; y++)
            for (size_t x = 0; x < w; x++)
                dst[(y * w + x) * c + ch] = rbuf[(y + kh - 1) * NW + (x + kw - 1)] * scale;
    }

    fftw_destroy_plan(pf); fftw_destroy_plan(pk); fftw_destroy_plan(pb);
    fftw_free(rbuf); fftw_free(kbuf); fftw_free(F); fftw_free(G);
    return true;
}

/* Is the transform cheaper than the dense loop here?
 *
 * Direct is w*h*c taps of kw*kh each. The transform is one kernel forward, then a forward and an
 * inverse per channel, over NH*NW points at roughly 5*N*log2(N) flops each, plus the pointwise
 * multiply. FFT_TAX is the one empirical number, and it was FITTED TO A MEASUREMENT rather than
 * guessed -- the first value here was guessed, at 6.0, and it was an order of magnitude out.
 *
 * The measurement, 512x512 with non-separable kernels, both paths forced:
 *
 *      k       direct      transform
 *      3x3      1.05 ms      ~2.7 ms
 *      5x5      2.60         ~2.7
 *      7x7      5.56         ~2.7
 *      9x9     10.61         ~2.7
 *     15x15    35.16         ~2.7
 *     21x21    78.90          2.63
 *     32x32   198.69          2.77
 *
 * The transform is flat in the kernel -- that is the whole point of it -- so the crossover is wherever
 * the dense loop passes ~2.7 ms, which is just under 5x5. Switching exactly there would trade a dense
 * loop for a transform of the same cost, so the tax is set to land the switch at 7x7, the first size
 * where the transform is a clear win (5.56 -> 2.7, about 2x). At 6.0 the model switched between 15x15
 * and 21x21 instead, leaving 35 ms on the table at 15x15 where 2.7 was available. */
#define FFT_TAX 0.6
static bool fft_is_cheaper(size_t w, size_t h, size_t c, size_t kw, size_t kh) {
    double direct = (double)w * (double)h * (double)c * (double)kw * (double)kh;
    size_t NH = fft_smooth(h + kh - 1), NW = fft_smooth(w + kw - 1);
    double N = (double)NH * (double)NW;
    double lg = log2(N > 2.0 ? N : 2.0);
    double fft = FFT_TAX * (2.0 * (double)c + 1.0) * N * lg;
    return fft < direct;
}
#endif /* USE_FFTW */

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
#ifdef USE_FFTW
    /* Not separable: the transform is the next thing to try, when the cost model says it wins. A
     * separable kernel is never sent here -- kw + kh taps beats any transform. */
    if (fft_is_cheaper(w, h, c, kw, kh) && convolve_fft(src, dst, w, h, c, k, kw, kh))
        return true;
#endif
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
                    /* A NULL support means a full rectangle -- the shape the separable path takes,
                     * used when its scratch allocation fails. */
                    if (sup && !sup[i * kw + j]) continue;
                    size_t sx = clampi((int64_t)x + (int64_t)j - cj, w);
                    double v = src[(sy * w + sx) * c + ch];
                    if (op == MORPH_DILATE) { if (v > acc) acc = v; }
                    else                    { if (v < acc) acc = v; }
                }
            }
            dst[(y * w + x) * c + ch] = acc;
        }
}

/* van Herk--Gil-Werman: a 1-D dilation or erosion in THREE comparisons per pixel, whatever the
 * structuring element's width.
 *
 * The separable form already reduced a k x k rectangle to 2k comparisons per pixel. This reduces each
 * axis to a constant, so the whole operation stops depending on the radius at all -- which is the
 * difference between morphology being usable at r = 20 and not. A benchmark put the crossover against
 * scipy at r ~ 4-6 and this is what closes it.
 *
 * THE IDEA. Cut the line into blocks of exactly k. Within each block compute a PREFIX maximum running
 * forwards and a SUFFIX maximum running backwards. Any window of width k straddles exactly two adjacent
 * blocks -- it cannot span three, because it is the same width as a block -- so its maximum is
 * max(suffix at the window's start, prefix at the window's end). Two lookups and one comparison,
 * regardless of k. The prefix and suffix arrays each cost one pass, so the total is three comparisons
 * per pixel amortised.
 *
 * It is EXACT, not approximate: max and min are associative and idempotent, so splitting a window at a
 * block boundary and recombining loses nothing. That matters here more than usual -- it means the fast
 * path must agree BIT-EXACTLY with the naive one, and a test asserts exactly that on random data at
 * several radii rather than merely checking the fast path against its own expectations.
 *
 * The caller passes a line already padded by r on each side, so boundary handling stays where it was
 * and this routine never has to know about it. */
static void morph_1d_vanherk(const double* in, double* out, size_t n, size_t k,
                             double* pre, double* suf, MorphOp op) {
    size_t m = n + k - 1;                 /* padded length; window i covers in[i .. i+k-1] */
    for (size_t b = 0; b < m; b += k) {
        size_t end = (b + k < m) ? b + k : m;   /* last block may be short */
        /* Prefix maxima, forwards within the block. */
        double acc = in[b];
        pre[b] = acc;
        for (size_t i = b + 1; i < end; i++) {
            double v = in[i];
            if (op == MORPH_DILATE) { if (v > acc) acc = v; }
            else                    { if (v < acc) acc = v; }
            pre[i] = acc;
        }
        /* Suffix maxima, backwards within the block. */
        acc = in[end - 1];
        suf[end - 1] = acc;
        for (size_t i = end - 1; i > b; i--) {
            double v = in[i - 1];
            if (op == MORPH_DILATE) { if (v > acc) acc = v; }
            else                    { if (v < acc) acc = v; }
            suf[i - 1] = acc;
        }
    }
    for (size_t i = 0; i < n; i++) {
        double a = suf[i], bb = pre[i + k - 1];
        out[i] = (op == MORPH_DILATE) ? (a > bb ? a : bb) : (a < bb ? a : bb);
    }
}

/* Separable max/min over a full rectangle: rows then columns. */
static void morph_separable(const double* src, double* dst, double* tmp,
                            size_t w, size_t h, size_t c,
                            size_t kh, size_t kw, MorphOp op) {
    size_t rh = kh / 2, rw = kw / 2;
    size_t maxline = (w > h ? w : h) + (kw > kh ? kw : kh);
    double* line = malloc(sizeof(double) * (maxline + 2));
    double* res  = malloc(sizeof(double) * (maxline + 2));
    double* pre  = malloc(sizeof(double) * (maxline + 2));
    double* suf  = malloc(sizeof(double) * (maxline + 2));
    if (!line || !res || !pre || !suf) {
        /* Without scratch there is no fast path; the caller's direct routine still answers. */
        free(line); free(res); free(pre); free(suf);
        morph_direct(src, dst, w, h, c, NULL, kh, kw, op);
        return;
    }

    /* Horizontal pass. The line is padded by rw on each side with the edge value, so the van Herk
     * routine sees a plain array and the replicate rule stays in one place. */
    for (size_t y = 0; y < h; y++)
        for (size_t ch = 0; ch < c; ch++) {
            for (size_t t = 0; t < rw; t++)          line[t] = src[(y * w + 0) * c + ch];
            for (size_t x = 0; x < w; x++)           line[rw + x] = src[(y * w + x) * c + ch];
            for (size_t t = 0; t < rw; t++)          line[rw + w + t] = src[(y * w + (w - 1)) * c + ch];
            morph_1d_vanherk(line, res, w, kw, pre, suf, op);
            for (size_t x = 0; x < w; x++) tmp[(y * w + x) * c + ch] = res[x];
        }

    /* Vertical pass over the horizontal result. */
    for (size_t x = 0; x < w; x++)
        for (size_t ch = 0; ch < c; ch++) {
            for (size_t t = 0; t < rh; t++)          line[t] = tmp[(0 * w + x) * c + ch];
            for (size_t y = 0; y < h; y++)           line[rh + y] = tmp[(y * w + x) * c + ch];
            for (size_t t = 0; t < rh; t++)          line[rh + h + t] = tmp[((h - 1) * w + x) * c + ch];
            morph_1d_vanherk(line, res, h, kh, pre, suf, op);
            for (size_t y = 0; y < h; y++) dst[(y * w + x) * c + ch] = res[y];
        }
    free(line); free(res); free(pre); free(suf);
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

/* ---- connected components -------------------------------------------------
 *
 * Two-pass union-find, the classic algorithm, and the two passes do genuinely different jobs.
 *
 * The FIRST pass walks in raster order and can only see neighbours it has already visited -- for
 * 8-connectivity that is W, NW, N, NE. A U-shaped region is the case that makes one pass
 * insufficient: the two arms get different labels because nothing has yet connected them, and the
 * base at the bottom is where they turn out to be the same component. Union-find records that
 * equivalence; the second pass applies it.
 *
 * The SECOND pass also RELABELS to 1..k in raster order of first appearance. Without that the labels
 * would be whatever the first pass happened to allocate, with gaps where two provisional labels were
 * later merged -- so `Max[labels]` would not be the component count, and no test could pin a label
 * pattern. Contiguous labels in scan order are what make both possible.
 *
 * CONNECTIVITY IS THE DISCRIMINATING PROPERTY, and it is exactly testable: two pixels touching only
 * at a corner are ONE component under 8-connectivity and TWO under 4. Every other property --
 * background staying 0, labels being contiguous, a single blob labelling 1 -- holds under either
 * rule, so a wrong default would pass all of them. CornerNeighbors -> True (8-connectivity) is the
 * default, matching Mathematica.
 *
 * RETURNS A PLAIN INTEGER MATRIX, NOT AN IMAGE, and this is a deliberate deviation. Mathematica
 * returns an Image; here that would be actively wrong, because Image type inference would call a
 * label array of 1..12 a "Byte" image and ImageData would then divide every label by 255. Labels are
 * indices, not brightnesses, and the one thing that must not happen to them is being scaled.
 */
static size_t cc_find(size_t* parent, size_t x) {
    while (parent[x] != x) { parent[x] = parent[parent[x]]; x = parent[x]; }
    return x;
}

static void cc_union(size_t* parent, size_t a, size_t b) {
    a = cc_find(parent, a); b = cc_find(parent, b);
    if (a == b) return;
    /* Union by lower index rather than by rank. Rank would be faster asymptotically; lower-index
     * keeps the representative deterministic, and the second pass relabels anyway so the extra
     * depth costs nothing measurable at image scale. */
    if (a < b) parent[b] = a; else parent[a] = b;
}

/* MorphologicalComponents[image] / [image, t] / with CornerNeighbors -> False */
static Expr* builtin_morphologicalcomponents(Expr* res) {
    size_t argc = res->data.function.arg_count;
    if (argc < 1 || argc > 3) return NULL;

    double thr = 0.0, im = 0.0;
    bool corners = true;
    size_t next_arg = 1;
    if (argc >= 2) {
        Expr* a1 = res->data.function.args[1];
        /* A threshold, or the option -- the option may appear in either position. */
        if (a1 && a1->type == EXPR_FUNCTION && a1->data.function.arg_count == 2
            && a1->data.function.head && a1->data.function.head->type == EXPR_SYMBOL
            && (a1->data.function.head->data.symbol.name == SYM_Rule
                || a1->data.function.head->data.symbol.name == SYM_RuleDelayed)) {
            /* handled below as an option */
        } else {
            if (!na_read_scalar(a1, &thr, &im) || im != 0.0) return NULL;
            next_arg = 2;
        }
    }
    for (size_t i = next_arg; i < argc; i++) {
        Expr* o = res->data.function.args[i];
        if (!o || o->type != EXPR_FUNCTION || o->data.function.arg_count != 2) return NULL;
        Expr* hd = o->data.function.head;
        if (!hd || hd->type != EXPR_SYMBOL) return NULL;
        if (hd->data.symbol.name != SYM_Rule && hd->data.symbol.name != SYM_RuleDelayed)
            return NULL;
        Expr* lhs = o->data.function.args[0];
        Expr* rhs = o->data.function.args[1];
        if (!lhs || lhs->type != EXPR_SYMBOL
            || strcmp(lhs->data.symbol.name, "CornerNeighbors") != 0) return NULL;
        if (!rhs || rhs->type != EXPR_SYMBOL) return NULL;
        if (strcmp(rhs->data.symbol.name, "True") == 0) corners = true;
        else if (strcmp(rhs->data.symbol.name, "False") == 0) corners = false;
        else return NULL;
    }

    size_t w = 0, h = 0; double* g = NULL;
    if (!img_grey_plane(res->data.function.args[0], &w, &h, &g)) return NULL;
    size_t n = w * h;

    size_t* lab = calloc(n, sizeof(size_t));
    size_t* parent = malloc(sizeof(size_t) * (n + 1));
    Expr* out = NULL;
    if (lab && parent) {
        size_t next = 1;
        parent[0] = 0;
        /* First pass: provisional labels, unioning already-seen neighbours. */
        for (size_t y = 0; y < h; y++) {
            for (size_t x = 0; x < w; x++) {
                size_t i = y * w + x;
                if (!(g[i] > thr)) continue;                 /* background */
                size_t best = 0;
                /* W, then NW, N, NE -- the visited half of the neighbourhood. */
                int64_t dxs[4] = { -1, -1, 0, 1 };
                int64_t dys[4] = {  0, -1, -1, -1 };
                for (int k = 0; k < 4; k++) {
                    if (!corners && (dxs[k] != 0 && dys[k] != 0)) continue;  /* 4-connectivity */
                    int64_t nx = (int64_t)x + dxs[k], ny = (int64_t)y + dys[k];
                    if (nx < 0 || (size_t)nx >= w || ny < 0) continue;
                    size_t j = (size_t)ny * w + (size_t)nx;
                    if (lab[j] == 0) continue;
                    if (best == 0) best = lab[j];
                    else cc_union(parent, best, lab[j]);
                }
                if (best == 0) { parent[next] = next; lab[i] = next; next++; }
                else lab[i] = cc_find(parent, best);
            }
        }
        /* Second pass: resolve to roots and RELABEL to 1..k in raster order of first appearance. */
        size_t* remap = calloc(next + 1, sizeof(size_t));
        if (remap) {
            size_t k = 0;
            for (size_t i = 0; i < n; i++) {
                if (lab[i] == 0) continue;
                size_t r = cc_find(parent, lab[i]);
                if (remap[r] == 0) remap[r] = ++k;
                lab[i] = remap[r];
            }
            /* A plain integer matrix -- see the note above on why not an Image. */
            Expr** rows = malloc(sizeof(Expr*) * h);
            if (rows) {
                bool ok = true;
                for (size_t y = 0; y < h; y++) rows[y] = NULL;
                for (size_t y = 0; y < h && ok; y++) {
                    Expr** cols = malloc(sizeof(Expr*) * w);
                    if (!cols) { ok = false; break; }
                    bool okc = true;
                    for (size_t x = 0; x < w; x++) {
                        cols[x] = expr_new_integer((int64_t)lab[y * w + x]);
                        if (!cols[x]) okc = false;
                    }
                    if (okc) rows[y] = expr_new_function(expr_new_symbol(SYM_List), cols, w);
                    else for (size_t x = 0; x < w; x++) expr_free(cols[x]);
                    free(cols);
                    if (!rows[y]) ok = false;
                }
                if (ok) out = expr_new_function(expr_new_symbol(SYM_List), rows, h);
                else for (size_t y = 0; y < h; y++) expr_free(rows[y]);
                free(rows);
            }
            free(remap);
        }
    }
    free(g); free(lab); free(parent);
    return out;
}

/* ---- rank filters ---------------------------------------------------------
 *
 * THE MEDIAN IS THE ONE OPERATOR HERE THAT IS NOT SEPARABLE, and that is worth stating plainly after
 * a week of leaning on separability. A sum, a maximum and a minimum all decompose: the sum over a
 * rectangle is the sum over rows of the sums over columns, and likewise for max and min, because all
 * three are ASSOCIATIVE and COMMUTATIVE reductions that ignore how the values are grouped. The median
 * is not: it depends on the RANK of a value within the whole window, and grouping destroys rank
 * information. The median of row-medians is a real filter, it is fast, and it is not the median --
 *
 *   {{1, 2, 9}, {3, 4, 5}, {6, 7, 8}}   true median of all nine = 5
 *                                        row medians {2, 4, 7}, their median = 4
 *
 * -- so the separable version is simply wrong, and a test pins the true value on exactly that window.
 * MeanFilter, by contrast, IS a convolution with a normalised box, and a test asserts the identity
 * against ImageConvolve rather than reimplementing it.
 *
 * WHY A MEDIAN AT ALL, when a Gaussian is cheaper: a median removes an isolated outlier EXACTLY,
 * where a Gaussian only attenuates it and smears it over the neighbourhood. On salt-and-pepper noise
 * that is the difference between clean and merely blurred, and it is the discriminating test here --
 * a single bright pixel in a constant field vanishes completely under the median and leaves a visible
 * bump under the mean.
 *
 * Insertion sort on the window rather than a histogram or a running median. For the radii people
 * actually use -- 1 to 3, so 9 to 49 values -- insertion sort on a small contiguous array beats both
 * on constant factors, and it is exact and obviously correct. A histogram median is the right choice
 * only for 8-bit data at large radii, and this works on doubles.
 */
static double window_median(double* buf, size_t n) {
    /* Insertion sort, then the middle element. For even n the LOWER middle is taken rather than
     * averaging the two: averaging would invent a value not present in the window, which for a
     * rank filter is exactly what a caller does not want -- the point of a median is that its
     * output is one of its inputs. Odd windows are the normal case anyway. */
    for (size_t i = 1; i < n; i++) {
        double v = buf[i];
        size_t j = i;
        while (j > 0 && buf[j - 1] > v) { buf[j] = buf[j - 1]; j--; }
        buf[j] = v;
    }
    return buf[(n - 1) / 2];
}

/* MedianFilter[image, r] -- median over a (2r+1) square neighbourhood. */
static Expr* builtin_medianfilter(Expr* res) {
    if (res->data.function.arg_count != 2) return NULL;
    double r = 0.0, im = 0.0;
    if (!na_read_scalar(res->data.function.args[1], &r, &im) || im != 0.0) return NULL;
    if (!(r >= 0.0) || r != floor(r) || r > 64.0) return NULL;
    size_t rr = (size_t)r, k = 2 * rr + 1;

    size_t w = 0, h = 0, c = 0; double* src = NULL;
    if (!image_load(res->data.function.args[0], &w, &h, &c, &src)) return NULL;
    double* dst = malloc(sizeof(double) * w * h * c);
    double* win = malloc(sizeof(double) * k * k);
    Expr* out = NULL;
    if (dst && win) {
        for (size_t y = 0; y < h; y++)
          for (size_t x = 0; x < w; x++)
            for (size_t ch = 0; ch < c; ch++) {
                size_t m = 0;
                for (size_t i = 0; i < k; i++) {
                    size_t sy = clampi((int64_t)y + (int64_t)i - (int64_t)rr, h);
                    for (size_t j = 0; j < k; j++) {
                        size_t sx = clampi((int64_t)x + (int64_t)j - (int64_t)rr, w);
                        win[m++] = src[(sy * w + sx) * c + ch];
                    }
                }
                dst[(y * w + x) * c + ch] = window_median(win, m);
            }
        out = image_build_real(dst, w, h, c);
    }
    free(src); free(dst); free(win);
    return out;
}

/* MeanFilter[image, r] -- the mean over a (2r+1) square, which IS a convolution with a normalised
 * box, so it is implemented as one rather than as a second averaging loop. Two implementations of
 * one identity is how the identity quietly stops holding -- the same reason GaussianFilter routes
 * through ImageConvolve. */
static Expr* builtin_meanfilter(Expr* res) {
    if (res->data.function.arg_count != 2) return NULL;
    double r = 0.0, im = 0.0;
    if (!na_read_scalar(res->data.function.args[1], &r, &im) || im != 0.0) return NULL;
    if (!(r >= 0.0) || r != floor(r) || r > 256.0) return NULL;
    size_t n = 2 * (size_t)r + 1;
    double* k = malloc(sizeof(double) * n * n);
    if (!k) return NULL;
    double inv = 1.0 / (double)(n * n);
    for (size_t i = 0; i < n * n; i++) k[i] = inv;

    size_t w = 0, h = 0, c = 0; double* src = NULL;
    if (!image_load(res->data.function.args[0], &w, &h, &c, &src)) { free(k); return NULL; }
    double* dst = malloc(sizeof(double) * w * h * c);
    Expr* out = NULL;
    if (dst) {
        convolve_dispatch(src, dst, w, h, c, k, n, n);
        out = image_build_real(dst, w, h, c);
    }
    free(src); free(dst); free(k);
    return out;
}

/* ---- distance transform ---------------------------------------------------
 *
 * EXACT Euclidean, by Felzenszwalb and Huttenlocher's algorithm, and the exactness is the reason for
 * choosing it over the obvious alternative. The classic two-pass CHAMFER transform propagates local
 * step costs and is approximate: it cannot represent sqrt(2) with integer steps, so a diagonal
 * distance comes out a few percent wrong. That error is invisible on a picture and fatal to a test --
 * it would force a tolerance where an equality is available.
 *
 * The exact version costs the same asymptotically. It computes, per axis,
 *
 *     D(x) = min over y of ( (x - y)^2 + f(y) )
 *
 * which is the LOWER ENVELOPE OF PARABOLAS: each y contributes a parabola of identical shape
 * translated to y and raised by f(y), and the transform is their pointwise minimum. Because the
 * parabolas all have the same curvature, any two intersect exactly once, so the envelope can be built
 * in one left-to-right sweep maintaining a stack of the parabolas still visible -- O(n) per row, no
 * sorting.
 *
 * SEPARABILITY IS EXACT HERE, unlike the median's, and for a reason worth naming: squared Euclidean
 * distance is a SUM over the axes, (dx^2 + dy^2), so minimising it decomposes into minimising per axis.
 * A median does not decompose because rank is not a sum. So the same word covers an exact
 * factorisation in one case and a wrong shortcut in the other, and which it is depends on whether the
 * quantity being reduced is additive.
 *
 * The test that separates exact from approximate is a 3-4-5 triangle: a single background pixel with a
 * foreground pixel three across and four down must read EXACTLY 5. A chamfer transform gives about
 * 5.03, and no picture would ever show the difference.
 */
#define DT_INF 1.0e30

/* One-dimensional squared-distance transform of f into d, using scratch v (indices) and z
 * (intersections). All four arrays are length n, z needs n + 1. */
static void dt_1d(const double* f, double* d, size_t n, size_t* v, double* z) {
    size_t k = 0;
    v[0] = 0;
    z[0] = -DT_INF;
    z[1] = DT_INF;
    for (size_t q = 1; q < n; q++) {
        /* Intersection of the parabola from q with the one currently on top of the stack. */
        double s = ((f[q] + (double)q * (double)q)
                    - (f[v[k]] + (double)v[k] * (double)v[k]))
                   / (2.0 * (double)q - 2.0 * (double)v[k]);
        /* If it lies left of where the top parabola became visible, that one is entirely hidden by
         * this one and is popped. The loop, not an if: one new parabola can hide several. */
        while (s <= z[k]) {
            if (k == 0) break;
            k--;
            s = ((f[q] + (double)q * (double)q)
                 - (f[v[k]] + (double)v[k] * (double)v[k]))
                / (2.0 * (double)q - 2.0 * (double)v[k]);
        }
        k++;
        v[k] = q;
        z[k] = s;
        z[k + 1] = DT_INF;
    }
    k = 0;
    for (size_t q = 0; q < n; q++) {
        while (z[k + 1] < (double)q) k++;
        double dx = (double)q - (double)v[k];
        d[q] = dx * dx + f[v[k]];
    }
}

/* DistanceTransform[image] / [image, t] -- for each pixel, the Euclidean distance to the nearest
 * BACKGROUND pixel, matching Mathematica and scipy: background pixels are 0, and the value rises
 * toward the interior of a blob. */
static Expr* builtin_distancetransform(Expr* res) {
    size_t argc = res->data.function.arg_count;
    if (argc != 1 && argc != 2) return NULL;
    double thr = 0.0, im = 0.0;
    if (argc == 2) {
        if (!na_read_scalar(res->data.function.args[1], &thr, &im) || im != 0.0) return NULL;
    }
    size_t w = 0, h = 0; double* g = NULL;
    if (!img_grey_plane(res->data.function.args[0], &w, &h, &g)) return NULL;
    size_t n = w * h;

    double* f = malloc(sizeof(double) * n);
    double* tmp = malloc(sizeof(double) * (w > h ? w : h));
    double* line = malloc(sizeof(double) * (w > h ? w : h));
    size_t* vv = malloc(sizeof(size_t) * (w > h ? w : h));
    double* zz = malloc(sizeof(double) * ((w > h ? w : h) + 1));
    Expr* out = NULL;
    if (f && tmp && line && vv && zz) {
        /* Seed: 0 at background, infinity at foreground. The transform then finds, for every pixel,
         * the nearest seed -- which is exactly "distance to the nearest background pixel". */
        for (size_t i = 0; i < n; i++) f[i] = (g[i] > thr) ? DT_INF : 0.0;

        for (size_t y = 0; y < h; y++) {
            for (size_t x = 0; x < w; x++) line[x] = f[y * w + x];
            dt_1d(line, tmp, w, vv, zz);
            for (size_t x = 0; x < w; x++) f[y * w + x] = tmp[x];
        }
        for (size_t x = 0; x < w; x++) {
            for (size_t y = 0; y < h; y++) line[y] = f[y * w + x];
            dt_1d(line, tmp, h, vv, zz);
            for (size_t y = 0; y < h; y++) f[y * w + x] = tmp[y];
        }
        /* The passes accumulate SQUARED distance, which is what makes the decomposition exact; the
         * square root is taken once at the end. Taking it per axis would be wrong, not merely
         * slower. */
        for (size_t i = 0; i < n; i++) f[i] = (f[i] >= DT_INF) ? DT_INF : sqrt(f[i]);
        out = image_build_real(f, w, h, 1);
    }
    free(g); free(f); free(tmp); free(line); free(vv); free(zz);
    return out;
}

/* ---- levels and tone adjustment -------------------------------------------
 *
 * ImageLevels returns DATA, not a plot. Mathematica's ImageHistogram draws a graphic and ImageLevels
 * gives the counts; in a computer algebra system the counts are the useful half, so that is what is
 * implemented and the plot is left to Histogram over the result.
 *
 * THE COUNTS SUM TO THE PIXEL COUNT. That is the property worth asserting, and it is exact: every
 * pixel lands in exactly one bin, so a total that disagrees means a bin boundary is wrong or a pixel
 * was dropped. Bit and Byte images use their natural levels -- 2 and 256 -- because those ARE the
 * distinct values; a Real image has no natural set, so it is binned into 256 over [0, 1] and the bin
 * count is stated rather than guessed at.
 *
 * IMAGEADJUST'S DEFAULT IS A FULL-RANGE STRETCH, which is unambiguous: subtract the minimum, divide by
 * the range, so the darkest pixel becomes exactly 0 and the brightest exactly 1. Two exact properties
 * follow and both are asserted -- the endpoints land on 0 and 1, and the operation is IDEMPOTENT,
 * since a second stretch of an already-stretched image is the identity. Idempotence is the one that
 * would catch an off-by-one in the range, because a slightly wrong divisor still looks like a
 * plausible contrast curve.
 *
 * A CONSTANT IMAGE HAS NO RANGE TO STRETCH, and dividing by zero is not the answer. It comes back
 * unchanged, which is the only defined choice: there is no contrast to expand, and mapping the single
 * value to 0 or to 1 would both be arbitrary.
 *
 * THE THREE-PARAMETER FORM IS OUR OWN DOCUMENTED CURVE, not a claim of bit-compatibility with
 * Mathematica, whose exact formula is not published in a form worth guessing at. Contrast pivots about
 * mid-grey, brightness is an offset, gamma is applied last on the unit interval; each step is stated
 * in the docstring so a caller can reproduce it rather than having to infer it.
 */
#define IMG_LEVEL_BINS 256

static Expr* builtin_imagelevels(Expr* res) {
    size_t argc = res->data.function.arg_count;
    if (argc != 1 && argc != 2) return NULL;

    size_t nbins = 0;
    if (argc == 2) {
        double b = 0.0, im = 0.0;
        if (!na_read_scalar(res->data.function.args[1], &b, &im) || im != 0.0) return NULL;
        if (!(b >= 2.0) || b != floor(b) || b > 65536.0) return NULL;
        nbins = (size_t)b;
    }

    Expr* img = res->data.function.args[0];
    ImgType t;
    bool is3 = image3d_info(img, NULL, NULL, NULL, NULL, &t);
    if (!is3 && !image_info(img, NULL, NULL, NULL, &t)) return NULL;

    size_t w = 0, h = 0, d = 1, c = 0; double* buf = NULL;
    if (is3) { if (!image3d_load(img, &w, &h, &d, &c, &buf)) return NULL; }
    else     { if (!image_load(img, &w, &h, &c, &buf)) return NULL; }
    size_t n = w * h * d * c;

    /* Natural levels where the type HAS them. A Bit image genuinely has two values and a Byte 256, so
     * binning them into anything else would invent structure. Only Real needs a choice. */
    if (nbins == 0) nbins = (t == IMG_BIT) ? 2 : IMG_LEVEL_BINS;

    double* count = calloc(nbins, sizeof(double));
    Expr* out = NULL;
    if (count) {
        for (size_t i = 0; i < n; i++) {
            double v = buf[i];
            if (v < 0.0) v = 0.0;
            if (v > 1.0) v = 1.0;
            size_t b = (size_t)(v * (double)(nbins - 1) + 0.5);
            if (b >= nbins) b = nbins - 1;
            count[b] += 1.0;
        }
        Expr** rows = malloc(sizeof(Expr*) * nbins);
        if (rows) {
            bool ok = true;
            for (size_t i = 0; i < nbins; i++) rows[i] = NULL;
            for (size_t i = 0; i < nbins && ok; i++) {
                Expr* pair[2];
                /* The level, on the same unit scale ImageData uses, so a level can be compared
                 * against a pixel value without rescaling. */
                pair[0] = expr_new_real((double)i / (double)(nbins - 1));
                pair[1] = expr_new_integer((int64_t)count[i]);
                if (pair[0] && pair[1])
                    rows[i] = expr_new_function(expr_new_symbol(SYM_List), pair, 2);
                else { expr_free(pair[0]); expr_free(pair[1]); }
                if (!rows[i]) ok = false;
            }
            if (ok) out = expr_new_function(expr_new_symbol(SYM_List), rows, nbins);
            else for (size_t i = 0; i < nbins; i++) expr_free(rows[i]);
            free(rows);
        }
        free(count);
    }
    free(buf);
    return out;
}

/* ImageAdjust[image] -- full-range stretch.
 * ImageAdjust[image, {c, b}] / [image, {c, b, g}] -- contrast, brightness, gamma. */
static Expr* builtin_imageadjust(Expr* res) {
    size_t argc = res->data.function.arg_count;
    if (argc != 1 && argc != 2) return NULL;

    double ct = 0.0, br = 0.0, gm = 1.0;
    bool params = false;
    if (argc == 2) {
        Expr* sp = res->data.function.args[1];
        if (!ker_is_list(sp)) return NULL;
        size_t k = sp->data.function.arg_count;
        if (k != 2 && k != 3) return NULL;
        double im = 0.0;
        if (!na_read_scalar(sp->data.function.args[0], &ct, &im) || im != 0.0) return NULL;
        if (!na_read_scalar(sp->data.function.args[1], &br, &im) || im != 0.0) return NULL;
        if (k == 3) {
            if (!na_read_scalar(sp->data.function.args[2], &gm, &im) || im != 0.0) return NULL;
            if (!(gm > 0.0)) return NULL;      /* a non-positive gamma is not a curve */
        }
        params = true;
    }

    Expr* img = res->data.function.args[0];
    bool is3 = image3d_info(img, NULL, NULL, NULL, NULL, NULL);
    size_t w = 0, h = 0, d = 1, c = 0; double* buf = NULL;
    if (is3) { if (!image3d_load(img, &w, &h, &d, &c, &buf)) return NULL; }
    else     { if (!image_load(img, &w, &h, &c, &buf)) return NULL; }
    size_t n = w * h * d * c;

    if (!params) {
        double lo = buf[0], hi = buf[0];
        for (size_t i = 1; i < n; i++) { if (buf[i] < lo) lo = buf[i]; if (buf[i] > hi) hi = buf[i]; }
        double span = hi - lo;
        /* No range: return unchanged rather than dividing by zero or picking an arbitrary end. */
        if (span > 0.0)
            for (size_t i = 0; i < n; i++) buf[i] = (buf[i] - lo) / span;
    } else {
        for (size_t i = 0; i < n; i++) {
            /* Contrast pivots about mid-grey so a contrast change does not also shift brightness;
             * brightness is then a plain offset; gamma applies last on the unit interval, which is
             * where it is defined. Clipping happens before gamma because a negative base has no real
             * power. */
            double v = (buf[i] - 0.5) * (1.0 + ct) + 0.5 + br;
            if (v < 0.0) v = 0.0;
            if (v > 1.0) v = 1.0;
            if (gm != 1.0) v = pow(v, 1.0 / gm);
            buf[i] = v;
        }
    }
    Expr* out = is3 ? image3d_build_real(buf, w, h, d, c)
                    : image_build_real(buf, w, h, c);
    free(buf);
    return out;
}

/* ---- correlation and template matching ------------------------------------
 *
 * CORRELATION IS CONVOLUTION WITH THE KERNEL REVERSED ON BOTH AXES, and implementing it that way rather
 * than as a second loop is the point: the two differ by exactly a reflection, so one of them should be
 * derived from the other or they will drift. A test asserts the identity
 * `ImageCorrelate[img, k] == ImageConvolve[img, Reverse[Reverse[k], 2]]`, which is the cleanest
 * statement of the relationship and would fail if either grew its own centring convention.
 *
 * The pair is also where the reflection becomes VISIBLE. On a symmetric kernel they agree exactly, so
 * every smoothing example looks identical; on a delta image with {{1, 2, 3}} convolution gives
 * {1, 2, 3} and correlation gives {3, 2, 1}. That is the same discriminator the convolution tests use,
 * now asserted from the other side.
 *
 * NORMALISED CROSS-CORRELATION is the template-matching method, and it earns its cost. Plain
 * correlation is maximised by BRIGHTNESS, not by similarity: a white patch beats a correct but darker
 * match, which makes raw correlation almost useless for finding a template. NCC subtracts the local
 * mean and divides by the local standard deviation, so it measures shape alone and is invariant to
 * both brightness offset and contrast scale.
 *
 * That buys an exact test. Where the template IS a crop of the image, the two windows are identical up
 * to nothing at all, so NCC there is exactly 1 -- and it is the global maximum. "The peak sits at the
 * crop's location" is a ground truth, not an accuracy figure, which is the same kind of assertion as
 * a k=1 classifier reproducing its training labels.
 *
 * A window of zero variance -- a flat patch -- has no shape to compare, so NCC is undefined there and
 * reports 0 rather than dividing by zero. Reporting 1 would be worse than wrong: a flat region would
 * then match every template perfectly.
 */
/* Correlation by REVERSING the kernel and convolving.
 *
 * The comment above says correlation and convolution differ by exactly a reflection, and this is what
 * acting on that costs: four lines instead of a second nested loop. Two things follow, and the second is
 * why the first version was wrong to duplicate the loop.
 *
 * It cannot drift from convolution, because it IS convolution -- the documented identity
 * `ImageCorrelate[img, k] == ImageConvolve[img, Reverse[Reverse[k], 2]]` holds by construction rather
 * than by two implementations happening to agree.
 *
 * And it inherits SEPARABILITY for free. A 5x5 box is rank 1, so the direct loop was paying 25
 * multiply-adds per pixel where the dispatcher pays 10 -- measured at 3.67 ms against scipy's 2.8, which
 * is what sent this back for a second look. Duplicating the loop did not merely repeat code; it silently
 * opted out of every optimisation the convolution path had accumulated. */
static void correlate_planes(const double* src, double* dst, size_t w, size_t h, size_t c,
                             const double* k, size_t kw, size_t kh) {
    double* rk = malloc(sizeof(double) * kw * kh);
    if (!rk) return;
    for (size_t i = 0; i < kh; i++)
        for (size_t j = 0; j < kw; j++)
            rk[i * kw + j] = k[(kh - 1 - i) * kw + (kw - 1 - j)];
    convolve_dispatch(src, dst, w, h, c, rk, kw, kh);
    free(rk);
}

/* Normalised cross-correlation of a greyscale plane against a template. */
/* Normalised cross-correlation by SUMMED-AREA TABLES.
 *
 * NCC at every position needs, per window, the cross term, the window sum and the window sum of
 * squares. Written directly that is three sweeps of the template over every pixel; the first version
 * here did two (one for the mean, one for the numerator and the deviation norm together) and cost
 * 21.9 ms against SciPy's 8.8 on 512x512 with a 32x32 template.
 *
 * Two identities remove the statistics from the inner loop entirely. With m = kw*kh,
 *
 *     sum (I - Ibar)(T - Tbar) = sum I*T - (sum I) * Tbar
 *     sum (I - Ibar)^2         = sum I^2 - (sum I)^2 / m
 *
 * so the only per-window quantities left are `sum I` and `sum I^2`, and a summed-area table answers
 * each in four lookups regardless of template size. The cross term is then a PLAIN correlation, which
 * means it can go through correlate_planes and inherit the separable path -- the same code that made
 * plain ImageCorrelate 1.94 ms. What was O(pixels * template) three times over becomes
 * O(pixels * template) once plus O(pixels).
 *
 * The tables are built through clampi rather than over a padded copy, so they encode exactly the
 * edge-replicated border the direct loop had, at no extra memory.
 *
 * ONE PROPERTY IS TRADED, and it is worth naming. The direct form computed sum(ds*dt) and sum(ds*ds)
 * from the same values in the same order, so a template matched against itself gave a peak of exactly
 * 1.0. Here the numerator and the variance come from different summations, so the peak is 1.0 within
 * a few ulp instead. The ARGMAX is unaffected -- it is an integer -- and that is the property template
 * matching actually rests on, so the tests assert the argmax exactly and the peak to 1e-12.
 */
static bool ncc_planes(const double* src, double* dst, size_t w, size_t h,
                       const double* t, size_t kw, size_t kh) {
    size_t m = kw * kh;
    int64_t ci = (int64_t)(kh / 2), cj = (int64_t)(kw / 2);

    /* The template's mean and deviation norm are position-independent: one pass, not w*h passes. */
    double tbar = 0.0;
    for (size_t i = 0; i < m; i++) tbar += t[i];
    tbar /= (double)m;
    double tnorm = 0.0;
    for (size_t i = 0; i < m; i++) { double dt = t[i] - tbar; tnorm += dt * dt; }
    tnorm = sqrt(tnorm);

    /* Tables are (PH+1) x (PW+1) over the window-aligned extent, with a zero first row and column so
     * the four-corner difference needs no boundary cases. */
    size_t PW = w + kw - 1, PH = h + kh - 1;
    size_t sw = PW + 1;
    double* s1 = calloc((PH + 1) * sw, sizeof(double));
    double* s2 = calloc((PH + 1) * sw, sizeof(double));
    double* cross = malloc(sizeof(double) * w * h);
    if (!s1 || !s2 || !cross) { free(s1); free(s2); free(cross); return false; }

    for (size_t y = 0; y < PH; y++) {
        size_t sy = clampi((int64_t)y - ci, h);
        for (size_t x = 0; x < PW; x++) {
            size_t sx = clampi((int64_t)x - cj, w);
            double v = src[sy * w + sx];
            s1[(y + 1) * sw + (x + 1)] = s1[y * sw + (x + 1)] + s1[(y + 1) * sw + x]
                                       - s1[y * sw + x] + v;
            s2[(y + 1) * sw + (x + 1)] = s2[y * sw + (x + 1)] + s2[(y + 1) * sw + x]
                                       - s2[y * sw + x] + v * v;
        }
    }

    /* The cross term, through the shared correlation so a separable template costs kw+kh rather
     * than kw*kh. Its border rule is the same clamp the tables encode. */
    correlate_planes(src, cross, w, h, 1, t, kw, kh);

    for (size_t y = 0; y < h; y++)
      for (size_t x = 0; x < w; x++) {
        size_t a = y * sw + x, b = y * sw + (x + kw);
        size_t cc = (y + kh) * sw + x, d = (y + kh) * sw + (x + kw);
        double sum1 = s1[d] - s1[b] - s1[cc] + s1[a];
        double sum2 = s2[d] - s2[b] - s2[cc] + s2[a];
        /* Cancellation can push a uniform window's variance a hair below zero; it is not negative,
         * it is zero, and sqrt of it would be NaN spreading through the result. */
        double var = sum2 - sum1 * sum1 / (double)m;
        double snorm = var > 0.0 ? sqrt(var) : 0.0;
        double num = cross[y * w + x] - sum1 * tbar;
        /* No variance on either side means no shape to compare: 0, not a division. */
        dst[y * w + x] = (snorm > 0.0 && tnorm > 0.0) ? num / (snorm * tnorm) : 0.0;
      }

    free(s1); free(s2); free(cross);
    return true;
}

/* ImageCorrelate[image, kernel] / [image, template, "NormalizedCrossCorrelation"] */
static Expr* builtin_imagecorrelate(Expr* res) {
    size_t argc = res->data.function.arg_count;
    if (argc != 2 && argc != 3) return NULL;

    bool ncc = false;
    if (argc == 3) {
        Expr* m = res->data.function.args[2];
        if (!m || m->type != EXPR_STRING) return NULL;
        if (strcmp(m->data.string, "NormalizedCrossCorrelation") == 0) ncc = true;
        else return NULL;
    }

    size_t kw = 0, kh = 0; double* k = NULL;
    if (!ker_load(res->data.function.args[1], &kh, &kw, &k)) return NULL;

    Expr* out = NULL;
    if (ncc) {
        /* NCC compares SHAPE, which is a property of brightness rather than of colour, so a colour
         * image is reduced to luminance first -- the same choice GradientFilter makes, and for the
         * same reason: combining per-channel scores needs an arbitrary rule. */
        size_t w = 0, h = 0; double* g = NULL;
        if (!img_grey_plane(res->data.function.args[0], &w, &h, &g)) { free(k); return NULL; }
        double* dst = malloc(sizeof(double) * w * h);
        if (dst && ncc_planes(g, dst, w, h, k, kw, kh))
            out = image_build_real(dst, w, h, 1);
        free(g); free(dst);
    } else {
        size_t w = 0, h = 0, c = 0; double* src = NULL;
        if (!image_load(res->data.function.args[0], &w, &h, &c, &src)) { free(k); return NULL; }
        double* dst = malloc(sizeof(double) * w * h * c);
        if (dst) {
            correlate_planes(src, dst, w, h, c, k, kw, kh);
            out = image_build_real(dst, w, h, c);
        }
        free(src); free(dst);
    }
    free(k);
    return out;
}

void imagefilter_init(void) {
    symtab_add_builtin("ImageCorrelate", builtin_imagecorrelate);
    symtab_get_def("ImageCorrelate")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("ImageCorrelate",
        "ImageCorrelate[image, kernel] correlates image with kernel: the kernel is NOT reflected, which "
        "is the only difference from ImageConvolve. The two are related exactly -- correlation equals "
        "convolution with the kernel reversed on both axes -- and they agree on any symmetric kernel, "
        "so the distinction only shows on an asymmetric one, where a delta with {{1,2,3}} gives "
        "{3,2,1} here and {1,2,3} convolved. "
        "ImageCorrelate[image, template, \"NormalizedCrossCorrelation\"] is template matching: it "
        "subtracts the local mean and divides by the local standard deviation, so it measures SHAPE and "
        "is invariant to brightness offset and contrast scale. Plain correlation is maximised by "
        "brightness rather than similarity -- a white patch beats a correct but darker match -- which is "
        "why raw correlation is a poor matcher. Where the template is a crop of the image the score is "
        "exactly 1 and is the global maximum. A flat window has no shape to compare and scores 0 rather "
        "than dividing by zero; scoring 1 would make every flat region match everything. Colour is "
        "reduced to luminance first for the NCC form.");
    symtab_add_builtin("ImageLevels", builtin_imagelevels);
    symtab_get_def("ImageLevels")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("ImageLevels",
        "ImageLevels[image] gives {{level, count}, ...}: the histogram as DATA, not a plot -- use "
        "Histogram over the result for a picture. ImageLevels[image, n] uses n bins. Levels are on "
        "the same unit scale as ImageData, so a level can be compared against a pixel value without "
        "rescaling. A \"Bit\" image uses its 2 natural levels and \"Byte\" its 256, because those "
        "ARE the distinct values; a \"Real\" image has no natural set and is binned into 256 over "
        "[0, 1]. The counts sum to the pixel count exactly, every pixel landing in one bin. Accepts "
        "volumes as well as planes.");

    symtab_add_builtin("ImageAdjust", builtin_imageadjust);
    symtab_get_def("ImageAdjust")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("ImageAdjust",
        "ImageAdjust[image] stretches to the full range: the darkest pixel becomes exactly 0 and the "
        "brightest exactly 1. It is IDEMPOTENT, a second stretch being the identity. A constant image "
        "has no range to stretch and comes back unchanged, since dividing by zero is not the answer "
        "and mapping the single value to either end would be arbitrary. "
        "ImageAdjust[image, {c, b}] and [image, {c, b, g}] apply contrast c, brightness b and gamma g "
        "by a curve stated here rather than inferred: v' = (v - 1/2)(1 + c) + 1/2 + b, clipped to "
        "[0, 1], then raised to the power 1/g. Contrast pivots about mid-grey so it does not also "
        "shift brightness; clipping precedes gamma because a negative base has no real power. This "
        "curve is Mathilda's documented choice, not a claim of bit-compatibility with Mathematica. "
        "Accepts volumes as well as planes.");
    symtab_add_builtin("DistanceTransform", builtin_distancetransform);
    symtab_get_def("DistanceTransform")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("DistanceTransform",
        "DistanceTransform[image] replaces each pixel by its EXACT Euclidean distance to the nearest "
        "background pixel; background pixels are 0, so the value rises toward the interior of a blob. "
        "DistanceTransform[image, t] takes pixels above t as foreground (default 0). "
        "Exact rather than the classic two-pass chamfer approximation, which cannot represent "
        "sqrt(2) with integer steps and so gets diagonal distances a few percent wrong -- invisible "
        "on a picture and fatal to a test. Uses Felzenszwalb and Huttenlocher's lower-envelope-of-"
        "parabolas method, O(n) per row with no sorting. Separability is EXACT here because squared "
        "Euclidean distance is a sum over the axes, so minimising it decomposes per axis; the square "
        "root is taken once at the end rather than per pass.");
    symtab_add_builtin("MedianFilter", builtin_medianfilter);
    symtab_get_def("MedianFilter")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("MedianFilter",
        "MedianFilter[image, r] replaces each pixel with the median over a (2r+1) x (2r+1) "
        "neighbourhood. Unlike a Gaussian it removes an isolated outlier EXACTLY rather than "
        "attenuating and smearing it, which is what makes it the filter for salt-and-pepper noise. "
        "It is also the one filter here that is NOT separable: a sum, a maximum and a minimum all "
        "decompose because they ignore grouping, but a median depends on a value's rank within the "
        "whole window, and grouping destroys rank -- the median of row medians of "
        "{{1,2,9},{3,4,5},{6,7,8}} is 4 where the true median is 5. For an even window the lower "
        "middle is taken rather than the average of the two, so the output is always one of the "
        "inputs.");

    symtab_add_builtin("MeanFilter", builtin_meanfilter);
    symtab_get_def("MeanFilter")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("MeanFilter",
        "MeanFilter[image, r] averages over a (2r+1) x (2r+1) neighbourhood. This IS a convolution "
        "with a normalised box, and it is implemented as one rather than as a separate averaging "
        "loop -- two implementations of one identity is how the identity quietly stops holding. "
        "Being a full rectangle the kernel is separable, so it costs kw + kh rather than kw * kh.");
    symtab_add_builtin("MorphologicalComponents", builtin_morphologicalcomponents);
    symtab_get_def("MorphologicalComponents")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("MorphologicalComponents",
        "MorphologicalComponents[image] labels the connected components of the foreground, giving "
        "an INTEGER MATRIX with background 0 and components numbered 1..k in raster order of first "
        "appearance. MorphologicalComponents[image, t] takes pixels above t as foreground "
        "(default 0, so nonzero is foreground). CornerNeighbors -> False uses 4-connectivity "
        "instead of the default 8. Two pixels touching only at a corner are ONE component under 8 "
        "and TWO under 4, which is the property that distinguishes the two rules -- every other "
        "property holds under either. "
        "A matrix rather than an Image, deliberately: Image type inference would call a label array "
        "of 1..12 a \"Byte\" image and ImageData would then divide every label by 255. Labels are "
        "indices, not brightnesses. Contiguous labels in scan order mean Max of the result is the "
        "component count.");
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
