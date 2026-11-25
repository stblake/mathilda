/* imagecolor.c -- ColorReplace, ColorQuantize and HistogramTransform.
 *
 * Three heads that act on an image's COLOURS rather than its geometry, and they share the one thing
 * that makes such operations awkward: a decision made per pixel needs a global view first. Replacing
 * a colour needs a distance rule, quantising needs a palette derived from every pixel, and
 * equalising needs the whole distribution. So each of these makes a pass to gather, then a pass to
 * write — which is why none of them fits the filter machinery in imagefilter.c.
 */

#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "expr.h"
#include "symtab.h"
#include "attr.h"
#include "image.h"
#include "imagecolor.h"
#include "sym_names.h"
#include "numeric.h"
#include "eval.h"

/* ------------------------------------------------------------------ helpers */

/* Rec. 601 luminance, the same weights ColorConvert uses. Kept identical on purpose: two heads
 * disagreeing about what "brightness" means is the sort of difference nobody notices until a
 * pipeline built from both stops round-tripping. */
static double luma(double r, double g, double b) { return 0.299 * r + 0.587 * g + 0.114 * b; }

/* A colour specification as three components.
 *
 * Accepts RGBColor[r, g, b], GrayLevel[v], a bare number, and {r, g, b}. A grey spec expands to
 * three equal components, so a grey image and a colour image can be compared in one code path. */
/* A bare numeric value (1/2, Pi/4) rather than a colour head. */
static bool head_is_numeric_value(const Expr* e)
{
    return e->type == EXPR_FUNCTION && e->data.function.head &&
           e->data.function.head->type == EXPR_SYMBOL &&
           e->data.function.head->data.symbol.name != SYM_List &&
           e->data.function.head->data.symbol.name != SYM_RGBColor &&
           e->data.function.head->data.symbol.name != SYM_GrayLevel &&
           expr_is_numeric_quantity(e);
}

static bool read_colour(const Expr* e, double rgb[3])
{
    if (!e) return false;
    if (e->type != EXPR_FUNCTION || head_is_numeric_value(e)) {
        double v;
        if (!image_real_value(e, &v)) return false;
        rgb[0] = rgb[1] = rgb[2] = v;
        return true;
    }
    if (e->type != EXPR_FUNCTION || !e->data.function.head ||
        e->data.function.head->type != EXPR_SYMBOL) return false;
    {
        const char* head = e->data.function.head->data.symbol.name;
        size_t n = e->data.function.arg_count;
        double v[3];
        size_t i;
        if ((head == SYM_GrayLevel || strcmp(head, "Graylevel") == 0) && n >= 1) {
            double g;
            if (!image_real_value(e->data.function.args[0], &g)) return false;
            rgb[0] = rgb[1] = rgb[2] = g;
            return true;
        }
        if ((head == SYM_RGBColor || head == SYM_List) && n >= 3) {
            for (i = 0; i < 3; i++) {
                if (!image_real_value(e->data.function.args[i], &v[i])) return false;
            }
            rgb[0] = v[0]; rgb[1] = v[1]; rgb[2] = v[2];
            return true;
        }
        /* RGBColor[grey] is not a thing, but GrayLevel-like one-argument colour heads are common
         * enough that a single component is read as grey rather than declined. */
        if (head == SYM_RGBColor && n == 1) {
            const Expr* a = e->data.function.args[0];
            if (a->type == EXPR_REAL) { rgb[0] = rgb[1] = rgb[2] = a->data.real; return true; }
        }
    }
    return false;
}

static bool is_grey_colour(const double c[3])
{
    return fabs(c[0] - c[1]) < 1e-12 && fabs(c[1] - c[2]) < 1e-12;
}

/* Rule or list of rules -> parallel arrays of colours. Returns the count, or 0 on anything that is
 * not a colour rule. */
static size_t read_rules(const Expr* e, double (*from)[3], double (*to)[3], size_t cap)
{
    const Expr* items[64];
    size_t n = 0, i;

    if (!e || e->type != EXPR_FUNCTION || !e->data.function.head ||
        e->data.function.head->type != EXPR_SYMBOL) return 0;
    {
        const char* head = e->data.function.head->data.symbol.name;
        if (head == SYM_Rule) {
            items[0] = e;
            n = 1;
        } else if (head == SYM_List) {
            n = e->data.function.arg_count;
            if (n == 0 || n > 64) return 0;
            for (i = 0; i < n; i++) items[i] = e->data.function.args[i];
        } else {
            return 0;
        }
    }
    if (n > cap) return 0;
    for (i = 0; i < n; i++) {
        const Expr* r = items[i];
        if (!r || r->type != EXPR_FUNCTION || !r->data.function.head ||
            r->data.function.head->type != EXPR_SYMBOL ||
            r->data.function.head->data.symbol.name != SYM_Rule ||
            r->data.function.arg_count != 2) return 0;
        if (!read_colour(r->data.function.args[0], from[i])) return 0;
        if (!read_colour(r->data.function.args[1], to[i])) return 0;
    }
    return n;
}

static bool as_double(const Expr* e, double* out)
{
    return image_real_value(e, out);   /* 1/2, Pi/4 as well as 0.5 (audit SetAlphaChannel#1) */
}

/* --------------------------------------------------------------- ColorReplace */

/* Mathematica's default tolerance. Stated as a constant because it is a JUDGEMENT: at 0 a
 * replacement matches only bit-identical colours, which after any filtering matches nothing at all. */
#define REPLACE_TOLERANCE 0.02

static Expr* builtin_colorreplace(Expr* res)
{
    double from[64][3], to[64][3];
    size_t nrules, i, k, n_px, argc;
    double tol = REPLACE_TOLERANCE;
    double* buf = NULL;
    double* out;
    size_t w = 0, h = 0, c = 0, oc;
    bool need_colour = false;
    Expr* r;

    if (res->type != EXPR_FUNCTION) return NULL;
    argc = res->data.function.arg_count;
    if (argc < 2 || argc > 3) return NULL;

    nrules = read_rules(res->data.function.args[1], from, to, 64);
    if (!nrules) return NULL;
    if (argc == 3) {
        if (!as_double(res->data.function.args[2], &tol) || tol < 0.0) {
            image_message("ColorReplace", "tol", "The tolerance specification %s should be a "
                          "non-negative real number.", res->data.function.args[2]);
            return NULL;
        }
    }
    if (!image_load(res->data.function.args[0], &w, &h, &c, &buf)) return NULL;

    /* A colour replacement on a grey image has to produce a COLOUR image, or the new colour would
     * be silently flattened to its luminance -- the caller asked for red and would get grey. */
    for (i = 0; i < nrules; i++)
        if (!is_grey_colour(to[i])) need_colour = true;
    /* Colour and alpha are counted separately: a grey+alpha image given a colour becomes RGB+alpha,
     * not luminance+alpha (audit ColorReplace#1). */
    size_t in_col = (c >= 3) ? 3 : 1;
    bool alpha = (c == 2 || c == 4);
    size_t out_col = (in_col == 1 && need_colour) ? 3 : in_col;
    oc = out_col + (alpha ? 1 : 0);

    n_px = w * h;
    out = (double*)malloc(n_px * oc * sizeof(double));
    if (!out) { free(buf); return NULL; }

    for (i = 0; i < n_px; i++) {
        const double* px = buf + i * c;
        double src[3];
        size_t hit = (size_t)-1;
        double best = 0.0;
        if (c >= 3) { src[0] = px[0]; src[1] = px[1]; src[2] = px[2]; }
        else        { src[0] = src[1] = src[2] = px[0]; }

        for (k = 0; k < nrules; k++) {
            double d0 = src[0] - from[k][0], d1 = src[1] - from[k][1], d2 = src[2] - from[k][2];
            double dist = sqrt(d0 * d0 + d1 * d1 + d2 * d2);
            /* NEAREST match within tolerance, not the first: with overlapping rules, "first wins"
             * makes the result depend on the order the caller happened to write them in. */
            if (dist <= tol && (hit == (size_t)-1 || dist < best)) { hit = k; best = dist; }
        }
        for (k = 0; k < out_col; k++) {
            double v;
            if (hit != (size_t)-1) v = (out_col == 3) ? to[hit][k] : luma(to[hit][0], to[hit][1], to[hit][2]);
            else                   v = (out_col == 3) ? src[k] : px[0];
            out[i * oc + k] = v;
        }
        /* An alpha channel passes through untouched: transparency is not a colour. */
        if (alpha) out[i * oc + out_col] = px[c - 1];
    }
    free(buf);
    r = image_build_real(out, w, h, oc);
    free(out);
    return r;
}

/* -------------------------------------------------------------- ColorQuantize */

/* A box of pixel indices, plus the colour range it spans. */
typedef struct { size_t start, len; } Box;

static double* qbuf;      /* the pixel colours being sorted, 3 doubles per entry */
static size_t* qidx;      /* index permutation */
static int qchan;         /* channel currently being sorted on */

static int qcmp(const void* a, const void* b)
{
    size_t ia = *(const size_t*)a, ib = *(const size_t*)b;
    double va = qbuf[ia * 3 + qchan], vb = qbuf[ib * 3 + qchan];
    if (va < vb) return -1;
    if (va > vb) return 1;
    /* Ties broken by index, so the palette does not depend on qsort's internal choices -- a
     * quantisation that differs run to run cannot be tested. */
    return ia < ib ? -1 : (ia > ib ? 1 : 0);
}

/* MEDIAN CUT rather than k-means. k-means gives a slightly better palette and needs a seed: the
 * result would then depend on the random stream, so the same image quantised twice could differ, and
 * no test could assert a palette. Median cut is deterministic, which for a documented, testable
 * builtin is worth more than a marginally lower error. */
/* ColorQuantize[image, {c1, c2, ...}]: each pixel replaced by its nearest palette colour
 * (Euclidean in RGB), alpha passed through. A grey image and an all-grey palette stay grey
 * (audit ColorQuantize#2). */
static Expr* colorquantize_palette(Expr* res)
{
    const Expr* pal = res->data.function.args[1];
    size_t np = pal->data.function.arg_count, i, k;
    if (np == 0 || np > 4096) return NULL;
    double (*pc)[3] = malloc(sizeof(double[3]) * np);
    if (!pc) return NULL;
    bool grey_pal = true;
    for (k = 0; k < np; k++) {
        if (!read_colour(pal->data.function.args[k], pc[k])) { free(pc); return NULL; }
        if (!is_grey_colour(pc[k])) grey_pal = false;
    }
    double* buf = NULL;
    size_t w = 0, h = 0, c = 0;
    if (!image_load(res->data.function.args[0], &w, &h, &c, &buf)) { free(pc); return NULL; }
    size_t in_col = c >= 3 ? 3 : 1;
    bool alpha = (c == 2 || c == 4);
    size_t out_col = (in_col == 1 && grey_pal) ? 1 : 3;
    size_t oc = out_col + (alpha ? 1 : 0);
    double* out = malloc(sizeof(double) * w * h * oc);
    Expr* r = NULL;
    if (out) {
        for (i = 0; i < w * h; i++) {
            const double* px = buf + i * c;
            double src[3];
            if (in_col == 3) { src[0] = px[0]; src[1] = px[1]; src[2] = px[2]; }
            else src[0] = src[1] = src[2] = px[0];
            size_t best = 0; double bd = 0.0;
            for (k = 0; k < np; k++) {
                double d0 = src[0] - pc[k][0], d1 = src[1] - pc[k][1], d2 = src[2] - pc[k][2];
                double d = d0 * d0 + d1 * d1 + d2 * d2;
                if (k == 0 || d < bd) { bd = d; best = k; }
            }
            for (k = 0; k < out_col; k++) out[i * oc + k] = pc[best][out_col == 1 ? 0 : k];
            if (alpha) out[i * oc + out_col] = px[c - 1];
        }
        r = image_build_real(out, w, h, oc);
    }
    free(out); free(buf); free(pc);
    return r;
}

static Expr* builtin_colorquantize(Expr* res)
{
    double* buf = NULL;
    size_t w = 0, h = 0, c = 0, n_px, i, j, k;
    long ncol;
    Box* boxes = NULL;
    size_t nboxes = 0;
    double* out = NULL;
    Expr* r = NULL;

    if (res->type != EXPR_FUNCTION) return NULL;
    if (res->data.function.arg_count == 3) {
        /* Dithering -> False is exactly what this computes -- every pixel mapped to its nearest
         * palette colour, no error diffusion -- so the option is accepted and dropped (audit
         * ColorQuantize#4). Dithering is not implemented, so any other value stays unevaluated. */
        const Expr* o = res->data.function.args[2];
        if (o->type == EXPR_FUNCTION && o->data.function.arg_count == 2
            && o->data.function.head->type == EXPR_SYMBOL
            && (o->data.function.head->data.symbol.name == SYM_Rule
                || o->data.function.head->data.symbol.name == SYM_RuleDelayed)
            && o->data.function.args[0]->type == EXPR_SYMBOL
            && strcmp(o->data.function.args[0]->data.symbol.name, "Dithering") == 0
            && o->data.function.args[1]->type == EXPR_SYMBOL
            && o->data.function.args[1]->data.symbol.name == SYM_False)
            return expr_new_function(expr_copy(res->data.function.head),
                (Expr*[]){ expr_copy(res->data.function.args[0]), expr_copy(res->data.function.args[1]) }, 2);
        return NULL;
    }
    if (res->data.function.arg_count != 2) return NULL;
    {
        const Expr* a = res->data.function.args[1];
        if (a && a->type == EXPR_FUNCTION && a->data.function.head->type == EXPR_SYMBOL
            && a->data.function.head->data.symbol.name == SYM_List)
            return colorquantize_palette(res);
        if (!a || a->type != EXPR_INTEGER) return NULL;
        ncol = (long)a->data.integer;
        if (ncol < 1 || ncol > 4096) return NULL;
    }
    if (!image_load(res->data.function.args[0], &w, &h, &c, &buf)) return NULL;
    n_px = w * h;
    if (!n_px) { free(buf); return NULL; }

    qbuf = (double*)malloc(n_px * 3 * sizeof(double));
    qidx = (size_t*)malloc(n_px * sizeof(size_t));
    boxes = (Box*)malloc((size_t)ncol * sizeof(Box));
    out = (double*)malloc(n_px * c * sizeof(double));
    if (!qbuf || !qidx || !boxes || !out) goto done;

    for (i = 0; i < n_px; i++) {
        const double* px = buf + i * c;
        if (c >= 3) { qbuf[i * 3] = px[0]; qbuf[i * 3 + 1] = px[1]; qbuf[i * 3 + 2] = px[2]; }
        else        { qbuf[i * 3] = qbuf[i * 3 + 1] = qbuf[i * 3 + 2] = px[0]; }
        qidx[i] = i;
    }
    boxes[0].start = 0; boxes[0].len = n_px; nboxes = 1;

    while ((long)nboxes < ncol) {
        /* Split the box with the widest single-channel spread. Widest spread rather than most
         * pixels: a large box of nearly identical colours does not need splitting, and a small one
         * spanning half the spectrum does. */
        size_t bi = 0, best_b = (size_t)-1;
        int best_ch = 0;
        double best_range = -1.0;
        for (bi = 0; bi < nboxes; bi++) {
            if (boxes[bi].len < 2) continue;
            for (k = 0; k < 3; k++) {
                double lo = 2.0, hi = -1.0;
                for (j = 0; j < boxes[bi].len; j++) {
                    double v = qbuf[qidx[boxes[bi].start + j] * 3 + k];
                    if (v < lo) lo = v;
                    if (v > hi) hi = v;
                }
                if (hi - lo > best_range) { best_range = hi - lo; best_b = bi; best_ch = (int)k; }
            }
        }
        if (best_b == (size_t)-1 || best_range <= 0.0) break;   /* every box is a single colour */

        qchan = best_ch;
        qsort(qidx + boxes[best_b].start, boxes[best_b].len, sizeof(size_t), qcmp);
        {
            /* Cut at a VALUE boundary: the position nearest the median at which the sorted channel
             * actually changes. Cutting at len/2 itself split a run of identical colours between two
             * boxes, so identical pixels came back as different colours and a box's mean invented a
             * colour the image does not hold (red and blue averaged to purple for an image that
             * already had only two colours). A boundary exists because best_range > 0. */
            size_t st = boxes[best_b].start, len = boxes[best_b].len, mid = len / 2, half = 0, dlt;
            for (dlt = 0; dlt < len && !half; dlt++) {
                size_t cand[2];
                int nc = 0, q;
                if (mid >= dlt) cand[nc++] = mid - dlt;
                if (dlt && mid + dlt < len) cand[nc++] = mid + dlt;
                for (q = 0; q < nc && !half; q++) {
                    size_t sp = cand[q];
                    if (sp >= 1 && sp < len
                        && qbuf[qidx[st + sp - 1] * 3 + qchan] < qbuf[qidx[st + sp] * 3 + qchan])
                        half = sp;
                }
            }
            if (!half) half = mid;      /* unreachable: the range is positive */
            Box lo = { boxes[best_b].start, half };
            Box hi = { boxes[best_b].start + half, boxes[best_b].len - half };
            boxes[best_b] = lo;
            boxes[nboxes++] = hi;
        }
    }

    /* Each box collapses to its mean colour. */
    for (i = 0; i < nboxes; i++) {
        double sum[3] = {0.0, 0.0, 0.0};
        for (j = 0; j < boxes[i].len; j++) {
            const double* q = qbuf + qidx[boxes[i].start + j] * 3;
            sum[0] += q[0]; sum[1] += q[1]; sum[2] += q[2];
        }
        if (!boxes[i].len) continue;
        sum[0] /= (double)boxes[i].len;
        sum[1] /= (double)boxes[i].len;
        sum[2] /= (double)boxes[i].len;
        for (j = 0; j < boxes[i].len; j++) {
            size_t p = qidx[boxes[i].start + j];
            if (c >= 3) {
                out[p * c] = sum[0]; out[p * c + 1] = sum[1]; out[p * c + 2] = sum[2];
                if (c == 4) out[p * c + 3] = buf[p * c + 3];
            } else {
                /* A grey pixel was copied into three identical channels, so their mean IS the
                 * grey value; recombining through the luminance weights (which do not sum to
                 * exactly 1.0 in binary) drifted 0.5 to 0.49999999999999994. */
                out[p * c] = sum[0];
                if (c == 2) out[p * c + 1] = buf[p * c + 1];
            }
        }
    }
    r = image_build_real(out, w, h, c);

done:
    free(buf); free(qbuf); free(qidx); free(boxes); free(out);
    qbuf = NULL; qidx = NULL;
    return r;
}

/* ---------------------------------------------------------- HistogramTransform */

#define HT_BINS 256

/* HistogramTransform[image] -- histogram equalisation.
 *
 * The mapping comes from the LUMINANCE distribution and is applied to every channel by a ratio, so
 * hue survives. Equalising each channel independently is the other obvious choice and it shifts
 * colour: a warm image comes back grey-ish, because equalising red and blue separately removes
 * exactly the imbalance that made it warm.
 */
static int ht_cmp(const void* a, const void* b)
{
    double x = *(const double*)a, y = *(const double*)b;
    return (x > y) - (x < y);
}

/* The target quantile function for HistogramTransform[image, ref] (a reference image's sorted
 * luminances) or HistogramTransform[image, dist] (Quantile[dist, p]), tabulated at the mid-rank of
 * every occupied bin. False when the reference is neither (audit HistogramTransform#1). */
static bool ht_target_map(const Expr* ref, const double* cdf, const size_t* hist, double* map)
{
    size_t rw = 0, rh = 0, rc = 0, i;
    double* rbuf = NULL;
    if (image_info(ref, NULL, NULL, NULL, NULL)) {
        if (!image_load(ref, &rw, &rh, &rc, &rbuf)) return false;
        size_t m = rw * rh;
        if (!m) { free(rbuf); return false; }
        double* lum = malloc(sizeof(double) * m);
        if (!lum) { free(rbuf); return false; }
        for (i = 0; i < m; i++) {
            const double* px = rbuf + i * rc;
            lum[i] = (rc >= 3) ? luma(px[0], px[1], px[2]) : px[0];
        }
        free(rbuf);
        qsort(lum, m, sizeof(double), ht_cmp);
        for (i = 0; i < HT_BINS; i++) {
            if (!hist[i]) continue;
            double pm = (i ? cdf[i - 1] : 0.0) * 0.5 + cdf[i] * 0.5;
            double q = ceil(pm * (double)m) - 1.0;
            size_t idx = q < 0.0 ? 0 : (q >= (double)m ? m - 1 : (size_t)q);
            map[i] = lum[idx];
        }
        free(lum);
        return true;
    }
    if (ref->type != EXPR_FUNCTION || ref->data.function.head->type != EXPR_SYMBOL) return false;
    for (i = 0; i < HT_BINS; i++) {
        if (!hist[i]) continue;
        double pm = (i ? cdf[i - 1] : 0.0) * 0.5 + cdf[i] * 0.5;
        Expr* call = expr_new_function(expr_new_symbol("Quantile"),
                         (Expr*[]){ expr_copy((Expr*)ref), expr_new_real(pm) }, 2);
        Expr* v = eval_and_free(expr_new_function(expr_new_symbol(SYM_N), (Expr*[]){ call }, 1));
        bool ok = image_real_value(v, &map[i]);
        expr_free(v);
        if (!ok) return false;
    }
    return true;
}

static Expr* builtin_histogramtransform(Expr* res)
{
    double* buf = NULL;
    size_t w = 0, h = 0, c = 0, n_px, i, k;
    size_t hist[HT_BINS];
    double cdf[HT_BINS];
    double map[HT_BINS];
    bool matched = false;
    double* out;
    Expr* r;

    if (res->type != EXPR_FUNCTION) return NULL;
    if (res->data.function.arg_count != 1 && res->data.function.arg_count != 2) return NULL;
    /* A volume equalises exactly as a plane does: the histogram of a volume is a histogram, so the
     * voxels are one flat buffer here and only the rebuild knows the rank (audit
     * HistogramTransform#3). */
    size_t d = 1;
    bool vol = image3d_info(res->data.function.args[0], NULL, NULL, NULL, NULL, NULL);
    if (vol ? !image3d_load(res->data.function.args[0], &w, &h, &d, &c, &buf)
            : !image_load(res->data.function.args[0], &w, &h, &c, &buf)) return NULL;
    n_px = w * h * d;
    if (!n_px) { free(buf); return NULL; }

    for (i = 0; i < HT_BINS; i++) hist[i] = 0;
    for (i = 0; i < n_px; i++) {
        const double* px = buf + i * c;
        double l = (c >= 3) ? luma(px[0], px[1], px[2]) : px[0];
        long b = (long)(l * (HT_BINS - 1) + 0.5);
        if (b < 0) b = 0;
        if (b >= HT_BINS) b = HT_BINS - 1;
        hist[b]++;
    }
    {
        size_t acc = 0;
        for (i = 0; i < HT_BINS; i++) {
            acc += hist[i];
            cdf[i] = (double)acc / (double)n_px;
        }
    }

    if (res->data.function.arg_count == 2) {
        if (!ht_target_map(res->data.function.args[1], cdf, hist, map)) { free(buf); return NULL; }
        matched = true;
    }

    out = (double*)malloc(n_px * c * sizeof(double));
    if (!out) { free(buf); return NULL; }
    for (i = 0; i < n_px; i++) {
        const double* px = buf + i * c;
        double l = (c >= 3) ? luma(px[0], px[1], px[2]) : px[0];
        long b = (long)(l * (HT_BINS - 1) + 0.5);
        double nl;
        if (b < 0) b = 0;
        if (b >= HT_BINS) b = HT_BINS - 1;
        nl = matched ? map[b] : cdf[b];
        if (c >= 3) {
            /* Scale the channels by the luminance's own change. A black pixel has no ratio to
             * scale, so it takes the new luminance in every channel -- grey, which is the only
             * hue-free answer available. */
            double f = (l > 1e-9) ? nl / l : 0.0;
            /* One factor for all three channels is what keeps the hue, so where it would push a
             * channel past 1 the FACTOR is capped, not each channel: clipping channels separately
             * turned an orange {1, 0.8, 0.2} yellow (audit HistogramTransform#2). The pixel then
             * reaches less than the target luminance, which is the price of keeping its colour. */
            double mx = px[0] > px[1] ? (px[0] > px[2] ? px[0] : px[2]) : (px[1] > px[2] ? px[1] : px[2]);
            if (mx > 0.0 && mx * f > 1.0) f = 1.0 / mx;
            for (k = 0; k < 3; k++) {
                double v = (l > 1e-9) ? px[k] * f : nl;
                out[i * c + k] = v < 0.0 ? 0.0 : (v > 1.0 ? 1.0 : v);
            }
            if (c == 4) out[i * c + 3] = px[3];
        } else {
            out[i * c] = nl;
            if (c == 2) out[i * c + 1] = px[1];
        }
    }
    free(buf);
    r = vol ? image3d_build_real(out, w, h, d, c) : image_build_real(out, w, h, c);
    free(out);
    return r;
}

void imagecolor_init(void)
{
    symtab_add_builtin("ColorReplace", builtin_colorreplace);
    symtab_get_def("ColorReplace")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("ColorReplace",
        "ColorReplace[image, old -> new] replaces every pixel within a tolerance of `old` by `new`; "
        "ColorReplace[image, {r1, r2, ...}] applies several rules and "
        "ColorReplace[image, rules, tol] sets the tolerance (default 0.02 -- at 0 only bit-identical "
        "colours match, which after any filtering is nothing at all). Distance is Euclidean in RGB, "
        "and where rules overlap the NEAREST wins rather than the first, so the answer does not "
        "depend on the order they were written. Colours may be RGBColor[r, g, b], GrayLevel[v], a "
        "number or {r, g, b}. Replacing a grey image's colour with a non-grey one produces a "
        "three-channel image, since flattening the new colour to its luminance would give grey when "
        "the caller asked for red. An alpha channel passes through: transparency is not a colour.");

    symtab_add_builtin("ColorQuantize", builtin_colorquantize);
    symtab_get_def("ColorQuantize")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("ColorQuantize",
        "ColorQuantize[image, {c1, c2, ...}] maps each pixel to the nearest palette colour. "
        "ColorQuantize[image, n] reduces the image to at most n colours by MEDIAN CUT: the box with "
        "the widest single-channel spread is split at its median until n boxes remain, and each "
        "collapses to its mean colour. Widest spread rather than most pixels, since a large box of "
        "nearly identical colours does not need splitting and a small one spanning half the "
        "spectrum does. Median cut rather than k-means because it is DETERMINISTIC -- a palette that "
        "depended on the random stream could not be tested or documented. The channel count is "
        "preserved and alpha passes through.");

    symtab_add_builtin("HistogramTransform", builtin_histogramtransform);
    symtab_get_def("HistogramTransform")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("HistogramTransform",
        "HistogramTransform[image, ref] matches the histogram of a reference image, and "
        "HistogramTransform[image, dist] that of a distribution. "
        "HistogramTransform[image] equalises the histogram, spreading the brightness distribution "
        "toward uniform over 256 bins by mapping each value through the cumulative distribution. "
        "The mapping is computed from the LUMINANCE and applied to every channel as a ratio, so hue "
        "survives; equalising each channel independently would shift colour, since it removes "
        "exactly the imbalance that makes an image warm or cool. A black pixel has no ratio to "
        "scale and takes the new luminance in every channel. Alpha passes through.");
}
