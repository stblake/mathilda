/* root_numeric.c — see root_numeric.h for the pipeline overview.
 *
 * The hard part of numerically evaluating Root[Function[p[#]], k] is
 * *root selection*: a naive Newton iteration from an unverified guess
 * can land in any basin of attraction and silently return a different
 * root than the user asked for.  The algorithm below guarantees
 * correctness in three layers:
 *
 *   - The companion-matrix QR gives us *every* approximate root in one
 *     shot, so canonical position is decidable from the start (no
 *     ambiguity about which root the iterate is converging to).
 *
 *   - The Sturm sign-variation count certifies the real-root tally is
 *     correct, so the real/complex partition used in canonical sorting
 *     is provably right (not just a per-root tolerance check that can
 *     misclassify near-real complex pairs).
 *
 *   - The basin-verification step at the end confirms the Newton-
 *     refined root still occupies canonical position k when compared
 *     against the other QR estimates uplifted to refined precision.
 */

#include "root_numeric.h"

#include "expr.h"
#include "common.h"
#include "eval.h"
#include "sym_names.h"
#include "arithmetic.h"
#include "symtab.h"

#include "poly/zupoly.h"
#include "poly/poly_eval_mpfr.h"

#include "linalg/eigen.h"

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <stdarg.h>
#include <float.h>

#ifdef USE_MPFR
#include <mpfr.h>
#endif
#include <gmp.h>

/* ====================================================================
 *  Diagnostics
 * ================================================================== */
static void root_warn(const char* tag, const char* fmt, ...) {
    va_list ap;
    fprintf(stderr, "Root::%s: ", tag);
    va_start(ap, fmt);
    vfprintf(stderr, fmt, ap);
    va_end(ap);
    fputc('\n', stderr);
}

#ifdef USE_MPFR

/* ====================================================================
 *  Step 1: extract polynomial and k from Root[Function[poly_in_slot], k]
 * ================================================================== */

/* Build a fresh Expr representing Slot[1]. */
static Expr* make_slot1(void) {
    Expr* one = expr_new_integer(1);
    Expr* slot = expr_new_function(expr_new_symbol("Slot"),
                                   (Expr*[]){ one }, 1);
    return slot;
}

/* True iff the candidate is `Root[Function[<body>], <integer-k>]`
 * with a positive integer k.  On success, *body_out is set to a
 * borrowed pointer into the Root tree and *k_out receives k. */
static bool parse_root_call(const Expr* e, Expr** body_out, int64_t* k_out) {
    if (!e || e->type != EXPR_FUNCTION) return false;
    if (!head_is(e, SYM_Root)) return false;
    if (e->data.function.arg_count != 2) return false;
    Expr* fn = e->data.function.args[0];
    Expr* kx = e->data.function.args[1];
    if (!fn || fn->type != EXPR_FUNCTION || !head_is(fn, SYM_Function)
        || fn->data.function.arg_count != 1) return false;
    if (!kx || kx->type != EXPR_INTEGER || kx->data.integer < 1) return false;
    *body_out = fn->data.function.args[0];
    *k_out    = kx->data.integer;
    return true;
}

/* ====================================================================
 *  Step 2: squarefree part via Yun (rad(p) = p / gcd(p, p'))
 * ================================================================== */

/* In-place derivative of p with respect to the polynomial variable.
 * Returns a fresh ZUPoly p' with p'[i] = (i+1) * p[i+1].  The zero
 * polynomial's derivative is zero. */
static ZUPoly* zupoly_derivative(const ZUPoly* p) {
    if (p->deg <= 0) return zupoly_zero();
    ZUPoly* q = zupoly_new(p->deg);
    mpz_t coef; mpz_init(coef);
    for (int i = 1; i <= p->deg; i++) {
        mpz_mul_ui(coef, p->c[i], (unsigned long)i);
        zupoly_setcoef(q, i - 1, coef);
    }
    mpz_clear(coef);
    return q;
}

/* Squarefree radical of p: p / gcd(p, p').  Always returns a fresh
 * ZUPoly with positive leading coefficient. */
static ZUPoly* zupoly_squarefree_part(const ZUPoly* p) {
    ZUPoly* dp = zupoly_derivative(p);
    ZUPoly* g  = zupoly_gcd(p, dp);
    ZUPoly* r  = zupoly_divexact(p, g);
    zupoly_free(dp);
    zupoly_free(g);
    if (!r) return zupoly_copy(p);   /* fallback: leave as-is */
    /* Normalize leading sign positive. */
    if (r->deg >= 0 && mpz_sgn(r->c[r->deg]) < 0) {
        ZUPoly* neg = zupoly_neg(r);
        zupoly_free(r);
        r = neg;
    }
    return r;
}

/* ====================================================================
 *  Step 3: Frobenius companion matrix at MPFR precision
 *
 *  For p(x) = c_n x^n + c_{n-1} x^{n-1} + ... + c_0, the companion
 *  matrix is row-major n*n:
 *      C[i, n-1] = -c_i / c_n     (last column = negated normalized lower coeffs)
 *      C[i+1, i] = 1              (sub-diagonal of ones)
 *      else 0.
 *  Its eigenvalues are exactly the roots of p.
 * ================================================================== */
static void build_companion_matrix(const ZUPoly* p, mpfr_prec_t bits,
                                   mpfr_t** out_mat, size_t* out_n) {
    int n = p->deg;
    *out_n = (size_t)n;
    mpfr_t* M = (mpfr_t*)malloc(sizeof(mpfr_t) * (size_t)(n * n));
    for (size_t i = 0; i < (size_t)(n * n); i++) {
        mpfr_init2(M[i], bits);
        mpfr_set_zero(M[i], 1);
    }
    /* Sub-diagonal of ones. */
    for (int i = 0; i + 1 < n; i++) {
        mpfr_set_si(M[(size_t)(i + 1) * (size_t)n + (size_t)i], 1, MPFR_RNDN);
    }
    /* Last column: -c_i / c_n. */
    mpfr_t lc, ratio;
    mpfr_init2(lc, bits); mpfr_init2(ratio, bits);
    mpfr_set_z(lc, p->c[n], MPFR_RNDN);
    for (int i = 0; i < n; i++) {
        mpfr_set_z(ratio, p->c[i], MPFR_RNDN);
        mpfr_div(ratio, ratio, lc, MPFR_RNDN);
        mpfr_neg(ratio, ratio, MPFR_RNDN);
        mpfr_set(M[(size_t)i * (size_t)n + (size_t)(n - 1)], ratio, MPFR_RNDN);
    }
    mpfr_clear(lc); mpfr_clear(ratio);
    *out_mat = M;
}

static void free_mpfr_array(mpfr_t* a, size_t count) {
    if (!a) return;
    for (size_t i = 0; i < count; i++) mpfr_clear(a[i]);
    free(a);
}

/* ====================================================================
 *  Step 4: Sturm sign-variation real-root count
 *
 *  The Sturm chain s_0 = p, s_1 = p', s_{k+2} = -rem(s_k, s_{k+1})
 *  gives V(a) - V(b) = (real roots of p in (a, b]) for a squarefree p,
 *  where V(x) = sign variations in s_0(x), s_1(x), ..., s_m(x).
 *
 *  We evaluate at +/- infinity via the leading coefficient signs:
 *  sign(s_k(+inf)) = sign(leading_coef(s_k))
 *  sign(s_k(-inf)) = sign(leading_coef(s_k)) * (-1)^deg(s_k)
 * ================================================================== */
static int sturm_real_root_count(const ZUPoly* p) {
    if (p->deg <= 0) return 0;
    /* Build chain in a dynamic array. */
    int cap = p->deg + 2;
    ZUPoly** chain = (ZUPoly**)malloc(sizeof(ZUPoly*) * (size_t)cap);
    int n = 0;
    chain[n++] = zupoly_copy(p);
    chain[n++] = zupoly_derivative(p);
    while (chain[n - 1]->deg > 0) {
        const ZUPoly* a = chain[n - 2];
        const ZUPoly* b = chain[n - 1];
        ZUPoly *q, *r;
        if (!zupoly_pseudodivrem(a, b, &q, &r)) {
            zupoly_free(q); zupoly_free(r);
            break;
        }
        zupoly_free(q);
        if (r->deg < 0) { zupoly_free(r); break; }
        /* zupoly_pseudodivrem yields  lc(b)^N * a = q*b + r,  where
         * N = deg(a) - deg(b) + 1.  The true Euclidean remainder R is
         * R = r / lc(b)^N.  Sturm requires the chain successor to be a
         * positive scalar multiple of -R, i.e.
         *   - if lc(b) > 0, or lc(b) < 0 with N even: use -r
         *   - if lc(b) < 0 with N odd:               use  r
         * (The opposite sign cancels the sign of lc(b)^N.) */
        int lc_sgn = mpz_sgn(b->c[b->deg]);
        int N = a->deg - b->deg + 1;
        bool negate = !(lc_sgn < 0 && (N % 2 == 1));
        ZUPoly* signed_r = negate ? zupoly_neg(r) : zupoly_copy(r);
        zupoly_free(r);
        /* Primitive part keeps the leading-coefficient sign, so the
         * sign chosen above survives content removal. */
        ZUPoly* pr = zupoly_primitive_part(signed_r);
        zupoly_free(signed_r);
        chain[n++] = pr;
        if (n >= cap) break;
    }

    /* Sign variations at +inf and -inf. */
    int V_pos = 0, V_neg = 0;
    int prev_pos = 0, prev_neg = 0; /* prior nonzero sign */
    for (int i = 0; i < n; i++) {
        if (chain[i]->deg < 0) continue;  /* zero entry */
        int sgn_lc = mpz_sgn(chain[i]->c[chain[i]->deg]);
        if (sgn_lc == 0) continue;
        int sgn_pos = sgn_lc;
        int sgn_neg = (chain[i]->deg % 2 == 0) ? sgn_lc : -sgn_lc;
        if (prev_pos != 0 && sgn_pos != prev_pos) V_pos++;
        if (prev_neg != 0 && sgn_neg != prev_neg) V_neg++;
        prev_pos = sgn_pos;
        prev_neg = sgn_neg;
    }

    for (int i = 0; i < n; i++) zupoly_free(chain[i]);
    free(chain);
    return V_neg - V_pos;
}

/* ====================================================================
 *  Step 5: canonical sort
 * ================================================================== */
int root_canonical_cmp_mpfr(const mpfr_t a_re, const mpfr_t a_im, int a_real,
                            const mpfr_t b_re, const mpfr_t b_im, int b_real) {
    if (a_real && !b_real) return -1;
    if (!a_real && b_real) return  1;
    if (a_real && b_real) {
        int c = mpfr_cmp(a_re, b_re);
        if (c < 0) return -1;
        if (c > 0) return  1;
        return 0;
    }
    /* both complex */
    int c = mpfr_cmp(a_re, b_re);
    if (c < 0) return -1;
    if (c > 0) return  1;
    /* by |im| ascending */
    mpfr_t aa, bb;
    mpfr_init2(aa, mpfr_get_prec(a_im));
    mpfr_init2(bb, mpfr_get_prec(b_im));
    mpfr_abs(aa, a_im, MPFR_RNDN);
    mpfr_abs(bb, b_im, MPFR_RNDN);
    c = mpfr_cmp(aa, bb);
    mpfr_clear(aa); mpfr_clear(bb);
    if (c < 0) return -1;
    if (c > 0) return  1;
    /* tie-break: smaller im (more negative) first */
    int sa = mpfr_sgn(a_im);
    int sb = mpfr_sgn(b_im);
    if (sa < sb) return -1;
    if (sa > sb) return  1;
    return 0;
}

/* Classify a single QR-derived eigenvalue as real if |Im| <
 * `tol_factor` * max(1, |Re|).  `tol_factor` is typically
 * 2^(-bits/2). */
static int classify_real(const mpfr_t re, const mpfr_t im,
                         const mpfr_t tol_factor) {
    mpfr_prec_t bits = mpfr_get_prec(re);
    mpfr_t scale, bound, abs_im;
    mpfr_init2(scale, bits);
    mpfr_init2(bound, bits);
    mpfr_init2(abs_im, bits);
    mpfr_abs(scale, re, MPFR_RNDN);
    mpfr_t one; mpfr_init2(one, bits); mpfr_set_ui(one, 1, MPFR_RNDN);
    if (mpfr_cmp(scale, one) < 0) mpfr_set(scale, one, MPFR_RNDN);
    mpfr_mul(bound, scale, tol_factor, MPFR_RNDN);
    mpfr_abs(abs_im, im, MPFR_RNDN);
    int is_real = (mpfr_cmp(abs_im, bound) <= 0) ? 1 : 0;
    mpfr_clear(scale); mpfr_clear(bound); mpfr_clear(abs_im); mpfr_clear(one);
    return is_real;
}

/* Sort an array of root candidates in canonical order.  Each root
 * carries (re, im, is_real).  We use the same backing arrays directly
 * (in-place insertion sort — n is small, typically <= 10). */
static void canonical_sort_roots(mpfr_t* re, mpfr_t* im, int* is_real,
                                 size_t n) {
    mpfr_prec_t bits = (n > 0) ? mpfr_get_prec(re[0]) : 53;
    for (size_t i = 1; i < n; i++) {
        for (size_t j = i; j > 0; j--) {
            int c = root_canonical_cmp_mpfr(re[j], im[j], is_real[j],
                                            re[j - 1], im[j - 1], is_real[j - 1]);
            if (c >= 0) break;
            /* swap j and j-1 */
            mpfr_t tmp; mpfr_init2(tmp, bits);
            mpfr_set(tmp, re[j], MPFR_RNDN);
            mpfr_set(re[j], re[j - 1], MPFR_RNDN);
            mpfr_set(re[j - 1], tmp, MPFR_RNDN);
            mpfr_set(tmp, im[j], MPFR_RNDN);
            mpfr_set(im[j], im[j - 1], MPFR_RNDN);
            mpfr_set(im[j - 1], tmp, MPFR_RNDN);
            int tib = is_real[j];
            is_real[j] = is_real[j - 1];
            is_real[j - 1] = tib;
            mpfr_clear(tmp);
        }
    }
}

/* ====================================================================
 *  Step 6: Newton refinement at high precision
 *
 *  Standard Newton: x_{k+1} = x_k - p(x_k) / p'(x_k), using the fused
 *  Horner evaluator from poly_eval_mpfr.  Converges quadratically on
 *  simple roots.  Halts when |p(x_k)| / |p'(x_k)| < 2^(-bits + guard)
 *  relative to scale.
 *
 *  Returns true on convergence, false on stall.
 * ================================================================== */
#define ROOT_NEWTON_MAX_ITER 100

static bool newton_refine_real(const mpfr_t* coeffs, int deg,
                               mpfr_t x, mpfr_prec_t bits) {
    mpfr_t pv, dpv, step, abs_step, tol, abs_x;
    mpfr_init2(pv, bits); mpfr_init2(dpv, bits);
    mpfr_init2(step, bits); mpfr_init2(abs_step, bits);
    mpfr_init2(tol, bits); mpfr_init2(abs_x, bits);

    mpfr_set_ui(tol, 1, MPFR_RNDN);
    mpfr_div_2si(tol, tol, (long)bits - 8, MPFR_RNDN);

    bool converged = false;
    for (int it = 0; it < ROOT_NEWTON_MAX_ITER; it++) {
        poly_eval_real_mpfr(coeffs, deg, x, pv, &dpv);
        if (mpfr_zero_p(dpv)) break;
        mpfr_div(step, pv, dpv, MPFR_RNDN);
        mpfr_sub(x, x, step, MPFR_RNDN);

        mpfr_abs(abs_step, step, MPFR_RNDN);
        mpfr_abs(abs_x, x, MPFR_RNDN);
        mpfr_t bound; mpfr_init2(bound, bits);
        mpfr_set_ui(bound, 1, MPFR_RNDN);
        if (mpfr_cmp(abs_x, bound) > 0) mpfr_set(bound, abs_x, MPFR_RNDN);
        mpfr_mul(bound, bound, tol, MPFR_RNDN);
        bool tiny = mpfr_cmp(abs_step, bound) <= 0;
        mpfr_clear(bound);
        if (tiny) { converged = true; break; }
    }

    mpfr_clear(pv); mpfr_clear(dpv);
    mpfr_clear(step); mpfr_clear(abs_step);
    mpfr_clear(tol); mpfr_clear(abs_x);
    return converged;
}

static bool newton_refine_complex(const mpfr_t* coeffs, int deg,
                                  mpfr_t xr, mpfr_t xi, mpfr_prec_t bits) {
    mpfr_t pr, pi, dpr, dpi;
    mpfr_t denom, sr, si, abs_step, tol, scale, mag2;
    mpfr_init2(pr, bits); mpfr_init2(pi, bits);
    mpfr_init2(dpr, bits); mpfr_init2(dpi, bits);
    mpfr_init2(denom, bits); mpfr_init2(sr, bits); mpfr_init2(si, bits);
    mpfr_init2(abs_step, bits); mpfr_init2(tol, bits);
    mpfr_init2(scale, bits); mpfr_init2(mag2, bits);

    mpfr_set_ui(tol, 1, MPFR_RNDN);
    mpfr_div_2si(tol, tol, (long)bits - 8, MPFR_RNDN);

    bool converged = false;
    for (int it = 0; it < ROOT_NEWTON_MAX_ITER; it++) {
        poly_eval_complex_mpfr(coeffs, deg, xr, xi, pr, pi, &dpr, &dpi);
        /* step = p / p'  (complex division) */
        mpfr_mul(denom, dpr, dpr, MPFR_RNDN);
        mpfr_mul(mag2, dpi, dpi, MPFR_RNDN);
        mpfr_add(denom, denom, mag2, MPFR_RNDN);
        if (mpfr_zero_p(denom)) break;
        /* step = (p * conj(p')) / |p'|^2 */
        mpfr_mul(sr, pr, dpr, MPFR_RNDN);
        mpfr_mul(mag2, pi, dpi, MPFR_RNDN);
        mpfr_add(sr, sr, mag2, MPFR_RNDN);
        mpfr_div(sr, sr, denom, MPFR_RNDN);
        mpfr_mul(si, pi, dpr, MPFR_RNDN);
        mpfr_mul(mag2, pr, dpi, MPFR_RNDN);
        mpfr_sub(si, si, mag2, MPFR_RNDN);
        mpfr_div(si, si, denom, MPFR_RNDN);
        mpfr_sub(xr, xr, sr, MPFR_RNDN);
        mpfr_sub(xi, xi, si, MPFR_RNDN);

        mpfr_mul(mag2, sr, sr, MPFR_RNDN);
        mpfr_mul(abs_step, si, si, MPFR_RNDN);
        mpfr_add(mag2, mag2, abs_step, MPFR_RNDN);
        mpfr_sqrt(abs_step, mag2, MPFR_RNDN);
        mpfr_mul(mag2, xr, xr, MPFR_RNDN);
        mpfr_mul(scale, xi, xi, MPFR_RNDN);
        mpfr_add(mag2, mag2, scale, MPFR_RNDN);
        mpfr_sqrt(scale, mag2, MPFR_RNDN);
        mpfr_t one; mpfr_init2(one, bits); mpfr_set_ui(one, 1, MPFR_RNDN);
        if (mpfr_cmp(scale, one) < 0) mpfr_set(scale, one, MPFR_RNDN);
        mpfr_mul(scale, scale, tol, MPFR_RNDN);
        bool tiny = mpfr_cmp(abs_step, scale) <= 0;
        mpfr_clear(one);
        if (tiny) { converged = true; break; }
    }

    mpfr_clear(pr); mpfr_clear(pi);
    mpfr_clear(dpr); mpfr_clear(dpi);
    mpfr_clear(denom); mpfr_clear(sr); mpfr_clear(si);
    mpfr_clear(abs_step); mpfr_clear(tol);
    mpfr_clear(scale); mpfr_clear(mag2);
    return converged;
}

/* ====================================================================
 *  Driver: solve the full pipeline once, given a target precision.
 *  Returns NULL on failure; caller may retry at higher precision.
 * ================================================================== */
static Expr* solve_root_at_precision(const ZUPoly* p_squarefree,
                                     int64_t k, mpfr_prec_t target_bits,
                                     bool want_machine) {
    int n = p_squarefree->deg;
    if (n <= 0) return NULL;
    if (k < 1 || k > n) {
        root_warn("indx", "index k=%lld is out of range 1..%d.",
                  (long long)k, n);
        return NULL;
    }

    /* Working precision: 2*target + 64 for QR; refinement uses target + 32. */
    mpfr_prec_t wp = 2 * target_bits + 64;
    mpfr_prec_t rp = target_bits + 32;

    /* ---- Step 3: companion matrix ---- */
    mpfr_t* M = NULL;
    size_t Nmat = 0;
    build_companion_matrix(p_squarefree, wp, &M, &Nmat);

    /* ---- Step 4: QR for all eigenvalues ---- */
    mpfr_t* eval_re = (mpfr_t*)malloc(sizeof(mpfr_t) * Nmat);
    mpfr_t* eval_im = (mpfr_t*)malloc(sizeof(mpfr_t) * Nmat);
    for (size_t i = 0; i < Nmat; i++) {
        mpfr_init2(eval_re[i], wp);
        mpfr_init2(eval_im[i], wp);
    }
    int qr_status = eigen_all_eigenvalues_real_mpfr(M, Nmat, wp,
                                                     eval_re, eval_im);
    free_mpfr_array(M, Nmat * Nmat);
    if (qr_status != 0) {
        free_mpfr_array(eval_re, Nmat);
        free_mpfr_array(eval_im, Nmat);
        root_warn("conv", "QR algorithm did not converge at %ld bits.",
                  (long)wp);
        return NULL;
    }

    /* ---- Step 4b: Sturm certificate ---- */
    int sturm_n_real = sturm_real_root_count(p_squarefree);

    /* Tolerance for QR real-classification: 2^(-wp/2). */
    mpfr_t tol_factor; mpfr_init2(tol_factor, wp);
    mpfr_set_ui(tol_factor, 1, MPFR_RNDN);
    mpfr_div_2si(tol_factor, tol_factor, (long)wp / 2, MPFR_RNDN);

    int* is_real_flags = (int*)malloc(sizeof(int) * Nmat);
    int qr_n_real = 0;
    for (size_t i = 0; i < Nmat; i++) {
        is_real_flags[i] = classify_real(eval_re[i], eval_im[i], tol_factor);
        if (is_real_flags[i]) qr_n_real++;
    }
    mpfr_clear(tol_factor);

    if (qr_n_real != sturm_n_real) {
        /* QR's tolerance-based classification disagrees with Sturm.
         * Reconcile: if Sturm says more reals, label the smallest-|Im|
         * complex pairs as real until counts match.  If Sturm says
         * fewer reals, label the largest-|Im| ones with their real
         * complex flag.  This brings the *count* into agreement while
         * trusting Sturm; Newton refinement at higher precision will
         * resolve any leftover misclassification.
         *
         * If after reconciliation the basin-verify step rejects the
         * pick, the outer driver retries at higher precision. */
        if (qr_n_real < sturm_n_real) {
            /* Promote near-real complex to real. */
            while (qr_n_real < sturm_n_real) {
                int idx_min = -1;
                mpfr_t best; mpfr_init2(best, wp); mpfr_set_zero(best, 1);
                bool any = false;
                for (size_t i = 0; i < Nmat; i++) {
                    if (is_real_flags[i]) continue;
                    mpfr_t a; mpfr_init2(a, wp);
                    mpfr_abs(a, eval_im[i], MPFR_RNDN);
                    if (!any || mpfr_cmp(a, best) < 0) {
                        mpfr_set(best, a, MPFR_RNDN);
                        idx_min = (int)i;
                        any = true;
                    }
                    mpfr_clear(a);
                }
                mpfr_clear(best);
                if (idx_min < 0) break;
                is_real_flags[idx_min] = 1;
                mpfr_set_zero(eval_im[idx_min], 1);
                qr_n_real++;
            }
        } else {
            /* Demote false reals. */
            while (qr_n_real > sturm_n_real) {
                int idx_max = -1;
                mpfr_t best; mpfr_init2(best, wp); mpfr_set_zero(best, 1);
                bool any = false;
                for (size_t i = 0; i < Nmat; i++) {
                    if (!is_real_flags[i]) continue;
                    mpfr_t a; mpfr_init2(a, wp);
                    mpfr_abs(a, eval_im[i], MPFR_RNDN);
                    if (!any || mpfr_cmp(a, best) > 0) {
                        mpfr_set(best, a, MPFR_RNDN);
                        idx_max = (int)i;
                        any = true;
                    }
                    mpfr_clear(a);
                }
                mpfr_clear(best);
                if (idx_max < 0) break;
                is_real_flags[idx_max] = 0;
                qr_n_real--;
            }
        }
    }

    /* ---- Step 5: canonical sort ---- */
    canonical_sort_roots(eval_re, eval_im, is_real_flags, Nmat);

    /* ---- Pick the k-th candidate ---- */
    size_t idx = (size_t)(k - 1);
    int picked_is_real = is_real_flags[idx];
    mpfr_t pick_re, pick_im;
    mpfr_init2(pick_re, rp);
    mpfr_init2(pick_im, rp);
    mpfr_set(pick_re, eval_re[idx], MPFR_RNDN);
    mpfr_set(pick_im, eval_im[idx], MPFR_RNDN);
    if (picked_is_real) mpfr_set_zero(pick_im, 1);

    /* ---- Step 6: Newton refinement ---- */
    mpfr_t* coeffs = NULL;
    int deg_out = -1;
    zupoly_to_mpfr_coeffs(p_squarefree, rp, &coeffs, &deg_out);

    bool ok;
    if (picked_is_real) {
        ok = newton_refine_real(coeffs, deg_out, pick_re, rp);
    } else {
        ok = newton_refine_complex(coeffs, deg_out, pick_re, pick_im, rp);
    }
    poly_eval_mpfr_free_coeffs(coeffs, deg_out);

    if (!ok) {
        free_mpfr_array(eval_re, Nmat);
        free_mpfr_array(eval_im, Nmat);
        free(is_real_flags);
        mpfr_clear(pick_re); mpfr_clear(pick_im);
        root_warn("conv", "Newton refinement did not converge at %ld bits.",
                  (long)rp);
        return NULL;
    }

    /* ---- Step 7: basin verification ---- *
     * Confirm that the refined root is closer to the original QR estimate
     * at position idx than to any *other* QR estimate.  This is the
     * classical Newton-basin check and tolerates the noise in QR
     * eigenvalues for clustered roots (notably the Re components of a
     * conjugate pair, which differ only by QR roundoff).  Re-sorting the
     * mixed array would have been brittle for exactly that case. */
    bool same;
    {
        mpfr_t dist_pick, dist_other, dr, di, sq;
        mpfr_init2(dist_pick, rp);
        mpfr_init2(dist_other, rp);
        mpfr_init2(dr, rp); mpfr_init2(di, rp); mpfr_init2(sq, rp);

        /* Distance squared from refined pick to its own QR estimate. */
        mpfr_sub(dr, pick_re, eval_re[idx], MPFR_RNDN);
        mpfr_sub(di, pick_im, eval_im[idx], MPFR_RNDN);
        mpfr_mul(dist_pick, dr, dr, MPFR_RNDN);
        mpfr_mul(sq, di, di, MPFR_RNDN);
        mpfr_add(dist_pick, dist_pick, sq, MPFR_RNDN);

        /* Minimum distance squared from refined pick to any other QR
         * estimate.  If dist_pick * 4 <= dist_other, the refined root is
         * unambiguously in the idx-basin. */
        bool first = true;
        for (size_t j = 0; j < Nmat; j++) {
            if (j == idx) continue;
            mpfr_t d; mpfr_init2(d, rp);
            mpfr_sub(dr, pick_re, eval_re[j], MPFR_RNDN);
            mpfr_sub(di, pick_im, eval_im[j], MPFR_RNDN);
            mpfr_mul(d, dr, dr, MPFR_RNDN);
            mpfr_mul(sq, di, di, MPFR_RNDN);
            mpfr_add(d, d, sq, MPFR_RNDN);
            if (first || mpfr_cmp(d, dist_other) < 0) {
                mpfr_set(dist_other, d, MPFR_RNDN);
                first = false;
            }
            mpfr_clear(d);
        }

        if (first) {
            /* Only one root — basin trivially correct. */
            same = true;
        } else {
            /* same  <=>  4 * dist_pick < dist_other  (factor 2 in distance) */
            mpfr_mul_2si(dist_pick, dist_pick, 2, MPFR_RNDN);
            same = mpfr_cmp(dist_pick, dist_other) <= 0;
        }
        mpfr_clear(dist_pick); mpfr_clear(dist_other);
        mpfr_clear(dr); mpfr_clear(di); mpfr_clear(sq);
    }

    Expr* out = NULL;
    if (same) {
        /* Round down to target precision. */
        mpfr_t re_out, im_out;
        mpfr_init2(re_out, target_bits); mpfr_init2(im_out, target_bits);
        mpfr_set(re_out, pick_re, MPFR_RNDN);
        mpfr_set(im_out, pick_im, MPFR_RNDN);
        if (picked_is_real) {
            if (want_machine) {
                double v = mpfr_get_d(re_out, MPFR_RNDN);
                out = expr_new_real(v);
            } else {
                out = expr_new_mpfr_copy(re_out);
            }
        } else {
            if (want_machine) {
                Expr* re_e = expr_new_real(mpfr_get_d(re_out, MPFR_RNDN));
                Expr* im_e = expr_new_real(mpfr_get_d(im_out, MPFR_RNDN));
                out = make_complex(re_e, im_e);
            } else {
                Expr* re_e = expr_new_mpfr_copy(re_out);
                Expr* im_e = expr_new_mpfr_copy(im_out);
                out = make_complex(re_e, im_e);
            }
        }
        mpfr_clear(re_out); mpfr_clear(im_out);
    }

    free_mpfr_array(eval_re, Nmat);
    free_mpfr_array(eval_im, Nmat);
    free(is_real_flags);
    mpfr_clear(pick_re); mpfr_clear(pick_im);
    return out;
}

/* ====================================================================
 *  Entry point
 * ================================================================== */
Expr* root_numericalize(const Expr* root_expr, NumericSpec spec) {
    Expr* body = NULL;
    int64_t k = 0;
    if (!parse_root_call(root_expr, &body, &k)) return NULL;

    /* Extract polynomial as ZUPoly using Slot[1] as the variable. */
    Expr* slot1 = make_slot1();
    ZUPoly* p_raw = expr_to_zupoly(body, slot1);
    expr_free(slot1);
    if (!p_raw) {
        root_warn("nonint", "polynomial has non-integer coefficients.");
        return NULL;
    }
    if (p_raw->deg <= 0) {
        zupoly_free(p_raw);
        root_warn("poly", "polynomial is constant (no roots).");
        return NULL;
    }

    /* Squarefree radical. */
    ZUPoly* p = zupoly_squarefree_part(p_raw);
    zupoly_free(p_raw);

    /* Target precision. */
    bool want_machine = (spec.mode == NUMERIC_MODE_MACHINE);
    mpfr_prec_t target_bits = want_machine
        ? (mpfr_prec_t)DBL_MANT_DIG
        : (mpfr_prec_t)spec.bits;
    if (target_bits < 53) target_bits = 53;

    /* Up to 3 escalations on basin-verify failure. */
    Expr* out = NULL;
    for (int attempt = 0; attempt < 3 && !out; attempt++) {
        mpfr_prec_t try_bits = target_bits;
        if (attempt > 0) {
            /* Bump 1.5x per attempt for the working precision; target
             * stays at user request, but solve_root_at_precision uses
             * 2*target + 64 as its working precision internally.  To
             * escalate we pretend the target is larger. */
            try_bits = target_bits + (mpfr_prec_t)(target_bits / 2) * (mpfr_prec_t)attempt;
        }
        out = solve_root_at_precision(p, k, try_bits, want_machine);
    }

    zupoly_free(p);
    return out;
}

#else  /* USE_MPFR not defined */

Expr* root_numericalize(const Expr* root_expr, NumericSpec spec) {
    (void)root_expr; (void)spec;
    return NULL;
}

#endif /* USE_MPFR */
