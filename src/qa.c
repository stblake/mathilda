/* qa.c — Q(α) algebraic-number type (Phase G1, Trager's algorithm).
 *
 * Implementation overview:
 *
 * - QAExt holds the minimal polynomial P_α as an mpq_t array.
 * - QANum holds an element of Q(α) as an mpq_t array indexed by
 *   α-degree, of length ext->deg.
 * - Multiplication uses a temp array of length 2*ext->deg - 1 (full
 *   polynomial product) followed by reduction modulo P_α.
 * - Inversion uses extended Euclidean over Q[y] between the QANum
 *   (lifted to its polynomial-in-y representation) and P_α.
 *
 * All polynomial work is done with plain mpq_t arrays; we deliberately
 * do NOT depend on the private QUPoly substrate in zupoly.c so that
 * Phase G1 stays a self-contained module. (zupoly.c's QUPoly is also
 * a candidate for promotion to a public substrate — see FACTOR_PLAN.md
 * §14 — but that refactor is orthogonal and deferred.) */

#include "qa.h"
#include <stdlib.h>
#include <string.h>

/* ====================== mpq_t-array helpers (private) =================== */

/* Allocate and zero-init an array of `n` mpq_t. */
static mpq_t* mpq_array_new_zero(size_t n) {
    mpq_t* a = (mpq_t*)malloc(sizeof(mpq_t) * n);
    for (size_t i = 0; i < n; i++) mpq_init(a[i]);
    return a;
}

/* Clear and free an array of `n` mpq_t. */
static void mpq_array_free(mpq_t* a, size_t n) {
    if (!a) return;
    for (size_t i = 0; i < n; i++) mpq_clear(a[i]);
    free(a);
}

/* Find the highest index i where a[i] != 0. Returns -1 if a is all zero. */
static long mpq_array_deg(const mpq_t* a, size_t n) {
    for (long i = (long)n - 1; i >= 0; i--) {
        if (mpq_sgn(a[i]) != 0) return i;
    }
    return -1;
}

/* ============================== QAExt =============================== */

QAExt* qaext_new(size_t deg) {
    QAExt* ext = (QAExt*)malloc(sizeof(QAExt));
    ext->deg = deg;
    ext->coef = mpq_array_new_zero(deg + 1);
    return ext;
}

void qaext_free(QAExt* ext) {
    if (!ext) return;
    mpq_array_free(ext->coef, ext->deg + 1);
    free(ext);
}

void qaext_set_coef_si(QAExt* ext, size_t i, long num, long den) {
    mpq_set_si(ext->coef[i], num, (unsigned long)den);
    mpq_canonicalize(ext->coef[i]);
}

void qaext_set_coef_mpq(QAExt* ext, size_t i, const mpq_t v) {
    mpq_set(ext->coef[i], v);
}

QAExt* qaext_sqrt_si(long c) {
    /* P_α(y) = y^2 - c */
    QAExt* ext = qaext_new(2);
    qaext_set_coef_si(ext, 0, -c, 1);
    qaext_set_coef_si(ext, 1, 0, 1);
    qaext_set_coef_si(ext, 2, 1, 1);
    return ext;
}

QAExt* qaext_root_si(long c, unsigned n) {
    /* P_α(y) = y^n - c */
    QAExt* ext = qaext_new(n);
    qaext_set_coef_si(ext, 0, -c, 1);
    for (unsigned i = 1; i < n; i++) qaext_set_coef_si(ext, i, 0, 1);
    qaext_set_coef_si(ext, n, 1, 1);
    return ext;
}

/* =========================== QANum construction ========================== */

QANum* qa_zero(const QAExt* ext) {
    QANum* a = (QANum*)malloc(sizeof(QANum));
    a->ext = ext;
    a->coef = mpq_array_new_zero(ext->deg);
    return a;
}

QANum* qa_one(const QAExt* ext) {
    QANum* a = qa_zero(ext);
    if (ext->deg > 0) mpq_set_si(a->coef[0], 1, 1);
    return a;
}

QANum* qa_alpha(const QAExt* ext) {
    QANum* a = qa_zero(ext);
    if (ext->deg > 1) mpq_set_si(a->coef[1], 1, 1);
    /* else: deg P_α = 1, so α is already determined as -P_α(0)/lc;
     * for the trivial-extension case we just return zero. */
    return a;
}

QANum* qa_from_mpq(const QAExt* ext, const mpq_t v) {
    QANum* a = qa_zero(ext);
    if (ext->deg > 0) mpq_set(a->coef[0], v);
    return a;
}

QANum* qa_from_si(const QAExt* ext, long num, long den) {
    QANum* a = qa_zero(ext);
    if (ext->deg > 0) {
        mpq_set_si(a->coef[0], num, (unsigned long)den);
        mpq_canonicalize(a->coef[0]);
    }
    return a;
}

QANum* qa_copy(const QANum* a) {
    QANum* r = qa_zero(a->ext);
    for (size_t i = 0; i < a->ext->deg; i++) mpq_set(r->coef[i], a->coef[i]);
    return r;
}

void qa_free(QANum* a) {
    if (!a) return;
    mpq_array_free(a->coef, a->ext->deg);
    free(a);
}

/* ============================== Predicates ============================== */

bool qa_is_zero(const QANum* a) {
    return mpq_array_deg(a->coef, a->ext->deg) < 0;
}

bool qa_is_one(const QANum* a) {
    if (a->ext->deg == 0) return false;
    if (mpq_cmp_si(a->coef[0], 1, 1) != 0) return false;
    for (size_t i = 1; i < a->ext->deg; i++) {
        if (mpq_sgn(a->coef[i]) != 0) return false;
    }
    return true;
}

bool qa_eq(const QANum* a, const QANum* b) {
    if (a->ext != b->ext) return false; /* extensions must match by identity */
    for (size_t i = 0; i < a->ext->deg; i++) {
        if (!mpq_equal(a->coef[i], b->coef[i])) return false;
    }
    return true;
}

/* =========================== Arithmetic =============================== */

QANum* qa_add(const QANum* a, const QANum* b) {
    QANum* r = qa_zero(a->ext);
    for (size_t i = 0; i < a->ext->deg; i++) {
        mpq_add(r->coef[i], a->coef[i], b->coef[i]);
    }
    return r;
}

QANum* qa_sub(const QANum* a, const QANum* b) {
    QANum* r = qa_zero(a->ext);
    for (size_t i = 0; i < a->ext->deg; i++) {
        mpq_sub(r->coef[i], a->coef[i], b->coef[i]);
    }
    return r;
}

QANum* qa_neg(const QANum* a) {
    QANum* r = qa_zero(a->ext);
    for (size_t i = 0; i < a->ext->deg; i++) {
        mpq_neg(r->coef[i], a->coef[i]);
    }
    return r;
}

QANum* qa_scale_mpq(const QANum* a, const mpq_t k) {
    QANum* r = qa_zero(a->ext);
    for (size_t i = 0; i < a->ext->deg; i++) {
        mpq_mul(r->coef[i], a->coef[i], k);
    }
    return r;
}

QANum* qa_scale_si(const QANum* a, long num, long den) {
    mpq_t k; mpq_init(k);
    mpq_set_si(k, num, (unsigned long)den);
    mpq_canonicalize(k);
    QANum* r = qa_scale_mpq(a, k);
    mpq_clear(k);
    return r;
}

/* Reduce a polynomial-in-α represented in `acc` (length `acc_n`) modulo
 * P_α. The result lives in acc[0..ext->deg-1]; acc[ext->deg..acc_n-1]
 * is left as garbage and not consulted by callers.
 *
 * Algorithm: standard polynomial long division. While acc has degree
 * d ≥ ext->deg, subtract (acc[d] / lc(P_α)) * y^(d - ext->deg) * P_α.
 * Each iteration zeroes out the current leading term and shifts d down
 * by 1 (or more if cancellation occurs). */
static void reduce_mod_palpha(mpq_t* acc, size_t acc_n, const QAExt* ext) {
    if (ext->deg == 0) return; /* trivial extension */
    long d = (long)acc_n - 1;
    /* Find true leading degree. */
    while (d >= 0 && mpq_sgn(acc[d]) == 0) d--;
    mpq_t coef, tmp;
    mpq_init(coef);
    mpq_init(tmp);
    while (d >= (long)ext->deg) {
        /* Quotient digit: q_d = acc[d] / lc(P_α). */
        mpq_div(coef, acc[d], ext->coef[ext->deg]);
        /* Subtract q_d * y^(d - ext->deg) * P_α from acc. The
         * y^(d - ext->deg) shift means P_α's coefficient i lands
         * at position (d - ext->deg) + i. */
        long shift = d - (long)ext->deg;
        for (size_t i = 0; i <= ext->deg; i++) {
            mpq_mul(tmp, coef, ext->coef[i]);
            mpq_sub(acc[shift + i], acc[shift + i], tmp);
        }
        /* acc[d] is now zero by construction; rescan downwards. */
        do { d--; } while (d >= 0 && mpq_sgn(acc[d]) == 0);
    }
    mpq_clear(coef);
    mpq_clear(tmp);
}

QANum* qa_mul(const QANum* a, const QANum* b) {
    const QAExt* ext = a->ext;
    if (ext->deg == 0) return qa_zero(ext);
    /* Compute the full polynomial product into a temp of length
     * 2*deg - 1, then reduce modulo P_α. */
    size_t n = 2 * ext->deg - 1;
    mpq_t* acc = mpq_array_new_zero(n);
    mpq_t tmp; mpq_init(tmp);
    for (size_t i = 0; i < ext->deg; i++) {
        if (mpq_sgn(a->coef[i]) == 0) continue;
        for (size_t j = 0; j < ext->deg; j++) {
            mpq_mul(tmp, a->coef[i], b->coef[j]);
            mpq_add(acc[i + j], acc[i + j], tmp);
        }
    }
    mpq_clear(tmp);
    reduce_mod_palpha(acc, n, ext);
    QANum* r = qa_zero(ext);
    for (size_t i = 0; i < ext->deg; i++) mpq_set(r->coef[i], acc[i]);
    mpq_array_free(acc, n);
    return r;
}

/* ============================== Inversion ============================== */
/* Extended Euclidean over Q[y] between a (lifted as a polynomial of
 * degree < ext->deg) and P_α. On exit:
 *   gcd is monic, stored in g_arr (of length g_len);
 *   s_arr satisfies s * a + t * P_α = g.
 *
 * If g is a unit (degree 0, i.e. a non-zero constant), then s_arr / g
 * is the inverse of a modulo P_α. If g has degree > 0, then a and P_α
 * share a non-trivial factor — possible only when P_α is reducible. */

/* Polynomial long division over Q[y]:
 *   a / b  =>  q (length cap) and r (in-place over a copy of a).
 * Caller passes a working buffer for the dividend; we mutate it. */
static void mpq_array_divrem(mpq_t* dividend, size_t dividend_n,
                              const mpq_t* divisor, size_t divisor_n,
                              mpq_t* q_out, size_t q_n) {
    long bd = (long)divisor_n - 1;
    while (bd >= 0 && mpq_sgn(divisor[bd]) == 0) bd--;
    if (bd < 0) return; /* divide-by-zero: caller should not reach here */
    long ad = (long)dividend_n - 1;
    while (ad >= 0 && mpq_sgn(dividend[ad]) == 0) ad--;

    /* Zero out q. */
    for (size_t i = 0; i < q_n; i++) mpq_set_ui(q_out[i], 0, 1);

    mpq_t coef, tmp;
    mpq_init(coef); mpq_init(tmp);
    while (ad >= bd) {
        long shift = ad - bd;
        mpq_div(coef, dividend[ad], divisor[bd]);
        if ((size_t)shift < q_n) mpq_set(q_out[shift], coef);
        for (long i = 0; i <= bd; i++) {
            mpq_mul(tmp, coef, divisor[i]);
            mpq_sub(dividend[shift + i], dividend[shift + i], tmp);
        }
        do { ad--; } while (ad >= 0 && mpq_sgn(dividend[ad]) == 0);
    }
    mpq_clear(coef); mpq_clear(tmp);
}

QANum* qa_inverse(const QANum* a) {
    if (qa_is_zero(a)) return NULL;
    const QAExt* ext = a->ext;
    size_t n = ext->deg + 1; /* working length: same as P_α */

    /* Initialise r0 = a (extended to length n), r1 = P_α.
     * s0 = 1, s1 = 0 (each as a polynomial in y, length n). */
    mpq_t* r0 = mpq_array_new_zero(n);
    mpq_t* r1 = mpq_array_new_zero(n);
    mpq_t* s0 = mpq_array_new_zero(n);
    mpq_t* s1 = mpq_array_new_zero(n);
    for (size_t i = 0; i < ext->deg; i++) mpq_set(r0[i], a->coef[i]);
    for (size_t i = 0; i <= ext->deg; i++) mpq_set(r1[i], ext->coef[i]);
    mpq_set_si(s0[0], 1, 1);

    mpq_t* q  = mpq_array_new_zero(n);
    mpq_t* qs = mpq_array_new_zero(n);
    mpq_t  tmp; mpq_init(tmp);

    while (mpq_array_deg(r1, n) >= 0) {
        /* (q, r2) = divrem(r0, r1). The remainder lives in-place in
         * a working buffer (we copy r0, then divrem mutates it into r2). */
        mpq_t* work = mpq_array_new_zero(n);
        for (size_t i = 0; i < n; i++) mpq_set(work[i], r0[i]);
        for (size_t i = 0; i < n; i++) mpq_set_ui(q[i], 0, 1);
        mpq_array_divrem(work, n, r1, n, q, n);
        /* `work` now holds r2 = r0 - q * r1. */

        /* s2 = s0 - q * s1.  Compute qs = q * s1 into a length-(2n-1)
         * scratch (truncated to n since deg s* < ext->deg ≤ n-1 and
         * deg q < ext->deg ≤ n-1, so deg(q*s1) ≤ 2*deg-2 < 2*(n-1)
         * but we never use slots ≥ n in s_*). */
        size_t conv_n = 2 * n - 1;
        mpq_t* prod = mpq_array_new_zero(conv_n);
        for (size_t i = 0; i < n; i++) {
            if (mpq_sgn(q[i]) == 0) continue;
            for (size_t j = 0; j < n; j++) {
                if (mpq_sgn(s1[j]) == 0) continue;
                mpq_mul(tmp, q[i], s1[j]);
                mpq_add(prod[i + j], prod[i + j], tmp);
            }
        }
        for (size_t i = 0; i < n; i++) {
            mpq_sub(qs[i], s0[i], prod[i]);
        }
        mpq_array_free(prod, conv_n);

        /* Rotate: r0,r1 <- r1,r2 ; s0,s1 <- s1,qs. */
        for (size_t i = 0; i < n; i++) mpq_set(r0[i], r1[i]);
        for (size_t i = 0; i < n; i++) mpq_set(r1[i], work[i]);
        for (size_t i = 0; i < n; i++) {
            mpq_swap(s0[i], s1[i]);
            mpq_set(s1[i], qs[i]);
        }
        mpq_array_free(work, n);
    }

    /* gcd is now r0. If its degree is 0, it's a non-zero constant unit;
     * inverse is s0 / r0[0]. Otherwise a is non-invertible. */
    QANum* result = NULL;
    long gd = mpq_array_deg(r0, n);
    if (gd == 0) {
        result = qa_zero(ext);
        mpq_t inv; mpq_init(inv);
        mpq_inv(inv, r0[0]);
        for (size_t i = 0; i < ext->deg; i++) {
            mpq_mul(result->coef[i], s0[i], inv);
        }
        mpq_clear(inv);
    }

    mpq_clear(tmp);
    mpq_array_free(r0, n);
    mpq_array_free(r1, n);
    mpq_array_free(s0, n);
    mpq_array_free(s1, n);
    mpq_array_free(q, n);
    mpq_array_free(qs, n);
    return result;
}

QANum* qa_div(const QANum* a, const QANum* b) {
    QANum* binv = qa_inverse(b);
    if (!binv) return NULL;
    QANum* r = qa_mul(a, binv);
    qa_free(binv);
    return r;
}
