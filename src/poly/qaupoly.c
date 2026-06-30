/* qaupoly.c — Q(α)[x] univariate polynomials (Phase G2, Trager).
 *
 * Layered on top of qa.{c,h}.  All polynomial operations are simple
 * generalisations of the integer / rational univariate operations:
 *
 *   - add/sub/neg: pairwise on coefficients.
 *   - mul: classical convolution with qa_mul + qa_add accumulator.
 *   - divrem: long division.  Q(α) is a field for irreducible P_α,
 *     so qa_inverse on the divisor's leading coefficient always
 *     succeeds; we return false defensively when it doesn't (which
 *     would indicate the caller passed a reducible P_α).
 *   - gcd: standard Euclidean PRS, normalised monic.
 *   - shift: f(x + c) via Horner expansion (build the result
 *     incrementally; no nested copies).
 *   - eval: Horner from the high coefficient down.
 *
 * We deliberately keep the arithmetic naïve (no Karatsuba, no PRS
 * pseudo-remainder optimisation) because Q(α) coefficients are
 * already heavyweight (mpq_t arrays of length deg P_α) and the
 * algorithms we plug into in G3+G4 are dominated by the resultant
 * computation, not the inner polynomial loop. */

#include "qaupoly.h"
#include <stdlib.h>
#include <string.h>

/* ============================ Construction ============================ */

QAUPoly* qaupoly_new(const QAExt* ext, int cap) {
    if (cap < 1) cap = 1;
    QAUPoly* p = (QAUPoly*)malloc(sizeof(QAUPoly));
    p->ext = ext;
    p->deg = -1;
    p->cap = cap;
    p->c = (QANum**)malloc(sizeof(QANum*) * (size_t)cap);
    for (int i = 0; i < cap; i++) p->c[i] = qa_zero(ext);
    return p;
}

QAUPoly* qaupoly_zero(const QAExt* ext) {
    return qaupoly_new(ext, 1);
}

QAUPoly* qaupoly_from_qa(const QANum* v) {
    QAUPoly* p = qaupoly_new(v->ext, 1);
    qa_free(p->c[0]);
    p->c[0] = qa_copy(v);
    p->deg = qa_is_zero(v) ? -1 : 0;
    return p;
}

QAUPoly* qaupoly_from_si(const QAExt* ext, long num, long den) {
    QAUPoly* p = qaupoly_new(ext, 1);
    qa_free(p->c[0]);
    p->c[0] = qa_from_si(ext, num, den);
    p->deg = qa_is_zero(p->c[0]) ? -1 : 0;
    return p;
}

QAUPoly* qaupoly_x(const QAExt* ext) {
    QAUPoly* p = qaupoly_new(ext, 2);
    qa_free(p->c[1]);
    p->c[1] = qa_one(ext);
    p->deg = 1;
    return p;
}

QAUPoly* qaupoly_copy(const QAUPoly* p) {
    int cap = (p->deg < 0) ? 1 : p->deg + 1;
    QAUPoly* r = qaupoly_new(p->ext, cap);
    for (int i = 0; i <= p->deg; i++) {
        qa_free(r->c[i]);
        r->c[i] = qa_copy(p->c[i]);
    }
    r->deg = p->deg;
    return r;
}

void qaupoly_free(QAUPoly* p) {
    if (!p) return;
    for (int i = 0; i < p->cap; i++) qa_free(p->c[i]);
    free(p->c);
    free(p);
}

/* ====================== Coefficient access ======================== */

static void qaupoly_grow(QAUPoly* p, int new_cap) {
    if (new_cap <= p->cap) return;
    p->c = (QANum**)realloc(p->c, sizeof(QANum*) * (size_t)new_cap);
    for (int i = p->cap; i < new_cap; i++) p->c[i] = qa_zero(p->ext);
    p->cap = new_cap;
}

void qaupoly_setcoef(QAUPoly* p, int i, const QANum* v) {
    if (i >= p->cap) qaupoly_grow(p, i + 1);
    qa_free(p->c[i]);
    p->c[i] = qa_copy(v);
    if (i > p->deg && !qa_is_zero(v)) p->deg = i;
}

const QANum* qaupoly_getcoef(const QAUPoly* p, int i) {
    if (i < 0 || i >= p->cap) {
        /* Out of range: return the c[0] slot (always exists, valid until
         * the polynomial is freed).  When p->deg < 0 c[0] is also zero,
         * so this is a defensible degenerate. */
        return p->c[0];
    }
    return p->c[i];
}

void qaupoly_normalize(QAUPoly* p) {
    int d = p->cap - 1;
    while (d >= 0 && qa_is_zero(p->c[d])) d--;
    p->deg = d;
}

/* ============================== Predicates ============================== */

bool qaupoly_is_zero(const QAUPoly* p) {
    return p->deg < 0;
}

bool qaupoly_eq(const QAUPoly* a, const QAUPoly* b) {
    if (a->deg != b->deg) return false;
    for (int i = 0; i <= a->deg; i++) {
        if (!qa_eq(a->c[i], b->c[i])) return false;
    }
    return true;
}

/* ============================== Arithmetic ============================== */

QAUPoly* qaupoly_add(const QAUPoly* a, const QAUPoly* b) {
    int max_deg = (a->deg > b->deg) ? a->deg : b->deg;
    if (max_deg < 0) return qaupoly_zero(a->ext);
    QAUPoly* r = qaupoly_new(a->ext, max_deg + 1);
    for (int i = 0; i <= max_deg; i++) {
        qa_free(r->c[i]);
        if (i <= a->deg && i <= b->deg) {
            r->c[i] = qa_add(a->c[i], b->c[i]);
        } else if (i <= a->deg) {
            r->c[i] = qa_copy(a->c[i]);
        } else {
            r->c[i] = qa_copy(b->c[i]);
        }
    }
    qaupoly_normalize(r);
    return r;
}

QAUPoly* qaupoly_sub(const QAUPoly* a, const QAUPoly* b) {
    int max_deg = (a->deg > b->deg) ? a->deg : b->deg;
    if (max_deg < 0) return qaupoly_zero(a->ext);
    QAUPoly* r = qaupoly_new(a->ext, max_deg + 1);
    for (int i = 0; i <= max_deg; i++) {
        qa_free(r->c[i]);
        if (i <= a->deg && i <= b->deg) {
            r->c[i] = qa_sub(a->c[i], b->c[i]);
        } else if (i <= a->deg) {
            r->c[i] = qa_copy(a->c[i]);
        } else {
            r->c[i] = qa_neg(b->c[i]);
        }
    }
    qaupoly_normalize(r);
    return r;
}

QAUPoly* qaupoly_neg(const QAUPoly* a) {
    if (a->deg < 0) return qaupoly_zero(a->ext);
    QAUPoly* r = qaupoly_new(a->ext, a->deg + 1);
    for (int i = 0; i <= a->deg; i++) {
        qa_free(r->c[i]);
        r->c[i] = qa_neg(a->c[i]);
    }
    r->deg = a->deg;
    return r;
}

QAUPoly* qaupoly_mul(const QAUPoly* a, const QAUPoly* b) {
    if (a->deg < 0 || b->deg < 0) return qaupoly_zero(a->ext);
    int rd = a->deg + b->deg;
    QAUPoly* r = qaupoly_new(a->ext, rd + 1);
    for (int i = 0; i <= a->deg; i++) {
        if (qa_is_zero(a->c[i])) continue;
        for (int j = 0; j <= b->deg; j++) {
            if (qa_is_zero(b->c[j])) continue;
            QANum* term = qa_mul(a->c[i], b->c[j]);
            QANum* sum  = qa_add(r->c[i + j], term);
            qa_free(term);
            qa_free(r->c[i + j]);
            r->c[i + j] = sum;
        }
    }
    qaupoly_normalize(r);
    return r;
}

QAUPoly* qaupoly_scale_qa(const QAUPoly* a, const QANum* k) {
    if (a->deg < 0 || qa_is_zero(k)) return qaupoly_zero(a->ext);
    QAUPoly* r = qaupoly_new(a->ext, a->deg + 1);
    for (int i = 0; i <= a->deg; i++) {
        qa_free(r->c[i]);
        r->c[i] = qa_mul(a->c[i], k);
    }
    qaupoly_normalize(r);
    return r;
}

QAUPoly* qaupoly_scale_si(const QAUPoly* a, long num, long den) {
    QANum* k = qa_from_si(a->ext, num, den);
    QAUPoly* r = qaupoly_scale_qa(a, k);
    qa_free(k);
    return r;
}

/* Long division.  Returns false if the divisor's leading coefficient
 * is not invertible in Q(α) — this should not happen for irreducible
 * P_α, but we guard defensively. */
bool qaupoly_divrem(const QAUPoly* a, const QAUPoly* b,
                    QAUPoly** q_out, QAUPoly** r_out) {
    if (b->deg < 0) return false;
    QANum* lc_inv = qa_inverse(b->c[b->deg]);
    if (!lc_inv) return false;

    if (a->deg < b->deg) {
        *q_out = qaupoly_zero(a->ext);
        *r_out = qaupoly_copy(a);
        qa_free(lc_inv);
        return true;
    }

    int q_deg = a->deg - b->deg;
    QAUPoly* q = qaupoly_new(a->ext, q_deg + 1);
    QAUPoly* r = qaupoly_copy(a);

    for (int i = q_deg; i >= 0; i--) {
        int rd = b->deg + i;
        if (rd > r->deg || qa_is_zero(r->c[rd])) continue;
        QANum* coef = qa_mul(r->c[rd], lc_inv);
        qa_free(q->c[i]);
        q->c[i] = qa_copy(coef);
        for (int j = 0; j <= b->deg; j++) {
            QANum* term = qa_mul(coef, b->c[j]);
            QANum* diff = qa_sub(r->c[j + i], term);
            qa_free(term);
            qa_free(r->c[j + i]);
            r->c[j + i] = diff;
        }
        qa_free(coef);
    }
    qaupoly_normalize(q);
    qaupoly_normalize(r);
    qa_free(lc_inv);
    *q_out = q;
    *r_out = r;
    return true;
}

QAUPoly* qaupoly_make_monic(const QAUPoly* p) {
    if (p->deg < 0) return NULL;
    QANum* lc_inv = qa_inverse(p->c[p->deg]);
    if (!lc_inv) return NULL;
    QAUPoly* r = qaupoly_scale_qa(p, lc_inv);
    qa_free(lc_inv);
    return r;
}

QAUPoly* qaupoly_gcd(const QAUPoly* a, const QAUPoly* b) {
    if (a->deg < 0 && b->deg < 0) return NULL;
    if (a->deg < 0) return qaupoly_make_monic(b);
    if (b->deg < 0) return qaupoly_make_monic(a);
    /* Monic Euclidean: r0=a, r1=b; loop r0,r1 := r1, r0 mod r1 until r1 = 0;
     * gcd is r0 (then make monic).
     *
     * Each remainder is normalised to monic immediately.  Over a field of
     * algebraic numbers Q(alpha) the bare Euclidean PRS suffers severe
     * coefficient swell: qaupoly_divrem multiplies the dividend through by
     * 1/lc(divisor), so a non-monic divisor inflates EVERY coefficient of the
     * remainder by the (large) rational height of that leading field element,
     * and the heights compound multiplicatively step over step (degree-12
     * cyclotomic inputs produced 50+ digit numerators, making the subsequent
     * PolynomialQuotient unusably slow).  Forcing each divisor monic keeps
     * 1/lc = 1, so no leading-coefficient inflation is injected and the
     * intermediate heights stay bounded.  Monic normalisation is GCD-
     * preserving (the GCD is defined up to a unit; the final make_monic
     * already canonicalises), so the result is identical -- only smaller. */
    QAUPoly* r0 = qaupoly_make_monic(a);
    QAUPoly* r1 = qaupoly_make_monic(b);
    if (!r0 || !r1) { qaupoly_free(r0); qaupoly_free(r1); return NULL; }
    while (!qaupoly_is_zero(r1)) {
        QAUPoly *q, *r2;
        if (!qaupoly_divrem(r0, r1, &q, &r2)) {
            qaupoly_free(r0);
            qaupoly_free(r1);
            return NULL;
        }
        qaupoly_free(q);
        qaupoly_free(r0);
        r0 = r1;
        /* Renormalise the new remainder to monic before it becomes the next
         * divisor (see rationale above). */
        if (!qaupoly_is_zero(r2)) {
            QAUPoly* m = qaupoly_make_monic(r2);
            if (m) { qaupoly_free(r2); r2 = m; }
        }
        r1 = r2;
    }
    qaupoly_free(r1);
    QAUPoly* g = qaupoly_make_monic(r0);
    qaupoly_free(r0);
    return g;
}

/* ============================ Substitution ============================ */

QANum* qaupoly_eval(const QAUPoly* p, const QANum* c) {
    if (p->deg < 0) return qa_zero(p->ext);
    /* Horner from the high coefficient down:
     *   acc = p[deg]
     *   for i = deg-1 .. 0: acc = acc * c + p[i] */
    QANum* acc = qa_copy(p->c[p->deg]);
    for (int i = p->deg - 1; i >= 0; i--) {
        QANum* mul = qa_mul(acc, c);
        QANum* sum = qa_add(mul, p->c[i]);
        qa_free(mul);
        qa_free(acc);
        acc = sum;
    }
    return acc;
}

QAUPoly* qaupoly_shift(const QAUPoly* p, const QANum* c) {
    /* Horner-style construction of p(x + c):
     *   acc(x) = p[deg]
     *   for i = deg-1 .. 0: acc(x) = acc(x) * (x + c) + p[i]
     * The result is the polynomial-in-x form. */
    if (p->deg < 0) return qaupoly_zero(p->ext);

    /* (x + c) as a length-2 polynomial. */
    QAUPoly* xpc = qaupoly_new(p->ext, 2);
    qa_free(xpc->c[0]); xpc->c[0] = qa_copy(c);
    qa_free(xpc->c[1]); xpc->c[1] = qa_one(p->ext);
    xpc->deg = 1;

    QAUPoly* acc = qaupoly_from_qa(p->c[p->deg]);
    for (int i = p->deg - 1; i >= 0; i--) {
        QAUPoly* prod = qaupoly_mul(acc, xpc);
        QAUPoly* term = qaupoly_from_qa(p->c[i]);
        QAUPoly* next = qaupoly_add(prod, term);
        qaupoly_free(prod);
        qaupoly_free(term);
        qaupoly_free(acc);
        acc = next;
    }
    qaupoly_free(xpc);
    return acc;
}
