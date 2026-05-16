/*
 * zupoly.c
 * --------
 * Implementation of univariate polynomials over Z with mpz_t
 * coefficients.  See zupoly.h for type/contract documentation.
 *
 * Implementation notes:
 *
 * 1. Coefficient storage.  We allocate `cap` mpz_t slots, all
 *    initialised at construction.  Slots beyond `deg` hold the value 0.
 *    This avoids the cost of mpz_init/mpz_clear on every coefficient
 *    update at the price of a small fixed memory overhead.
 *
 * 2. Capacity growth.  When an update needs more slots than `cap`, we
 *    double `cap` (or grow to fit, whichever is larger) and mpz_init
 *    the new slots to zero.  Old slots survive the realloc unchanged.
 *
 * 3. Normalization.  After arithmetic we walk from `cap-1` down to
 *    locate the highest non-zero slot.  This keeps the `deg` field
 *    accurate so that subsequent operations see the correct shape.
 *
 * 4. GCD algorithm.  We use the subresultant pseudo-remainder
 *    sequence (Brown-Collins).  This avoids both the rational-
 *    coefficient growth of plain Euclid over Z and the simple
 *    primitive-prs algorithm's quadratic-in-coefficient-size cost.
 *    For modest degrees (the typical case in this codebase) the
 *    algorithm runs in a small multiple of the input size.
 *
 * 5. Pseudo-division.  Implements the textbook recipe:
 *      lc(b)^(deg(a)-deg(b)+1) * a = q*b + r,  deg(r) < deg(b).
 *    Required for the subresultant chain (where rational coefficients
 *    are unacceptable) and for non-monic exact division.
 */

#include <stdio.h>     /* must precede zupoly.h's <gmp.h> for gmp_fprintf */
#include <stdlib.h>
#include <string.h>

#include "zupoly.h"
#include "expr.h"
#include "expand.h"
#include "eval.h"
#include "poly.h"

/* ====================================================================== */
/*  Internal helpers                                                      */
/* ====================================================================== */

/* Grow `p->cap` to at least `new_cap`, mpz_init'ing the new slots to 0.
 * Existing coefficients survive unchanged. */
static void zupoly_reserve(ZUPoly* p, int new_cap) {
    if (new_cap <= p->cap) return;
    int target = (p->cap > 0) ? p->cap : 1;
    while (target < new_cap) target *= 2;
    p->c = (mpz_t*)realloc(p->c, sizeof(mpz_t) * (size_t)target);
    for (int i = p->cap; i < target; i++) {
        mpz_init(p->c[i]);
    }
    p->cap = target;
}

/* ====================================================================== */
/*  Construction / destruction                                            */
/* ====================================================================== */

ZUPoly* zupoly_new(int cap) {
    if (cap < 1) cap = 1;
    ZUPoly* p = (ZUPoly*)malloc(sizeof(ZUPoly));
    p->deg = -1;
    p->cap = cap;
    p->c = (mpz_t*)malloc(sizeof(mpz_t) * (size_t)cap);
    for (int i = 0; i < cap; i++) {
        mpz_init(p->c[i]);
    }
    return p;
}

ZUPoly* zupoly_zero(void) { return zupoly_new(1); }

ZUPoly* zupoly_from_int(int64_t n) {
    ZUPoly* p = zupoly_new(1);
    if (n != 0) {
        mpz_set_si(p->c[0], (long)n);
        p->deg = 0;
    }
    return p;
}

ZUPoly* zupoly_copy(const ZUPoly* p) {
    int cap = (p->deg < 0) ? 1 : p->deg + 1;
    ZUPoly* q = zupoly_new(cap);
    if (p->deg >= 0) {
        for (int i = 0; i <= p->deg; i++) mpz_set(q->c[i], p->c[i]);
        q->deg = p->deg;
    }
    return q;
}

void zupoly_free(ZUPoly* p) {
    if (!p) return;
    for (int i = 0; i < p->cap; i++) mpz_clear(p->c[i]);
    free(p->c);
    free(p);
}

/* ====================================================================== */
/*  Coefficient access                                                    */
/* ====================================================================== */

void zupoly_setcoef(ZUPoly* p, int i, const mpz_t v) {
    if (i < 0) return;
    zupoly_reserve(p, i + 1);
    mpz_set(p->c[i], v);
    /* Update deg.  If we wrote a non-zero at i > deg, deg becomes i.
     * If we wrote zero at the leading position, walk down. */
    if (mpz_sgn(v) != 0) {
        if (i > p->deg) p->deg = i;
    } else if (i == p->deg) {
        while (p->deg >= 0 && mpz_sgn(p->c[p->deg]) == 0) p->deg--;
    }
}

void zupoly_setcoef_si(ZUPoly* p, int i, int64_t v) {
    mpz_t tmp; mpz_init(tmp); mpz_set_si(tmp, (long)v);
    zupoly_setcoef(p, i, tmp);
    mpz_clear(tmp);
}

const mpz_t* zupoly_getcoef(const ZUPoly* p, int i) {
    if (i < 0 || i > p->deg) return NULL;
    return (const mpz_t*)&p->c[i];
}

void zupoly_normalize(ZUPoly* p) {
    int d = p->cap - 1;
    while (d >= 0 && mpz_sgn(p->c[d]) == 0) d--;
    p->deg = d;
}

/* ====================================================================== */
/*  Predicates                                                            */
/* ====================================================================== */

bool zupoly_is_zero(const ZUPoly* p) { return p->deg < 0; }

bool zupoly_eq(const ZUPoly* a, const ZUPoly* b) {
    if (a->deg != b->deg) return false;
    for (int i = 0; i <= a->deg; i++) {
        if (mpz_cmp(a->c[i], b->c[i]) != 0) return false;
    }
    return true;
}

int zupoly_cmp(const ZUPoly* a, const ZUPoly* b) {
    if (a->deg != b->deg) return (a->deg < b->deg) ? -1 : 1;
    for (int i = a->deg; i >= 0; i--) {
        int c = mpz_cmp(a->c[i], b->c[i]);
        if (c != 0) return c;
    }
    return 0;
}

/* ====================================================================== */
/*  Arithmetic                                                            */
/* ====================================================================== */

ZUPoly* zupoly_add(const ZUPoly* a, const ZUPoly* b) {
    int max_deg = (a->deg > b->deg) ? a->deg : b->deg;
    if (max_deg < 0) return zupoly_zero();
    ZUPoly* r = zupoly_new(max_deg + 1);
    for (int i = 0; i <= max_deg; i++) {
        if (i <= a->deg) mpz_set(r->c[i], a->c[i]);
        if (i <= b->deg) mpz_add(r->c[i], r->c[i], b->c[i]);
    }
    zupoly_normalize(r);
    return r;
}

ZUPoly* zupoly_sub(const ZUPoly* a, const ZUPoly* b) {
    int max_deg = (a->deg > b->deg) ? a->deg : b->deg;
    if (max_deg < 0) return zupoly_zero();
    ZUPoly* r = zupoly_new(max_deg + 1);
    for (int i = 0; i <= max_deg; i++) {
        if (i <= a->deg) mpz_set(r->c[i], a->c[i]);
        if (i <= b->deg) mpz_sub(r->c[i], r->c[i], b->c[i]);
    }
    zupoly_normalize(r);
    return r;
}

ZUPoly* zupoly_mul(const ZUPoly* a, const ZUPoly* b) {
    if (a->deg < 0 || b->deg < 0) return zupoly_zero();
    int rd = a->deg + b->deg;
    ZUPoly* r = zupoly_new(rd + 1);
    mpz_t tmp; mpz_init(tmp);
    for (int i = 0; i <= a->deg; i++) {
        if (mpz_sgn(a->c[i]) == 0) continue;
        for (int j = 0; j <= b->deg; j++) {
            mpz_mul(tmp, a->c[i], b->c[j]);
            mpz_add(r->c[i + j], r->c[i + j], tmp);
        }
    }
    mpz_clear(tmp);
    zupoly_normalize(r);
    return r;
}

ZUPoly* zupoly_neg(const ZUPoly* a) {
    ZUPoly* r = zupoly_copy(a);
    for (int i = 0; i <= r->deg; i++) mpz_neg(r->c[i], r->c[i]);
    return r;
}

ZUPoly* zupoly_scale(const ZUPoly* a, const mpz_t s) {
    if (mpz_sgn(s) == 0 || a->deg < 0) return zupoly_zero();
    ZUPoly* r = zupoly_new(a->deg + 1);
    for (int i = 0; i <= a->deg; i++) mpz_mul(r->c[i], a->c[i], s);
    r->deg = a->deg;
    return r;
}

ZUPoly* zupoly_scale_si(const ZUPoly* a, int64_t s) {
    mpz_t tmp; mpz_init_set_si(tmp, (long)s);
    ZUPoly* r = zupoly_scale(a, tmp);
    mpz_clear(tmp);
    return r;
}

bool zupoly_divrem_monic(const ZUPoly* a, const ZUPoly* b,
                         ZUPoly** q_out, ZUPoly** r_out) {
    if (b->deg < 0) return false;     /* divide by zero */
    /* General case: we don't require b to be monic, but the
     * coefficient arithmetic stays in Z only when the leading
     * coefficient of the running remainder is divisible by lc(b)
     * at every step.  For non-monic b that fails in general --
     * use zupoly_pseudodivrem when arbitrary b is needed.
     *
     * Implementation: long division.  Quotient has degree
     * deg(a) - deg(b) at most. */
    if (a->deg < b->deg) {
        *q_out = zupoly_zero();
        *r_out = zupoly_copy(a);
        return true;
    }
    int q_deg = a->deg - b->deg;
    ZUPoly* q = zupoly_new(q_deg + 1);
    ZUPoly* r = zupoly_copy(a);
    mpz_t coef, tmp; mpz_init(coef); mpz_init(tmp);
    for (int i = q_deg; i >= 0; i--) {
        int rd = b->deg + i;
        if (rd > r->deg || mpz_sgn(r->c[rd]) == 0) continue;
        /* If lc(b) does not divide r->c[rd], the integer division
         * truncates and the algorithm produces an incorrect quotient.
         * In monic case (lc(b) == 1) divexact is trivially correct. */
        if (!mpz_divisible_p(r->c[rd], b->c[b->deg])) {
            mpz_clear(coef); mpz_clear(tmp);
            zupoly_free(q); zupoly_free(r);
            *q_out = NULL; *r_out = NULL;
            return false;
        }
        mpz_divexact(coef, r->c[rd], b->c[b->deg]);
        mpz_set(q->c[i], coef);
        for (int j = 0; j <= b->deg; j++) {
            mpz_mul(tmp, coef, b->c[j]);
            mpz_sub(r->c[j + i], r->c[j + i], tmp);
        }
    }
    mpz_clear(coef); mpz_clear(tmp);
    zupoly_normalize(q);
    zupoly_normalize(r);
    *q_out = q;
    *r_out = r;
    return true;
}

bool zupoly_pseudodivrem(const ZUPoly* a, const ZUPoly* b,
                         ZUPoly** q_out, ZUPoly** r_out) {
    if (b->deg < 0) return false;
    if (a->deg < b->deg) {
        *q_out = zupoly_zero();
        *r_out = zupoly_copy(a);
        return true;
    }
    /* Pre-multiply a by lc(b)^(deg(a) - deg(b) + 1), then run plain
     * long division (which now is exact because all leading
     * coefficients of intermediate remainders are divisible by lc(b)
     * thanks to the prefactor). */
    int e = a->deg - b->deg + 1;
    mpz_t pf; mpz_init(pf);
    mpz_pow_ui(pf, b->c[b->deg], (unsigned long)e);
    ZUPoly* a_scaled = zupoly_scale(a, pf);
    mpz_clear(pf);

    int q_deg = a->deg - b->deg;
    ZUPoly* q = zupoly_new(q_deg + 1);
    ZUPoly* r = a_scaled;  /* take ownership */
    mpz_t coef, tmp; mpz_init(coef); mpz_init(tmp);
    for (int i = q_deg; i >= 0; i--) {
        int rd = b->deg + i;
        if (rd > r->deg || mpz_sgn(r->c[rd]) == 0) continue;
        /* By construction lc(b) divides r->c[rd]; assert this would be
         * a programming-error guard.  In production we just trust the
         * recipe. */
        mpz_divexact(coef, r->c[rd], b->c[b->deg]);
        mpz_set(q->c[i], coef);
        for (int j = 0; j <= b->deg; j++) {
            mpz_mul(tmp, coef, b->c[j]);
            mpz_sub(r->c[j + i], r->c[j + i], tmp);
        }
    }
    mpz_clear(coef); mpz_clear(tmp);
    zupoly_normalize(q);
    zupoly_normalize(r);
    *q_out = q;
    *r_out = r;
    return true;
}

ZUPoly* zupoly_divexact(const ZUPoly* a, const ZUPoly* b) {
    if (b->deg < 0) return NULL;
    /* Trivial cases. */
    if (a->deg < 0) return zupoly_zero();
    if (a->deg < b->deg) return NULL;  /* not exact */
    ZUPoly *q, *r;
    if (!zupoly_divrem_monic(a, b, &q, &r)) {
        /* lc(b) didn't divide some intermediate; not exact over Z. */
        return NULL;
    }
    if (!zupoly_is_zero(r)) {
        zupoly_free(q); zupoly_free(r);
        return NULL;
    }
    zupoly_free(r);
    return q;
}

/* Subresultant pseudo-remainder sequence -- Brown-Collins.  Returns
 * an associate of gcd(a, b) over Z.  The returned polynomial is
 * primitive with positive leading coefficient.
 *
 * Reference: Knuth TAOCP vol 2 §4.6.1, or Geddes-Czapor-Labahn ch. 7. */
ZUPoly* zupoly_gcd(const ZUPoly* a, const ZUPoly* b) {
    /* Edge cases. */
    if (zupoly_is_zero(a) && zupoly_is_zero(b)) return zupoly_zero();
    if (zupoly_is_zero(a)) return zupoly_primitive_part(b);
    if (zupoly_is_zero(b)) return zupoly_primitive_part(a);

    /* Order so deg(F) >= deg(G). */
    ZUPoly* F;
    ZUPoly* G;
    if (a->deg >= b->deg) { F = zupoly_copy(a); G = zupoly_copy(b); }
    else                  { F = zupoly_copy(b); G = zupoly_copy(a); }

    /* Strip contents up front so subsequent psr stays small. */
    mpz_t cF, cG, dGCD;
    mpz_init(cF); mpz_init(cG); mpz_init(dGCD);
    zupoly_content(F, cF);
    zupoly_content(G, cG);
    mpz_gcd(dGCD, cF, cG);
    {
        ZUPoly* tmp = zupoly_primitive_part(F); zupoly_free(F); F = tmp;
        tmp = zupoly_primitive_part(G); zupoly_free(G); G = tmp;
    }

    mpz_t g, h; mpz_init_set_si(g, 1); mpz_init_set_si(h, 1);
    /* Subresultant PRS loop. */
    while (!zupoly_is_zero(G)) {
        int delta = F->deg - G->deg;
        ZUPoly *q, *r;
        if (!zupoly_pseudodivrem(F, G, &q, &r)) {
            /* Should never happen for non-zero G. */
            zupoly_free(F); zupoly_free(G); zupoly_free(q); zupoly_free(r);
            mpz_clear(cF); mpz_clear(cG); mpz_clear(dGCD);
            mpz_clear(g); mpz_clear(h);
            return zupoly_zero();
        }
        zupoly_free(q);
        /* Replace F <- G, then G <- r / (g * h^delta). */
        zupoly_free(F); F = G;
        if (zupoly_is_zero(r)) { G = r; break; }

        /* Divisor: g * h^delta. */
        mpz_t denom, hpow; mpz_init(denom); mpz_init(hpow);
        mpz_pow_ui(hpow, h, (unsigned long)delta);
        mpz_mul(denom, g, hpow);
        ZUPoly* r_red = zupoly_new(r->deg + 1);
        for (int i = 0; i <= r->deg; i++) {
            mpz_divexact(r_red->c[i], r->c[i], denom);
        }
        r_red->deg = r->deg;
        zupoly_free(r);
        G = r_red;

        /* Update g and h.
         *   g_new = lc(F)
         *   h_new = h^(1-delta) * g_new^delta
         * computed as h * (g_new / h)^delta when delta >= 1. */
        mpz_set(g, F->c[F->deg]);
        if (delta == 0) {
            /* h unchanged */
        } else if (delta == 1) {
            mpz_set(h, g);
        } else {
            mpz_t gh; mpz_init(gh);
            mpz_pow_ui(gh, g, (unsigned long)delta);
            mpz_t hd; mpz_init(hd);
            mpz_pow_ui(hd, h, (unsigned long)(delta - 1));
            mpz_divexact(h, gh, hd);
            mpz_clear(gh); mpz_clear(hd);
        }
        mpz_clear(denom); mpz_clear(hpow);
    }

    /* F is an associate of gcd(a, b).  Normalise to primitive,
     * positive-lead, then multiply by the integer-content gcd. */
    ZUPoly* pp = zupoly_primitive_part(F);
    if (pp->deg >= 0 && mpz_sgn(pp->c[pp->deg]) < 0) {
        for (int i = 0; i <= pp->deg; i++) mpz_neg(pp->c[i], pp->c[i]);
    }
    ZUPoly* result = zupoly_scale(pp, dGCD);
    zupoly_free(pp);

    zupoly_free(F); zupoly_free(G);
    mpz_clear(cF); mpz_clear(cG); mpz_clear(dGCD);
    mpz_clear(g); mpz_clear(h);
    return result;
}

void zupoly_content(const ZUPoly* p, mpz_t out) {
    mpz_set_ui(out, 0);
    for (int i = 0; i <= p->deg; i++) {
        mpz_gcd(out, out, p->c[i]);
    }
    mpz_abs(out, out);
}

ZUPoly* zupoly_primitive_part(const ZUPoly* p) {
    if (zupoly_is_zero(p)) return zupoly_zero();
    mpz_t c; mpz_init(c);
    zupoly_content(p, c);
    if (mpz_cmp_ui(c, 1) == 0) {
        mpz_clear(c);
        return zupoly_copy(p);
    }
    ZUPoly* r = zupoly_new(p->deg + 1);
    for (int i = 0; i <= p->deg; i++) mpz_divexact(r->c[i], p->c[i], c);
    r->deg = p->deg;
    mpz_clear(c);
    return r;
}

/* ====================================================================== */
/*  Substitution / evaluation                                             */
/* ====================================================================== */

void zupoly_eval(const ZUPoly* p, const mpz_t alpha, mpz_t out) {
    /* Horner. */
    mpz_set_ui(out, 0);
    for (int i = p->deg; i >= 0; i--) {
        mpz_mul(out, out, alpha);
        mpz_add(out, out, p->c[i]);
    }
}

void zupoly_eval_si(const ZUPoly* p, int64_t alpha, mpz_t out) {
    mpz_t a; mpz_init_set_si(a, (long)alpha);
    zupoly_eval(p, a, out);
    mpz_clear(a);
}

ZUPoly* zupoly_shift_si(const ZUPoly* p, int64_t alpha) {
    /* p(x + alpha): Horner expansion in (x + a).
     *   q(x) = c_d
     *   for k from d-1 down to 0:
     *       q(x) := q(x) * (x + a) + c_k
     * Each step multiplies the running polynomial by (x + a) in place.
     * Multiplication by (x + a) yields:
     *     new_q[i] = q[i-1] + a * q[i]    for 1 <= i <= deg
     *     new_q[deg+1] = q[deg]           (the new leading term)
     *     new_q[0] = a * q[0]
     * The slot order during the in-place update matters: we must
     * iterate from HIGH to LOW so that q[i-1] is still the old value
     * when we read it.  An earlier version of this routine iterated
     * low-to-high and silently dropped the q[i-1] contribution into
     * coefficient i, producing p(x+a) = p(x+a)/2 patterns. */
    if (p->deg < 0) return zupoly_zero();
    ZUPoly* q = zupoly_new(p->deg + 1);
    mpz_set(q->c[0], p->c[p->deg]);
    q->deg = 0;
    mpz_t a; mpz_init_set_si(a, (long)alpha);
    mpz_t tmp; mpz_init(tmp);
    for (int k = p->deg - 1; k >= 0; k--) {
        zupoly_reserve(q, q->deg + 2);
        /* The new top slot (q->c[q->deg + 1]) starts as zero.  Walk
         * indices from top down so each read sees the original
         * pre-multiply value. */
        for (int i = q->deg + 1; i >= 1; i--) {
            mpz_mul(tmp, q->c[i], a);
            mpz_add(tmp, tmp, q->c[i - 1]);
            mpz_set(q->c[i], tmp);
        }
        mpz_mul(q->c[0], q->c[0], a);
        q->deg += 1;
        mpz_add(q->c[0], q->c[0], p->c[k]);
        zupoly_normalize(q);
    }
    mpz_clear(a); mpz_clear(tmp);
    return q;
}

/* ====================================================================== */
/*  Conversion to / from Mathilda Expr                                     */
/* ====================================================================== */

ZUPoly* expr_to_zupoly(const struct Expr* e, const struct Expr* var) {
    /* Use the existing poly machinery to extract integer coefficients
     * c[0..deg] from `e` viewed as a polynomial in `var`.  Any non-
     * integer coefficient causes a NULL return. */
    Expr* expanded = expr_expand((Expr*)e);
    int deg = get_degree_poly(expanded, (Expr*)var);
    if (deg < 0) {
        expr_free(expanded);
        return zupoly_zero();  /* the expression is the zero polynomial */
    }
    ZUPoly* p = zupoly_new(deg + 1);
    for (int i = 0; i <= deg; i++) {
        Expr* c = get_coeff(expanded, (Expr*)var, i);
        if (c->type == EXPR_INTEGER) {
            mpz_set_si(p->c[i], (long)c->data.integer);
        } else if (c->type == EXPR_BIGINT) {
            mpz_set(p->c[i], c->data.bigint);
        } else if (is_zero_poly(c)) {
            /* leave as 0 */
        } else {
            /* Non-integer coefficient -- bail out. */
            expr_free(c);
            zupoly_free(p);
            expr_free(expanded);
            return NULL;
        }
        expr_free(c);
    }
    zupoly_normalize(p);
    expr_free(expanded);
    return p;
}

struct Expr* zupoly_to_expr(const ZUPoly* p, const struct Expr* var) {
    if (zupoly_is_zero(p)) return expr_new_integer(0);

    /* Build sum from nonzero terms. */
    int n_terms = 0;
    for (int i = 0; i <= p->deg; i++) if (mpz_sgn(p->c[i]) != 0) n_terms++;

    Expr** terms = (Expr**)malloc(sizeof(Expr*) * (size_t)n_terms);
    int idx = 0;
    for (int i = 0; i <= p->deg; i++) {
        if (mpz_sgn(p->c[i]) == 0) continue;
        Expr* coef_expr;
        if (mpz_fits_slong_p(p->c[i])) {
            coef_expr = expr_new_integer((int64_t)mpz_get_si(p->c[i]));
        } else {
            coef_expr = expr_new_bigint_from_mpz(p->c[i]);
        }
        Expr* term;
        if (i == 0) {
            term = coef_expr;
        } else {
            Expr* var_part;
            if (i == 1) var_part = expr_copy((Expr*)var);
            else {
                Expr* args[2] = { expr_copy((Expr*)var), expr_new_integer(i) };
                var_part = expr_new_function(expr_new_symbol("Power"), args, 2);
                var_part = eval_and_free(var_part);
            }
            if (mpz_cmp_ui(p->c[i], 1) == 0) {
                expr_free(coef_expr);
                term = var_part;
            } else if (mpz_cmp_si(p->c[i], -1) == 0) {
                expr_free(coef_expr);
                Expr* args[2] = { expr_new_integer(-1), var_part };
                term = expr_new_function(expr_new_symbol("Times"), args, 2);
                term = eval_and_free(term);
            } else {
                Expr* args[2] = { coef_expr, var_part };
                term = expr_new_function(expr_new_symbol("Times"), args, 2);
                term = eval_and_free(term);
            }
        }
        terms[idx++] = term;
    }
    Expr* result;
    if (n_terms == 1) {
        result = terms[0];
    } else {
        result = expr_new_function(expr_new_symbol("Plus"), terms, n_terms);
        result = eval_and_free(result);
    }
    free(terms);
    return result;
}

/* ====================================================================== */
/*  QUPoly: rational-coefficient univariate polynomial (private)          */
/*                                                                        */
/*  Used solely as the internal substrate for zupoly_diophantine.  We     */
/*  do not expose this type because the public API of zupoly should       */
/*  remain integer-typed; the rational machinery is an implementation     */
/*  detail of the Diophantine solver.                                     */
/*                                                                        */
/*  Storage parallels ZUPoly: dense array of mpq_t, all initialised to    */
/*  zero, deg tracking the highest non-zero coefficient.  Slots beyond    */
/*  deg are valid mpq_t holding zero, available for in-place updates.    */
/* ====================================================================== */

typedef struct {
    int     deg;
    int     cap;
    mpq_t*  c;
} QUPoly;

static QUPoly* qupoly_new(int cap) {
    if (cap < 1) cap = 1;
    QUPoly* p = (QUPoly*)malloc(sizeof(QUPoly));
    p->deg = -1;
    p->cap = cap;
    p->c = (mpq_t*)malloc(sizeof(mpq_t) * (size_t)cap);
    for (int i = 0; i < cap; i++) mpq_init(p->c[i]);
    return p;
}

static void qupoly_free(QUPoly* p) {
    if (!p) return;
    for (int i = 0; i < p->cap; i++) mpq_clear(p->c[i]);
    free(p->c);
    free(p);
}

static void qupoly_normalize(QUPoly* p) {
    int d = p->cap - 1;
    while (d >= 0 && mpq_sgn(p->c[d]) == 0) d--;
    p->deg = d;
}

static QUPoly* qupoly_copy(const QUPoly* a) {
    int cap = (a->deg < 0) ? 1 : a->deg + 1;
    QUPoly* r = qupoly_new(cap);
    for (int i = 0; i <= a->deg; i++) mpq_set(r->c[i], a->c[i]);
    r->deg = a->deg;
    return r;
}

/* Construct from a ZUPoly (each integer coefficient is a rational
 * with denominator 1). */
static QUPoly* qupoly_from_zupoly(const ZUPoly* p) {
    int cap = (p->deg < 0) ? 1 : p->deg + 1;
    QUPoly* r = qupoly_new(cap);
    for (int i = 0; i <= p->deg; i++) {
        mpq_set_z(r->c[i], p->c[i]);
    }
    r->deg = p->deg;
    return r;
}

/* If every coefficient of p is integer-valued (denominator 1 after
 * canonicalisation), allocate and return a ZUPoly copy.  Otherwise
 * return NULL. */
static ZUPoly* qupoly_to_zupoly_if_integer(const QUPoly* p) {
    /* mpq_t coefficients are kept canonical by mpq_canonicalize after
     * arithmetic; we rely on mpz_cmp_ui(denominator, 1) here. */
    int cap = (p->deg < 0) ? 1 : p->deg + 1;
    ZUPoly* r = zupoly_new(cap);
    for (int i = 0; i <= p->deg; i++) {
        mpz_ptr den = mpq_denref(p->c[i]);
        if (mpz_cmp_ui(den, 1) != 0) {
            zupoly_free(r);
            return NULL;
        }
        mpz_set(r->c[i], mpq_numref(p->c[i]));
    }
    r->deg = p->deg;
    return r;
}

/* Arithmetic operations on QUPoly (only what we need for EE). */

static QUPoly* qupoly_sub(const QUPoly* a, const QUPoly* b) {
    int max_deg = (a->deg > b->deg) ? a->deg : b->deg;
    if (max_deg < 0) return qupoly_new(1);
    QUPoly* r = qupoly_new(max_deg + 1);
    for (int i = 0; i <= max_deg; i++) {
        if (i <= a->deg) mpq_set(r->c[i], a->c[i]);
        if (i <= b->deg) mpq_sub(r->c[i], r->c[i], b->c[i]);
    }
    qupoly_normalize(r);
    return r;
}

static QUPoly* qupoly_mul(const QUPoly* a, const QUPoly* b) {
    if (a->deg < 0 || b->deg < 0) return qupoly_new(1);
    int rd = a->deg + b->deg;
    QUPoly* r = qupoly_new(rd + 1);
    mpq_t tmp; mpq_init(tmp);
    for (int i = 0; i <= a->deg; i++) {
        if (mpq_sgn(a->c[i]) == 0) continue;
        for (int j = 0; j <= b->deg; j++) {
            mpq_mul(tmp, a->c[i], b->c[j]);
            mpq_add(r->c[i + j], r->c[i + j], tmp);
        }
    }
    mpq_clear(tmp);
    qupoly_normalize(r);
    return r;
}

/* Long division over Q[x].  Returns true unless b is zero. */
static bool qupoly_divrem(const QUPoly* a, const QUPoly* b,
                          QUPoly** q_out, QUPoly** r_out) {
    if (b->deg < 0) return false;
    if (a->deg < b->deg) {
        *q_out = qupoly_new(1);
        *r_out = qupoly_copy(a);
        return true;
    }
    int q_deg = a->deg - b->deg;
    QUPoly* q = qupoly_new(q_deg + 1);
    QUPoly* r = qupoly_copy(a);
    mpq_t coef, tmp; mpq_init(coef); mpq_init(tmp);
    for (int i = q_deg; i >= 0; i--) {
        int rd = b->deg + i;
        if (rd > r->deg || mpq_sgn(r->c[rd]) == 0) continue;
        mpq_div(coef, r->c[rd], b->c[b->deg]);
        mpq_set(q->c[i], coef);
        for (int j = 0; j <= b->deg; j++) {
            mpq_mul(tmp, coef, b->c[j]);
            mpq_sub(r->c[j + i], r->c[j + i], tmp);
        }
    }
    mpq_clear(coef); mpq_clear(tmp);
    qupoly_normalize(q);
    qupoly_normalize(r);
    *q_out = q;
    *r_out = r;
    return true;
}

static bool qupoly_is_zero(const QUPoly* p) { return p->deg < 0; }

/* Extended Euclidean over Q[x].  On exit:
 *   *g_out = gcd(a, b), monic.
 *   *s_out, *t_out: cofactors with s*a + t*b = g.
 *
 * Always succeeds for non-zero inputs (Q[x] is Euclidean). */
static void qupoly_xgcd(const QUPoly* a, const QUPoly* b,
                        QUPoly** g_out, QUPoly** s_out, QUPoly** t_out) {
    /* Standard cofactor recurrence:
     *   r_0 = a,  r_1 = b
     *   s_0 = 1,  s_1 = 0
     *   t_0 = 0,  t_1 = 1
     *   r_{i+1} = r_{i-1} - q_i * r_i
     *   s_{i+1} = s_{i-1} - q_i * s_i
     *   t_{i+1} = t_{i-1} - q_i * t_i
     * Loop until r_{i+1} = 0; then (r_i, s_i, t_i) is the answer. */
    QUPoly* r0 = qupoly_copy(a);
    QUPoly* r1 = qupoly_copy(b);
    QUPoly* s0 = qupoly_new(1); mpq_set_ui(s0->c[0], 1, 1); s0->deg = 0;
    QUPoly* s1 = qupoly_new(1);
    QUPoly* t0 = qupoly_new(1);
    QUPoly* t1 = qupoly_new(1); mpq_set_ui(t1->c[0], 1, 1); t1->deg = 0;

    while (!qupoly_is_zero(r1)) {
        QUPoly *q, *r2;
        qupoly_divrem(r0, r1, &q, &r2);

        /* s_next = s0 - q * s1; same shape for t. */
        QUPoly* qs1 = qupoly_mul(q, s1);
        QUPoly* s2 = qupoly_sub(s0, qs1);
        qupoly_free(qs1);

        QUPoly* qt1 = qupoly_mul(q, t1);
        QUPoly* t2 = qupoly_sub(t0, qt1);
        qupoly_free(qt1);

        qupoly_free(r0); r0 = r1; r1 = r2;
        qupoly_free(s0); s0 = s1; s1 = s2;
        qupoly_free(t0); t0 = t1; t1 = t2;
        qupoly_free(q);
    }

    /* Normalise: gcd is r0; make it monic by scaling all three. */
    if (!qupoly_is_zero(r0) && mpq_cmp_ui(r0->c[r0->deg], 1, 1) != 0) {
        mpq_t inv; mpq_init(inv);
        mpq_inv(inv, r0->c[r0->deg]);
        for (int i = 0; i <= r0->deg; i++) mpq_mul(r0->c[i], r0->c[i], inv);
        for (int i = 0; i <= s0->deg; i++) mpq_mul(s0->c[i], s0->c[i], inv);
        for (int i = 0; i <= t0->deg; i++) mpq_mul(t0->c[i], t0->c[i], inv);
        mpq_clear(inv);
    }

    qupoly_free(r1); qupoly_free(s1); qupoly_free(t1);
    *g_out = r0;
    *s_out = s0;
    *t_out = t0;
}

/* ====================================================================== */
/*  zupoly_diophantine                                                    */
/* ====================================================================== */

bool zupoly_diophantine(const ZUPoly* u, const ZUPoly* v, const ZUPoly* e,
                        ZUPoly** du_out, ZUPoly** dv_out) {
    *du_out = NULL;
    *dv_out = NULL;

    /* Monic preconditions. */
    if (u->deg < 0 || v->deg < 0) return false;
    if (mpz_cmp_ui(u->c[u->deg], 1) != 0) return false;
    if (mpz_cmp_ui(v->c[v->deg], 1) != 0) return false;

    /* Run extended Euclidean over Q[x] to obtain s, t with
     *    s * u + t * v = gcd(u, v).
     * For coprime monic u, v the gcd is 1. */
    QUPoly* qu = qupoly_from_zupoly(u);
    QUPoly* qv = qupoly_from_zupoly(v);
    QUPoly *qg, *qs, *qt;
    qupoly_xgcd(qu, qv, &qg, &qs, &qt);

    /* Coprime check: gcd must be the unit constant 1. */
    if (qg->deg != 0 || mpq_cmp_ui(qg->c[0], 1, 1) != 0) {
        qupoly_free(qu); qupoly_free(qv);
        qupoly_free(qg); qupoly_free(qs); qupoly_free(qt);
        return false;
    }
    qupoly_free(qg);

    /* delta_u = (t * e) mod u in Q[x]. */
    QUPoly* qe = qupoly_from_zupoly(e);
    QUPoly* qte = qupoly_mul(qt, qe);
    QUPoly *q_quot, *qdu;
    qupoly_divrem(qte, qu, &q_quot, &qdu);
    qupoly_free(q_quot);
    qupoly_free(qte);

    /* delta_v = (e - delta_u * v) / u  (exact in Q[x]). */
    QUPoly* qduv = qupoly_mul(qdu, qv);
    QUPoly* qe_minus = qupoly_sub(qe, qduv);
    qupoly_free(qduv);
    QUPoly* qdv, *qrem;
    qupoly_divrem(qe_minus, qu, &qdv, &qrem);
    /* The remainder should be exactly zero; if it isn't, something
     * went wrong upstream (e.g. e wasn't congruent to delta_u * v
     * mod u, which contradicts the Diophantine derivation). */
    if (!qupoly_is_zero(qrem)) {
        qupoly_free(qrem); qupoly_free(qdv); qupoly_free(qe_minus);
        qupoly_free(qdu); qupoly_free(qe);
        qupoly_free(qu); qupoly_free(qv); qupoly_free(qs); qupoly_free(qt);
        return false;
    }
    qupoly_free(qrem); qupoly_free(qe_minus);

    /* Convert delta_u, delta_v to ZUPoly if both are integer. */
    ZUPoly* du = qupoly_to_zupoly_if_integer(qdu);
    ZUPoly* dv = qupoly_to_zupoly_if_integer(qdv);
    qupoly_free(qdu); qupoly_free(qdv);
    qupoly_free(qe);
    qupoly_free(qu); qupoly_free(qv); qupoly_free(qs); qupoly_free(qt);

    if (!du || !dv) {
        if (du) zupoly_free(du);
        if (dv) zupoly_free(dv);
        return false;
    }
    *du_out = du;
    *dv_out = dv;
    return true;
}

/* ====================================================================== */
/*  Debug                                                                 */
/* ====================================================================== */

void zupoly_print(const ZUPoly* p, const char* var) {
    if (zupoly_is_zero(p)) { fprintf(stderr, "0"); return; }
    bool first = true;
    for (int i = 0; i <= p->deg; i++) {
        if (mpz_sgn(p->c[i]) == 0) continue;
        if (!first) fprintf(stderr, " + ");
        first = false;
        if (i == 0) gmp_fprintf(stderr, "%Zd", p->c[i]);
        else if (mpz_cmp_ui(p->c[i], 1) == 0) {
            if (i == 1) fprintf(stderr, "%s", var);
            else        fprintf(stderr, "%s^%d", var, i);
        } else if (mpz_cmp_si(p->c[i], -1) == 0) {
            if (i == 1) fprintf(stderr, "-%s", var);
            else        fprintf(stderr, "-%s^%d", var, i);
        } else {
            if (i == 1) gmp_fprintf(stderr, "%Zd*%s", p->c[i], var);
            else        gmp_fprintf(stderr, "%Zd*%s^%d", p->c[i], var, i);
        }
    }
}
