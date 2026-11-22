/* radcanon.c -- the canonical form of an exact product of rational radicals.
 *
 * The problem.  A product of a Gaussian-rational coefficient and powers
 * r_i^(e_i) of positive rationals r_i with rational, non-integer exponents
 * has many spellings of one value: -1/3 Sqrt[3/2] = -1/Sqrt[6],
 * 2 Sqrt[3/2] = Sqrt[6], Sqrt[6]/3 = Sqrt[2/3], 6^(1/3)/2 = 3^(1/3)/2^(2/3),
 * (2/3)^(5/3) = 2/3 (2/3)^(2/3).  The evaluator needs one of them, chosen
 * by an algorithm rather than by a list of shapes, so that every product
 * of this kind lands on it.
 *
 * The form.  Write the value over pairwise-coprime integers b_k > 1:
 *     c * prod_k b_k^(t_k),   t_k rational.
 * Then
 *   1. each t_k splits as k + f with k = trunc(t_k) (towards zero) and
 *      0 < |f| < 1: b^k joins the rational coefficient, b^f stays a radical
 *      (2^(3/2) = 2 Sqrt[2], 2^(-3/2) = 1/(2 Sqrt[2]): the radical keeps
 *      the sign of t, so an integer is never moved across the fraction bar);
 *   2. radicals with the same |f| merge into one power: those with f > 0
 *      form the numerator P, those with f < 0 the denominator N, giving
 *      P^|f| when N = 1, N^-|f| when P = 1, and (P/N)^|f| otherwise
 *      (Sqrt[2] Sqrt[3] = Sqrt[6], 3^(1/2) 2^(-1/2) = Sqrt[3/2],
 *      2^(-1/2) 3^(-1/2) = 1/Sqrt[6]).  Radicals with different |f| stay
 *      apart (2^(1/3) 3^(2/3) is not merged into 18^(1/3)).
 * Every step is an identity for positive bases (b^(k+f) = b^k b^f, and
 * x^f y^f = (x y)^f for x, y > 0), so the form is exact.  It is unique
 * given the basis, and the basis is unique too: it is the gcd-free basis
 * (factor refinement) of all the integers involved, which needs gcds only,
 * never integer factorisation, so a 60-digit coefficient against a
 * 60-digit radicand costs a handful of gcds.  Merging by |f| then recombines
 * basis elements that the refinement separated needlessly, which is why
 * the result does not depend on how far the refinement went: two coprime
 * pieces with equal exponents give the same power whether or not they were
 * split.  Power applies the same algorithm to one radical after factoring
 * its base, which is what makes 12^(1/3) = 2^(2/3) 3^(1/3); within a
 * product the inputs are already Power's output, so their bases are
 * square-free in that sense and no factoring is needed.
 *
 * A complex coefficient takes part through its content: z = c u with
 * c = gcd(numerators)/lcm(denominators) > 0 and u a primitive Gaussian
 * integer, and only c meets the radicals ((3 + 6 I)/Sqrt[3] =
 * (1 + 2 I) Sqrt[3], (1 + 3 I)/2 Sqrt[2] = (1 + 3 I)/Sqrt[2]).
 *
 * A negative base.  (-n)^e = (-1)^e n^e on the principal branch for every
 * n > 0, so one negative power (-n)^e with 0 < e < 1 joins the basis
 * through n^e, and its sign (-1)^e goes back onto the merged radical whose
 * exponent magnitude is e: (-a/b)^e, with that group's P/N as a/b, is
 * (-1)^e (a/b)^e because a/b > 0.  When no group has |f| = e the sign
 * stays as (-1)^e: (-6)^(1/3)/2 = (-3)^(1/3)/2^(2/3),
 * (-2)^(1/3) 3^(1/3) = (-6)^(1/3), (-2)^(1/3)/4 = (-1)^(1/3)/(2 2^(2/3)).
 * Two negative powers never take part together (their signs would have to
 * meet: (-2)^(1/3) (-3)^(1/3) is not (-6)^(1/3) times anything simple),
 * which the callers enforce.
 *
 * This is the reference's automatic form on all of the ~180 products
 * checked, with one exception: the reference leaves
 * 2^(-2/3) 3^(-1/3) 5^(2/3) unmerged while merging 2^(-2/3) 5^(2/3) on its
 * own, which is not a rule but an artefact of its pairwise search; here
 * the two 2/3 powers merge, (5/2)^(2/3)/3^(1/3). */
#include "radcanon.h"
#include "sym_names.h"
#include "arithmetic.h"   /* make_complex */
#include <stdlib.h>
#include <string.h>

void rc_atoms_init(RcAtoms* at) {
    at->a = NULL;
    at->n = at->cap = 0;
}

void rc_atoms_clear(RcAtoms* at) {
    for (size_t i = 0; i < at->n; i++) {
        mpz_clear(at->a[i].v);
        mpq_clear(at->a[i].e);
    }
    free(at->a);
    at->a = NULL;
    at->n = at->cap = 0;
}

void rc_atoms_add(RcAtoms* at, const mpz_t v, const mpq_t e) {
    if (mpz_cmp_ui(v, 1) <= 0 || mpq_sgn(e) == 0) return;
    if (at->n == at->cap) {
        at->cap = at->cap ? 2 * at->cap : 8;
        at->a = realloc(at->a, at->cap * sizeof(RcAtom));
    }
    mpz_init_set(at->a[at->n].v, v);
    mpq_init(at->a[at->n].e);
    mpq_set(at->a[at->n].e, e);
    at->n++;
}

static Expr* rc_int(const mpz_t z) {
    return expr_bigint_normalize(expr_new_bigint_from_mpz(z));
}

/* num/den as an Integer or a Rational (den > 0, gcd 1). */
static Expr* rc_ratio(const mpz_t num, const mpz_t den) {
    if (mpz_cmp_ui(den, 1) == 0) return rc_int(num);
    return expr_new_function(expr_new_symbol(SYM_Rational),
                             (Expr*[]){ rc_int(num), rc_int(den) }, 2);
}

static Expr* rc_mpq(const mpq_t q) {
    return rc_ratio(mpq_numref(q), mpq_denref(q));
}

/* The gcd-free basis of the values: pairwise coprime, each value a product
 * of powers of basis elements.  Refines any two elements sharing a factor
 * g into g, x/g, y/g until no two share one (the standard factor
 * refinement); terminates because the product of the elements never grows
 * and each step either shrinks it or removes a duplicate. */
static size_t rc_basis(const RcAtoms* at, mpz_t** out) {
    size_t nb = 0, cap = at->n * 2 + 4;
    mpz_t* B = malloc(cap * sizeof(mpz_t));
    for (size_t i = 0; i < at->n; i++) {
        bool dup = false;
        for (size_t k = 0; k < nb && !dup; k++) dup = mpz_cmp(B[k], at->a[i].v) == 0;
        if (!dup) mpz_init_set(B[nb++], at->a[i].v);
    }
    mpz_t g, x, y;
    mpz_inits(g, x, y, NULL);
    bool again = true;
    while (again) {
        again = false;
        for (size_t i = 0; i < nb && !again; i++)
            for (size_t j = i + 1; j < nb && !again; j++) {
                mpz_gcd(g, B[i], B[j]);
                if (mpz_cmp_ui(g, 1) == 0) continue;
                /* strip every power of g at once (mpz_remove), so a value
                 * like 10^100000 against 2 refines in a few steps, not
                 * 100000; what still shares a smaller factor with g is
                 * refined on a later pass */
                mpz_remove(x, B[i], g);
                mpz_remove(y, B[j], g);
                mpz_set(B[i], g);
                mpz_set(B[j], x);
                if (nb == cap) { cap *= 2; B = realloc(B, cap * sizeof(mpz_t)); }
                mpz_init_set(B[nb++], y);
                /* drop 1s and duplicates */
                size_t w = 0;
                for (size_t r = 0; r < nb; r++) {
                    bool drop = mpz_cmp_ui(B[r], 1) == 0;
                    for (size_t q = 0; q < w && !drop; q++) drop = mpz_cmp(B[q], B[r]) == 0;
                    if (drop) { mpz_clear(B[r]); continue; }
                    if (w != r) { mpz_init_set(B[w], B[r]); mpz_clear(B[r]); }
                    w++;
                }
                nb = w;
                again = true;
            }
    }
    mpz_clears(g, x, y, NULL);
    *out = B;
    return nb;
}

/* One radical b^f of the result, before merging. */
typedef struct {
    mpz_t b;
    mpq_t f;      /* 0 < |f| < 1 */
    mpq_t af;     /* |f| */
} RcRad;

static int rc_rad_cmp(const void* x, const void* y) {
    const RcRad* a = x;
    const RcRad* b = y;
    int c = mpq_cmp(a->af, b->af);
    if (c) return c;
    return mpz_cmp(a->b, b->b);
}

size_t rc_canonical(const RcAtoms* at, const mpq_t coef_re, const mpq_t coef_im,
                    const mpq_t* neg_e, Expr*** out) {
    /* the content c > 0 of the coefficient, and its unit u = coef / c */
    mpz_t cn, cd, t;
    mpz_inits(cn, cd, t, NULL);
    mpq_t ure, uim;
    mpq_inits(ure, uim, NULL);
    if (mpq_sgn(coef_im) == 0) {
        mpz_abs(cn, mpq_numref(coef_re));
        mpz_set(cd, mpq_denref(coef_re));
        mpq_set_si(ure, mpq_sgn(coef_re), 1);
    } else {
        mpz_gcd(cn, mpq_numref(coef_re), mpq_numref(coef_im));   /* gcd(0, b) = |b| */
        mpz_lcm(cd, mpq_denref(coef_re), mpq_denref(coef_im));
        mpq_t c;
        mpq_init(c);
        mpq_set_num(c, cn);
        mpq_set_den(c, cd);
        mpq_canonicalize(c);
        mpq_div(ure, coef_re, c);
        mpq_div(uim, coef_im, c);
        mpq_clear(c);
    }

    /* atoms of the content and the given powers, over one coprime basis */
    RcAtoms all;
    rc_atoms_init(&all);
    {
        mpq_t one, mone;
        mpq_inits(one, mone, NULL);
        mpq_set_si(one, 1, 1);
        mpq_set_si(mone, -1, 1);
        rc_atoms_add(&all, cn, one);
        rc_atoms_add(&all, cd, mone);
        mpq_clears(one, mone, NULL);
    }
    for (size_t i = 0; i < at->n; i++) rc_atoms_add(&all, at->a[i].v, at->a[i].e);

    mpz_t* B;
    size_t nb = rc_basis(&all, &B);
    mpq_t* T = malloc((nb ? nb : 1) * sizeof(mpq_t));
    for (size_t k = 0; k < nb; k++) mpq_init(T[k]);
    {
        mpz_t r;
        mpz_init(r);
        mpq_t c;
        mpq_init(c);
        for (size_t i = 0; i < all.n; i++) {
            mpz_set(r, all.a[i].v);
            for (size_t k = 0; k < nb && mpz_cmp_ui(r, 1) != 0; k++) {
                unsigned long m = mpz_remove(r, r, B[k]);   /* multiplicity */
                if (m) {
                    mpq_set_ui(c, m, 1);
                    mpq_mul(c, c, all.a[i].e);
                    mpq_add(T[k], T[k], c);
                }
            }
        }
        mpz_clear(r);
        mpq_clear(c);
    }

    /* integer parts into the coefficient, fractional parts as radicals */
    mpz_set_ui(cn, 1);
    mpz_set_ui(cd, 1);
    RcRad* R = malloc((nb ? nb : 1) * sizeof(RcRad));
    size_t nr = 0;
    for (size_t k = 0; k < nb; k++) {
        if (mpq_sgn(T[k]) == 0) continue;
        mpz_tdiv_q(t, mpq_numref(T[k]), mpq_denref(T[k]));     /* trunc(T) */
        if (mpz_sgn(t) != 0) {
            /* |trunc(T)| fits a machine word: it is bounded by the atoms'
             * multiplicities times their exponents */
            mpz_t p, at_;
            mpz_inits(p, at_, NULL);
            mpz_abs(at_, t);
            mpz_pow_ui(p, B[k], mpz_get_ui(at_));
            mpz_clear(at_);
            if (mpz_sgn(t) > 0) mpz_mul(cn, cn, p);
            else mpz_mul(cd, cd, p);
            mpz_clear(p);
        }
        mpq_t f;
        mpq_init(f);
        mpq_set_z(f, t);
        mpq_sub(f, T[k], f);
        if (mpq_sgn(f) != 0) {
            mpz_init_set(R[nr].b, B[k]);
            mpq_init(R[nr].f);
            mpq_set(R[nr].f, f);
            mpq_init(R[nr].af);
            mpq_abs(R[nr].af, f);
            nr++;
        }
        mpq_clear(f);
    }
    if (nr > 1) qsort(R, nr, sizeof(RcRad), rc_rad_cmp);

    size_t cap = nr + 3, no = 0;
    Expr** F = malloc(cap * sizeof(Expr*));

    /* the number: (cn/cd) u */
    {
        mpq_t c, re, im;
        mpq_inits(c, re, im, NULL);
        mpq_set_num(c, cn);
        mpq_set_den(c, cd);
        mpq_canonicalize(c);
        mpq_mul(re, c, ure);
        mpq_mul(im, c, uim);
        if (mpq_sgn(im) != 0)
            F[no++] = make_complex(rc_mpq(re), rc_mpq(im));
        else if (!(mpz_cmp_ui(mpq_numref(re), 1) == 0 && mpz_cmp_ui(mpq_denref(re), 1) == 0))
            F[no++] = rc_mpq(re);
        mpq_clears(c, re, im, NULL);
    }

    /* the merged radicals, one per |f| */
    bool sign_placed = neg_e == NULL;
    mpz_t P, N;
    mpz_inits(P, N, NULL);
    for (size_t i = 0; i < nr; ) {
        size_t j = i;
        mpz_set_ui(P, 1);
        mpz_set_ui(N, 1);
        while (j < nr && mpq_equal(R[j].af, R[i].af)) {
            if (mpq_sgn(R[j].f) > 0) mpz_mul(P, P, R[j].b);
            else mpz_mul(N, N, R[j].b);
            j++;
        }
        Expr* base;
        Expr* ex;
        if (!sign_placed && mpq_equal(R[i].af, *neg_e)) {
            /* (-1)^e (P/N)^e = (-P/N)^e */
            mpz_neg(P, P);
            base = rc_ratio(P, N);
            ex = rc_mpq(R[i].af);
            sign_placed = true;
        } else if (mpz_cmp_ui(N, 1) == 0) {
            base = rc_int(P);
            ex = rc_mpq(R[i].af);
        } else if (mpz_cmp_ui(P, 1) == 0) {
            mpq_t nf;
            mpq_init(nf);
            mpq_neg(nf, R[i].af);
            base = rc_int(N);
            ex = rc_mpq(nf);
            mpq_clear(nf);
        } else {
            base = rc_ratio(P, N);
            ex = rc_mpq(R[i].af);
        }
        if (no == cap) { cap *= 2; F = realloc(F, cap * sizeof(Expr*)); }
        F[no++] = expr_new_function(expr_new_symbol(SYM_Power), (Expr*[]){ base, ex }, 2);
        i = j;
    }
    if (!sign_placed) {
        if (no == cap) { cap *= 2; F = realloc(F, cap * sizeof(Expr*)); }
        F[no++] = expr_new_function(expr_new_symbol(SYM_Power),
                                    (Expr*[]){ expr_new_integer(-1), rc_mpq(*neg_e) }, 2);
    }

    mpz_clears(P, N, NULL);
    for (size_t i = 0; i < nr; i++) {
        mpz_clear(R[i].b);
        mpq_clear(R[i].f);
        mpq_clear(R[i].af);
    }
    free(R);
    for (size_t k = 0; k < nb; k++) { mpz_clear(B[k]); mpq_clear(T[k]); }
    free(B);
    free(T);
    rc_atoms_clear(&all);
    mpq_clears(ure, uim, NULL);
    mpz_clears(cn, cd, t, NULL);
    *out = F;
    return no;
}
