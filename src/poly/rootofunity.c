/* rootofunity.c — cyclotomic-polynomial construction + root-of-unity
 * recognition for the Q(α) extension substrate.  See rootofunity.h. */

#include "rootofunity.h"

#include "expr.h"
#include "arithmetic.h"   /* is_rational */
#include "sym_names.h"    /* SYM_Power, SYM_Rational, SYM_Complex, SYM_I */

#include <gmp.h>
#include <stdlib.h>

/* ====================================================================== */
/*  Integer polynomials (dense, coef[i] = coefficient of y^i)             */
/* ====================================================================== */

typedef struct { size_t deg; mpz_t* c; } IPoly;

static IPoly ipoly_alloc(size_t deg) {
    IPoly p;
    p.deg = deg;
    p.c = (mpz_t*)malloc(sizeof(mpz_t) * (deg + 1));
    for (size_t i = 0; i <= deg; i++) mpz_init(p.c[i]);
    return p;
}

static void ipoly_free(IPoly* p) {
    if (!p->c) return;
    for (size_t i = 0; i <= p->deg; i++) mpz_clear(p->c[i]);
    free(p->c);
    p->c = NULL;
}

static IPoly ipoly_copy(const IPoly* a) {
    IPoly r = ipoly_alloc(a->deg);
    for (size_t i = 0; i <= a->deg; i++) mpz_set(r.c[i], a->c[i]);
    return r;
}

/* Exact division of A by a MONIC divisor B (B->c[B->deg] == 1).  The
 * quotient is integer and exact by construction for the cyclotomic
 * recursion; any remainder is discarded (it is zero). */
static IPoly ipoly_div_monic(const IPoly* A, const IPoly* B) {
    size_t a = A->deg, b = B->deg;
    IPoly R = ipoly_copy(A);
    size_t qd = (a >= b) ? (a - b) : 0;
    IPoly Q = ipoly_alloc(qd);
    mpz_t t;
    mpz_init(t);
    if (a >= b) {
        for (size_t i = a - b + 1; i-- > 0; ) {
            mpz_set(Q.c[i], R.c[i + b]);          /* B monic ⇒ divide by 1 */
            for (size_t j = 0; j <= b; j++) {
                mpz_mul(t, Q.c[i], B->c[j]);
                mpz_sub(R.c[i + j], R.c[i + j], t);
            }
        }
    }
    mpz_clear(t);
    ipoly_free(&R);
    return Q;
}

/* Φ_n(y) by the exact divisor recursion.  Returns a fresh IPoly (caller
 * frees with ipoly_free).  n ≥ 1. */
static IPoly cyclotomic_poly(unsigned long n) {
    if (n == 1) {
        IPoly r = ipoly_alloc(1);          /* y − 1 */
        mpz_set_si(r.c[0], -1);
        mpz_set_si(r.c[1], 1);
        return r;
    }
    IPoly num = ipoly_alloc(n);            /* y^n − 1 */
    mpz_set_si(num.c[0], -1);
    mpz_set_si(num.c[n], 1);
    for (unsigned long d = 1; d < n; d++) {
        if (n % d) continue;
        IPoly phid = cyclotomic_poly(d);
        IPoly q = ipoly_div_monic(&num, &phid);
        ipoly_free(&phid);
        ipoly_free(&num);
        num = q;
    }
    return num;                            /* now equals Φ_n */
}

QAExt* qaext_cyclotomic(unsigned long n) {
    if (n == 0) return NULL;
    IPoly phi = cyclotomic_poly(n);
    QAExt* ext = qaext_new(phi.deg);
    mpq_t v;
    mpq_init(v);
    for (size_t i = 0; i <= phi.deg; i++) {
        mpq_set_z(v, phi.c[i]);            /* integer coefficient → Q */
        mpq_canonicalize(v);
        qaext_set_coef_mpq(ext, i, v);
    }
    mpq_clear(v);
    ipoly_free(&phi);
    return ext;
}

/* ====================================================================== */
/*  Root-of-unity recognition                                             */
/* ====================================================================== */

bool expr_is_root_of_unity_pow(const Expr* e, int64_t* p_out, int64_t* q_out) {
    if (!e) return false;

    /* I  →  (-1)^(1/2);   Complex[0, ±1]. */
    if (e->type == EXPR_SYMBOL && e->data.symbol == SYM_I) {
        *p_out = 1; *q_out = 2; return true;
    }
    if (e->type == EXPR_FUNCTION
        && e->data.function.head
        && e->data.function.head->type == EXPR_SYMBOL
        && e->data.function.head->data.symbol == SYM_Complex
        && e->data.function.arg_count == 2) {
        Expr* re = e->data.function.args[0];
        Expr* im = e->data.function.args[1];
        if (re->type == EXPR_INTEGER && re->data.integer == 0
            && im->type == EXPR_INTEGER) {
            if (im->data.integer == 1)  { *p_out = 1; *q_out = 2; return true; }
            if (im->data.integer == -1) { *p_out = 3; *q_out = 2; return true; }
        }
        return false;
    }

    /* Power[-1, p/q]  (the canonical form for every other root of unity). */
    if (e->type == EXPR_FUNCTION
        && e->data.function.head
        && e->data.function.head->type == EXPR_SYMBOL
        && e->data.function.head->data.symbol == SYM_Power
        && e->data.function.arg_count == 2) {
        Expr* base = e->data.function.args[0];
        Expr* exp  = e->data.function.args[1];
        if (base->type == EXPR_INTEGER && base->data.integer == -1) {
            int64_t p, q;
            if (is_rational(exp, &p, &q) && q >= 1) {
                *p_out = p; *q_out = q; return true;
            }
        }
    }
    return false;
}
