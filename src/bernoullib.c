/* Mathilda -- Bernoulli numbers and polynomials.
 *
 *   BernoulliB[n]      Bernoulli number      B_n
 *   BernoulliB[n, x]   Bernoulli polynomial  B_n(x)
 *
 * The Bernoulli polynomials are the coefficients of the generating function
 *   t e^(x t) / (e^t - 1) = Sum_{n>=0} B_n(x) t^n / n!,
 * and the Bernoulli numbers are B_n = B_n(0). For odd n the numbers vanish
 * except B_1 = -1/2.
 *
 * Evaluation is layered so each argument shape takes the cheapest exact or
 * numeric route:
 *
 *   exact non-negative integer n   ->  exact rational B_n (cached recurrence)
 *   inexact integer-valued n       ->  the rational, numericalised (Real/MPFR)
 *   BernoulliB[n, x]               ->  Sum_j C(n,j) B_{n-j} x^j, built with
 *                                      exact rational coefficients then
 *                                      evaluated (exact x stays exact,
 *                                      inexact x evaluates numerically)
 *   everything else                ->  stays symbolic (return NULL)
 *
 * Bernoulli numbers B_n are computed by the recurrence
 *   B_0 = 1,   B_m = -1/(m+1) Sum_{k=0}^{m-1} C(m+1, k) B_k   (m >= 1),
 * using exact GMP rationals in a lazily-grown, process-lifetime cache. (The
 * same construction lives in src/zeta.c and src/polygamma.c; it is replicated
 * here to keep this module self-contained.)
 *
 * Attributes: Listable, Protected.
 */
#include "bernoullib.h"

#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <gmp.h>

#include "attr.h"
#include "eval.h"          /* eval_and_free */
#include "symtab.h"

#ifdef USE_MPFR
#include <mpfr.h>
#endif

/* Largest index for which BernoulliB[n] builds the exact number, and largest
 * degree for which BernoulliB[n, x] expands the polynomial. The recurrence is
 * O(n^2) in growing big-rationals, so these guard against accidental runaway;
 * beyond them the call stays symbolic. */
#define BERNOULLI_NUMBER_CAP 5000
#define BERNOULLI_POLY_CAP   1000

/* ------------------------------------------------------------------ */
/* Bernoulli numbers B_k (exact, lazily cached, process-lifetime)      */
/* ------------------------------------------------------------------ */

static mpq_t* g_bern = NULL;
static size_t g_bern_len = 0;

static void bern_ensure(size_t upto) {
    if (upto + 1 <= g_bern_len) return;
    size_t newlen = upto + 1;
    mpq_t* grown = (mpq_t*)realloc(g_bern, newlen * sizeof(mpq_t));
    if (!grown) return; /* out of memory: leave cache as-is, callers degrade */
    g_bern = grown;

    mpq_t sum, term, factor;
    mpq_inits(sum, term, factor, (mpq_ptr)0);
    mpz_t binz;
    mpz_init(binz);

    for (size_t m = g_bern_len; m < newlen; m++) {
        mpq_init(g_bern[m]);
        if (m == 0) { mpq_set_ui(g_bern[m], 1, 1); continue; }
        mpq_set_ui(sum, 0, 1);
        for (size_t k = 0; k < m; k++) {
            mpz_bin_uiui(binz, (unsigned long)(m + 1), (unsigned long)k);
            mpq_set_z(factor, binz);
            mpq_mul(term, factor, g_bern[k]);
            mpq_add(sum, sum, term);
        }
        mpq_set_si(term, -1, (unsigned long)(m + 1));
        mpq_canonicalize(term);
        mpq_mul(g_bern[m], sum, term);
        mpq_canonicalize(g_bern[m]);
    }
    g_bern_len = newlen;
    mpz_clear(binz);
    mpq_clears(sum, term, factor, (mpq_ptr)0);
}

/* B_idx into `out`. Odd indices above 1 are exactly zero. */
static void bern_get_q(mpq_t out, size_t idx) {
    if (idx > 1 && (idx & 1u)) { mpq_set_ui(out, 0, 1); return; }
    bern_ensure(idx);
    if (idx < g_bern_len) mpq_set(out, g_bern[idx]);
    else mpq_set_ui(out, 0, 1);
}

/* Build a canonical Integer/BigInt/Rational from an mpq via Times[num,
 * Power[den, -1]] so the evaluator normalises arbitrary-size components. */
static Expr* bern_expr_from_mpq(const mpq_t q) {
    mpz_t num, den;
    mpz_init(num); mpz_init(den);
    mpq_get_num(num, q);
    mpq_get_den(den, q);
    Expr* en = expr_bigint_normalize(expr_new_bigint_from_mpz(num));
    Expr* out;
    if (mpz_cmp_ui(den, 1) == 0) {
        out = en;
    } else {
        Expr* ed = expr_bigint_normalize(expr_new_bigint_from_mpz(den));
        Expr* inv = expr_new_function(expr_new_symbol("Power"),
                        (Expr*[]){ ed, expr_new_integer(-1) }, 2);
        out = eval_and_free(expr_new_function(expr_new_symbol("Times"),
                        (Expr*[]){ en, inv }, 2));
    }
    mpz_clear(num); mpz_clear(den);
    return out;
}

/* ------------------------------------------------------------------ */
/* First-argument index extraction                                     */
/* ------------------------------------------------------------------ */

typedef enum { BIDX_NONE = 0, BIDX_EXACT, BIDX_REAL, BIDX_MPFR } bidx_kind;

/* Recognise a non-negative integer first argument. Sets *idx to its value and,
 * for an MPFR operand, *prec to its precision. Returns:
 *   BIDX_EXACT  exact Integer / BigInt
 *   BIDX_REAL   machine Real holding an integer value
 *   BIDX_MPFR   MPFR holding an integer value (numeric, arbitrary precision)
 *   BIDX_NONE   negative, non-integral, out of range, or non-numeric. */
static bidx_kind bern_index(const Expr* e, unsigned long* idx, long* prec) {
    if (!e) return BIDX_NONE;
    if (e->type == EXPR_INTEGER) {
        if (e->data.integer < 0) return BIDX_NONE;
        *idx = (unsigned long)e->data.integer;
        return BIDX_EXACT;
    }
    if (e->type == EXPR_BIGINT) {
        if (mpz_sgn(e->data.bigint) < 0 || !mpz_fits_ulong_p(e->data.bigint))
            return BIDX_NONE;
        *idx = mpz_get_ui(e->data.bigint);
        return BIDX_EXACT;
    }
    if (e->type == EXPR_REAL) {
        double v = e->data.real;
        if (!(v >= 0.0) || floor(v) != v) return BIDX_NONE;
        *idx = (unsigned long)v;
        return BIDX_REAL;
    }
#ifdef USE_MPFR
    if (e->type == EXPR_MPFR) {
        if (mpfr_sgn(e->data.mpfr) < 0 || !mpfr_integer_p(e->data.mpfr) ||
            !mpfr_fits_ulong_p(e->data.mpfr, MPFR_RNDN))
            return BIDX_NONE;
        *idx  = mpfr_get_ui(e->data.mpfr, MPFR_RNDN);
        *prec = (long)mpfr_get_prec(e->data.mpfr);
        return BIDX_MPFR;
    }
#endif
    (void)prec;
    return BIDX_NONE;
}

/* ------------------------------------------------------------------ */
/* BernoulliB[n]  (the Bernoulli number)                               */
/* ------------------------------------------------------------------ */

static Expr* bern_one_arg(Expr* arg) {
    unsigned long n = 0;
    long prec = 53;
    bidx_kind k = bern_index(arg, &n, &prec);
    if (k == BIDX_NONE) return NULL;          /* symbolic / negative / non-integer */
    if (n > BERNOULLI_NUMBER_CAP) return NULL; /* too large: stay symbolic */

    mpq_t q;
    mpq_init(q);
    bern_get_q(q, (size_t)n);

    Expr* out;
    if (k == BIDX_EXACT) {
        out = bern_expr_from_mpq(q);
    } else if (k == BIDX_REAL) {
        out = expr_new_real(mpq_get_d(q));
    } else { /* BIDX_MPFR */
#ifdef USE_MPFR
        out = expr_new_mpfr_bits((mpfr_prec_t)prec);
        mpfr_set_q(out->data.mpfr, q, MPFR_RNDN);
#else
        out = bern_expr_from_mpq(q);
#endif
    }
    mpq_clear(q);
    return out;
}

/* ------------------------------------------------------------------ */
/* BernoulliB[n, x]  (the Bernoulli polynomial)                        */
/* ------------------------------------------------------------------ */

/* B_n(x) = Sum_{j=0}^{n} C(n,j) B_{n-j} x^j. Builds the symbolic polynomial
 * and evaluates it once. Coefficients C(n,j) B_{n-j} are exact GMP rationals,
 * so large n keep full precision; zero coefficients (B_{n-j} = 0 for odd
 * n-j > 1) are skipped. */
static Expr* bernoulli_polynomial(size_t n, Expr* x) {
    Expr** terms = (Expr**)malloc((n + 1) * sizeof(Expr*));
    if (!terms) return NULL;
    size_t count = 0;

    mpq_t coeff, bq, binq;
    mpq_inits(coeff, bq, binq, (mpq_ptr)0);
    mpz_t binz;
    mpz_init(binz);

    for (size_t j = 0; j <= n; j++) {
        bern_get_q(bq, n - j);                 /* B_{n-j} */
        if (mpq_sgn(bq) == 0) continue;        /* vanishing coefficient */
        mpz_bin_uiui(binz, (unsigned long)n, (unsigned long)j); /* C(n, j) */
        mpq_set_z(binq, binz);
        mpq_mul(coeff, binq, bq);
        mpq_canonicalize(coeff);

        Expr* c = bern_expr_from_mpq(coeff);
        Expr* term;
        if (j == 0) {
            term = c;                          /* x^0 = 1 */
        } else {
            Expr* xj = expr_new_function(expr_new_symbol("Power"),
                          (Expr*[]){ expr_copy(x), expr_new_integer((int64_t)j) }, 2);
            term = expr_new_function(expr_new_symbol("Times"),
                          (Expr*[]){ c, xj }, 2);
        }
        terms[count++] = term;
    }

    mpz_clear(binz);
    mpq_clears(coeff, bq, binq, (mpq_ptr)0);

    Expr* poly;
    if (count == 0) {
        poly = expr_new_integer(0);            /* unreachable: x^n term is nonzero */
    } else if (count == 1) {
        poly = terms[0];
    } else {
        poly = expr_new_function(expr_new_symbol("Plus"), terms, count);
    }
    free(terms);
    return eval_and_free(poly);
}

static Expr* bern_two_arg(Expr* nexpr, Expr* x) {
    unsigned long n = 0;
    long prec = 53;
    bidx_kind k = bern_index(nexpr, &n, &prec);
    if (k == BIDX_NONE) return NULL;           /* symbolic / negative / non-integer n */
    if (n > BERNOULLI_POLY_CAP) return NULL;    /* too large: stay symbolic */
    return bernoulli_polynomial((size_t)n, x);
}

/* ------------------------------------------------------------------ */
/* Builtin entry point                                                 */
/* ------------------------------------------------------------------ */

/* Print a Mathematica-compatible argt diagnostic for a wrong argument count
 * and return NULL so the evaluator leaves the call unevaluated. */
static Expr* bern_emit_argt(size_t argc) {
    fprintf(stderr,
            "BernoulliB::argt: BernoulliB called with %zu argument%s; "
            "1 or 2 arguments are expected.\n",
            argc, argc == 1 ? "" : "s");
    return NULL;
}

Expr* builtin_bernoullib(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    Expr** args = res->data.function.args;

    if (argc == 1) return bern_one_arg(args[0]);
    if (argc == 2) return bern_two_arg(args[0], args[1]);
    return bern_emit_argt(argc);
}

/* ------------------------------------------------------------------ */
/* Registration                                                        */
/* ------------------------------------------------------------------ */

void bernoullib_init(void) {
    symtab_add_builtin("BernoulliB", builtin_bernoullib);
    symtab_get_def("BernoulliB")->attributes |= (ATTR_LISTABLE | ATTR_PROTECTED);
    /* Docstring lives in info.c (info_init). */
}
