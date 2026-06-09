/* Mathilda -- Euler numbers and polynomials.
 *
 *   EulerE[n]      Euler number      E_n
 *   EulerE[n, x]   Euler polynomial  E_n(x)
 *
 * The Euler polynomials are the coefficients of the generating function
 *   2 e^(x t) / (e^t + 1) = Sum_{n>=0} E_n(x) t^n / n!,
 * and the Euler numbers are E_n = 2^n E_n(1/2). For odd n the numbers vanish;
 * E_0 = 1, E_2 = -1, E_4 = 5, E_6 = -61, ...
 *
 * Evaluation is layered so each argument shape takes the cheapest exact or
 * numeric route:
 *
 *   exact non-negative integer n   ->  exact integer E_n (cached recurrence)
 *   inexact integer-valued n       ->  the integer, numericalised (Real/MPFR)
 *   EulerE[n, x]                   ->  the degree-n polynomial in monomial
 *                                      form with exact rational coefficients,
 *                                      then evaluated (exact x stays exact,
 *                                      inexact x evaluates numerically)
 *   EulerE[n, 1/2], symbolic n     ->  2^-n EulerE[n]
 *   everything else                ->  stays symbolic (return NULL)
 *
 * Euler numbers E_n are computed by the recurrence
 *   E_0 = 1,   E_{2m} = -Sum_{k=0}^{m-1} C(2m, 2k) E_{2k}   (m >= 1),
 * with odd-index numbers identically zero, using exact GMP integers in a
 * lazily-grown, process-lifetime cache.
 *
 * The Euler polynomial coefficients are obtained from the Taylor expansion
 * about x = 1/2.  Writing E_n(x) = Sum_i c_i x^i, the algebra collapses (the
 * powers of two in C(n,j) E_{n-j}/2^{n-j} (x-1/2)^j combine to a single
 * 2^{n-i}) to the all-integer inner sum
 *   S_i = Sum_{j=i}^{n} (-1)^{j-i} C(n,j) C(j,i) E_{n-j},   c_i = S_i / 2^{n-i}.
 *
 * Attributes: Listable, Protected.
 */
#include "eulere.h"

#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <gmp.h>

#include "attr.h"
#include "eval.h"          /* eval_and_free */
#include "sym_names.h"     /* SYM_Rational */
#include "symtab.h"

#ifdef USE_MPFR
#include <mpfr.h>
#endif

/* Largest index for which EulerE[n] builds the exact number, and largest
 * degree for which EulerE[n, x] expands the polynomial. The recurrence is
 * O(n^2) in growing big-integers, so these guard against accidental runaway;
 * beyond them the call stays symbolic. */
#define EULER_NUMBER_CAP 5000
#define EULER_POLY_CAP   1000

/* ------------------------------------------------------------------ */
/* Euler numbers E_k (exact integers, lazily cached, process-lifetime) */
/* ------------------------------------------------------------------ */

static mpz_t* g_euler = NULL;
static size_t g_euler_len = 0;

static void euler_ensure(size_t upto) {
    if (upto + 1 <= g_euler_len) return;
    size_t newlen = upto + 1;
    mpz_t* grown = (mpz_t*)realloc(g_euler, newlen * sizeof(mpz_t));
    if (!grown) return; /* out of memory: leave cache as-is, callers degrade */
    g_euler = grown;

    mpz_t sum, term, binz;
    mpz_inits(sum, term, binz, (mpz_ptr)0);

    for (size_t m = g_euler_len; m < newlen; m++) {
        mpz_init(g_euler[m]);
        if (m == 0) { mpz_set_ui(g_euler[m], 1); continue; }
        if (m & 1u) { mpz_set_ui(g_euler[m], 0); continue; } /* odd -> 0 */
        /* E_m = -Sum_{k even, k < m} C(m, k) E_k. */
        mpz_set_ui(sum, 0);
        for (size_t k = 0; k < m; k += 2) {
            mpz_bin_uiui(binz, (unsigned long)m, (unsigned long)k);
            mpz_mul(term, binz, g_euler[k]);
            mpz_add(sum, sum, term);
        }
        mpz_neg(g_euler[m], sum);
    }
    g_euler_len = newlen;
    mpz_clears(sum, term, binz, (mpz_ptr)0);
}

/* E_idx into `out`. Odd indices are exactly zero. */
static void euler_get_z(mpz_t out, size_t idx) {
    if (idx & 1u) { mpz_set_ui(out, 0); return; }
    euler_ensure(idx);
    if (idx < g_euler_len) mpz_set(out, g_euler[idx]);
    else mpz_set_ui(out, 0);
}

/* Build a canonical Integer/BigInt from an mpz. */
static Expr* euler_expr_from_mpz(const mpz_t z) {
    return expr_bigint_normalize(expr_new_bigint_from_mpz(z));
}

/* Build a canonical Integer/BigInt/Rational from an mpq via Times[num,
 * Power[den, -1]] so the evaluator normalises arbitrary-size components. */
static Expr* euler_expr_from_mpq(const mpq_t q) {
    mpz_t num, den;
    mpz_init(num); mpz_init(den);
    mpq_get_num(num, q);
    mpq_get_den(den, q);
    Expr* en = euler_expr_from_mpz(num);
    Expr* out;
    if (mpz_cmp_ui(den, 1) == 0) {
        out = en;
    } else {
        Expr* ed = euler_expr_from_mpz(den);
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

typedef enum { EIDX_NONE = 0, EIDX_EXACT, EIDX_REAL, EIDX_MPFR } eidx_kind;

/* Recognise a non-negative integer first argument. Sets *idx to its value and,
 * for an MPFR operand, *prec to its precision. Returns:
 *   EIDX_EXACT  exact Integer / BigInt
 *   EIDX_REAL   machine Real holding an integer value
 *   EIDX_MPFR   MPFR holding an integer value (numeric, arbitrary precision)
 *   EIDX_NONE   negative, non-integral, out of range, or non-numeric. */
static eidx_kind euler_index(const Expr* e, unsigned long* idx, long* prec) {
    if (!e) return EIDX_NONE;
    if (e->type == EXPR_INTEGER) {
        if (e->data.integer < 0) return EIDX_NONE;
        *idx = (unsigned long)e->data.integer;
        return EIDX_EXACT;
    }
    if (e->type == EXPR_BIGINT) {
        if (mpz_sgn(e->data.bigint) < 0 || !mpz_fits_ulong_p(e->data.bigint))
            return EIDX_NONE;
        *idx = mpz_get_ui(e->data.bigint);
        return EIDX_EXACT;
    }
    if (e->type == EXPR_REAL) {
        double v = e->data.real;
        if (!(v >= 0.0) || floor(v) != v) return EIDX_NONE;
        *idx = (unsigned long)v;
        return EIDX_REAL;
    }
#ifdef USE_MPFR
    if (e->type == EXPR_MPFR) {
        if (mpfr_sgn(e->data.mpfr) < 0 || !mpfr_integer_p(e->data.mpfr) ||
            !mpfr_fits_ulong_p(e->data.mpfr, MPFR_RNDN))
            return EIDX_NONE;
        *idx  = mpfr_get_ui(e->data.mpfr, MPFR_RNDN);
        *prec = (long)mpfr_get_prec(e->data.mpfr);
        return EIDX_MPFR;
    }
#endif
    (void)prec;
    return EIDX_NONE;
}

/* True when e is the exact rational 1/2 (Rational[1, 2]). */
static bool euler_is_half(const Expr* e) {
    if (!e || e->type != EXPR_FUNCTION) return false;
    Expr* h = e->data.function.head;
    if (!h || h->type != EXPR_SYMBOL || h->data.symbol != SYM_Rational)
        return false;
    if (e->data.function.arg_count != 2) return false;
    Expr* a = e->data.function.args[0];
    Expr* b = e->data.function.args[1];
    return a && b && a->type == EXPR_INTEGER && a->data.integer == 1 &&
           b->type == EXPR_INTEGER && b->data.integer == 2;
}

/* ------------------------------------------------------------------ */
/* EulerE[n]  (the Euler number)                                       */
/* ------------------------------------------------------------------ */

static Expr* euler_one_arg(Expr* arg) {
    unsigned long n = 0;
    long prec = 53;
    eidx_kind k = euler_index(arg, &n, &prec);
    if (k == EIDX_NONE) return NULL;          /* symbolic / negative / non-integer */
    if (n > EULER_NUMBER_CAP) return NULL;     /* too large: stay symbolic */

    mpz_t z;
    mpz_init(z);
    euler_get_z(z, (size_t)n);

    Expr* out;
    if (k == EIDX_EXACT) {
        out = euler_expr_from_mpz(z);
    } else if (k == EIDX_REAL) {
        out = expr_new_real(mpz_get_d(z));
    } else { /* EIDX_MPFR */
#ifdef USE_MPFR
        out = expr_new_mpfr_bits((mpfr_prec_t)prec);
        mpfr_set_z(out->data.mpfr, z, MPFR_RNDN);
#else
        out = euler_expr_from_mpz(z);
#endif
    }
    mpz_clear(z);
    return out;
}

/* ------------------------------------------------------------------ */
/* EulerE[n, x]  (the Euler polynomial)                                */
/* ------------------------------------------------------------------ */

/* E_n(x) = Sum_i c_i x^i with c_i = S_i / 2^{n-i} and the integer inner sum
 * S_i = Sum_{j=i}^{n} (-1)^{j-i} C(n,j) C(j,i) E_{n-j}. Builds the symbolic
 * polynomial in monomial form and evaluates it once. Coefficients are exact
 * GMP rationals, so large n keep full precision; zero coefficients are
 * skipped. */
static Expr* euler_polynomial(size_t n, Expr* x) {
    Expr** terms = (Expr**)malloc((n + 1) * sizeof(Expr*));
    if (!terms) return NULL;
    size_t count = 0;

    mpz_t S, em, cnj, cji, prod;
    mpz_inits(S, em, cnj, cji, prod, (mpz_ptr)0);
    mpq_t coeff, two_pow;
    mpq_inits(coeff, two_pow, (mpq_ptr)0);

    for (size_t i = 0; i <= n; i++) {
        mpz_set_ui(S, 0);
        for (size_t j = i; j <= n; j++) {
            euler_get_z(em, n - j);            /* E_{n-j} (0 for odd n-j) */
            if (mpz_sgn(em) == 0) continue;
            mpz_bin_uiui(cnj, (unsigned long)n, (unsigned long)j); /* C(n, j) */
            mpz_bin_uiui(cji, (unsigned long)j, (unsigned long)i); /* C(j, i) */
            mpz_mul(prod, cnj, cji);
            mpz_mul(prod, prod, em);
            if ((j - i) & 1u) mpz_sub(S, S, prod);
            else              mpz_add(S, S, prod);
        }
        if (mpz_sgn(S) == 0) continue;          /* vanishing coefficient */

        /* coeff = S / 2^{n-i}. */
        mpq_set_z(coeff, S);
        mpz_set_ui(em, 1);
        mpz_mul_2exp(em, em, (mp_bitcnt_t)(n - i)); /* 2^{n-i} */
        mpq_set_z(two_pow, em);
        mpq_div(coeff, coeff, two_pow);
        mpq_canonicalize(coeff);

        Expr* c = euler_expr_from_mpq(coeff);
        Expr* term;
        if (i == 0) {
            term = c;                          /* x^0 = 1 */
        } else {
            Expr* xi = expr_new_function(expr_new_symbol("Power"),
                          (Expr*[]){ expr_copy(x), expr_new_integer((int64_t)i) }, 2);
            term = expr_new_function(expr_new_symbol("Times"),
                          (Expr*[]){ c, xi }, 2);
        }
        terms[count++] = term;
    }

    mpz_clears(S, em, cnj, cji, prod, (mpz_ptr)0);
    mpq_clears(coeff, two_pow, (mpq_ptr)0);

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

/* EulerE[n, 1/2] -> 2^-n EulerE[n] for symbolic n. */
static Expr* euler_half_rule(Expr* nexpr) {
    Expr* base   = expr_new_integer(2);
    Expr* negn   = expr_new_function(expr_new_symbol("Times"),
                       (Expr*[]){ expr_new_integer(-1), expr_copy(nexpr) }, 2);
    Expr* twopow = expr_new_function(expr_new_symbol("Power"),
                       (Expr*[]){ base, negn }, 2);
    Expr* en     = expr_new_function(expr_new_symbol("EulerE"),
                       (Expr*[]){ expr_copy(nexpr) }, 1);
    return eval_and_free(expr_new_function(expr_new_symbol("Times"),
                       (Expr*[]){ twopow, en }, 2));
}

static Expr* euler_two_arg(Expr* nexpr, Expr* x) {
    unsigned long n = 0;
    long prec = 53;
    eidx_kind k = euler_index(nexpr, &n, &prec);
    if (k == EIDX_NONE) {                       /* symbolic / negative / non-integer n */
        if (euler_is_half(x)) return euler_half_rule(nexpr);
        return NULL;
    }
    if (n > EULER_POLY_CAP) return NULL;         /* too large: stay symbolic */
    return euler_polynomial((size_t)n, x);
}

/* ------------------------------------------------------------------ */
/* Builtin entry point                                                 */
/* ------------------------------------------------------------------ */

/* Print a Mathematica-compatible argt diagnostic for a wrong argument count
 * and return NULL so the evaluator leaves the call unevaluated. */
static Expr* euler_emit_argt(size_t argc) {
    fprintf(stderr,
            "EulerE::argt: EulerE called with %zu argument%s; "
            "1 or 2 arguments are expected.\n",
            argc, argc == 1 ? "" : "s");
    return NULL;
}

Expr* builtin_eulere(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    Expr** args = res->data.function.args;

    if (argc == 1) return euler_one_arg(args[0]);
    if (argc == 2) return euler_two_arg(args[0], args[1]);
    return euler_emit_argt(argc);
}

/* ------------------------------------------------------------------ */
/* Registration                                                        */
/* ------------------------------------------------------------------ */

void eulere_init(void) {
    symtab_add_builtin("EulerE", builtin_eulere);
    symtab_get_def("EulerE")->attributes |= (ATTR_LISTABLE | ATTR_PROTECTED);
    /* Docstring lives in info.c (info_init). */
}
