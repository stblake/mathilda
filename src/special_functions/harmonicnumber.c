/* Mathilda -- HarmonicNumber: generalized (order-r) harmonic numbers.
 *
 *   HarmonicNumber[n]     H_n     = Sum_{i=1}^n 1/i
 *   HarmonicNumber[n, r]  H_n^(r) = Sum_{i=1}^n 1/i^r
 *
 * Rather than carry bespoke numeric kernels, HarmonicNumber reduces to the
 * primitives the system already provides and lets the evaluator finish the job:
 *
 *   n -> Infinity             ->  Zeta[r]
 *   r a non-positive integer  ->  Faulhaber polynomial (built from BernoulliB)
 *   n a non-negative integer,
 *     r a positive integer    ->  exact rational by GMP binary splitting
 *                                 (n <= HN_EXPAND_CAP; ~10 ms at n = 10^5)
 *     r inexact and numeric   ->  explicit sum for n <= HN_SUM_CAP, else the
 *                                 analytic reduction below
 *     r anything else         ->  explicit Plus  1 + 2^-r + ... + n^-r  for
 *                                 n <= HN_SUM_CAP, else symbolic
 *   inexact argument          ->  N[ Zeta[r] - Zeta[r, n+1] ]   (r != 1)
 *                                 N[ EulerGamma + PolyGamma[0, n+1] ]  (r == 1)
 *                                 (a machine-real n with r == 1 goes straight
 *                                  to the ND kernel's helper sf_machine_harmonic)
 *   otherwise                 ->  symbolic (return NULL)
 *
 * The analytic identity  H_n^(r) = Zeta[r] - Zeta[r, n+1]  (and its r == 1
 * digamma special case) carries arbitrary precision and complex arguments
 * straight through Zeta / PolyGamma.
 *
 * Memory: builtin_harmonicnumber takes ownership of `res` but must not free it
 * (the evaluator does).  Every Expr* built here is owned and either handed to
 * expr_new_function (which adopts it) or released via eval_and_free.
 *
 * Attributes: Listable, NumericFunction, Protected.
 */
#include "harmonicnumber.h"

#include "sym_names.h"
#include "expr.h"
#include "eval.h"          /* evaluate, eval_and_free */
#include "symtab.h"
#include "attr.h"
#include "arithmetic.h"    /* is_infinity_sym */
#include "numeric.h"       /* numeric_min_inexact_bits, numeric_bits_to_digits */
#include "message.h"   /* mth_message: Quiet/Check funnel */
#include "sf_machine.h"    /* sf_machine_harmonic: the ND kernel's helper */

#include <stdbool.h>
#include <string.h>
#include <stdint.h>
#include <stdio.h>
#include <math.h>
#include <gmp.h>

/* Largest non-negative integer n for which HarmonicNumber[n, r] with a positive
 * integer order r evaluates to its exact rational.  Beyond the cap the result
 * stays symbolic, and N[] then takes the analytic route on the numericized
 * argument (N[HarmonicNumber[10^9], 100] is sub-millisecond). */
#define HN_EXPAND_CAP 100000L

/* Bound on the size of that exact rational, in bits of the unreduced
 * denominator r * log2(n!).  HarmonicNumber[10^5] needs ~1.5e6 bits; the bound
 * only keeps a pathological order (HarmonicNumber[10^5, 10^6], a ~10^11-bit
 * number) from exhausting memory -- such an input stays symbolic. */
#define HN_EXACT_MAX_BITS 268435456.0   /* 2^28 */

/* Largest n for which an order r that is not an integer (symbolic, rational,
 * complex, or inexact) expands to the explicit sum 1 + 2^-r + ... + n^-r --
 * the threshold the reference applies (HarmonicNumber[4, r] expands,
 * HarmonicNumber[5, r] does not).  Past it an exact r stays symbolic and an
 * inexact numeric r takes the analytic route, instead of building and
 * evaluating an n-term sum (HarmonicNumber[10^5, 1/2] used to hang). */
#define HN_SUM_CAP 4

/* ------------------------------------------------------------------ */
/*  Tiny builders                                                      */
/* ------------------------------------------------------------------ */

static Expr* hn_int(int64_t v) { return expr_new_integer(v); }

/* head[a] from an owned arg. */
static Expr* hn_fn1(const char* head, Expr* a) {
    return expr_new_function(expr_new_symbol(head), (Expr*[]){ a }, 1);
}
/* head[a, b] from owned args. */
static Expr* hn_fn2(const char* head, Expr* a, Expr* b) {
    return expr_new_function(expr_new_symbol(head), (Expr*[]){ a, b }, 2);
}

/* ------------------------------------------------------------------ */
/*  Predicates                                                         */
/* ------------------------------------------------------------------ */

/* Print `HarmonicNumber::argt: HarmonicNumber called with N arguments; 1 or 2
 * arguments are expected.` (Mathematica uses `argt` for the 1-OR-2 case). */
static Expr* hn_emit_argt(size_t argc) {
    mth_message("HarmonicNumber", "argt", "HarmonicNumber called with %zu argument%s; 1 or 2 arguments are expected.", argc, argc == 1 ? "" : "s");
    return NULL;
}

/* True if e contains an inexact numeric leaf (Real or MPFR) anywhere. */
static bool hn_contains_inexact(const Expr* e) {
    if (!e) return false;
    if (e->type == EXPR_REAL) return true;
#ifdef USE_MPFR
    if (e->type == EXPR_MPFR) return true;
#endif
    if (e->type == EXPR_FUNCTION) {
        if (hn_contains_inexact(e->data.function.head)) return true;
        for (size_t i = 0; i < e->data.function.arg_count; i++)
            if (hn_contains_inexact(e->data.function.args[i])) return true;
    }
    return false;
}

/* True if n can be given a numeric value: a number leaf, a Rational/Complex of
 * numbers, a Constant symbol (Pi, E, EulerGamma, ...), or any expression that
 * already carries an inexact leaf.  Guards the numeric path against generic
 * free symbols (HarmonicNumber[x, 2.5] must stay symbolic in x). */
static bool hn_numericizable(const Expr* e) {
    if (!e) return false;
    switch (e->type) {
        case EXPR_INTEGER:
        case EXPR_BIGINT:
        case EXPR_REAL:
            return true;
#ifdef USE_MPFR
        case EXPR_MPFR:
            return true;
#endif
        case EXPR_SYMBOL: {
            SymbolDef* d = symtab_lookup(e->data.symbol.name);
            return d && (d->attributes & ATTR_CONSTANT);
        }
        case EXPR_FUNCTION: {
            const Expr* h = e->data.function.head;
            if (h->type == EXPR_SYMBOL &&
                (h->data.symbol.name == SYM_Rational ||
                 h->data.symbol.name == SYM_Complex)) {
                for (size_t i = 0; i < e->data.function.arg_count; i++)
                    if (!hn_numericizable(e->data.function.args[i])) return false;
                return true;
            }
            return hn_contains_inexact(e);
        }
        default:
            return false;
    }
}

/* If e is an exact non-negative integer (Integer >= 0, or a non-negative BigInt
 * within machine range), store it in *out and return true. */
static bool hn_nonneg_int(const Expr* e, int64_t* out) {
    if (e->type == EXPR_INTEGER) {
        if (e->data.integer < 0) return false;
        *out = e->data.integer;
        return true;
    }
    if (e->type == EXPR_BIGINT) {
        if (mpz_sgn(e->data.bigint) < 0) return false;
        if (!mpz_fits_slong_p(e->data.bigint)) return false;
        *out = (int64_t)mpz_get_si(e->data.bigint);
        return true;
    }
    return false;
}

/* If e is an exact integer (Integer, or BigInt within machine range), store it
 * in *out and return true. */
static bool hn_exact_int(const Expr* e, int64_t* out) {
    if (e->type == EXPR_INTEGER) { *out = e->data.integer; return true; }
    if (e->type == EXPR_BIGINT) {
        if (!mpz_fits_slong_p(e->data.bigint)) return false;
        *out = (int64_t)mpz_get_si(e->data.bigint);
        return true;
    }
    return false;
}

/* True if r is numerically the value 1 (exact Integer 1, Real 1.0, MPFR 1). */
static bool hn_is_one(const Expr* r) {
    if (r->type == EXPR_INTEGER) return r->data.integer == 1;
    if (r->type == EXPR_REAL)    return r->data.real == 1.0;
#ifdef USE_MPFR
    if (r->type == EXPR_MPFR)    return mpfr_cmp_ui(r->data.mpfr, 1) == 0;
#endif
    return false;
}

/* ------------------------------------------------------------------ */
/*  Case builders                                                      */
/* ------------------------------------------------------------------ */

/* n a non-negative integer: HarmonicNumber[n, r] = Sum_{i=1}^n i^-r.
 * Builds Plus[ Power[1,-r], ..., Power[n,-r] ] and evaluates it.  n == 0 gives
 * the empty sum 0.  `r` is borrowed (copied per term). */
static Expr* hn_finite_sum(int64_t n, const Expr* r) {
    if (n == 0) return hn_int(0);

    Expr** terms = malloc(sizeof(Expr*) * (size_t)n);
    for (int64_t i = 1; i <= n; i++) {
        /* exponent = -r (kept symbolic; the evaluator folds Power[i, -r]) */
        Expr* negr = hn_fn2("Times", hn_int(-1), expr_copy((Expr*)r));
        terms[i - 1] = hn_fn2("Power", hn_int(i), negr);
    }
    Expr* plus = expr_new_function(expr_new_symbol("Plus"), terms, (size_t)n);
    free(terms);
    return eval_and_free(plus);
}

/* Binary splitting for the exact partial sum
 *     Sum_{i=a}^{b-1} i^-r  =  P / Q.
 * Halves combine as P = P_L Q_R + P_R Q_L, Q = Q_L Q_R, so the work is a few
 * balanced big multiplications per level -- O(M(N) log N) for an N-bit result
 * -- instead of the n sequential rational additions (each O(N) and growing) of
 * a left-to-right Plus, which made HarmonicNumber[10^4] take ~0.8 s and
 * HarmonicNumber[10^5] effectively hang.  Q is the full product Prod i^r; one
 * gcd at the end reduces the result.  (Reducing at intermediate levels too was
 * measured slower at every threshold tried -- the gcds cost more than the
 * smaller operands save.) */
#define HN_BS_LEAF   16

static void hn_bsplit(unsigned long a, unsigned long b, unsigned long r,
                      mpz_t P, mpz_t Q) {
    if (b - a <= HN_BS_LEAF) {
        mpz_t t;
        mpz_init(t);
        mpz_set_ui(P, 0);
        mpz_set_ui(Q, 1);
        for (unsigned long i = a; i < b; i++) {
            mpz_ui_pow_ui(t, i, r);          /* P/Q + 1/t = (P t + Q) / (Q t) */
            mpz_mul(P, P, t);
            mpz_add(P, P, Q);
            mpz_mul(Q, Q, t);
        }
        mpz_clear(t);
        return;
    }
    unsigned long m = a + (b - a) / 2;
    mpz_t P2, Q2;
    mpz_init(P2); mpz_init(Q2);
    hn_bsplit(a, m, r, P, Q);
    hn_bsplit(m, b, r, P2, Q2);
    mpz_mul(P, P, Q2);
    mpz_addmul(P, P2, Q);
    mpz_mul(Q, Q, Q2);
    mpz_clear(P2); mpz_clear(Q2);
}

/* n a non-negative integer, r a positive integer: the exact rational H_n^(r).
 * Returns NULL (stay symbolic) when it would exceed HN_EXACT_MAX_BITS. */
static Expr* hn_exact_rational(int64_t n, int64_t r) {
    if (n == 0) return hn_int(0);
    if (n == 1) return hn_int(1);
    /* r * log2(n!) bounds the unreduced denominator. */
    double bits = (double)r * lgamma((double)n + 1.0) / log(2.0);
    if (!(bits <= HN_EXACT_MAX_BITS)) return NULL;
    mpz_t P, Q;
    mpz_init(P); mpz_init(Q);
    hn_bsplit(1UL, (unsigned long)n + 1UL, (unsigned long)r, P, Q);
    Expr* out = make_rational_mpz(P, Q);   /* the final gcd reduction */
    mpz_clear(P); mpz_clear(Q);
    return out;
}

/* r a non-positive integer (-m, m >= 0): the Faulhaber closed form
 *   Sum_{i=1}^n i^m = n^m + 1/(m+1) Sum_{j=0}^m C(m+1,j) BernoulliB[j] n^(m+1-j),
 * a polynomial in n (with the B_1 = +1/2 summation convention realised by the
 * extra n^m term).  m == 0 gives n.  `n` is borrowed. */
static Expr* hn_faulhaber(const Expr* n, int64_t m) {
    if (m == 0) return expr_copy((Expr*)n);

    /* inner = Sum_{j=0}^m C(m+1,j) * BernoulliB[j] * n^(m+1-j) */
    size_t nt = (size_t)m + 1;
    Expr** terms = malloc(sizeof(Expr*) * nt);
    for (int64_t j = 0; j <= m; j++) {
        Expr* binom = hn_fn2("Binomial", hn_int(m + 1), hn_int(j));
        Expr* bern  = hn_fn1("BernoulliB", hn_int(j));
        Expr* powr  = hn_fn2("Power", expr_copy((Expr*)n), hn_int(m + 1 - j));
        Expr* coef  = hn_fn2("Times", binom, bern);
        terms[j]    = hn_fn2("Times", coef, powr);
    }
    Expr* inner = expr_new_function(expr_new_symbol("Plus"), terms, nt);
    free(terms);

    /* P = n^m + inner/(m+1) */
    Expr* scaled = hn_fn2("Times", hn_fn2("Power", hn_int(m + 1), hn_int(-1)), inner);
    Expr* nm     = hn_fn2("Power", expr_copy((Expr*)n), hn_int(m));
    Expr* sum    = hn_fn2("Plus", nm, scaled);

    /* Expand to the canonical polynomial form. */
    return eval_and_free(hn_fn1("Expand", sum));
}

/* Inexact / numericizable argument: reduce to the analytic form and N it.
 *   r == 1 : N[ EulerGamma + PolyGamma[0, n+1] ]
 *   else   : N[ Zeta[r] - Zeta[r, n+1] ]
 * Precision is taken from the inexact inputs (machine unless an MPFR argument
 * raises it).  `n`, `r` are borrowed. */
static Expr* hn_numeric(const Expr* n, const Expr* r) {
    /* Force numericization (so Constant n such as E resolve too).  Use the
     * input precision: machine -> N[reduction]; MPFR -> N[reduction, digits]. */
    long bits = 0;
#ifdef USE_MPFR
    long bn = numeric_min_inexact_bits(n);
    long br = numeric_min_inexact_bits(r);
    if (bn > 0) bits = bn;
    if (br > 0 && (bits == 0 || br < bits)) bits = br;
#endif
    long digits = 0;
    if (bits > 53) {
        digits = (long)numeric_bits_to_digits(bits);
        if (digits < 1) digits = 1;
    }

    /* An exact integer n (reached with an inexact order, as in
     * HarmonicNumber[10^5, 2.5]) is numericized first: left exact,
     * PolyGamma[0, n+1] and Zeta[r, n+1] would evaluate to their exact
     * closed forms -- an n-term sum -- before the outer N ever ran. */
    Expr* nv = expr_copy((Expr*)n);
    if (n->type == EXPR_INTEGER || n->type == EXPR_BIGINT)
        nv = eval_and_free(digits > 0 ? hn_fn2("N", nv, hn_int(digits))
                                      : hn_fn1("N", nv));
    Expr* np1 = hn_fn2("Plus", nv, hn_int(1));

    Expr* reduction;
    if (hn_is_one(r)) {
        Expr* pg = hn_fn2("PolyGamma", hn_int(0), np1);
        reduction = hn_fn2("Plus", expr_new_symbol("EulerGamma"), pg);
    } else {
        Expr* zr  = hn_fn1("Zeta", expr_copy((Expr*)r));
        Expr* zra = hn_fn2("Zeta", expr_copy((Expr*)r), np1);
        Expr* neg = hn_fn2("Times", hn_int(-1), zra);
        reduction = hn_fn2("Plus", zr, neg);
    }

    Expr* ncall = digits > 0 ? hn_fn2("N", reduction, hn_int(digits))
                             : hn_fn1("N", reduction);
    return eval_and_free(ncall);
}

/* ------------------------------------------------------------------ */
/*  Entry point                                                        */
/* ------------------------------------------------------------------ */

Expr* builtin_harmonicnumber(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    if (argc < 1 || argc > 2) return hn_emit_argt(argc);

    Expr* n = res->data.function.args[0];
    /* The order r defaults to 1.  Hold an owned copy so every case can
     * expr_copy it freely; release it once before returning. */
    Expr* r = (argc == 2) ? expr_copy(res->data.function.args[1]) : hn_int(1);

    Expr* result = NULL;
    int64_t ni, ri;

    if (is_infinity_sym(n)) {
        /* HarmonicNumber[Infinity, r] = Zeta[r]. */
        result = eval_and_free(hn_fn1("Zeta", expr_copy(r)));
    } else if (n->type == EXPR_REAL && hn_is_one(r)) {
        /* Machine-real n, order 1: the ND kernel's own helper, so the scalar
         * answer bit-matches the packed / NDArray / Compile surfaces.  A pole
         * (a negative integer n) declines to the analytic route. */
        double v;
        if (sf_machine_harmonic(n->data.real, &v))
            result = expr_new_real(v);
        else
            result = hn_numeric(n, r);
    } else if (hn_exact_int(r, &ri) && ri <= 0) {
        /* r a non-positive integer: Faulhaber polynomial in n -- |r|+1 terms
         * whatever n is.  For an integer n with fewer terms than that, the
         * plain sum of integers 1^-r + ... + n^-r is the cheaper route. */
        if (hn_nonneg_int(n, &ni) && ni <= 1 - ri)
            result = hn_finite_sum(ni, r);
        else
            result = hn_faulhaber(n, -ri);
    } else if (hn_nonneg_int(n, &ni)) {
        if (hn_exact_int(r, &ri)) {
            /* r a positive integer: the exact rational, within the cap. */
            if (ni <= HN_EXPAND_CAP) result = hn_exact_rational(ni, ri);
        } else if (ni <= HN_SUM_CAP) {
            /* Any other order: the short explicit sum. */
            result = hn_finite_sum(ni, r);
        } else if (hn_contains_inexact(r) && hn_numericizable(r)) {
            /* Inexact numeric order: analytic reduction. */
            result = hn_numeric(n, r);
        }
    } else if ((hn_contains_inexact(n) || hn_contains_inexact(r))
               && hn_numericizable(n)) {
        /* Inexact, numericizable argument: analytic reduction, numericized. */
        result = hn_numeric(n, r);
    }
    /* else: everything else stays symbolic (result == NULL). */

    expr_free(r);
    return result;
}

/* ------------------------------------------------------------------ */
/*  Registration                                                       */
/* ------------------------------------------------------------------ */

void harmonicnumber_init(void) {
    symtab_add_builtin("HarmonicNumber", builtin_harmonicnumber);
    symtab_get_def("HarmonicNumber")->attributes |=
        (ATTR_LISTABLE | ATTR_NUMERICFUNCTION | ATTR_PROTECTED);
    /* Docstring lives in info.c (info_init). */
}
