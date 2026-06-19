/* Mathilda -- HarmonicNumber: generalized (order-r) harmonic numbers.
 *
 *   HarmonicNumber[n]     H_n     = Sum_{i=1}^n 1/i
 *   HarmonicNumber[n, r]  H_n^(r) = Sum_{i=1}^n 1/i^r
 *
 * Rather than carry bespoke numeric kernels, HarmonicNumber reduces to the
 * primitives the system already provides and lets the evaluator finish the job:
 *
 *   n a non-negative integer  ->  explicit finite sum  Sum_{i=1}^n i^-r
 *                                 (combines to an exact rational for integer r,
 *                                  stays an explicit Plus for symbolic/complex r)
 *   n -> Infinity             ->  Zeta[r]
 *   r a non-positive integer  ->  Faulhaber polynomial (built from BernoulliB)
 *   inexact argument          ->  N[ Zeta[r] - Zeta[r, n+1] ]   (r != 1)
 *                                 N[ EulerGamma + PolyGamma[0, n+1] ]  (r == 1)
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

#include <stdbool.h>
#include <string.h>
#include <stdint.h>
#include <stdio.h>
#include <math.h>
#include <gmp.h>

/* Largest non-negative integer n for which HarmonicNumber[n, r] expands to the
 * explicit finite sum Sum_{i=1}^n i^-r.  An exact integer argument is never
 * numerically contaminated, so beyond the cap the result simply stays symbolic
 * (there is no cheaper numeric fallback to take). */
#define HN_EXPAND_CAP 100000L

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
    fprintf(stderr,
            "HarmonicNumber::argt: HarmonicNumber called with %zu argument%s; "
            "1 or 2 arguments are expected.\n",
            argc, argc == 1 ? "" : "s");
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
            SymbolDef* d = symtab_lookup(e->data.symbol);
            return d && (d->attributes & ATTR_CONSTANT);
        }
        case EXPR_FUNCTION: {
            const Expr* h = e->data.function.head;
            if (h->type == EXPR_SYMBOL &&
                (strcmp(h->data.symbol, "Rational") == 0 ||
                 strcmp(h->data.symbol, "Complex") == 0)) {
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
    Expr* np1 = hn_fn2("Plus", expr_copy((Expr*)n), hn_int(1));

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

    /* Force numericization (so Constant n such as E resolve too).  Use the
     * input precision: machine -> N[reduction]; MPFR -> N[reduction, digits]. */
    long bn = numeric_min_inexact_bits(n);
    long br = numeric_min_inexact_bits(r);
    long bits = 0;
    if (bn > 0) bits = bn;
    if (br > 0 && (bits == 0 || br < bits)) bits = br;

    Expr* ncall;
    if (bits > 53) {
        long digits = (long)numeric_bits_to_digits(bits);
        if (digits < 1) digits = 1;
        ncall = hn_fn2("N", reduction, hn_int(digits));
    } else {
        ncall = hn_fn1("N", reduction);
    }
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
    } else if (hn_nonneg_int(n, &ni) && ni <= HN_EXPAND_CAP) {
        /* n a non-negative integer (within the cap): explicit finite sum. */
        result = hn_finite_sum(ni, r);
    } else if (hn_exact_int(r, &ri) && ri <= 0) {
        /* r a non-positive integer: Faulhaber polynomial in n. */
        result = hn_faulhaber(n, -ri);
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
