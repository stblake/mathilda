/* Mathilda -- Lucas numbers and Lucas polynomials.
 *
 * Surface forms
 * -------------
 *   LucasL[n]      L_n, the nth Lucas number.
 *   LucasL[n, x]   L_n(x), the nth Lucas polynomial.
 *
 * Evaluation strategy
 * -------------------
 * The builtin follows the same two-tier philosophy as Fibonacci (and the
 * rest of the system): exact arithmetic is computed directly in C, while
 * inexact / symbolic-constant reductions are expressed as expression trees
 * and handed back to the evaluator (and, for numeric requests, to
 * `numericalize`).
 *
 *   * Exact integer order n:
 *       - LucasL[n]      : GMP fast-doubling of the Fibonacci pair
 *                          (F_m, F_{m+1}), then L_m = 2 F_{m+1} - F_m,
 *                          O(log n) big-integer math. Negative orders via
 *                          L_{-m} = (-1)^m L_m.
 *       - LucasL[n, x]   : the recurrence L_k = x L_{k-1} + L_{k-2} with
 *                          L_0 = 2, L_1 = x, evaluated (Expand-ed) at each
 *                          step so the partial result stays a canonical
 *                          (and, for numeric x, fully reduced) expression.
 *
 *   * Inexact / non-integer order (Real, MPFR, or Complex with an inexact
 *     part) -- the generalized closed forms, built symbolically and then
 *     numericalized at the precision carried by the inputs:
 *
 *         LucasL[n]    = phi^n + Cos[Pi n] phi^-n,  phi = GoldenRatio.
 *         LucasL[n, x] = beta^n + Cos[Pi n] beta^-n,
 *                        beta = (x + Sqrt[x^2 + 4]) / 2.
 *
 *     Because GoldenRatio / Pi / Cos / Sqrt / Power already carry numeric
 *     paths (machine and MPFR), `numericalize` drives the whole reduction;
 *     this also yields complex results for complex arguments for free.
 *
 *   * Anything else (purely symbolic order, exact non-integer order with no
 *     N applied) -> NULL, leaving the call unevaluated.
 *
 * Memory: the builtin honours the ownership contract -- it never frees `res`,
 * returns a fresh Expr* on success or NULL otherwise, and clears every
 * GMP / temporary tree it allocates.
 */

#include "lucas.h"
#include "eval.h"
#include "symtab.h"
#include "numeric.h"

#include <gmp.h>
#include <stdbool.h>
#include <stdint.h>
#include <string.h>

/* ---------------------------------------------------------------------- */
/* Tiny expression builders                                                */
/* ---------------------------------------------------------------------- */

static Expr* mk_int(int64_t v) { return expr_new_integer(v); }
static Expr* mk_sym(const char* s) { return expr_new_symbol(s); }

static Expr* mk_fn1(const char* name, Expr* a) {
    Expr* args[1] = { a };
    return expr_new_function(mk_sym(name), args, 1);
}

static Expr* mk_fn2(const char* name, Expr* a, Expr* b) {
    Expr* args[2] = { a, b };
    return expr_new_function(mk_sym(name), args, 2);
}

/* ---------------------------------------------------------------------- */
/* Numeric classification                                                  */
/* ---------------------------------------------------------------------- */

/* True iff `e` carries an inexact numeric value anywhere relevant: a
 * machine Real, an MPFR real, or a Complex[...] whose parts do. Exact
 * leaves (Integer, BigInt, Rational) are deliberately excluded -- WL keeps
 * those symbolic until N is applied. */
static bool contains_inexact(const Expr* e) {
    if (!e) return false;
    if (e->type == EXPR_REAL) return true;
#ifdef USE_MPFR
    if (e->type == EXPR_MPFR) return true;
#endif
    if (e->type == EXPR_FUNCTION &&
        e->data.function.head->type == EXPR_SYMBOL &&
        strcmp(e->data.function.head->data.symbol.name, "Complex") == 0) {
        for (size_t i = 0; i < e->data.function.arg_count; i++) {
            if (contains_inexact(e->data.function.args[i])) return true;
        }
    }
    return false;
}

/* Precision spec matching the inputs: MPFR (at the combined bit width) when
 * any operand carries MPFR precision, otherwise machine double. `b` may be
 * NULL for the one-argument form. */
static NumericSpec lucas_spec(const Expr* a, const Expr* b) {
#ifdef USE_MPFR
    if (numeric_any_mpfr(a, b)) {
        NumericSpec s;
        s.mode = NUMERIC_MODE_MPFR;
        s.bits = numeric_combined_bits(a, b, 0);
        return s;
    }
#endif
    (void)a; (void)b;
    return numeric_machine_spec();
}

/* ---------------------------------------------------------------------- */
/* Exact integer Lucas numbers (fast doubling)                             */
/* ---------------------------------------------------------------------- */

/* Compute the Fibonacci pair (F(m), F(m+1)) for m >= 0 using the
 * fast-doubling identities
 *     F(2k)   = F(k) [2 F(k+1) - F(k)]
 *     F(2k+1) = F(k+1)^2 + F(k)^2
 * scanning the bits of m from most- to least-significant. `fm` and `fm1`
 * must be live mpz_t (they receive F(m) and F(m+1)); `m` is read-only. */
static void fib_pair_mpz(mpz_t fm, mpz_t fm1, const mpz_t m) {
    mpz_t a, b, c, d, t;
    mpz_inits(a, b, c, d, t, NULL);
    mpz_set_ui(a, 0);   /* F(0) */
    mpz_set_ui(b, 1);   /* F(1) */

    size_t bits = mpz_sizeinbase(m, 2);   /* >= 1 even for m == 0 */
    for (size_t i = bits; i-- > 0; ) {
        /* c = F(2k) = a (2b - a) */
        mpz_mul_2exp(t, b, 1);
        mpz_sub(t, t, a);
        mpz_mul(c, a, t);
        /* d = F(2k+1) = a^2 + b^2 */
        mpz_mul(t, a, a);
        mpz_mul(d, b, b);
        mpz_add(d, d, t);
        if (mpz_tstbit(m, i)) {
            mpz_set(a, d);          /* F(2k+1) */
            mpz_add(b, c, d);       /* F(2k+2) */
        } else {
            mpz_set(a, c);          /* F(2k)   */
            mpz_set(b, d);          /* F(2k+1) */
        }
    }
    mpz_set(fm, a);     /* F(m)   */
    mpz_set(fm1, b);    /* F(m+1) */
    mpz_clears(a, b, c, d, t, NULL);
}

/* out <- L(m) for m >= 0, via L(m) = 2 F(m+1) - F(m). `out` must be a live
 * mpz_t; `m` is read-only. */
static void lucas_mpz(mpz_t out, const mpz_t m) {
    mpz_t fm, fm1;
    mpz_inits(fm, fm1, NULL);
    fib_pair_mpz(fm, fm1, m);
    mpz_mul_2exp(out, fm1, 1);   /* 2 F(m+1) */
    mpz_sub(out, out, fm);       /* - F(m)   */
    mpz_clears(fm, fm1, NULL);
}

/* LucasL[n] for an exact integer-like `n` (Integer or BigInt). */
static Expr* lucas_number(const Expr* n) {
    mpz_t nv, m, l;
    expr_to_mpz(n, nv);              /* inits nv */
    mpz_init(m);
    mpz_abs(m, nv);
    mpz_init(l);
    lucas_mpz(l, m);
    /* L_{-m} = (-1)^m L_m: negate when n < 0 and m is odd. */
    if (mpz_sgn(nv) < 0 && mpz_odd_p(m)) mpz_neg(l, l);
    Expr* r = expr_bigint_normalize(expr_new_bigint_from_mpz(l));
    mpz_clears(nv, m, l, NULL);
    return r;
}

/* ---------------------------------------------------------------------- */
/* Lucas polynomials (exact integer order)                                 */
/* ---------------------------------------------------------------------- */

/* L_m(x) for m >= 0 via the recurrence L_k = x L_{k-1} + L_{k-2}, with
 * L_0 = 2 and L_1 = x. Each step is Expand-ed so the running value stays a
 * flat, canonical polynomial in `x` (and collapses to a number when `x` is
 * numeric). `x` is read-only. */
static Expr* lucas_poly(long m, const Expr* x) {
    if (m == 0) return mk_int(2);
    if (m == 1) return expr_copy((Expr*)x);

    Expr* lm2 = mk_int(2);              /* L_0 */
    Expr* lm1 = expr_copy((Expr*)x);    /* L_1 */
    for (long k = 2; k <= m; k++) {
        Expr* term = mk_fn2("Times", expr_copy((Expr*)x), expr_copy(lm1));
        Expr* sum  = mk_fn2("Plus", term, expr_copy(lm2));
        Expr* lk   = eval_and_free(mk_fn1("Expand", sum));
        expr_free(lm2);
        lm2 = lm1;
        lm1 = lk;
    }
    expr_free(lm2);
    return lm1;   /* L_m(x) */
}

/* LucasL[n, x] for an exact integer-like `n`, returning the polynomial in
 * `x` (or a number if `x` is numeric). NULL if the order does not fit in a
 * machine long (degree too large to materialise). */
static Expr* lucas_polynomial(const Expr* n, const Expr* x) {
    mpz_t nv;
    expr_to_mpz(n, nv);             /* inits nv */
    if (!mpz_fits_slong_p(nv)) { mpz_clear(nv); return NULL; }
    long ni = mpz_get_si(nv);
    mpz_clear(nv);

    long m = ni < 0 ? -ni : ni;
    Expr* poly = lucas_poly(m, x);
    /* L_{-m}(x) = (-1)^m L_m(x): negate when n < 0 and m is odd. */
    if (ni < 0 && (m % 2) != 0) {
        poly = eval_and_free(mk_fn2("Times", mk_int(-1), poly));
    }
    return poly;
}

/* ---------------------------------------------------------------------- */
/* Generalized closed forms (inexact / non-integer order)                  */
/* ---------------------------------------------------------------------- */

/* Build base^n + Cos[Pi n] base^-n as an expression tree, taking ownership
 * of the freshly built `base` and a private copy of the numeric order `n`.
 * Caller owns the result. */
static Expr* closed_form(Expr* base, const Expr* n) {
    Expr* pos    = mk_fn2("Power", expr_copy(base), expr_copy((Expr*)n));
    Expr* neg    = mk_fn2("Power", base,
                          mk_fn2("Times", mk_int(-1), expr_copy((Expr*)n)));
    Expr* cosarg = mk_fn1("Cos", mk_fn2("Times", mk_sym("Pi"),
                                        expr_copy((Expr*)n)));
    return mk_fn2("Plus", pos, mk_fn2("Times", cosarg, neg));
}

/* LucasL[n] for inexact / complex order: phi^n + Cos[Pi n] phi^-n. */
static Expr* lucas_number_numeric(const Expr* n) {
    Expr* form = closed_form(mk_sym("GoldenRatio"), n);
    Expr* r = numericalize(form, lucas_spec(n, NULL));
    expr_free(form);
    return r;
}

/* LucasL[n, x] for inexact order with numeric x:
 *     beta = (x + Sqrt[x^2 + 4]) / 2,
 *     beta^n + Cos[Pi n] beta^-n. */
static Expr* lucas_polynomial_numeric(const Expr* n, const Expr* x) {
    Expr* x2p4 = mk_fn2("Plus", mk_fn2("Power", expr_copy((Expr*)x), mk_int(2)),
                        mk_int(4));
    Expr* sqrtd = mk_fn1("Sqrt", x2p4);                 /* Sqrt[x^2 + 4] */
    Expr* beta  = mk_fn2("Times",
                         mk_fn2("Plus", expr_copy((Expr*)x), sqrtd),
                         mk_fn2("Power", mk_int(2), mk_int(-1)));
    Expr* form  = closed_form(beta, n);
    Expr* r = numericalize(form, lucas_spec(n, x));
    expr_free(form);
    return r;
}

/* ---------------------------------------------------------------------- */
/* Builtin entry point                                                     */
/* ---------------------------------------------------------------------- */

Expr* builtin_lucasl(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;

    if (argc == 1) {
        Expr* n = res->data.function.args[0];
        if (expr_is_integer_like(n))  return lucas_number(n);
        if (contains_inexact(n))      return lucas_number_numeric(n);
        return NULL;                  /* symbolic / exact non-integer */
    }

    if (argc == 2) {
        Expr* n = res->data.function.args[0];
        Expr* x = res->data.function.args[1];
        if (expr_is_integer_like(n))  return lucas_polynomial(n, x);
        if (contains_inexact(n) && expr_is_numeric_like(x))
            return lucas_polynomial_numeric(n, x);
        return NULL;                  /* symbolic order, or symbolic x */
    }

    return NULL;
}

/* ---------------------------------------------------------------------- */
/* Registration                                                            */
/* ---------------------------------------------------------------------- */

void lucas_init(void) {
    symtab_add_builtin("LucasL", builtin_lucasl);
}
