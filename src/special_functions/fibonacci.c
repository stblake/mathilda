/* Mathilda -- Fibonacci numbers and Fibonacci polynomials.
 *
 * Surface forms
 * -------------
 *   Fibonacci[n]      F_n, the nth Fibonacci number.
 *   Fibonacci[n, x]   F_n(x), the nth Fibonacci polynomial.
 *
 * Evaluation strategy
 * -------------------
 * The builtin follows the same two-tier philosophy as the rest of the
 * system: exact arithmetic is computed directly in C, while inexact /
 * symbolic-constant reductions are expressed as expression trees and handed
 * back to the evaluator (and, for numeric requests, to `numericalize`).
 *
 *   * Exact integer order n:
 *       - Fibonacci[n]      : GMP fast-doubling, O(log n) big-integer math.
 *                             Negative orders via F_{-m} = (-1)^{m+1} F_m.
 *       - Fibonacci[n, x]   : the recurrence F_k = x F_{k-1} + F_{k-2},
 *                             evaluated at each step so the partial result
 *                             stays a canonical (and, for numeric x, fully
 *                             reduced) expression.
 *
 *   * Inexact / non-integer order (Real, MPFR, or Complex with an inexact
 *     part) -- the generalized closed forms, built symbolically and then
 *     numericalized at the precision carried by the inputs:
 *
 *         Fibonacci[n]    = (phi^n - Cos[Pi n] phi^-n) / Sqrt[5],
 *                           phi = GoldenRatio.
 *         Fibonacci[n, x] = (beta^n - Cos[Pi n] beta^-n) / Sqrt[x^2 + 4],
 *                           beta = (x + Sqrt[x^2 + 4]) / 2.
 *
 *     Because GoldenRatio / Pi / Cos / Sqrt / Power already carry numeric
 *     paths (machine and MPFR), `numericalize` drives the whole reduction;
 *     this also yields complex results for complex arguments for free.
 *
 *   * Fibonacci[n, x] with an exact non-integer order n (e.g. a Rational)
 *     and an exact numeric x evaluates the SAME closed form exactly, but
 *     only keeps the result when it collapses to a number -- so
 *     Fibonacci[1/2, 0] -> 1/2 while Fibonacci[1/2, x] (symbolic x) and a
 *     non-collapsing exact x stay unevaluated, matching the one-argument
 *     convention that exact non-integer orders stay symbolic until they
 *     reduce to a value.
 *
 *   * Anything else (purely symbolic order, one-argument exact non-integer
 *     order with no N applied) -> NULL, leaving the call unevaluated.
 *
 * Memory: the builtin honours the ownership contract -- it never frees `res`,
 * returns a fresh Expr* on success or NULL otherwise, and clears every
 * GMP / temporary tree it allocates.
 */

#include "fibonacci.h"
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
static NumericSpec fib_spec(const Expr* a, const Expr* b) {
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
/* Exact integer Fibonacci numbers (fast doubling)                         */
/* ---------------------------------------------------------------------- */

/* out <- F_m for m >= 0, using the fast-doubling identities
 *     F(2k)   = F(k) [2 F(k+1) - F(k)]
 *     F(2k+1) = F(k+1)^2 + F(k)^2
 * scanning the bits of m from most- to least-significant. `out` must be a
 * live mpz_t; `m` is read-only. */
static void fib_mpz(mpz_t out, const mpz_t m) {
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
    mpz_set(out, a);
    mpz_clears(a, b, c, d, t, NULL);
}

/* Fibonacci[n] for an exact integer-like `n` (Integer or BigInt). */
static Expr* fib_number(const Expr* n) {
    mpz_t nv, m, f;
    expr_to_mpz(n, nv);              /* inits nv */
    mpz_init(m);
    mpz_abs(m, nv);
    mpz_init(f);
    fib_mpz(f, m);
    /* F_{-m} = (-1)^{m+1} F_m: negate when n < 0 and m is even. */
    if (mpz_sgn(nv) < 0 && mpz_even_p(m)) mpz_neg(f, f);
    Expr* r = expr_bigint_normalize(expr_new_bigint_from_mpz(f));
    mpz_clears(nv, m, f, NULL);
    return r;
}

/* ---------------------------------------------------------------------- */
/* Fibonacci polynomials (exact integer order)                             */
/* ---------------------------------------------------------------------- */

/* F_m(x) for m >= 0 via the recurrence F_k = x F_{k-1} + F_{k-2}. Each step
 * is Expand-ed so the running value stays a flat, canonical polynomial in `x`
 * (and collapses to a number when `x` is numeric). `x` is read-only. */
static Expr* fib_poly(long m, const Expr* x) {
    if (m == 0) return mk_int(0);
    if (m == 1) return mk_int(1);

    Expr* fm2 = mk_int(0);   /* F_0 */
    Expr* fm1 = mk_int(1);   /* F_1 */
    for (long k = 2; k <= m; k++) {
        Expr* term = mk_fn2("Times", expr_copy((Expr*)x), expr_copy(fm1));
        Expr* sum  = mk_fn2("Plus", term, expr_copy(fm2));
        Expr* fk   = eval_and_free(mk_fn1("Expand", sum));
        expr_free(fm2);
        fm2 = fm1;
        fm1 = fk;
    }
    expr_free(fm2);
    return fm1;   /* F_m(x) */
}

/* Fibonacci[n, x] for an exact integer-like `n`, returning the polynomial in
 * `x` (or a number if `x` is numeric). NULL if the order does not fit in a
 * machine long (degree too large to materialise). */
static Expr* fib_polynomial(const Expr* n, const Expr* x) {
    mpz_t nv;
    expr_to_mpz(n, nv);             /* inits nv */
    if (!mpz_fits_slong_p(nv)) { mpz_clear(nv); return NULL; }
    long ni = mpz_get_si(nv);
    mpz_clear(nv);

    long m = ni < 0 ? -ni : ni;
    Expr* poly = fib_poly(m, x);
    /* F_{-m}(x) = (-1)^{m+1} F_m(x): negate when n < 0 and m is even. */
    if (ni < 0 && (m % 2) == 0) {
        poly = eval_and_free(mk_fn2("Times", mk_int(-1), poly));
    }
    return poly;
}

/* ---------------------------------------------------------------------- */
/* Generalized closed forms (inexact / non-integer order)                  */
/* ---------------------------------------------------------------------- */

/* Build (base^n - Cos[Pi n] base^-n) / denom as an expression tree, taking
 * ownership of the freshly built `base`, `denom` and a private copy of the
 * numeric order `n`. Caller owns the result. */
static Expr* closed_form(Expr* base, const Expr* n, Expr* denom) {
    Expr* pos    = mk_fn2("Power", expr_copy(base), expr_copy((Expr*)n));
    Expr* neg    = mk_fn2("Power", base,
                          mk_fn2("Times", mk_int(-1), expr_copy((Expr*)n)));
    Expr* cosarg = mk_fn1("Cos", mk_fn2("Times", mk_sym("Pi"),
                                        expr_copy((Expr*)n)));
    Expr* numer  = mk_fn2("Plus", pos,
                          mk_fn2("Times", mk_int(-1),
                                 mk_fn2("Times", cosarg, neg)));
    return mk_fn2("Times", numer, mk_fn2("Power", denom, mk_int(-1)));
}

/* Fibonacci[n] for inexact / complex order: (phi^n - Cos[Pi n] phi^-n)/Sqrt[5]. */
static Expr* fib_number_numeric(const Expr* n) {
    Expr* form = closed_form(mk_sym("GoldenRatio"), n, mk_fn1("Sqrt", mk_int(5)));
    Expr* r = numericalize(form, fib_spec(n, NULL));
    expr_free(form);
    return r;
}

/* Fibonacci[n, x] via the generalized closed form, for a non-integer numeric
 * order n and a numeric x:
 *     beta = (x + Sqrt[x^2 + 4]) / 2,
 *     (beta^n - Cos[Pi n] beta^-n) / Sqrt[x^2 + 4].
 *
 * When any input is inexact the form is numericalized at the inputs'
 * precision (e.g. Fibonacci[1/2, 3.2] -> 0.4948...). When every input is
 * exact the form is evaluated exactly and the result is kept only if it
 * collapses to a number -- so Fibonacci[1/2, 0] -> 1/2, while exact inputs
 * that would leave a radical stay unevaluated (NULL), mirroring the
 * one-argument rule that exact non-integer orders stay symbolic. */
static Expr* fib_polynomial_closed(const Expr* n, const Expr* x, bool inexact) {
    Expr* x2p4 = mk_fn2("Plus", mk_fn2("Power", expr_copy((Expr*)x), mk_int(2)),
                        mk_int(4));
    Expr* sqrtd = mk_fn1("Sqrt", x2p4);                 /* Sqrt[x^2 + 4] */
    Expr* beta  = mk_fn2("Times",
                         mk_fn2("Plus", expr_copy((Expr*)x), expr_copy(sqrtd)),
                         mk_fn2("Power", mk_int(2), mk_int(-1)));
    Expr* form  = closed_form(beta, n, sqrtd);          /* adopts sqrtd as denom */

    if (inexact) {
        Expr* r = numericalize(form, fib_spec(n, x));
        expr_free(form);
        return r;
    }

    /* Exact path: evaluate, but accept only a numeric result. */
    Expr* r = eval_and_free(form);
    if (expr_is_numeric_like(r)) return r;
    expr_free(r);
    return NULL;
}

/* ---------------------------------------------------------------------- */
/* Builtin entry point                                                     */
/* ---------------------------------------------------------------------- */

Expr* builtin_fibonacci(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;

    if (argc == 1) {
        Expr* n = res->data.function.args[0];
        if (expr_is_integer_like(n))  return fib_number(n);
        if (contains_inexact(n))      return fib_number_numeric(n);
        return NULL;                  /* symbolic / exact non-integer */
    }

    if (argc == 2) {
        Expr* n = res->data.function.args[0];
        Expr* x = res->data.function.args[1];
        if (expr_is_integer_like(n))  return fib_polynomial(n, x);
        /* Non-integer numeric order with numeric x -> generalized closed
         * form (numeric when any input is inexact, exact-then-collapse
         * otherwise). Symbolic order or symbolic x stays unevaluated. */
        if (expr_is_numeric_like(n) && expr_is_numeric_like(x))
            return fib_polynomial_closed(n, x, contains_inexact(n) ||
                                               contains_inexact(x));
        return NULL;                  /* symbolic order, or symbolic x */
    }

    return NULL;
}

/* ---------------------------------------------------------------------- */
/* Registration                                                            */
/* ---------------------------------------------------------------------- */

void fibonacci_init(void) {
    symtab_add_builtin("Fibonacci", builtin_fibonacci);
}
