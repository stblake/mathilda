/* beta.c -- the Euler beta function and its incomplete / generalized forms.
 *
 *   Beta[a, b]         = Gamma(a) Gamma(b) / Gamma(a+b)
 *   Beta[z, a, b]      = Int_0^z t^(a-1) (1-t)^(b-1) dt          (incomplete)
 *   Beta[z0, z1, a, b] = Beta[z1, a, b] - Beta[z0, a, b]         (generalized)
 *
 * Rather than re-deriving the transcendental machinery, Beta is assembled from
 * functions that already exist and are well tested:
 *
 *   - Beta[a, b] is reduced to a ratio of Gamma calls and handed back to the
 *     evaluator. The existing Gamma builtin already closes exact integers,
 *     half-integers, rationals (-> Sqrt[Pi] forms), machine / arbitrary-
 *     precision reals, and complex arguments, so Beta inherits all of that for
 *     free. Non-positive-integer poles (where a, b, or a+b hits a gamma pole)
 *     are detected up front: a surviving pole gives ComplexInfinity, while a
 *     cancelling pair of poles reduces by the finite limit of the gamma ratio.
 *
 *   - Beta[z, a, b] is reduced through Hypergeometric2F1 via
 *       B_z(a, b) = z^a / a * 2F1(a, 1-b; a+1; z),
 *     which terminates to an exact closed form when b is a positive integer
 *     (1-b a non-positive integer) and evaluates numerically (real or complex)
 *     otherwise.
 *
 * Memory: builtin_beta takes ownership of res. It returns a freshly built tree
 * (the evaluator frees res) or NULL to leave the call unevaluated (the
 * evaluator keeps res). It never frees res itself. Every intermediate tree is
 * consumed by eval_and_free. */

#include "beta.h"

#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>

#include <gmp.h>

#include "arithmetic.h"   /* is_rational, make_rational */
#include "attr.h"
#include "eval.h"          /* eval_and_free */
#include "expr.h"          /* expr_is_integer_like, expr_to_mpz, expr_is_numeric_like */
#include "symtab.h"

/* ------------------------------------------------------------------ */
/* Small symbolic builders                                            */
/* ------------------------------------------------------------------ */

static Expr* mk_sym(const char* s) { return expr_new_symbol(s); }

static Expr* fn1(const char* head, Expr* a) {
    return expr_new_function(mk_sym(head), (Expr*[]){ a }, 1);
}
static Expr* fn2(const char* head, Expr* a, Expr* b) {
    return expr_new_function(mk_sym(head), (Expr*[]){ a, b }, 2);
}
static Expr* fn4(const char* head, Expr* a, Expr* b, Expr* c, Expr* d) {
    return expr_new_function(mk_sym(head), (Expr*[]){ a, b, c, d }, 4);
}

/* Normalized integer leaf from a GMP value (demotes to int64 when it fits). */
static Expr* mpz_to_expr(const mpz_t z) {
    return expr_bigint_normalize(expr_new_bigint_from_mpz(z));
}

/* Gamma[copy of x]. */
static Expr* gamma_of(const Expr* x) { return fn1("Gamma", expr_copy((Expr*)x)); }

/* ------------------------------------------------------------------ */
/* Integer predicates                                                 */
/* ------------------------------------------------------------------ */

/* True when e is an exact integer equal to v. */
static bool int_eq(const Expr* e, long v) {
    if (!expr_is_integer_like(e)) return false;
    mpz_t z; mpz_init(z); expr_to_mpz(e, z);
    bool r = (mpz_cmp_si(z, v) == 0);
    mpz_clear(z);
    return r;
}

/* True when e is an exact positive integer. */
static bool pos_int(const Expr* e) {
    if (!expr_is_integer_like(e)) return false;
    mpz_t z; mpz_init(z); expr_to_mpz(e, z);
    bool r = (mpz_sgn(z) > 0);
    mpz_clear(z);
    return r;
}

/* ------------------------------------------------------------------ */
/* Beta[a, b]: poles and the gamma-ratio reduction                    */
/* ------------------------------------------------------------------ */

/* Classify the pole structure of Beta[a, b] on the integer lattice.
 *
 *   Beta[a, b] = Gamma(a) Gamma(b) / Gamma(a+b).
 * Gamma has simple poles at the non-positive integers. Counting pole order
 *   p = [a in Z<=0] + [b in Z<=0] - [a+b in Z<=0]
 * gives three regimes:
 *   p  > 0 : a surviving pole            -> ComplexInfinity
 *   p == 0 : a cancelling pair of poles  -> finite limit of the gamma ratio
 *   p  < 0 : only the denominator is a pole; the ratio -> 0, which the Gamma
 *            route already produces, so this is left to the caller (returns 0).
 *
 * For the cancelling case exactly one of a, b is a non-positive integer (if
 * both were, p would be >= 1); call it -m, and the other argument B is then a
 * positive integer. With a+b = -k the limit is
 *   Beta = (-1)^B * k! (B-1)! / m!.
 *
 * Returns 1 and sets *out when the result is determined (ComplexInfinity or the
 * finite value); returns 0 (leaving *out untouched) when a, b are not a lattice
 * special and the caller should continue with the numeric reduction. */
static int beta_lattice_case(const Expr* a, const Expr* b, Expr** out) {
    bool ai = expr_is_integer_like(a);
    bool bi = expr_is_integer_like(b);

    mpz_t A, B;
    mpz_init(A); mpz_init(B);
    if (ai) expr_to_mpz(a, A);
    if (bi) expr_to_mpz(b, B);

    bool a_nonpos = ai && mpz_sgn(A) <= 0;
    bool b_nonpos = bi && mpz_sgn(B) <= 0;

    if (!a_nonpos && !b_nonpos) {
        mpz_clear(A); mpz_clear(B);
        return 0;                       /* no numerator pole: continue */
    }

    /* a+b is a non-positive integer only when both a and b are integers. */
    bool s_nonpos = false;
    mpz_t S; mpz_init(S);
    if (ai && bi) { mpz_add(S, A, B); s_nonpos = mpz_sgn(S) <= 0; }

    int p = (a_nonpos ? 1 : 0) + (b_nonpos ? 1 : 0) - (s_nonpos ? 1 : 0);

    if (p > 0) {
        *out = mk_sym("ComplexInfinity");
        mpz_clear(A); mpz_clear(B); mpz_clear(S);
        return 1;
    }

    /* p == 0: cancelling poles. Identify the non-positive arg (-m) and the
     * other positive-integer arg B; k = -(a+b). */
    mpz_t m, k, bo;
    mpz_init(m); mpz_init(k); mpz_init(bo);
    if (a_nonpos) { mpz_neg(m, A); mpz_set(bo, B); }
    else          { mpz_neg(m, B); mpz_set(bo, A); }
    mpz_neg(k, S);                      /* k = -(a+b) */

    /* Beta = (-1)^B * k! (B-1)! / m!. */
    Expr* sign = expr_new_integer(mpz_odd_p(bo) ? -1 : 1);
    Expr* kE   = mpz_to_expr(k);
    mpz_t bo1; mpz_init(bo1); mpz_sub_ui(bo1, bo, 1);
    Expr* bo1E = mpz_to_expr(bo1);
    Expr* mE   = mpz_to_expr(m);

    Expr* form = expr_new_function(mk_sym("Times"),
        (Expr*[]){ sign,
                   fn1("Factorial", kE),
                   fn1("Factorial", bo1E),
                   fn2("Power", fn1("Factorial", mE), expr_new_integer(-1)) }, 4);
    *out = eval_and_free(form);

    mpz_clear(bo1);
    mpz_clear(m); mpz_clear(k); mpz_clear(bo);
    mpz_clear(A); mpz_clear(B); mpz_clear(S);
    return 1;
}

/* Beta[a, b] = Gamma(a) Gamma(b) / Gamma(a+b). */
static Expr* beta_two_arg(Expr* a, Expr* b) {
    Expr* lattice = NULL;
    if (beta_lattice_case(a, b, &lattice)) return lattice;

    /* Reduce only when both arguments are numeric quantities (exact or
     * inexact); a symbolic argument keeps Beta[a, b] unevaluated, matching
     * the Wolfram Language. */
    if (!expr_is_numeric_like(a) || !expr_is_numeric_like(b)) return NULL;

    /* When one argument is a positive integer n, the gamma ratio collapses to
     *   Beta[n, b] = (n-1)! / Pochhammer[b, n],
     * a finite product that the evaluator reduces to an exact value even when
     * the other argument is rational (e.g. Beta[3, 1/3] = 27/14). The plain
     * Gamma route would instead leave an unsimplified ratio of unrelated gamma
     * values such as Gamma[1/3]/Gamma[10/3]. */
    if (pos_int(a)) {
        Expr* nm1  = fn2("Plus", expr_copy(a), expr_new_integer(-1));
        Expr* poch = fn2("Pochhammer", expr_copy(b), expr_copy(a));
        return eval_and_free(fn2("Times", fn1("Factorial", nm1),
                                 fn2("Power", poch, expr_new_integer(-1))));
    }
    if (pos_int(b)) {
        Expr* nm1  = fn2("Plus", expr_copy(b), expr_new_integer(-1));
        Expr* poch = fn2("Pochhammer", expr_copy(a), expr_copy(b));
        return eval_and_free(fn2("Times", fn1("Factorial", nm1),
                                 fn2("Power", poch, expr_new_integer(-1))));
    }

    Expr* ab  = fn2("Plus", expr_copy(a), expr_copy(b));
    Expr* form = expr_new_function(mk_sym("Times"),
        (Expr*[]){ gamma_of(a), gamma_of(b),
                   fn2("Power", fn1("Gamma", ab), expr_new_integer(-1)) }, 3);
    return eval_and_free(form);
}

/* ------------------------------------------------------------------ */
/* Beta[z, a, b]: incomplete beta via 2F1                             */
/* ------------------------------------------------------------------ */

/* B_z(a, b) = z^a / a * Hypergeometric2F1[a, 1-b, a+1, z]. */
static Expr* beta_incomplete(Expr* z, Expr* a, Expr* b) {
    /* 1-b and a+1. */
    Expr* one_minus_b = fn2("Subtract", expr_new_integer(1), expr_copy(b));
    Expr* a_plus_1    = fn2("Plus", expr_copy(a), expr_new_integer(1));
    Expr* hyp = fn4("Hypergeometric2F1",
                    expr_copy(a), one_minus_b, a_plus_1, expr_copy(z));

    Expr* form = expr_new_function(mk_sym("Times"),
        (Expr*[]){ fn2("Power", expr_copy(z), expr_copy(a)),  /* z^a   */
                   fn2("Power", expr_copy(a), expr_new_integer(-1)), /* 1/a */
                   hyp }, 3);
    return eval_and_free(form);
}

static Expr* beta_three_arg(Expr* z, Expr* a, Expr* b) {
    if (int_eq(z, 0)) return expr_new_integer(0);            /* B_0(a,b) = 0 */
    if (int_eq(z, 1)) {                                      /* B_1(a,b) = B(a,b) */
        /* Hand back the two-argument form for the evaluator to reduce (it
         * stays Beta[a, b] when a, b are symbolic, matching the Wolfram
         * Language). */
        return expr_new_function(mk_sym("Beta"),
                   (Expr*[]){ expr_copy(a), expr_copy(b) }, 2);
    }

    /* Reduce when the call is fully numeric, or when b is a positive integer
     * (the 2F1 series then terminates to an exact closed form even for
     * symbolic z and a). Otherwise stay symbolic. */
    bool numeric = expr_is_numeric_like(z) &&
                   expr_is_numeric_like(a) && expr_is_numeric_like(b);
    if (!numeric && !pos_int(b)) return NULL;

    return beta_incomplete(z, a, b);
}

/* ------------------------------------------------------------------ */
/* Beta[z0, z1, a, b]: generalized incomplete beta                    */
/* ------------------------------------------------------------------ */

static Expr* beta_four_arg(Expr* z0, Expr* z1, Expr* a, Expr* b) {
    /* Keep the four-argument form symbolic unless every argument is numeric;
     * the difference Beta[z1,a,b] - Beta[z0,a,b] is then a concrete number, and
     * symbolic differentiation is handled directly in deriv.c. */
    if (!expr_is_numeric_like(z0) || !expr_is_numeric_like(z1) ||
        !expr_is_numeric_like(a)  || !expr_is_numeric_like(b)) {
        return NULL;
    }
    Expr* hi = expr_new_function(mk_sym("Beta"),
                   (Expr*[]){ expr_copy(z1), expr_copy(a), expr_copy(b) }, 3);
    Expr* lo = expr_new_function(mk_sym("Beta"),
                   (Expr*[]){ expr_copy(z0), expr_copy(a), expr_copy(b) }, 3);
    return eval_and_free(fn2("Subtract", hi, lo));
}

/* ------------------------------------------------------------------ */
/* Builtin entry point                                                */
/* ------------------------------------------------------------------ */

/* Wolfram-compatible argb diagnostic for a wrong argument count. */
static Expr* beta_emit_argb(size_t argc) {
    fprintf(stderr,
            "Beta::argb: Beta called with %zu argument%s; "
            "between 2 and 4 arguments are expected.\n",
            argc, argc == 1 ? "" : "s");
    return NULL;
}

Expr* builtin_beta(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    Expr** args = res->data.function.args;

    switch (argc) {
        case 2: return beta_two_arg(args[0], args[1]);
        case 3: return beta_three_arg(args[0], args[1], args[2]);
        case 4: return beta_four_arg(args[0], args[1], args[2], args[3]);
        default: return beta_emit_argb(argc);
    }
}

/* ------------------------------------------------------------------ */
/* Registration                                                       */
/* ------------------------------------------------------------------ */

void beta_init(void) {
    symtab_add_builtin("Beta", builtin_beta);
    symtab_get_def("Beta")->attributes |=
        (ATTR_LISTABLE | ATTR_NUMERICFUNCTION | ATTR_PROTECTED | ATTR_READPROTECTED);
    /* Docstring lives in info.c (info_init). */
}
