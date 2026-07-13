/* coprimeq.c -- CoprimeQ[] (incl. Gaussian-integer coprimality).
 * Split from numbertheory.c; see numbertheory.h and
 * numbertheory_internal.h for the subsystem layout. */

#include "numbertheory.h"
#include "numbertheory_internal.h"
#include "arithmetic.h"
#include "eval.h"
#include "sym_names.h"
#include "internal.h"
#include "print.h"
#include "symtab.h"
#include "attr.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <gmp.h>

/* ---- CoprimeQ ----------------------------------------------------------- */

/* Extract a Gaussian-integer expression into two initialized mpz_t parts.
 * Accepts a rational integer (machine int or BigInt; imaginary part 0) or an
 * exact Gaussian integer Complex[a, b] with integer-like parts.  Returns true
 * on success (re/im are set); false for anything else (rationals, reals,
 * symbols, ...).  re and im must already be mpz_init'd by the caller. */
static bool coprimeq_to_gaussian(Expr* e, mpz_t re, mpz_t im) {
    if (expr_is_integer_like(e)) {
        expr_to_mpz(e, re);
        mpz_set_ui(im, 0);
        return true;
    }
    Expr *r, *i;
    if (is_complex(e, &r, &i) && expr_is_integer_like(r) && expr_is_integer_like(i)) {
        expr_to_mpz(r, re);
        expr_to_mpz(i, im);
        return true;
    }
    return false;
}

/* Round-to-nearest integer division for the Gaussian Euclidean step:
 * q = round(num / den) with den > 0, computed as floor((2 num + den)/(2 den)).
 * All operands are mpz_t; q must be init'd by the caller. */
static void coprimeq_round_div(mpz_t q, const mpz_t num, const mpz_t den) {
    mpz_t tnum, tden;
    mpz_init(tnum);
    mpz_init(tden);
    mpz_mul_2exp(tnum, num, 1);     /* 2*num */
    mpz_add(tnum, tnum, den);       /* 2*num + den */
    mpz_mul_2exp(tden, den, 1);     /* 2*den */
    mpz_fdiv_q(q, tnum, tden);      /* floor division (den > 0) */
    mpz_clears(tnum, tden, NULL);
}

/* True iff two Gaussian integers a = a1 + a2 i and b = b1 + b2 i are coprime
 * in Z[i], i.e. their Gaussian GCD is a unit (norm 1).  Uses the Euclidean
 * algorithm with round-to-nearest division; the remainder norm strictly
 * decreases, guaranteeing termination.  Inputs are consumed by value (copied
 * internally), so the caller's mpz_t are left untouched. */
static bool gaussian_pair_coprime(const mpz_t a1, const mpz_t a2,
                                  const mpz_t b1, const mpz_t b2) {
    mpz_t x1, x2, y1, y2, q1, q2, num, t, nb, r1, r2;
    mpz_inits(x1, x2, y1, y2, q1, q2, num, t, nb, r1, r2, NULL);
    mpz_set(x1, a1); mpz_set(x2, a2);
    mpz_set(y1, b1); mpz_set(y2, b2);

    /* Loop while y != 0. */
    while (mpz_sgn(y1) != 0 || mpz_sgn(y2) != 0) {
        /* norm(y) = y1^2 + y2^2. */
        mpz_mul(nb, y1, y1);
        mpz_mul(t, y2, y2);
        mpz_add(nb, nb, t);

        /* q = round(x * conj(y) / norm(y)).
         * x*conj(y) = (x1 y1 + x2 y2) + (x2 y1 - x1 y2) i. */
        mpz_mul(num, x1, y1);
        mpz_mul(t, x2, y2);
        mpz_add(num, num, t);
        coprimeq_round_div(q1, num, nb);

        mpz_mul(num, x2, y1);
        mpz_mul(t, x1, y2);
        mpz_sub(num, num, t);
        coprimeq_round_div(q2, num, nb);

        /* r = x - q*y; q*y = (q1 y1 - q2 y2) + (q1 y2 + q2 y1) i. */
        mpz_mul(r1, q1, y1);
        mpz_mul(t, q2, y2);
        mpz_sub(r1, r1, t);
        mpz_sub(r1, x1, r1);

        mpz_mul(r2, q1, y2);
        mpz_mul(t, q2, y1);
        mpz_add(r2, r2, t);
        mpz_sub(r2, x2, r2);

        /* x <- y; y <- r. */
        mpz_set(x1, y1); mpz_set(x2, y2);
        mpz_set(y1, r1); mpz_set(y2, r2);
    }

    /* x now holds the GCD; coprime iff norm(x) == 1 (x is a unit). */
    mpz_mul(t, x1, x1);
    mpz_mul(nb, x2, x2);
    mpz_add(t, t, nb);
    bool coprime = (mpz_cmp_ui(t, 1) == 0);

    mpz_clears(x1, x2, y1, y2, q1, q2, num, t, nb, r1, r2, NULL);
    return coprime;
}

/* CoprimeQ[n1, n2, ...] yields True when the arguments are pairwise relatively
 * prime, and False otherwise.  Over the ordinary integers a pair is coprime
 * when GCD == 1; with GaussianIntegers -> True (or when any argument is an
 * exact Gaussian integer) coprimality is tested in Z[i] via the Gaussian
 * Euclidean algorithm.  Machine ints and BigInts are handled uniformly through
 * GMP.  As a *Q predicate CoprimeQ always returns a Boolean: anything that is
 * not a manifestly coprime integer or Gaussian integer (rationals, reals,
 * symbols, malformed options) makes the result False.  CoprimeQ[] is False and
 * CoprimeQ[n] (a single argument, no pairs) is True.  Listable threading over
 * list arguments is performed by the evaluator before this builtin runs. */
Expr* builtin_coprimeq(Expr* res) {
    if (res->type != EXPR_FUNCTION) return expr_new_symbol(SYM_False);
    size_t argc = res->data.function.arg_count;
    Expr** args = res->data.function.args;

    /* Scan arguments: collect the numeric operands and pick up an optional
     * GaussianIntegers -> True | False rule (which, because CoprimeQ is
     * Orderless, may appear at any position).  Any other Rule, or a
     * GaussianIntegers option with a non-Boolean value, is a malformed option
     * and forces a False result. */
    bool gaussian = false;
    Expr** nums = NULL;
    size_t nnums = 0;
    if (argc > 0) nums = (Expr**)malloc(argc * sizeof(Expr*));
    for (size_t i = 0; i < argc; i++) {
        Expr* a = args[i];
        if (a->type == EXPR_FUNCTION &&
            a->data.function.head->type == EXPR_SYMBOL &&
            a->data.function.head->data.symbol.name == SYM_Rule &&
            a->data.function.arg_count == 2) {
            Expr* key = a->data.function.args[0];
            Expr* val = a->data.function.args[1];
            if (key->type == EXPR_SYMBOL && key->data.symbol.name == SYM_GaussianIntegers &&
                val->type == EXPR_SYMBOL &&
                (val->data.symbol.name == SYM_True || val->data.symbol.name == SYM_False)) {
                gaussian = (val->data.symbol.name == SYM_True);
                continue;
            }
            free(nums);
            return expr_new_symbol(SYM_False);  /* malformed / unknown option */
        }
        nums[nnums++] = a;
    }

    /* Per the documented contract: no operands -> False, one operand -> True. */
    if (nnums == 0) { free(nums); return expr_new_symbol(SYM_False); }
    if (nnums == 1) { free(nums); return expr_new_symbol(SYM_True); }

    /* If any operand is an exact Gaussian integer (non-real), switch to Z[i]. */
    if (!gaussian) {
        for (size_t i = 0; i < nnums; i++) {
            if (is_gaussian_integer(nums[i])) { gaussian = true; break; }
        }
    }

    bool result = true;

    if (gaussian) {
        /* Extract every operand as a Gaussian integer up front; any operand
         * that is not integer-like or a Gaussian integer is not manifestly
         * coprime -> False. */
        mpz_t* re = (mpz_t*)malloc(nnums * sizeof(mpz_t));
        mpz_t* im = (mpz_t*)malloc(nnums * sizeof(mpz_t));
        bool ok = true;
        size_t got = 0;
        for (size_t i = 0; i < nnums; i++) {
            mpz_init(re[i]);
            mpz_init(im[i]);
            got = i + 1;
            if (!coprimeq_to_gaussian(nums[i], re[i], im[i])) { ok = false; break; }
        }
        if (!ok) {
            result = false;
        } else {
            for (size_t i = 0; i < nnums && result; i++) {
                for (size_t j = i + 1; j < nnums && result; j++) {
                    if (!gaussian_pair_coprime(re[i], im[i], re[j], im[j])) result = false;
                }
            }
        }
        for (size_t i = 0; i < got; i++) mpz_clears(re[i], im[i], NULL);
        free(re);
        free(im);
    } else {
        /* Ordinary integers.  Every operand must be integer-like. */
        for (size_t i = 0; i < nnums; i++) {
            if (!expr_is_integer_like(nums[i])) { result = false; break; }
        }
        if (result) {
            mpz_t* v = (mpz_t*)malloc(nnums * sizeof(mpz_t));
            for (size_t i = 0; i < nnums; i++) {
                mpz_init(v[i]);
                expr_to_mpz(nums[i], v[i]);
            }
            mpz_t g;
            mpz_init(g);
            for (size_t i = 0; i < nnums && result; i++) {
                for (size_t j = i + 1; j < nnums && result; j++) {
                    mpz_gcd(g, v[i], v[j]);
                    if (mpz_cmp_ui(g, 1) != 0) result = false;
                }
            }
            mpz_clear(g);
            for (size_t i = 0; i < nnums; i++) mpz_clear(v[i]);
            free(v);
        }
    }

    free(nums);
    return expr_new_symbol(result ? SYM_True : SYM_False);
}
