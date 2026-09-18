/*
 * test_dsolve_m55_stress.c — anti-overfit stress families for M55
 * (DSolve`SpecialFunctionForm: the generalised power-potential recogniser).
 *
 *   Single power     y'' + A x^m y == 0
 *                        -> Sqrt[x] Z_{1/(m+2)}(kappa x^((m+2)/2)),  Z = BesselJ/Y.
 *   Two-term         y'' + (alpha x^(2c) + beta x^(c-1)) y == 0
 *                        -> Coulomb/Whittaker -> x^((1-d)/2) WhittakerM[k,+-mu,z]
 *                           (d = c+1), emitted as the verifiable Hypergeometric1F1 form.
 *   y'-carrying      y'' + A x^m y' + B x^(m-1) y == 0
 *                        -> the Liouville normal form is the two-term family; solved
 *                           via the normal-form pre-pass.
 *
 * Every family carries a SYMBOLIC exponent (the corpus shape) so the recogniser's
 * symbolic path is exercised (a numeric exponent would take the pre-existing
 * numeric pure-power Bessel row instead).  Solved through the PINNED method, then
 * verified by NUMERIC back-substitution on the ORIGINAL residual at instantiated
 * parameters (the solutions carry symbolic-exponent Bessel / 1F1 that zero_test
 * cannot discharge); Head === List is checked first so a decline cannot pass
 * vacuously.  The `n`/`nn`-named exponent families additionally guard the
 * ds_residual_numeric_zero / zero_test precision-ladder fix.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>

#include "core.h"
#include "eval.h"
#include "expr.h"
#include "parse.h"
#include "print.h"
#include "symtab.h"
#include "test_utils.h"

static char* eval_str(const char* input) {
    Expr* p = parse_expression(input);
    ASSERT(p != NULL);
    Expr* e = evaluate(p);
    expr_free(p);
    char* s = expr_to_string(e);
    expr_free(e);
    return s;
}
static bool lang_true(const char* input) {
    char* s = eval_str(input);
    bool ok = (strcmp(s, "True") == 0);
    if (!ok) fprintf(stderr, "  expected True: %s  =>  %s\n", input, s);
    free(s);
    return ok;
}
#define ASSERT_TRUE(input) ASSERT_MSG(lang_true(input), "expected True: %s", (input))

/* Solve `eqn` via the PINNED DSolve`SpecialFunctionForm, require Head === List,
 * then require the ORIGINAL residual `resid` to vanish numerically under `subs`. */
static void res_ok(const char* eqn, const char* resid, const char* subs) {
    char buf[2600];
    snprintf(buf, sizeof(buf),
        "With[{sol = DSolve`SpecialFunctionForm[%s, y, x]}, "
        "Head[sol] === List && Length[sol] >= 1 && "
        "Abs[N[(%s) /. sol[[1]] /. {%s}, 20]] < 10^-6]",
        eqn, resid, subs);
    ASSERT_TRUE(buf);
}

/* Family A: single power  y'' + a x^k y == 0  -> Bessel.  a>0 and a<0, k integer
 * and fractional. */
static void t_m55_stress_single_power(void) {
    const char* subs[] = {
        "a -> 1,  k -> 2,    C[1] -> 6/5, C[2] -> 7/10, x -> 3/5",
        "a -> -1, k -> 3,    C[1] -> 9/10, C[2] -> 11/10, x -> 4/5",
        "a -> 2,  k -> 1/2,  C[1] -> 6/5, C[2] -> 7/10, x -> 7/10",
        "a -> 1,  k -> -1/2, C[1] -> 6/5, C[2] -> 7/10, x -> 3/5",
        "a -> 3,  k -> 7/3,  C[1] -> 6/5, C[2] -> 7/10, x -> 1/2",
        "a -> -2, k -> 5/2,  C[1] -> 6/5, C[2] -> 7/10, x -> 4/5",
    };
    for (int i = 0; i < 6; i++)
        res_ok("y''[x] + a x^k y[x] == 0", "y''[x] + a x^k y[x]", subs[i]);
}

/* Family B: two-term  y'' + (a x^(2k) + b x^(k-1)) y == 0  -> Whittaker/1F1. */
static void t_m55_stress_two_term(void) {
    const char* subs[] = {
        "a -> -1, b -> 2,  k -> 3,   C[1] -> 6/5, C[2] -> 7/10, x -> 7/10",
        "a -> 1,  b -> 1,  k -> 2,   C[1] -> 6/5, C[2] -> 7/10, x -> 1/2",
        "a -> 2,  b -> -3, k -> 5/2, C[1] -> 6/5, C[2] -> 7/10, x -> 4/5",
        "a -> -2, b -> 1,  k -> 1/2, C[1] -> 6/5, C[2] -> 7/10, x -> 3/5",
        "a -> 1,  b -> 3,  k -> 4,   C[1] -> 6/5, C[2] -> 7/10, x -> 1/2",
    };
    for (int i = 0; i < 5; i++)
        res_ok("y''[x] + (a x^(2 k) + b x^(k - 1)) y[x] == 0",
               "y''[x] + (a x^(2 k) + b x^(k - 1)) y[x]", subs[i]);
}

/* Family C: y'-carrying  y'' + a x^k y' + b x^(k-1) y == 0  (normal-form pre-pass). */
static void t_m55_stress_yprime(void) {
    const char* subs[] = {
        "a -> 1, b -> 2,  k -> 3,   C[1] -> 6/5, C[2] -> 7/10, x -> 3/5",
        "a -> 2, b -> 1,  k -> 2,   C[1] -> 6/5, C[2] -> 7/10, x -> 1/2",
        "a -> 1, b -> -1, k -> 1/2, C[1] -> 6/5, C[2] -> 7/10, x -> 7/10",
        "a -> 3, b -> 2,  k -> 5/2, C[1] -> 6/5, C[2] -> 7/10, x -> 1/2",
    };
    for (int i = 0; i < 4; i++)
        res_ok("y''[x] + a x^k y'[x] + b x^(k - 1) y[x] == 0",
               "y''[x] + a x^k y'[x] + b x^(k - 1) y[x]", subs[i]);
}

/* Family D: the exponent symbol named `n` / `nn` must NOT hang (guards the
 * ds_residual_numeric_zero / zero_test precision-ladder fix). */
static void t_m55_stress_n_name_guard(void) {
    res_ok("y''[x] + x^n y[x] == 0", "y''[x] + x^n y[x]",
           "n -> 7/3, C[1] -> 6/5, C[2] -> 7/10, x -> 3/5");
    res_ok("y''[x] + (a x^(2 n) + b x^(n - 1)) y[x] == 0",
           "y''[x] + (a x^(2 n) + b x^(n - 1)) y[x]",
           "a -> 1, b -> 1, n -> 2, C[1] -> 6/5, C[2] -> 7/10, x -> 1/2");
    res_ok("y''[x] + (a x^(2 nn) + b x^(nn - 1)) y[x] == 0",
           "y''[x] + (a x^(2 nn) + b x^(nn - 1)) y[x]",
           "a -> -1, b -> 2, nn -> 3, C[1] -> 6/5, C[2] -> 7/10, x -> 7/10");
}

int main(void) {
    symtab_init();
    core_init();
    test_load_init_m();

    TEST(t_m55_stress_single_power);
    TEST(t_m55_stress_two_term);
    TEST(t_m55_stress_yprime);
    TEST(t_m55_stress_n_name_guard);

    printf("All DSolve M55 stress tests passed.\n");
    return 0;
}
