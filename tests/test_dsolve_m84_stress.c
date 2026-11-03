/*
 * test_dsolve_m84_stress.c — anti-overfit stress families for M84 (§2.2.58–59).
 *
 * M84 generalised DSolve`ChangeOfVariable, which had rationalised only trig
 * functions of the BARE variable (t = Cos/Sin/Tan[x]):
 *
 *   F1  exponential coefficients, t = E^(k x) with k the base frequency of every
 *       E^(c x) present (symbolic k included).  The transformed equation keeps a
 *       Y' term, so this also pins the dsolve_cv_active exemption of the
 *       SpecialFunctionForm normal-form pre-pass (top-level-only otherwise) and the
 *       denominator-cleared hand-off (a nested DSolve skips the linear normaliser,
 *       so the Euler form t^2 Y'' + a Y was unreachable from Y'' + a/t^2 Y).
 *   F2  scaled trig arguments (Cot[2 x], Cot[b x], Cot[x/2]) via s = k x plus
 *       TrigExpand, and the t = Sin[x]^2 substitution.
 *   F3  hyperbolic coefficients through TrigToExp, including the t = E^(2 k x)
 *       retry an even coefficient set (Tanh, Coth) needs.
 *   F4  the IVP branch fix in dsolve_fit_constants: of two fitted constants that
 *       both meet y(x0) == y0, the one whose body satisfies the ODE AT x0 wins
 *       (y' == 2 Sqrt[y], y[a] == b returned (x - a - Sqrt[b])^2 through v0.6147).
 *   F5  controls: the HEAD bare-x Legendre path, and a non-linear exponent that
 *       must decline fast.
 *
 * Every solve asserts Head === List first, so a declining method cannot pass
 * vacuously; every answer is checked by its numeric residual in the ORIGINAL
 * equation with every parameter and constant instantiated.
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

/* m84ok[eq]: DSolve answers (a List) and the residual of its first solution, with
 * C[1], C[2] and the parameters a, b, c, k, n instantiated, is below 1e-8 at three
 * points.  $-prefixed locals: plain names would collide with the parameters. */
static void define_helpers(void) {
    free(eval_str(
        "m84res[$eq_, $sol_] := Module[{$r = ($eq[[1]] - $eq[[2]]) /. $sol[[1]] /. "
        "{C[1] -> 0.8, C[2] -> 1.1} /. {a -> 0.7, b -> 1.3, c -> 0.9, k -> 0.6, n -> 0.35}}, "
        "Max[Table[Abs[N[$r /. x -> $p]], {$p, {0.4, 0.9, 1.3}}]]]"));
    free(eval_str(
        "m84ok[$eq_] := Module[{$s = TimeConstrained[DSolve[$eq, y, x], 20, $Aborted]}, "
        "Head[$s] === List && m84res[$eq, $s] < 10^-8]"));
}

static void check_all(const char* const* eqs, size_t n) {
    char q[600];
    for (size_t i = 0; i < n; i++) {
        snprintf(q, sizeof q, "m84ok[%s]", eqs[i]);
        ASSERT_TRUE(q);
    }
}

/* F1 — t = E^(k x): Bessel / Whittaker in t, numeric and symbolic k. */
static void t_m84_exponential(void) {
    static const char* const E[] = {
        "y''[x] + c E^(k x) y[x] == 0",                       /* symbolic k */
        "y''[x] + a y'[x] + (b + c E^(2 x)) y[x] == 0",       /* Y' term kept */
        "y''[x] + (b + c E^(x/2)) y[x] == 0",                 /* k = 1/2 */
        "y''[x] - (a^2 - c E^(3 x)) y[x] == 0",
        "y''[x] + (a + b E^x + c E^(2 x)) y[x] == 0",         /* Whittaker in E^x */
    };
    check_all(E, sizeof E / sizeof E[0]);
}

/* F2 — scaled trig arguments and t = Sin[x]^2. */
static void t_m84_scaled_trig(void) {
    static const char* const E[] = {
        "y''[x] - Cot[2 x] y'[x] + 2 y[x] == 0",
        "y''[x] + a Cot[b x] y'[x] + c y[x] == 0",            /* symbolic base k = b */
        "y''[x] + Tan[x] y'[x] + a Cot[x]^2 y[x] == 0",       /* Euler after t = Sin[x] */
        "y''[x] + 2 Cot[2 x] y'[x] + 4 n (n + 1) y[x] == 0",  /* Legendre in Cos[2 x] */
        "y''[x] + Cot[x/2] y'[x]/2 + n (n + 1)/4 y[x] == 0",  /* k = 1/2 */
    };
    check_all(E, sizeof E / sizeof E[0]);
}

/* F3 — hyperbolic coefficients (TrigToExp, then t = E^(k x) or E^(2 k x)). */
static void t_m84_hyperbolic(void) {
    static const char* const E[] = {
        "y''[x] + a Tanh[x] y'[x] + b y[x] == 0",
        "y''[x] + a Coth[x] y'[x] + b y[x] == 0",
        "y''[x] + 2 a Tanh[2 x] y'[x] + b y[x] == 0",
    };
    check_all(E, sizeof E / sizeof E[0]);
}

/* F4 — the IVP branch holds AT the initial point.  The residual is sampled just
 * either side of x0, where the wrong root (x - a - Sqrt[b])^2 has y' < 0. */
static void t_m84_ivp_local_branch(void) {
    ASSERT_TRUE("Module[{$s = DSolve[{y'[x] == 2 Sqrt[y[x]], y[a] == b}, y, x], $r}, "
                "$r = (y'[x] - 2 Sqrt[y[x]]) /. $s[[1]] /. {a -> 1/2, b -> 3/10}; "
                "Head[$s] === List && Max[Abs[N[$r /. x -> 1/2 + #]] & /@ {-1/100, 1/100}] < 10^-10]");
    ASSERT_TRUE("Module[{$s = DSolve[{y'[x] == 2 Sqrt[y[x]], y[a] == b}, y, x]}, "
                "Head[$s] === List && Simplify[(y[a] /. $s[[1]]) - b] === 0]");
    ASSERT_TRUE("DSolve[{y'[x] == 2 Sqrt[y[x]], y[1] == 4}, y, x] === {{y -> Function[{x}, (1 + x)^2]}}");
    ASSERT_TRUE("DSolve[{y'[x] == 3 y[x]^(2/3), y[0] == 8}, y, x] === {{y -> Function[{x}, (2 + x)^3]}}");
}

/* F5 — controls. */
static void t_m84_controls(void) {
    /* the HEAD bare-x trig path (Legendre in Cos[x]) is untouched */
    ASSERT_TRUE("m84ok[y''[x] + Cot[x] y'[x] + n (n + 1) y[x] == 0]");
    /* an exponent that is not linear in x has no base frequency: a fast decline */
    ASSERT_TRUE("Head[TimeConstrained[DSolve`ChangeOfVariable[y''[x] + E^(x^2) y[x] == 0, y, x], "
                "5, $Aborted]] === DSolve`ChangeOfVariable");
    /* mixed trig + exponential coefficients: neither family applies */
    ASSERT_TRUE("Head[TimeConstrained[DSolve`ChangeOfVariable[y''[x] + (Sin[x] + E^x) y[x] == 0, "
                "y, x], 10, $Aborted]] === DSolve`ChangeOfVariable");
}

/* F6 — two v0.6154 regressions caught by the campaign's final per-record A/B.
 * (a) PolynomialGCD's generalized-variable abstraction left x-free RADICAL
 * constants as coefficients, so the gcd went to the number-field engine over
 * Q((-3)^(1/3), (-1)^(2/3), 2^(1/3)): this IVP took 12.9 s (0.04 s before).
 * (b) Integrate split off an x-free constant carrying a sum, which then survived
 * unnormalised into SolvableForY's induced-ODE solution, so it was no longer
 * solved for the function and the record declined. */
static void t_m84_v6154_followups(void) {
    ASSERT_TRUE("Module[{$t = AbsoluteTime[], $s}, "
                "$s = DSolve[{(x^2 + 3 x) y'[x] == y[x]^3 + 2 y[x], y[1] == 1}, y, x]; "
                "Head[$s] === List && AbsoluteTime[] - $t < 3]");
    ASSERT_TRUE("Head[TimeConstrained[DSolve[(-y[x] + y'[x] x)^2 == y'[x]^2 - 2 y[x] y'[x]/x + 1, "
                "y, x], 30, $Aborted]] === List");
    /* the unit convention of the generalized gcd: a numeric factor is dropped */
    ASSERT_TRUE("PolynomialGCD[Sqrt[2] Sqrt[a] x, Sqrt[2] Sqrt[a] x^2] === Sqrt[a] x");
    ASSERT_TRUE("PolynomialGCD[(x + Sqrt[2]) Sqrt[a], (x + Sqrt[2]) x Sqrt[a]] === Sqrt[a] (Sqrt[2] + x)");
}

int main(void) {
    alarm(400);                          /* ~25 solves of 0.1–3 s each */
    setvbuf(stdout, NULL, _IOLBF, 0);    /* the shared alarm() loses buffered output */
    symtab_init();
    core_init();
    test_load_init_m();
    define_helpers();

    TEST(t_m84_exponential);
    TEST(t_m84_scaled_trig);
    TEST(t_m84_hyperbolic);
    TEST(t_m84_ivp_local_branch);
    TEST(t_m84_controls);
    TEST(t_m84_v6154_followups);

    printf("All DSolve M84 stress tests passed.\n");
    return 0;
}
