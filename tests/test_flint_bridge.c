/*
 * test_flint_bridge.c
 * -------------------
 * Extensive unit tests for the FLINT-backed Expr<->polynomial bridge (M1):
 * the rational multivariate GCD path src/poly/flint_bridge.c.
 *
 * Two avenues are exercised:
 *   (1) the C entry point flint_multivariate_gcd() directly — granular, and
 *       able to assert the NULL "bail to classical path" contract;
 *   (2) the scaffolding builtin Flint`GCD[...] through the evaluator — full
 *       integration (dispatch, Expr round-trip, re-evaluation/canonicalisation).
 *
 * GCD is defined up to a unit; FLINT returns the *monic* representative (lead
 * coefficient 1 under lex, smaller variable name = greater), so the expected
 * surface forms below are that canonical monic form. We compare the printed
 * form after evaluating the bridge's (deliberately unsimplified) output tree.
 *
 * When the build has no FLINT (USE_FLINT undefined) the bridge is all stubs;
 * main() detects this via flint_bridge_available() and skips cleanly (exit 0),
 * matching the graceful-degrade policy.
 */
#include "expr.h"
#include "eval.h"
#include "core.h"
#include "symtab.h"
#include "parse.h"
#include "print.h"
#include "flint_bridge.h"
#include "test_utils.h"

#include <stdio.h>
#include <string.h>
#include <stdlib.h>

/* ------------------------------------------------------------------ */
/*  Helpers                                                            */
/* ------------------------------------------------------------------ */

/* Parse + evaluate a source string to a canonical owned Expr. */
static Expr* eval_str(const char* src) {
    Expr* parsed = parse_expression(src);
    ASSERT(parsed != NULL);
    Expr* e = evaluate(parsed);
    expr_free(parsed);
    return e;
}

/*
 * Direct path: flint_multivariate_gcd(eval(a), eval(b)) must evaluate to
 * `expected` (printed form).
 */
static void check_gcd(const char* a, const char* b, const char* expected) {
    Expr* ea = eval_str(a);
    Expr* eb = eval_str(b);
    Expr* g  = flint_multivariate_gcd(ea, eb);
    if (g == NULL) {
        fprintf(stderr, "FAIL: flint_multivariate_gcd returned NULL for "
                        "gcd(%s, %s) (expected %s)\n", a, b, expected);
        exit(1);
    }
    Expr* ge = evaluate(g);
    char* s = expr_to_string(ge);
    if (strcmp(s, expected) != 0) {
        fprintf(stderr, "FAIL: gcd(%s, %s)\n  expected: %s\n  got:      %s\n",
                a, b, expected, s);
        exit(1);
    }
    free(s);
    expr_free(g);
    expr_free(ge);
    expr_free(ea);
    expr_free(eb);
}

/* GCD is commutative: assert gcd(a,b) == gcd(b,a) == expected. */
static void check_gcd_sym(const char* a, const char* b, const char* expected) {
    check_gcd(a, b, expected);
    check_gcd(b, a, expected);
}

/* Direct path: the bridge must decline (NULL) — input outside Q[x_1..x_n]. */
static void check_gcd_null(const char* a, const char* b) {
    Expr* ea = eval_str(a);
    Expr* eb = eval_str(b);
    Expr* g  = flint_multivariate_gcd(ea, eb);
    if (g != NULL) {
        char* s = expr_to_string(g);
        fprintf(stderr, "FAIL: gcd(%s, %s) should bail (NULL) but returned %s\n",
                a, b, s);
        free(s);
        exit(1);
    }
    expr_free(ea);
    expr_free(eb);
}

/* Direct number-field path: flint_numberfield_gcd(eval a, eval b) -> expected. */
static void check_gcd_nf(const char* a, const char* b, const char* expected) {
    Expr* ea = eval_str(a);
    Expr* eb = eval_str(b);
    Expr* g  = flint_numberfield_gcd(ea, eb);
    if (g == NULL) {
        fprintf(stderr, "FAIL: flint_numberfield_gcd returned NULL for "
                        "gcd(%s, %s) (expected %s)\n", a, b, expected);
        exit(1);
    }
    Expr* ge = evaluate(g);
    char* s = expr_to_string(ge);
    if (strcmp(s, expected) != 0) {
        fprintf(stderr, "FAIL: nf gcd(%s, %s)\n  expected: %s\n  got:      %s\n",
                a, b, expected, s);
        exit(1);
    }
    free(s);
    expr_free(g);
    expr_free(ge);
    expr_free(ea);
    expr_free(eb);
}

/* The number-field path must decline (NULL): no radical, tower, or multivariate. */
static void check_gcd_nf_null(const char* a, const char* b) {
    Expr* ea = eval_str(a);
    Expr* eb = eval_str(b);
    Expr* g  = flint_numberfield_gcd(ea, eb);
    if (g != NULL) {
        char* s = expr_to_string(g);
        fprintf(stderr, "FAIL: nf gcd(%s, %s) should bail (NULL) but returned %s\n",
                a, b, s);
        free(s);
        exit(1);
    }
    expr_free(ea);
    expr_free(eb);
}

/* Builtin path: Flint`GCD[...] through the evaluator. */
static void check_builtin(const char* input, const char* expected) {
    Expr* e = eval_str(input);
    char* s = expr_to_string(e);
    if (strcmp(s, expected) != 0) {
        fprintf(stderr, "FAIL: %s\n  expected: %s\n  got:      %s\n",
                input, expected, s);
        exit(1);
    }
    free(s);
    expr_free(e);
}

/* ------------------------------------------------------------------ */
/*  Univariate over Q                                                  */
/* ------------------------------------------------------------------ */

static void test_univariate(void) {
    check_gcd_sym("x^2 - 1", "x^2 + 2 x + 1", "1 + x");
    check_gcd_sym("x^2 - 1", "x - 1", "-1 + x");
    check_gcd_sym("x^3 - 1", "x^2 - 1", "-1 + x");
    check_gcd_sym("(x-1)*(x-2)*(x-3)", "(x-2)*(x-3)*(x-4)", "6 - 5 x + x^2");
    /* identical operands -> the (monic) polynomial itself */
    check_gcd("x^2 + 3 x + 2", "x^2 + 3 x + 2", "2 + 3 x + x^2");
    /* high degree: x^5 - 1 divides x^10 - 1, so the gcd is x^5 - 1 */
    check_gcd_sym("x^10 - 1", "x^5 - 1", "-1 + x^5");
    /* coprime cyclotomic pieces share only Phi_1 = x - 1 */
    check_gcd_sym("x^6 - 1", "x^4 - 1", "-1 + x^2");
}

/* ------------------------------------------------------------------ */
/*  Multivariate over Q                                                */
/* ------------------------------------------------------------------ */

static void test_multivariate(void) {
    check_gcd_sym("(x+y)*(x-2 y)", "(x+y)*(x+3 y)", "x + y");
    check_gcd_sym("x^3 y - x y^3", "x^2 y^2 - y^4", "x^2 y - y^3");
    /* common factor is a product of two binomials */
    check_gcd_sym("(x+y)*(x-y)*(x+1)", "(x+y)*(x-y)*(x+2)",
                  "x^2 - y^2");
    /* three variables, shared linear factor */
    check_gcd_sym("(x+y+z)*(x-z)", "(x+y+z)*(y+z)", "x + y + z");
    /* coprime multivariate -> 1 */
    check_gcd_sym("x + y", "x - y", "1");
}

/* ------------------------------------------------------------------ */
/*  Monic normalisation: canonical sign / leading coefficient         */
/* ------------------------------------------------------------------ */

static void test_monic_normalisation(void) {
    /* gcd(x^2 - y^2, x - y): lead term x has coeff 1 already -> x - y. */
    check_gcd_sym("x^2 - y^2", "x - y", "x - y");
    /* Content is stripped and result made monic: gcd(2x+2, 4x+4) = x + 1. */
    check_gcd_sym("2 x + 2", "4 x + 4", "1 + x");
    /* Negative leading coefficient is normalised away: gcd(-x-1, -2x-2). */
    check_gcd_sym("-x - 1", "-2 x - 2", "1 + x");
    /* Scaling both inputs by a rational leaves the monic gcd unchanged. */
    check_gcd_sym("(2/3) x^2 - (2/3)", "(5/7) x - 5/7", "-1 + x");
}

/* ------------------------------------------------------------------ */
/*  Rational coefficients                                              */
/* ------------------------------------------------------------------ */

static void test_rational_coeffs(void) {
    check_gcd_sym("x^2 - 1/4", "x - 1/2", "-1/2 + x");
    check_gcd_sym("x^2 - 1/9", "x + 1/3", "1/3 + x");
    check_gcd_sym("x^2 - 1/4", "x^2 + x + 1/4", "1/2 + x");
    /* large denominators */
    check_gcd_sym("x - 1/1000000000000", "x^2 - 1/1000000000000000000000000",
                  "-1/1000000000000 + x");
}

/* ------------------------------------------------------------------ */
/*  Big-integer coefficients (beyond int64)                           */
/* ------------------------------------------------------------------ */

static void test_bigint_coeffs(void) {
    /* 10^19 > 2^63: forces the EXPR_BIGINT conversion path both ways. */
    check_gcd_sym("(x + 10000000000000000000)*(x - 1)",
                  "(x + 10000000000000000000)*(x + 2)",
                  "10000000000000000000 + x");
    check_gcd_sym("x^2 - 100000000000000000000000000000000000000",
                  "x - 10000000000000000000",
                  "-10000000000000000000 + x");
}

/* ------------------------------------------------------------------ */
/*  Zero / identity edge cases                                        */
/* ------------------------------------------------------------------ */

static void test_zero_and_identity(void) {
    /* gcd(0, p) = p (monic) */
    check_gcd("0", "x^2 - 1", "-1 + x^2");
    check_gcd("x^2 - 1", "0", "-1 + x^2");
    /* gcd(p, 1) = 1 */
    check_gcd_sym("x^2 - 1", "1", "1");
    /* gcd where one divides the other */
    check_gcd_sym("x^2 - 1", "(x-1)*(x+1)*(x^2+1)", "-1 + x^2");
}

/* ------------------------------------------------------------------ */
/*  Determinism: result independent of variable discovery order       */
/* ------------------------------------------------------------------ */

static void test_determinism(void) {
    /* Same polynomials, written with variables introduced in opposite order. */
    check_gcd("y^2 + 2 x y + x^2", "y + x", "x + y");
    check_gcd("x^2 + 2 x y + y^2", "x + y", "x + y");
    /* Repeated calls are stable. */
    for (int i = 0; i < 5; i++)
        check_gcd("(a+b)*(a-b)", "(a+b)*(a+2 b)", "a + b");
}

/* ------------------------------------------------------------------ */
/*  Bail-out contract: NULL for inputs outside Q[x_1..x_n]            */
/* ------------------------------------------------------------------ */

static void test_bail_out_of_scope(void) {
    check_gcd_null("Sqrt[2] x", "x");          /* algebraic generator (M3) */
    check_gcd_null("x + Sqrt[3]", "x");        /* irrational constant term */
    check_gcd_null("x - 1.5", "x - 1");        /* inexact real coefficient */
    check_gcd_null("x^2 - 2.0", "x");          /* inexact real coefficient */
    check_gcd_null("1/x", "x");                /* negative power (Laurent)  */
    check_gcd_null("x^n", "x");                /* symbolic exponent         */
    check_gcd_null("6", "4");                  /* pure numeric -> classical */
    check_gcd_null("2/3", "4/9");              /* pure rational             */
    check_gcd_null("Sin[x]", "x");             /* transcendental head       */
}

/* ------------------------------------------------------------------ */
/*  Number field Q(sqrt d) (M2)                                        */
/* ------------------------------------------------------------------ */

static void test_numberfield(void) {
    /* The rigorous split the rational/free-variable view cannot make. */
    check_gcd_nf("x^2 - 2", "x - Sqrt[2]", "-Sqrt[2] + x");
    check_gcd_nf("x - Sqrt[2]", "x^2 - 2", "-Sqrt[2] + x");   /* commutative */
    check_gcd_nf("x^2 - 2", "x + Sqrt[2]", "Sqrt[2] + x");
    check_gcd_nf("x^2 - 3", "x - Sqrt[3]", "-Sqrt[3] + x");
    check_gcd_nf("x^2 - 5", "x + Sqrt[5]", "Sqrt[5] + x");
    /* non-square, non-prime radicand */
    check_gcd_nf("x^2 - 6", "x - Sqrt[6]", "-Sqrt[6] + x");
    /* shared quadratic factor survives a coprime cofactor */
    check_gcd_nf("x^2 - 2", "(x - Sqrt[2])*(x - 5)", "-Sqrt[2] + x");
    /* mixed rational + alpha coefficient reconstruction (c0 + c1 Sqrt[d]) */
    check_gcd_nf("(x - 1 - Sqrt[2])*(x - 1)", "(x - 1 - Sqrt[2])*(x - 2)",
                 "-1 - Sqrt[2] + x");
    check_gcd_nf("(x - 3 - 2 Sqrt[2])*(x + 1)", "(x - 3 - 2 Sqrt[2])*(x - 1)",
                 "-3 - 2 Sqrt[2] + x");
    /* leading coefficient in the field is normalised away (monic) */
    check_gcd_nf("x^2 - 2", "Sqrt[2]*(x - Sqrt[2])", "-Sqrt[2] + x");
    check_gcd_nf("Sqrt[2] x", "x", "x");
    /* conjugate factors are coprime over the field */
    check_gcd_nf("x - Sqrt[2]", "x + Sqrt[2]", "1");

    /* Bail-outs -> NULL (handed to other paths / later milestones). */
    check_gcd_nf_null("x^2 - 1", "x - 1");        /* no radical: rational path */
    check_gcd_nf_null("x - Sqrt[2]", "x - Sqrt[3]"); /* tower (two radicals)   */
    check_gcd_nf_null("x - Sqrt[2]", "y - Sqrt[2]"); /* two variables          */
    check_gcd_nf_null("x - Sqrt[x]", "x");        /* radicand not an integer   */
}

/* ------------------------------------------------------------------ */
/*  Cyclotomic field Q(zeta_n) (M2)                                    */
/* ------------------------------------------------------------------ */

static void check_gcd_cyc(const char* a, const char* b, const char* expected) {
    Expr* ea = eval_str(a);
    Expr* eb = eval_str(b);
    Expr* g  = flint_cyclotomic_gcd(ea, eb);
    if (g == NULL) {
        fprintf(stderr, "FAIL: flint_cyclotomic_gcd returned NULL for "
                        "gcd(%s, %s) (expected %s)\n", a, b, expected);
        exit(1);
    }
    Expr* ge = evaluate(g);
    char* s = expr_to_string(ge);
    if (strcmp(s, expected) != 0) {
        fprintf(stderr, "FAIL: cyc gcd(%s, %s)\n  expected: %s\n  got:      %s\n",
                a, b, expected, s);
        exit(1);
    }
    free(s);
    expr_free(g);
    expr_free(ge);
    expr_free(ea);
    expr_free(eb);
}

static void check_gcd_cyc_null(const char* a, const char* b) {
    Expr* ea = eval_str(a);
    Expr* eb = eval_str(b);
    Expr* g  = flint_cyclotomic_gcd(ea, eb);
    if (g != NULL) {
        char* s = expr_to_string(g);
        fprintf(stderr, "FAIL: cyc gcd(%s, %s) should bail (NULL) but returned %s\n",
                a, b, s);
        free(s);
        exit(1);
    }
    expr_free(ea);
    expr_free(eb);
}

static void test_cyclotomic(void) {
    /* The result is returned in the canonical generator basis (-1)^(1/Q);
     * exact forms below were verified equal to the input root-of-unity via
     * PossibleZeroQ (see the semantic cross-checks at the end). */

    /* Q(zeta_3): x^2 + x + 1 = (x - zeta_3)(x - zeta_3^2). zeta_3 = (-1)^(2/3);
     * the field generator is (-1)^(1/3) = zeta_6, and zeta_3 = zeta_6 - 1. */
    check_gcd_cyc("x^2 + x + 1", "x - (-1)^(2/3)", "1 - (-1)^(1/3) + x");
    /* Phi_6 = x^2 - x + 1, root (-1)^(1/3) directly */
    check_gcd_cyc("x^2 - x + 1", "x - (-1)^(1/3)", "-(-1)^(1/3) + x");
    check_gcd_cyc("(x - (-1)^(1/3))*(x - 1)", "(x - (-1)^(1/3))*(x - 2)",
                  "-(-1)^(1/3) + x");
    /* Q(zeta_5): Phi_5 split by a primitive 5th root (-1)^(2/5) */
    check_gcd_cyc("x^4 + x^3 + x^2 + x + 1", "x - (-1)^(2/5)", "-(-1)^(2/5) + x");
    /* shared root-of-unity linear factor survives a coprime cofactor */
    check_gcd_cyc("(x - (-1)^(2/5))*(x - 2)", "(x - (-1)^(2/5))*(x - 3)",
                  "-(-1)^(2/5) + x");
    /* coprime cases -> 1 (each carries a root of unity, so the path engages) */
    check_gcd_cyc("x - (-1)^(1/3)", "x - 1", "1");
    check_gcd_cyc("x - (-1)^(1/3)", "x - (-1)^(2/3)", "1");

    /* Semantic cross-checks, independent of the output basis. */
    check_builtin("PossibleZeroQ[Flint`GCD[x^2 + x + 1, x - (-1)^(2/3)] "
                  "- (x - (-1)^(2/3))]", "True");
    check_builtin("PossibleZeroQ[Flint`GCD[x^4 + x^3 + x^2 + x + 1, x - (-1)^(2/5)] "
                  "- (x - (-1)^(2/5))]", "True");

    /* Bail-outs -> NULL. */
    check_gcd_cyc_null("x^2 - 1", "x - 1");                 /* no root of unity */
    check_gcd_cyc_null("x - Sqrt[2]", "x - (-1)^(1/3)");    /* radical: a tower */
    check_gcd_cyc_null("x - (-1)^(1/3)", "y - (-1)^(1/3)"); /* two variables    */
    check_gcd_cyc_null("x - I", "x - 1");                   /* Complex[0,1] deferred */
}

/* ------------------------------------------------------------------ */
/*  Radical tower Q(sqrt d_1, ..., sqrt d_r) (M2)                       */
/* ------------------------------------------------------------------ */

static void check_gcd_tower(const char* a, const char* b, const char* expected) {
    Expr* ea = eval_str(a);
    Expr* eb = eval_str(b);
    Expr* g  = flint_tower_gcd(ea, eb);
    if (g == NULL) {
        fprintf(stderr, "FAIL: flint_tower_gcd returned NULL for "
                        "gcd(%s, %s) (expected %s)\n", a, b, expected);
        exit(1);
    }
    Expr* ge = evaluate(g);
    char* s = expr_to_string(ge);
    if (strcmp(s, expected) != 0) {
        fprintf(stderr, "FAIL: tower gcd(%s, %s)\n  expected: %s\n  got:      %s\n",
                a, b, expected, s);
        exit(1);
    }
    free(s);
    expr_free(g);
    expr_free(ge);
    expr_free(ea);
    expr_free(eb);
}

static void check_gcd_tower_null(const char* a, const char* b) {
    Expr* ea = eval_str(a);
    Expr* eb = eval_str(b);
    Expr* g  = flint_tower_gcd(ea, eb);
    if (g != NULL) {
        char* s = expr_to_string(g);
        fprintf(stderr, "FAIL: tower gcd(%s, %s) should bail (NULL) but returned %s\n",
                a, b, s);
        free(s);
        exit(1);
    }
    expr_free(ea);
    expr_free(eb);
}

static void test_tower(void) {
    /* Q(sqrt 2, sqrt 3): shared linear factor with a single-radical coefficient */
    check_gcd_tower("(x - Sqrt[2])*(x - Sqrt[3])", "(x - Sqrt[2])*(x + Sqrt[3])",
                    "-Sqrt[2] + x");
    /* mixed coefficient sqrt2 + sqrt3 reconstructed in the product basis */
    check_gcd_tower("(x - Sqrt[2] - Sqrt[3])*(x - 1)",
                    "(x - Sqrt[2] - Sqrt[3])*(x - 2)", "-Sqrt[2] - Sqrt[3] + x");
    /* a shared quadratic over the tower: cross term sqrt2*sqrt3 = sqrt6 */
    check_gcd_tower("(x - Sqrt[2])*(x - Sqrt[3])*(x - 1)",
                    "(x - Sqrt[2])*(x - Sqrt[3])*(x - 2)",
                    "Sqrt[6] + (-Sqrt[2] - Sqrt[3]) x + x^2");
    /* one input has the radical only inside a quadratic */
    check_gcd_tower("x^2 - 5", "(x - Sqrt[5])*(x - Sqrt[2])", "-Sqrt[5] + x");
    /* coprime over the tower -> 1 */
    check_gcd_tower("x - Sqrt[2]", "x - Sqrt[3]", "1");
    /* three independent radicals: Q(sqrt2, sqrt3, sqrt5), degree 8 */
    check_gcd_tower("(x - Sqrt[2] - Sqrt[3] - Sqrt[5])*(x - 1)",
                    "(x - Sqrt[2] - Sqrt[3] - Sqrt[5])*(x + 1)",
                    "-Sqrt[2] - Sqrt[3] - Sqrt[5] + x");

    /* Semantic cross-checks. */
    check_builtin("PossibleZeroQ[Flint`GCD[(x - Sqrt[2] - Sqrt[3])*(x - 1), "
                  "(x - Sqrt[2] - Sqrt[3])*(x - 2)] - (x - Sqrt[2] - Sqrt[3])]", "True");

    /* Bail-outs -> NULL. */
    check_gcd_tower_null("x^2 - 2", "x - Sqrt[2]");        /* one radical: nf path */
    check_gcd_tower_null("x - (-1)^(1/3)", "x - 1");       /* root of unity        */
    check_gcd_tower_null("x - Sqrt[2] - Sqrt[3]", "y - Sqrt[2]"); /* two variables */
}

/* ------------------------------------------------------------------ */
/*  Parametric radical Q(t_1..t_p)(sqrt k), k a symbol (M3)            */
/* ------------------------------------------------------------------ */

static void check_gcd_param(const char* a, const char* b, const char* expected) {
    Expr* ea = eval_str(a);
    Expr* eb = eval_str(b);
    Expr* g  = flint_parametric_sqrt_gcd(ea, eb);
    if (g == NULL) {
        fprintf(stderr, "FAIL: flint_parametric_sqrt_gcd NULL for gcd(%s, %s) "
                        "(expected %s)\n", a, b, expected);
        exit(1);
    }
    char* s = expr_to_string(g);
    if (strcmp(s, expected) != 0) {
        fprintf(stderr, "FAIL: param gcd(%s, %s)\n  expected: %s\n  got:      %s\n",
                a, b, expected, s);
        exit(1);
    }
    free(s);
    expr_free(g);
    expr_free(ea);
    expr_free(eb);
}

static void check_gcd_param_null(const char* a, const char* b) {
    Expr* ea = eval_str(a);
    Expr* eb = eval_str(b);
    Expr* g  = flint_parametric_sqrt_gcd(ea, eb);
    if (g != NULL) {
        char* s = expr_to_string(g);
        fprintf(stderr, "FAIL: param gcd(%s, %s) should bail (NULL) but got %s\n",
                a, b, s);
        free(s);
        exit(1);
    }
    expr_free(ea);
    expr_free(eb);
}

/* Assert `head[op]` terminates and preserves value: the reduced form equals the
 * input at a fixed rational sample point (deterministic — PossibleZeroQ's
 * random sampling misfires on these huge k^(9/2)-coefficient expressions even
 * though the exact difference is 0). Substitutes every symbol that can occur
 * (a,b,k,u,x,y) with a generic non-pole rational and checks |N[diff]| ~ 0. */
static void check_reduce_preserving_head(const char* head, const char* op) {
    char buf[9000];
    snprintf(buf, sizeof buf,
             "Abs[N[(%s[%s] - (%s)) /. "
             "{a -> 2, b -> 3, k -> 5, u -> 7/10, v -> 9/10, w -> 13/10, x -> 7/10, y -> 11/10, z -> 3/10}, 30]] < 10^-20",
             head, op, op);
    Expr* e = eval_str(buf);
    char* s = expr_to_string(e);
    if (strcmp(s, "True") != 0) {
        fprintf(stderr, "FAIL: %s value-preservation\n  op: %.100s...\n  got: %s\n",
                head, op, s);
        exit(1);
    }
    free(s); expr_free(e);
}

/* Used for the large parametric Q(a,b,k)(Sqrt k) operands captured from the
 * Goursat Sqrt[k] descent — before the FLINT parametric radical path these
 * canonic = Cancel[Together[...]] calls hung. */
static void check_reduces_preserving(const char* op) {
    check_reduce_preserving_head("Cancel", op);
    check_reduce_preserving_head("Together", op);
}

/* Auto-captured Cancel/Together operands from the Goursat Sqrt[k] descent
 * (Integrate[(k x^2-1)/((a k x+b)(b x+a) Sqrt[x(1-x)(1-k x)]), x]).
 * Each is a rational function over the parametric radical field Q(a,b,k)(Sqrt k)
 * that the descent's canonic = Cancel[Together[...]] must reduce fast (these
 * hung before the FLINT parametric path was wired in). We assert termination +
 * value-preservation (the reduced form equals the input). */
static void test_goursat_descent_operands(void) {
    check_reduces_preserving("(-1/Sqrt[k] - u/Sqrt[k])/(1 - u)");
    check_reduces_preserving("(-1 + 2 u^2 - u^4 - 2 Sqrt[k] + 2 u^4 Sqrt[k] - k + 2 u^2 k - u^4 k)/k");
    check_reduces_preserving("(-1 + 2 u - u^2 - 2 Sqrt[k] + 2 u^2 Sqrt[k] - k + 2 u k - u^2 k)/(-1 + 2 Sqrt[k] - k)");
    check_reduces_preserving("0 + ((-1 - 1/k - 2/Sqrt[k]) u^0)/(-1 - 1/k + 2/Sqrt[k]) + ((2 + 2/k) u^1)/(-1 - 1/k + 2/Sqrt[k]) + ((-1 - 1/k + 2/Sqrt[k]) u^2)/(-1 - 1/k + 2/Sqrt[k])");
    check_reduces_preserving("((-1/Sqrt[k] - u/Sqrt[k])/(1 - u) - (-1/Sqrt[k] - u/Sqrt[k])^2/(1 - u)^2 + (k (-1/Sqrt[k] - u/Sqrt[k])^3)/(1 - u)^3 - (k (-1/Sqrt[k] - u/Sqrt[k])^2)/(1 - u)^2) (1 - u)^4");
    check_reduces_preserving("(8 a b k + 4 b^2 k + 8 a b k^(3/2) + 4 b^2 k^(3/2) + 4 a^2 k^2 + 4 a^2 k^(5/2) + u (8 a b k + 4 b^2 k - 8 a b k^(3/2) - 4 b^2 k^(3/2) + 4 a^2 k^2 - 4 a^2 k^(5/2)))/(-b^4 - 2 b^4 Sqrt[k] + 2 a^2 b^2 k - b^4 k + 4 a^2 b^2 k^(3/2) - a^4 k^2 + 2 a^2 b^2 k^2 - 2 a^4 k^(5/2) - a^4 k^3 + u (2 b^4 + 12 a^2 b^2 k + 16 a b^3 k + 2 b^4 k + 2 a^4 k^2 + 16 a^3 b k^2 + 12 a^2 b^2 k^2 + 2 a^4 k^3) + u^2 (-b^4 + 2 b^4 Sqrt[k] + 2 a^2 b^2 k - b^4 k - 4 a^2 b^2 k^(3/2) - a^4 k^2 + 2 a^2 b^2 k^2 + 2 a^4 k^(5/2) - a^4 k^3))");
    check_reduces_preserving("(8 a b k + 4 b^2 k + 8 a b k^(3/2) + 4 b^2 k^(3/2) + 4 a^2 k^2 + 4 a^2 k^(5/2) + u (8 a b k + 4 b^2 k - 8 a b k^(3/2) - 4 b^2 k^(3/2) + 4 a^2 k^2 - 4 a^2 k^(5/2)))/(v (-b^4 - 2 b^4 Sqrt[k] + 2 a^2 b^2 k - b^4 k + 4 a^2 b^2 k^(3/2) - a^4 k^2 + 2 a^2 b^2 k^2 - 2 a^4 k^(5/2) - a^4 k^3 + u (2 b^4 + 12 a^2 b^2 k + 16 a b^3 k + 2 b^4 k + 2 a^4 k^2 + 16 a^3 b k^2 + 12 a^2 b^2 k^2 + 2 a^4 k^3) + u^2 (-b^4 + 2 b^4 Sqrt[k] + 2 a^2 b^2 k - b^4 k - 4 a^2 b^2 k^(3/2) - a^4 k^2 + 2 a^2 b^2 k^2 + 2 a^4 k^(5/2) - a^4 k^3)))");
    check_reduces_preserving("(0 + (8 a b k + 4 b^2 k + 8 a b k^(3/2) + 4 b^2 k^(3/2) + 4 a^2 k^2 + 4 a^2 k^(5/2)) u^0 + (8 a b k + 4 b^2 k - 8 a b k^(3/2) - 4 b^2 k^(3/2) + 4 a^2 k^2 - 4 a^2 k^(5/2)) u^1)/(0 + (-b^4 - 2 b^4 Sqrt[k] + 2 a^2 b^2 k - b^4 k + 4 a^2 b^2 k^(3/2) - a^4 k^2 + 2 a^2 b^2 k^2 - 2 a^4 k^(5/2) - a^4 k^3) u^0 + (2 b^4 + 12 a^2 b^2 k + 16 a b^3 k + 2 b^4 k + 2 a^4 k^2 + 16 a^3 b k^2 + 12 a^2 b^2 k^2 + 2 a^4 k^3) u^1 + (-b^4 + 2 b^4 Sqrt[k] + 2 a^2 b^2 k - b^4 k - 4 a^2 b^2 k^(3/2) - a^4 k^2 + 2 a^2 b^2 k^2 + 2 a^4 k^(5/2) - a^4 k^3) u^2)");
    check_reduces_preserving("(8 a b k + 8 u^2 a b k + 4 b^2 k + 4 u^2 b^2 k + 8 a b k^(3/2) - 8 u^2 a b k^(3/2) + 4 b^2 k^(3/2) - 4 u^2 b^2 k^(3/2) + 4 a^2 k^2 + 4 u^2 a^2 k^2 + 4 a^2 k^(5/2) - 4 u^2 a^2 k^(5/2))/(-b^4 + 2 u^2 b^4 - u^4 b^4 - 2 b^4 Sqrt[k] + 2 u^4 b^4 Sqrt[k] + 2 a^2 b^2 k + 12 u^2 a^2 b^2 k + 2 u^4 a^2 b^2 k + 16 u^2 a b^3 k - b^4 k + 2 u^2 b^4 k - u^4 b^4 k + 4 a^2 b^2 k^(3/2) - 4 u^4 a^2 b^2 k^(3/2) - a^4 k^2 + 2 u^2 a^4 k^2 - u^4 a^4 k^2 + 16 u^2 a^3 b k^2 + 2 a^2 b^2 k^2 + 12 u^2 a^2 b^2 k^2 + 2 u^4 a^2 b^2 k^2 - 2 a^4 k^(5/2) + 2 u^4 a^4 k^(5/2) - a^4 k^3 + 2 u^2 a^4 k^3 - u^4 a^4 k^3)");
    check_reduces_preserving("((8 u a b k + 8 u^3 a b k + 4 u b^2 k + 4 u^3 b^2 k + 8 u a b k^(3/2) - 8 u^3 a b k^(3/2) + 4 u b^2 k^(3/2) - 4 u^3 b^2 k^(3/2) + 4 u a^2 k^2 + 4 u^3 a^2 k^2 + 4 u a^2 k^(5/2) - 4 u^3 a^2 k^(5/2))/(-b^4 + 2 u^2 b^4 - u^4 b^4 - 2 b^4 Sqrt[k] + 2 u^4 b^4 Sqrt[k] + 2 a^2 b^2 k + 12 u^2 a^2 b^2 k + 2 u^4 a^2 b^2 k + 16 u^2 a b^3 k - b^4 k + 2 u^2 b^4 k - u^4 b^4 k + 4 a^2 b^2 k^(3/2) - 4 u^4 a^2 b^2 k^(3/2) - a^4 k^2 + 2 u^2 a^4 k^2 - u^4 a^4 k^2 + 16 u^2 a^3 b k^2 + 2 a^2 b^2 k^2 + 12 u^2 a^2 b^2 k^2 + 2 u^4 a^2 b^2 k^2 - 2 a^4 k^(5/2) + 2 u^4 a^4 k^(5/2) - a^4 k^3 + 2 u^2 a^4 k^3 - u^4 a^4 k^3))/u");
    check_reduces_preserving("(8 a b k + 4 b^2 k + 8 a b k^(3/2) + 4 b^2 k^(3/2) + 4 a^2 k^2 + 4 a^2 k^(5/2) + u (8 a b k + 4 b^2 k - 8 a b k^(3/2) - 4 b^2 k^(3/2) + 4 a^2 k^2 - 4 a^2 k^(5/2)))/(Sqrt[(-1 + 2 u - u^2 - 2 Sqrt[k] + 2 u^2 Sqrt[k] - k + 2 u k - u^2 k)/(-1 + 2 Sqrt[k] - k)] (-b^4 - 2 b^4 Sqrt[k] + 2 a^2 b^2 k - b^4 k + 4 a^2 b^2 k^(3/2) - a^4 k^2 + 2 a^2 b^2 k^2 - 2 a^4 k^(5/2) - a^4 k^3 + u (2 b^4 + 12 a^2 b^2 k + 16 a b^3 k + 2 b^4 k + 2 a^4 k^2 + 16 a^3 b k^2 + 12 a^2 b^2 k^2 + 2 a^4 k^3) + u^2 (-b^4 + 2 b^4 Sqrt[k] + 2 a^2 b^2 k - b^4 k - 4 a^2 b^2 k^(3/2) - a^4 k^2 + 2 a^2 b^2 k^2 + 2 a^4 k^(5/2) - a^4 k^3)))");
    check_reduces_preserving("(-2 a b - b^2 - a^2 k + (a^2 k^3 (-1/Sqrt[k] - u/Sqrt[k])^4)/(1 - u)^4 - 2 (a^2 k^3 (-1/Sqrt[k] - u/Sqrt[k])^3)/(1 - u)^3 + 2 (a^2 k^2 (-1/Sqrt[k] - u/Sqrt[k]))/(1 - u) - 4 (a b k^2 (-1/Sqrt[k] - u/Sqrt[k])^3)/(1 - u)^3 + 2 (a b k^2 (-1/Sqrt[k] - u/Sqrt[k])^4)/(1 - u)^4 + 4 (a b k (-1/Sqrt[k] - u/Sqrt[k]))/(1 - u) + (b^2 k^2 (-1/Sqrt[k] - u/Sqrt[k])^4)/(1 - u)^4 - 2 (b^2 k^2 (-1/Sqrt[k] - u/Sqrt[k])^3)/(1 - u)^3 + 2 (b^2 k (-1/Sqrt[k] - u/Sqrt[k]))/(1 - u))/(2 a^2 b^2 + 2 a b^3 + 2 a^3 b k + 2 a^2 b^2 k - 2 (a^4 k^2 (-1/Sqrt[k] - u/Sqrt[k])^2)/(1 - u)^2 - 2 (a^4 k^3 (-1/Sqrt[k] - u/Sqrt[k])^2)/(1 - u)^2 + 2 (a^4 k^3 (-1/Sqrt[k] - u/Sqrt[k])^3)/(1 - u)^3 + 2 (a^4 k^2 (-1/Sqrt[k] - u/Sqrt[k]))/(1 - u) - 4 (a^3 b k^2 (-1/Sqrt[k] - u/Sqrt[k])^2)/(1 - u)^2 + 2 (a^3 b k^3 (-1/Sqrt[k] - u/Sqrt[k])^4)/(1 - u)^4 - 4 (a^2 b^2 k^2 (-1/Sqrt[k] - u/Sqrt[k])^3)/(1 - u)^3 - 4 (a^2 b^2 k (-1/Sqrt[k] - u/Sqrt[k]))/(1 - u) + 2 (a^2 b^2 k^2 (-1/Sqrt[k] - u/Sqrt[k])^4)/(1 - u)^4 + 2 (a^2 b^2 k^3 (-1/Sqrt[k] - u/Sqrt[k])^4)/(1 - u)^4 - 4 (a b^3 k (-1/Sqrt[k] - u/Sqrt[k])^2)/(1 - u)^2 + 2 (a b^3 k^2 (-1/Sqrt[k] - u/Sqrt[k])^4)/(1 - u)^4 - 2 (b^4 (-1/Sqrt[k] - u/Sqrt[k])^2)/(1 - u)^2 + 2 (b^4 (-1/Sqrt[k] - u/Sqrt[k]))/(1 - u) - 2 (b^4 k (-1/Sqrt[k] - u/Sqrt[k])^2)/(1 - u)^2 + 2 (b^4 k (-1/Sqrt[k] - u/Sqrt[k])^3)/(1 - u)^3)");
    check_reduces_preserving("(a^2 (-32 u^4 k^2 + u (4 k^2 + 4 k^(5/2)) + u^2 (-16 k^2 - 16 k^(5/2)) + u^3 (28 k^2 + 20 k^(5/2)) + u^5 (28 k^2 - 20 k^(5/2)) + u^6 (-16 k^2 + 16 k^(5/2)) + u^7 (4 k^2 - 4 k^(5/2))) + a b (-64 u^4 k + u (8 k + 8 k^(3/2)) + u^2 (-32 k - 32 k^(3/2)) + u^3 (56 k + 40 k^(3/2)) + u^5 (56 k - 40 k^(3/2)) + u^6 (-32 k + 32 k^(3/2)) + u^7 (8 k - 8 k^(3/2))) + b^2 (-32 u^4 k + u (4 k + 4 k^(3/2)) + u^2 (-16 k - 16 k^(3/2)) + u^3 (28 k + 20 k^(3/2)) + u^5 (28 k - 20 k^(3/2)) + u^6 (-16 k + 16 k^(3/2)) + u^7 (4 k - 4 k^(3/2))))/(a^4 (-k^2 - 2 k^(5/2) - k^3 + u^4 (10 k^2 + 10 k^3) + u (4 k^2 + 8 k^(5/2) + 4 k^3) + u^2 (-4 k^2 - 12 k^(5/2) - 4 k^3) + u^3 (-4 k^2 + 8 k^(5/2) - 4 k^3) + u^5 (-4 k^2 - 8 k^(5/2) - 4 k^3) + u^6 (-4 k^2 + 12 k^(5/2) - 4 k^3) + u^7 (4 k^2 - 8 k^(5/2) + 4 k^3) + u^8 (-k^2 + 2 k^(5/2) - k^3)) + a^3 b (16 u^2 k^2 - 64 u^3 k^2 + 96 u^4 k^2 - 64 u^5 k^2 + 16 u^6 k^2) + a^2 b^2 (2 k + 4 k^(3/2) + 2 k^2 + u^4 (76 k + 76 k^2) + u (-8 k - 16 k^(3/2) - 8 k^2) + u^2 (24 k + 24 k^(3/2) + 24 k^2) + u^3 (-56 k - 16 k^(3/2) - 56 k^2) + u^5 (-56 k + 16 k^(3/2) - 56 k^2) + u^6 (24 k - 24 k^(3/2) + 24 k^2) + u^7 (-8 k + 16 k^(3/2) - 8 k^2) + u^8 (2 k - 4 k^(3/2) + 2 k^2)) + a b^3 (16 u^2 k - 64 u^3 k + 96 u^4 k - 64 u^5 k + 16 u^6 k) + b^4 (-1 - 2 Sqrt[k] - k + u^4 (10 + 10 k) + u (4 + 8 Sqrt[k] + 4 k) + u^2 (-4 - 12 Sqrt[k] - 4 k) + u^3 (-4 + 8 Sqrt[k] - 4 k) + u^5 (-4 - 8 Sqrt[k] - 4 k) + u^6 (-4 + 12 Sqrt[k] - 4 k) + u^7 (4 - 8 Sqrt[k] + 4 k) + u^8 (-1 + 2 Sqrt[k] - k)))");
    check_reduces_preserving("(2 (-(-1/(-1 + 2 Sqrt[k] - k) - 2 Sqrt[k]/(-1 + 2 Sqrt[k] - k) - k/(-1 + 2 Sqrt[k] - k)) Sqrt[-1/(-1 + 2 Sqrt[k] - k) - k/(-1 + 2 Sqrt[k] - k) + 2 Sqrt[k]/(-1 + 2 Sqrt[k] - k)] + u (2/(-1 + 2 Sqrt[k] - k) + 2 k/(-1 + 2 Sqrt[k] - k)) - u^2 Sqrt[-1/(-1 + 2 Sqrt[k] - k) - k/(-1 + 2 Sqrt[k] - k) + 2 Sqrt[k]/(-1 + 2 Sqrt[k] - k)]) (8 a b k + 4 b^2 k + 8 a b k^(3/2) + 4 b^2 k^(3/2) + 4 a^2 k^2 + 4 a^2 k^(5/2) + ((u^2 + 1/(-1 + 2 Sqrt[k] - k) + k/(-1 + 2 Sqrt[k] - k) + 2 Sqrt[k]/(-1 + 2 Sqrt[k] - k)) (8 a b k + 4 b^2 k - 8 a b k^(3/2) - 4 b^2 k^(3/2) + 4 a^2 k^2 - 4 a^2 k^(5/2)))/(2/(-1 + 2 Sqrt[k] - k) + 2 k/(-1 + 2 Sqrt[k] - k) - 2 u Sqrt[-1/(-1 + 2 Sqrt[k] - k) - k/(-1 + 2 Sqrt[k] - k) + 2 Sqrt[k]/(-1 + 2 Sqrt[k] - k)])))/((u + (Sqrt[-1/(-1 + 2 Sqrt[k] - k) - k/(-1 + 2 Sqrt[k] - k) + 2 Sqrt[k]/(-1 + 2 Sqrt[k] - k)] (u^2 + 1/(-1 + 2 Sqrt[k] - k) + k/(-1 + 2 Sqrt[k] - k) + 2 Sqrt[k]/(-1 + 2 Sqrt[k] - k)))/(2/(-1 + 2 Sqrt[k] - k) + 2 k/(-1 + 2 Sqrt[k] - k) - 2 u Sqrt[-1/(-1 + 2 Sqrt[k] - k) - k/(-1 + 2 Sqrt[k] - k) + 2 Sqrt[k]/(-1 + 2 Sqrt[k] - k)])) (2/(-1 + 2 Sqrt[k] - k) + 2 k/(-1 + 2 Sqrt[k] - k) - 2 u Sqrt[-1/(-1 + 2 Sqrt[k] - k) - k/(-1 + 2 Sqrt[k] - k) + 2 Sqrt[k]/(-1 + 2 Sqrt[k] - k)])^2 (-b^4 - 2 b^4 Sqrt[k] + 2 a^2 b^2 k - b^4 k + 4 a^2 b^2 k^(3/2) - a^4 k^2 + 2 a^2 b^2 k^2 - 2 a^4 k^(5/2) - a^4 k^3 + ((u^2 + 1/(-1 + 2 Sqrt[k] - k) + k/(-1 + 2 Sqrt[k] - k) + 2 Sqrt[k]/(-1 + 2 Sqrt[k] - k)) (2 b^4 + 12 a^2 b^2 k + 16 a b^3 k + 2 b^4 k + 2 a^4 k^2 + 16 a^3 b k^2 + 12 a^2 b^2 k^2 + 2 a^4 k^3))/(2/(-1 + 2 Sqrt[k] - k) + 2 k/(-1 + 2 Sqrt[k] - k) - 2 u Sqrt[-1/(-1 + 2 Sqrt[k] - k) - k/(-1 + 2 Sqrt[k] - k) + 2 Sqrt[k]/(-1 + 2 Sqrt[k] - k)]) + ((u^2 + 1/(-1 + 2 Sqrt[k] - k) + k/(-1 + 2 Sqrt[k] - k) + 2 Sqrt[k]/(-1 + 2 Sqrt[k] - k))^2 (-b^4 + 2 b^4 Sqrt[k] + 2 a^2 b^2 k - b^4 k - 4 a^2 b^2 k^(3/2) - a^4 k^2 + 2 a^2 b^2 k^2 + 2 a^4 k^(5/2) - a^4 k^3))/(2/(-1 + 2 Sqrt[k] - k) + 2 k/(-1 + 2 Sqrt[k] - k) - 2 u Sqrt[-1/(-1 + 2 Sqrt[k] - k) - k/(-1 + 2 Sqrt[k] - k) + 2 Sqrt[k]/(-1 + 2 Sqrt[k] - k)])^2))");
}

static void test_parametric(void) {
    /* Q(k)(sqrt k): x^2 - k = (x - sqrt k)(x + sqrt k). The GCD is monic in
     * FLINT's variable order, so it is a unit multiple of x - sqrt k. */
    check_gcd_param("x^2 - k", "x - Sqrt[k]", "Sqrt[k] - x");
    check_gcd_param("x^2 - k", "x + Sqrt[k]", "Sqrt[k] + x");
    /* With parameters a, b: the shared factor (x - sqrt k) is extracted; the
     * cofactors (a x + b) and (x + k) are coprime. */
    check_gcd_param("(x - Sqrt[k])*(a x + b)", "(x - Sqrt[k])*(x + k)", "Sqrt[k] - x");
    /* coprime over the field -> unit (1) */
    check_gcd_param("x - Sqrt[k]", "x + 1", "1");
    /* shared repeated radical factor (x - sqrt k)^2 (verified semantically,
     * since the raw unit sign of the GCD is FLINT-order dependent). */
    check_builtin("PossibleZeroQ[Flint`GCD[(x - Sqrt[k])^2*(x + a), "
                  "(x - Sqrt[k])^2*(x + b)] - (x - Sqrt[k])^2]", "True");

    /* Semantic cross-checks through the consumer (canonical surface form). */
    check_builtin("Cancel[(x^2 - k)/(x - Sqrt[k])]", "Sqrt[k] + x");
    check_builtin("PossibleZeroQ[PolynomialGCD[(x - Sqrt[k])*(a x^2 + b x + k), "
                  "(x - Sqrt[k])*(x^3 + a k)] - (Sqrt[k] - x)]", "True");

    /* Bail-outs -> NULL (handled by other paths, or out of scope). */
    check_gcd_param_null("x^2 - 2", "x - Sqrt[2]");   /* integer radicand: nf path */
    check_gcd_param_null("x^2 - 1", "x - 1");         /* no radical                */
    check_gcd_param_null("x - Sqrt[k]", "x - Sqrt[m]"); /* two distinct radicands  */
    check_gcd_param_null("x - k^(1/3)", "x - 1");     /* cube root: higher degree  */
}

/* ------------------------------------------------------------------ */
/*  Builtin path through the evaluator (Flint`GCD[...])               */
/* ------------------------------------------------------------------ */
/*  Parametric radical resultant (the Rothstein-Trager bottleneck)     */
/* ------------------------------------------------------------------ */

static void test_parametric_resultant(void) {
    /* Res_x(x - sqrt k, x + sqrt k) = 2 sqrt k. */
    check_builtin("Resultant[x - Sqrt[k], x + Sqrt[k], x]", "2 Sqrt[k]");
    /* Shared root sqrt k -> resultant 0. */
    check_builtin("Resultant[x^2 - k, x - Sqrt[k], x]", "0");
    check_builtin("Resultant[(x - Sqrt[k])*(x + a), x - Sqrt[k], x]", "0");
    /* The k^(9/2)-coefficient case that timed out (>90 s) under the classical
     * subresultant PRS is now instant via the sqrt(k) -> S collapse; assert it
     * agrees with the classical resultant evaluated at a concrete k. */
    check_builtin("PossibleZeroQ[(Resultant[(a x^2 - k)*(x - Sqrt[k]), "
                  "(b x + k)*(x + Sqrt[k]), x] /. {a -> 2, b -> 3, k -> 5}) "
                  "- Resultant[(2 x^2 - 5)*(x - Sqrt[5]), (3 x + 5)*(x + Sqrt[5]), x]]",
                  "True");
}

/* ------------------------------------------------------------------ */
/*  Parametric radical factoring (the Apart/Factor blocker)            */
/* ------------------------------------------------------------------ */

static void test_parametric_factor(void) {
    /* Splits over the field Q(a,b,sqrt k) — collapse to Q(a,b,S), factor over Q. */
    check_builtin("Factor[u^2 - 2 Sqrt[k] u + k]", "(Sqrt[k] - u)^2");
    check_builtin("Factor[u^2 - (a + b) Sqrt[k] u + a b k]",
                  "(a Sqrt[k] - u) (b Sqrt[k] - u)");
    /* No Sqrt[symbol] present -> parametric hook bails, stays over Q. */
    check_builtin("Factor[u^2 - k]", "-k + u^2");
    check_builtin("Factor[x^2 - 1]", "(-1 + x) (1 + x)");
    /* FactorSquareFree over the field. */
    check_builtin("FactorSquareFree[u^2 - 2 Sqrt[k] u + k]", "(Sqrt[k] - u)^2");
    /* Correctness: the factorisation multiplies back to the input. */
    check_builtin("PossibleZeroQ[Factor[a u^2 - a k + b Sqrt[k] u^2 - b k Sqrt[k]] "
                  "- (a u^2 - a k + b Sqrt[k] u^2 - b k Sqrt[k])]", "True");
}

/* ------------------------------------------------------------------ */
/*  Field xgcd / divrem over Q(params)(sqrt k) via gr_poly<fmpz_mpoly_q> */
/* ------------------------------------------------------------------ */

static void test_parametric_field_ops(void) {
    /* PolynomialExtendedGCD over Q(p,q,sqrt k): Bézout identity u a + v b = g.
     * The classical Euclidean xgcd over the radical coefficients blows up; this
     * routes through gr_poly over the rational-function ring fmpz_mpoly_q. */
    check_builtin("With[{a = x^2 + p Sqrt[k] x + q k, b = x - Sqrt[k]}, "
                  "PossibleZeroQ[Module[{r = PolynomialExtendedGCD[a, b, x]}, "
                  "r[[2,1]] a + r[[2,2]] b - r[[1]]]]]", "True");
    /* Two degree-2 factors with k^(3/2) coefficients (the Goursat Apart case). */
    check_builtin("With[{a = x^2 + (p + q Sqrt[k]) x + (p q k + Sqrt[k]), "
                  "b = x^2 - (p - q Sqrt[k]) x + p k}, "
                  "PossibleZeroQ[Module[{r = PolynomialExtendedGCD[a, b, x]}, "
                  "r[[2,1]] a + r[[2,2]] b - r[[1]]]]]", "True");
    /* Quotient/Remainder over the field: a == q b + r. */
    check_builtin("With[{a = x^3 + p Sqrt[k] x^2 + q k, b = x^2 + Sqrt[k] x - k}, "
                  "PossibleZeroQ[Expand[PolynomialQuotient[a, b, x] b "
                  "+ PolynomialRemainder[a, b, x] - a]]]", "True");
    /* Plain-rational xgcd unchanged (hook bails, classical path). */
    check_builtin("PolynomialExtendedGCD[x^2 - 1, x - 1, x]", "{-1 + x, {0, 1}}");

    /* Field-normalize: Cancel/Together of a Plus of fractions over Q(a,b,k)(sqrt k)
     * — the case where extract-num/den saw den=1 and dropped to the slow QA path
     * (the LogToReal A^2+B^2 bottleneck: 3.5s -> 0.0004s). Value preserved. */
    check_builtin("PossibleZeroQ[Cancel[b^4/(b^2 + 2 a b Sqrt[k] + a^2 k)^2 "
                  "+ a^4 k^2/(b^2 + 2 a b Sqrt[k] + a^2 k)^2, Extension -> Automatic] "
                  "- (b^4 + a^4 k^2)/(b^2 + 2 a b Sqrt[k] + a^2 k)^2]", "True");
    check_builtin("PossibleZeroQ[Together[1/(x - Sqrt[k]) + 1/(x + a Sqrt[k]), "
                  "Extension -> Automatic] - (1/(x - Sqrt[k]) + 1/(x + a Sqrt[k]))]", "True");
}

/* ------------------------------------------------------------------ */

static void test_builtin_path(void) {
    check_builtin("Flint`GCD[(x+y)*(x-2 y), (x+y)*(x+3 y)]", "x + y");
    check_builtin("Flint`GCD[x^2 - 1, x^2 + 2 x + 1]", "1 + x");
    check_builtin("Flint`GCD[x^2 - 1/4, x - 1/2]", "-1/2 + x");
    check_builtin("Flint`GCD[x^3 y - x y^3, x^2 y^2 - y^4]", "x^2 y - y^3");
    check_builtin("Flint`GCD[x^2 - 2, x - 1]", "1");
    /* number-field path is reached through the same builtin */
    check_builtin("Flint`GCD[x^2 - 2, x - Sqrt[2]]", "-Sqrt[2] + x");
    check_builtin("Flint`GCD[Sqrt[2] x, x]", "x");
    /* out-of-scope: stays unevaluated (head printed back, deferring to other
     * paths) — never silently wrong */
    check_builtin("Flint`GCD[6, 4]", "Flint`GCD[6, 4]");
    /* listability is not claimed; wrong arity stays unevaluated */
    check_builtin("Flint`GCD[x^2 - 1]", "Flint`GCD[-1 + x^2]");
}

/* ------------------------------------------------------------------ */
/*  Consumers: PolynomialGCD / Cancel / Together / Apart now reduce     */
/*  over the algebraic-extension fields (the payoff of wiring the       */
/*  bridge into builtin_polynomialgcd and Cancel's divide-back).        */
/* ------------------------------------------------------------------ */

static void test_consumers(void) {
    /* PolynomialGCD is FLINT-backed for the extension cases. */
    check_builtin("PolynomialGCD[x^2 - 2, x - Sqrt[2]]", "-Sqrt[2] + x");
    check_builtin("PolynomialGCD[x^2 + x + 1, x - (-1)^(2/3)]", "1 - (-1)^(1/3) + x");

    /* Cancel now reduces over Q(sqrt d) — previously returned unreduced. */
    check_builtin("Cancel[(x^2 - 2)/(x - Sqrt[2])]", "Sqrt[2] + x");
    check_builtin("Cancel[(x^2 - 3)/(x + Sqrt[3])]", "-Sqrt[3] + x");
    /* ... over Q(zeta_n) ... */
    check_builtin("Cancel[(x^2 + x + 1)/(x - (-1)^(2/3))]", "(-1)^(1/3) + x");
    /* ... and over a radical tower. */
    check_builtin("Cancel[((x - Sqrt[2] - Sqrt[3])*(x - 1))/"
                  "((x - Sqrt[2] - Sqrt[3])*(x - 2))]", "(-1 + x)/(-2 + x)");

    /* Together over Q(sqrt 2). */
    check_builtin("Together[1/(x - Sqrt[2]) + 1/(x + Sqrt[2])]", "(2 x)/(-2 + x^2)");

    /* Semantic check independent of surface form. */
    check_builtin("PossibleZeroQ[Cancel[(x^2 - 2)/(x - Sqrt[2])] - (x + Sqrt[2])]", "True");

    /* Plain-rational cases are unchanged (no FLINT extension path taken). */
    check_builtin("Cancel[(x^2 - 1)/(x - 1)]", "1 + x");
    check_builtin("Cancel[(2 x + 2)/(4 x + 4)]", "1/2");
    check_builtin("PolynomialGCD[x^2 - 1, x^2 + 2 x + 1]", "1 + x");
}

/* ------------------------------------------------------------------ */
/*  Attributes                                                         */
/* ------------------------------------------------------------------ */

static void test_attributes(void) {
    check_builtin("MemberQ[Attributes[Flint`GCD], Protected]", "True");
}

/* ------------------------------------------------------------------ */

int main(void) {
    symtab_init();
    core_init();

    if (!flint_bridge_available()) {
        printf("FLINT not compiled in (USE_FLINT off); skipping bridge tests.\n");
        return 0;
    }

    TEST(test_univariate);
    TEST(test_multivariate);
    TEST(test_monic_normalisation);
    TEST(test_rational_coeffs);
    TEST(test_bigint_coeffs);
    TEST(test_zero_and_identity);
    TEST(test_determinism);
    TEST(test_bail_out_of_scope);
    TEST(test_numberfield);
    TEST(test_cyclotomic);
    TEST(test_tower);
    TEST(test_parametric);
    TEST(test_parametric_resultant);
    TEST(test_parametric_factor);
    TEST(test_parametric_field_ops);
    TEST(test_goursat_descent_operands);
    TEST(test_builtin_path);
    TEST(test_consumers);
    TEST(test_attributes);

    printf("All FLINT bridge tests passed!\n");
    return 0;
}
