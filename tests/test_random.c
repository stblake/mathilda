#include "expr.h"
#include "eval.h"
#include "core.h"
#include "symtab.h"
#include "test_utils.h"
#include "parse.h"
#include "print.h"
#include <string.h>
#include <stdlib.h>
#include <math.h>

/* Helper: parse, evaluate, return real value. Returns NAN on error. */
static double eval_to_real(const char* input) {
    Expr* e = parse_expression(input);
    Expr* r = evaluate(e);
    double val = NAN;
    if (r->type == EXPR_REAL) val = r->data.real;
    else if (r->type == EXPR_INTEGER) val = (double)r->data.integer;
    expr_free(e);
    expr_free(r);
    return val;
}

/* Helper: parse, evaluate, return integer value. Returns -9999 on error. */
static int64_t eval_to_int(const char* input) {
    Expr* e = parse_expression(input);
    Expr* r = evaluate(e);
    int64_t val = -9999;
    if (r->type == EXPR_INTEGER) val = r->data.integer;
    expr_free(e);
    expr_free(r);
    return val;
}

/* Helper: parse, evaluate, return string (caller must free). */
static char* eval_to_str(const char* input) {
    Expr* e = parse_expression(input);
    Expr* r = evaluate(e);
    char* s = expr_to_string(r);
    expr_free(e);
    expr_free(r);
    return s;
}

/* ---- RandomInteger[] ---- */
void test_randominteger_no_args(void) {
    /* RandomInteger[] should return 0 or 1 */
    for (int i = 0; i < 100; i++) {
        int64_t v = eval_to_int("RandomInteger[]");
        ASSERT_MSG(v == 0 || v == 1,
                   "RandomInteger[] returned %lld, expected 0 or 1",
                   (long long)v);
    }
}

/* ---- RandomInteger[imax] ---- */
void test_randominteger_single_arg(void) {
    /* RandomInteger[0] should always return 0 */
    for (int i = 0; i < 20; i++) {
        int64_t v = eval_to_int("RandomInteger[0]");
        ASSERT_MSG(v == 0, "RandomInteger[0] returned %lld", (long long)v);
    }

    /* RandomInteger[10] should return values in [0, 10] */
    for (int i = 0; i < 200; i++) {
        int64_t v = eval_to_int("RandomInteger[10]");
        ASSERT_MSG(v >= 0 && v <= 10,
                   "RandomInteger[10] returned %lld", (long long)v);
    }

    /* RandomInteger[1] should return 0 or 1 */
    int seen0 = 0, seen1 = 0;
    for (int i = 0; i < 200; i++) {
        int64_t v = eval_to_int("RandomInteger[1]");
        ASSERT_MSG(v == 0 || v == 1,
                   "RandomInteger[1] returned %lld", (long long)v);
        if (v == 0) seen0 = 1;
        if (v == 1) seen1 = 1;
    }
    ASSERT_MSG(seen0 && seen1, "RandomInteger[1] never produced both 0 and 1");
}

/* ---- RandomInteger[{imin, imax}] ---- */
void test_randominteger_range(void) {
    /* RandomInteger[{5, 5}] should always be 5 */
    for (int i = 0; i < 20; i++) {
        int64_t v = eval_to_int("RandomInteger[{5, 5}]");
        ASSERT_MSG(v == 5, "RandomInteger[{5,5}] returned %lld", (long long)v);
    }

    /* RandomInteger[{-3, 3}] should be in [-3, 3] */
    for (int i = 0; i < 200; i++) {
        int64_t v = eval_to_int("RandomInteger[{-3, 3}]");
        ASSERT_MSG(v >= -3 && v <= 3,
                   "RandomInteger[{-3,3}] returned %lld", (long long)v);
    }

    /* RandomInteger[{10, 20}] should be in [10, 20] */
    for (int i = 0; i < 200; i++) {
        int64_t v = eval_to_int("RandomInteger[{10, 20}]");
        ASSERT_MSG(v >= 10 && v <= 20,
                   "RandomInteger[{10,20}] returned %lld", (long long)v);
    }
}

/* ---- RandomInteger[range, n] - flat list ---- */
void test_randominteger_flat_list(void) {
    /* RandomInteger[{0, 9}, 5] should give a list of 5 elements */
    char* s = eval_to_str("Length[RandomInteger[{0, 9}, 5]]");
    ASSERT_STR_EQ(s, "5");
    free(s);

    /* All elements should be integers in [0, 9] */
    /* Use Table-style check: Min and Max of the list */
    s = eval_to_str("Min[RandomInteger[{0, 9}, 1000]] >= 0");
    ASSERT_STR_EQ(s, "True");
    free(s);

    s = eval_to_str("Max[RandomInteger[{0, 9}, 1000]] <= 9");
    ASSERT_STR_EQ(s, "True");
    free(s);

    /* RandomInteger[5, 0] should give an empty list */
    s = eval_to_str("RandomInteger[5, 0]");
    ASSERT_STR_EQ(s, "{}");
    free(s);
}

/* ---- RandomInteger[range, {n1, n2}] - multi-dimensional ---- */
void test_randominteger_multidim(void) {
    /* RandomInteger[{0, 1}, {3, 4}] should give a 3x4 matrix */
    char* s = eval_to_str("Dimensions[RandomInteger[{0, 1}, {3, 4}]]");
    ASSERT_STR_EQ(s, "{3, 4}");
    free(s);

    /* RandomInteger[{0, 1}, {2, 3, 4}] should give a 2x3x4 tensor */
    s = eval_to_str("Dimensions[RandomInteger[{0, 1}, {2, 3, 4}]]");
    ASSERT_STR_EQ(s, "{2, 3, 4}");
    free(s);

    /* All values in range */
    s = eval_to_str("Min[Flatten[RandomInteger[{1, 6}, {10, 10}]]] >= 1");
    ASSERT_STR_EQ(s, "True");
    free(s);

    s = eval_to_str("Max[Flatten[RandomInteger[{1, 6}, {10, 10}]]] <= 6");
    ASSERT_STR_EQ(s, "True");
    free(s);
}

/* ---- SeedRandom reproducibility ---- */
void test_seedrandom(void) {
    /* SeedRandom[42] should make RandomInteger deterministic */
    assert_eval_eq("SeedRandom[42]", "Null", 0);

    /* Generate a sequence */
    int64_t seq1[10];
    for (int i = 0; i < 10; i++) {
        seq1[i] = eval_to_int("RandomInteger[{0, 1000}]");
    }

    /* Re-seed with same value */
    assert_eval_eq("SeedRandom[42]", "Null", 0);

    /* Should get same sequence */
    for (int i = 0; i < 10; i++) {
        int64_t v = eval_to_int("RandomInteger[{0, 1000}]");
        ASSERT_MSG(v == seq1[i],
                   "SeedRandom reproducibility failed at index %d: %lld != %lld",
                   i, (long long)v, (long long)seq1[i]);
    }
}

/* ---- SeedRandom[] reseeds from entropy ---- */
void test_seedrandom_no_args(void) {
    assert_eval_eq("SeedRandom[]", "Null", 0);
}

/* ---- Symbolic arguments remain unevaluated ---- */
void test_randominteger_symbolic(void) {
    assert_eval_eq("RandomInteger[x]", "RandomInteger[x]", 0);
    assert_eval_eq("RandomInteger[{a, b}]", "RandomInteger[{a, b}]", 0);
    assert_eval_eq("RandomInteger[5, x]", "RandomInteger[5, x]", 0);
}

/* ---- SeedRandom returns Null ---- */
void test_seedrandom_returns_null(void) {
    assert_eval_eq("SeedRandom[123]", "Null", 0);
}

/* ---- Bignum range ---- */
void test_randominteger_bignum(void) {
    /* RandomInteger with a large range that exceeds int64 */
    /* SeedRandom for determinism, then just check it doesn't crash
       and returns an integer */
    Expr* e = parse_expression("SeedRandom[99]");
    Expr* r = evaluate(e);
    expr_free(e);
    expr_free(r);

    /* 10^20 is a bignum */
    e = parse_expression("IntegerQ[RandomInteger[10^20]]");
    r = evaluate(e);
    char* s = expr_to_string(r);
    ASSERT_STR_EQ(s, "True");
    free(s);
    expr_free(e);
    expr_free(r);

    /* Range with bignums: RandomInteger[{-10^20, 10^20}] */
    e = parse_expression("IntegerQ[RandomInteger[{-10^20, 10^20}]]");
    r = evaluate(e);
    s = expr_to_string(r);
    ASSERT_STR_EQ(s, "True");
    free(s);
    expr_free(e);
    expr_free(r);
}

/* ---- Negative range: RandomInteger[{-10, -5}] ---- */
void test_randominteger_negative_range(void) {
    for (int i = 0; i < 100; i++) {
        int64_t v = eval_to_int("RandomInteger[{-10, -5}]");
        ASSERT_MSG(v >= -10 && v <= -5,
                   "RandomInteger[{-10,-5}] returned %lld", (long long)v);
    }
}

/* ---- Coverage check: all values in a small range appear ---- */
void test_randominteger_coverage(void) {
    /* With SeedRandom and enough trials, RandomInteger[{1,3}]
       should produce 1, 2, and 3 */
    int seen[4] = {0};
    for (int i = 0; i < 200; i++) {
        int64_t v = eval_to_int("RandomInteger[{1, 3}]");
        ASSERT_MSG(v >= 1 && v <= 3,
                   "RandomInteger[{1,3}] returned %lld", (long long)v);
        seen[v] = 1;
    }
    ASSERT_MSG(seen[1] && seen[2] && seen[3],
               "Not all values 1,2,3 were produced by RandomInteger[{1,3}]");
}

/* ---- Attributes check ---- */
void test_randominteger_attributes(void) {
    /* RandomInteger should be Protected */
    char* s = eval_to_str("MemberQ[Attributes[RandomInteger], Protected]");
    ASSERT_STR_EQ(s, "True");
    free(s);

    s = eval_to_str("MemberQ[Attributes[SeedRandom], Protected]");
    ASSERT_STR_EQ(s, "True");
    free(s);
}

/* ---- Information/docstring ---- */
void test_randominteger_info(void) {
    /* Verify ?RandomInteger returns something (not unevaluated) */
    Expr* e = parse_expression("Information[RandomInteger]");
    Expr* r = evaluate(e);
    /* Information returns a string */
    ASSERT(r->type == EXPR_STRING);
    expr_free(e);
    expr_free(r);
}

/* ---- Edge: RandomInteger[range, {n}] is a flat list ---- */
void test_randominteger_single_dim_list(void) {
    char* s = eval_to_str("Dimensions[RandomInteger[{0, 1}, {5}]]");
    ASSERT_STR_EQ(s, "{5}");
    free(s);
}

/* ---- RandomReal[] ---- */
void test_randomreal_no_args(void) {
    /* RandomReal[] should return a real in [0, 1) */
    for (int i = 0; i < 100; i++) {
        double v = eval_to_real("RandomReal[]");
        ASSERT_MSG(!isnan(v), "RandomReal[] did not return a number");
        ASSERT_MSG(v >= 0.0 && v < 1.0,
                   "RandomReal[] returned %g, expected in [0, 1)", v);
    }
}

/* ---- RandomReal[xmax] ---- */
void test_randomreal_single_arg(void) {
    /* RandomReal[0] should always return 0.0 */
    for (int i = 0; i < 20; i++) {
        double v = eval_to_real("RandomReal[0]");
        ASSERT_MSG(v == 0.0, "RandomReal[0] returned %g", v);
    }

    /* RandomReal[10] should return values in [0, 10) */
    for (int i = 0; i < 200; i++) {
        double v = eval_to_real("RandomReal[10]");
        ASSERT_MSG(v >= 0.0 && v < 10.0,
                   "RandomReal[10] returned %g", v);
    }

    /* RandomReal[1.5] should return values in [0, 1.5) */
    for (int i = 0; i < 100; i++) {
        double v = eval_to_real("RandomReal[1.5]");
        ASSERT_MSG(v >= 0.0 && v < 1.5,
                   "RandomReal[1.5] returned %g", v);
    }
}

/* ---- RandomReal[{xmin, xmax}] ---- */
void test_randomreal_range(void) {
    /* RandomReal[{5.0, 5.0}] should always be 5.0 */
    for (int i = 0; i < 20; i++) {
        double v = eval_to_real("RandomReal[{5.0, 5.0}]");
        ASSERT_MSG(v == 5.0, "RandomReal[{5.0,5.0}] returned %g", v);
    }

    /* RandomReal[{-3, 3}] should be in [-3, 3) */
    for (int i = 0; i < 200; i++) {
        double v = eval_to_real("RandomReal[{-3, 3}]");
        ASSERT_MSG(v >= -3.0 && v < 3.0,
                   "RandomReal[{-3,3}] returned %g", v);
    }

    /* RandomReal[{10, 20}] should be in [10, 20) */
    for (int i = 0; i < 200; i++) {
        double v = eval_to_real("RandomReal[{10, 20}]");
        ASSERT_MSG(v >= 10.0 && v < 20.0,
                   "RandomReal[{10,20}] returned %g", v);
    }
}

/* ---- RandomReal with rational range ---- */
void test_randomreal_rational_range(void) {
    /* RandomReal[{1/4, 3/4}] should be in [0.25, 0.75) */
    for (int i = 0; i < 100; i++) {
        double v = eval_to_real("RandomReal[{1/4, 3/4}]");
        ASSERT_MSG(v >= 0.25 && v < 0.75,
                   "RandomReal[{1/4,3/4}] returned %g", v);
    }
}

/* ---- RandomReal[range, n] - flat list ---- */
void test_randomreal_flat_list(void) {
    /* RandomReal[{0, 1}, 5] should give a list of 5 elements */
    char* s = eval_to_str("Length[RandomReal[{0, 1}, 5]]");
    ASSERT_STR_EQ(s, "5");
    free(s);

    /* RandomReal[5, 0] should give an empty list */
    s = eval_to_str("RandomReal[5, 0]");
    ASSERT_STR_EQ(s, "{}");
    free(s);

    /* All elements should be reals (Head == Real) */
    s = eval_to_str("SeedRandom[99]; Head[RandomReal[{0, 1}, 3][[1]]]");
    ASSERT_STR_EQ(s, "Real");
    free(s);
}

/* ---- RandomReal[range, {n1, n2}] - multi-dimensional ---- */
void test_randomreal_multidim(void) {
    /* RandomReal[{0, 1}, {3, 4}] should give a 3x4 matrix */
    char* s = eval_to_str("Dimensions[RandomReal[{0, 1}, {3, 4}]]");
    ASSERT_STR_EQ(s, "{3, 4}");
    free(s);

    /* RandomReal[{0, 1}, {2, 3, 4}] should give a 2x3x4 tensor */
    s = eval_to_str("Dimensions[RandomReal[{0, 1}, {2, 3, 4}]]");
    ASSERT_STR_EQ(s, "{2, 3, 4}");
    free(s);
}

/* ---- SeedRandom reproducibility with RandomReal ---- */
void test_randomreal_seedrandom(void) {
    /* SeedRandom[42] should make RandomReal deterministic */
    assert_eval_eq("SeedRandom[42]", "Null", 0);

    /* Generate a sequence */
    double seq1[10];
    for (int i = 0; i < 10; i++) {
        seq1[i] = eval_to_real("RandomReal[]");
    }

    /* Re-seed with same value */
    assert_eval_eq("SeedRandom[42]", "Null", 0);

    /* Should get same sequence */
    for (int i = 0; i < 10; i++) {
        double v = eval_to_real("RandomReal[]");
        ASSERT_MSG(v == seq1[i],
                   "SeedRandom reproducibility failed at index %d: %g != %g",
                   i, v, seq1[i]);
    }
}

/* ---- Symbolic arguments remain unevaluated ---- */
void test_randomreal_symbolic(void) {
    assert_eval_eq("RandomReal[x]", "RandomReal[x]", 0);
    assert_eval_eq("RandomReal[{a, b}]", "RandomReal[{a, b}]", 0);
    assert_eval_eq("RandomReal[5, x]", "RandomReal[5, x]", 0);
}

/* ---- Attributes check ---- */
void test_randomreal_attributes(void) {
    char* s = eval_to_str("MemberQ[Attributes[RandomReal], Protected]");
    ASSERT_STR_EQ(s, "True");
    free(s);
}

/* ---- Information/docstring ---- */
void test_randomreal_info(void) {
    Expr* e = parse_expression("Information[RandomReal]");
    Expr* r = evaluate(e);
    ASSERT(r->type == EXPR_STRING);
    expr_free(e);
    expr_free(r);
}

/* ---- RandomReal[{xmin, xmax}] with negative range ---- */
void test_randomreal_negative_range(void) {
    for (int i = 0; i < 100; i++) {
        double v = eval_to_real("RandomReal[{-10, -5}]");
        ASSERT_MSG(v >= -10.0 && v < -5.0,
                   "RandomReal[{-10,-5}] returned %g", v);
    }
}

/* ---- Coverage: values spread across range ---- */
void test_randomreal_coverage(void) {
    /* With enough trials, RandomReal[{0,1}] should produce values
       in both [0, 0.5) and [0.5, 1) */
    int seen_low = 0, seen_high = 0;
    for (int i = 0; i < 200; i++) {
        double v = eval_to_real("RandomReal[]");
        if (v < 0.5) seen_low = 1;
        if (v >= 0.5) seen_high = 1;
    }
    ASSERT_MSG(seen_low && seen_high,
               "RandomReal[] did not produce values in both halves of [0,1)");
}

/* ---- RandomReal[{xmin, xmax}] with single dim list ---- */
void test_randomreal_single_dim_list(void) {
    char* s = eval_to_str("Dimensions[RandomReal[{0, 1}, {5}]]");
    ASSERT_STR_EQ(s, "{5}");
    free(s);
}

/* ---- RandomReal always returns EXPR_REAL type ---- */
void test_randomreal_type(void) {
    /* Even with integer range args, result should be Real */
    char* s = eval_to_str("Head[RandomReal[{0, 10}]]");
    ASSERT_STR_EQ(s, "Real");
    free(s);

    s = eval_to_str("Head[RandomReal[]]");
    ASSERT_STR_EQ(s, "Real");
    free(s);
}

/* ==== RandomComplex tests ==== */

/* Helper: parse, evaluate, extract Re and Im parts as doubles.
 * Returns true if the result is Complex[re, im] with real parts. */
static bool eval_to_complex(const char* input, double* re, double* im) {
    Expr* e = parse_expression(input);
    Expr* r = evaluate(e);
    /* Check if it's Complex[re, im] */
    if (r->type == EXPR_FUNCTION &&
        r->data.function.head->type == EXPR_SYMBOL &&
        strcmp(r->data.function.head->data.symbol.name, "Complex") == 0 &&
        r->data.function.arg_count == 2) {
        Expr* re_expr = r->data.function.args[0];
        Expr* im_expr = r->data.function.args[1];
        if (re_expr->type == EXPR_REAL) *re = re_expr->data.real;
        else if (re_expr->type == EXPR_INTEGER) *re = (double)re_expr->data.integer;
        else { expr_free(e); expr_free(r); return false; }
        if (im_expr->type == EXPR_REAL) *im = im_expr->data.real;
        else if (im_expr->type == EXPR_INTEGER) *im = (double)im_expr->data.integer;
        else { expr_free(e); expr_free(r); return false; }
        expr_free(e);
        expr_free(r);
        return true;
    }
    expr_free(e);
    expr_free(r);
    return false;
}

void test_randomcomplex_no_args(void) {
    /* RandomComplex[] should return Complex with re,im in [0, 1) */
    for (int i = 0; i < 100; i++) {
        double re, im;
        bool ok = eval_to_complex("RandomComplex[]", &re, &im);
        ASSERT_MSG(ok, "RandomComplex[] did not return a Complex");
        ASSERT_MSG(re >= 0.0 && re < 1.0,
                   "RandomComplex[] re=%g out of [0,1)", re);
        ASSERT_MSG(im >= 0.0 && im < 1.0,
                   "RandomComplex[] im=%g out of [0,1)", im);
    }
}

void test_randomcomplex_single_complex_arg(void) {
    /* RandomComplex[2 + 3 I] -> re in [0,2), im in [0,3) */
    for (int i = 0; i < 100; i++) {
        double re, im;
        bool ok = eval_to_complex("RandomComplex[2 + 3 I]", &re, &im);
        ASSERT_MSG(ok, "RandomComplex[2+3I] did not return a Complex");
        ASSERT_MSG(re >= 0.0 && re < 2.0,
                   "RandomComplex[2+3I] re=%g out of [0,2)", re);
        ASSERT_MSG(im >= 0.0 && im < 3.0,
                   "RandomComplex[2+3I] im=%g out of [0,3)", im);
    }
}

void test_randomcomplex_range(void) {
    /* RandomComplex[{1 + 2 I, 3 + 5 I}] -> re in [1,3), im in [2,5) */
    for (int i = 0; i < 100; i++) {
        double re, im;
        bool ok = eval_to_complex("RandomComplex[{1 + 2 I, 3 + 5 I}]", &re, &im);
        ASSERT_MSG(ok, "RandomComplex[{1+2I,3+5I}] did not return a Complex");
        ASSERT_MSG(re >= 1.0 && re < 3.0,
                   "RandomComplex range re=%g out of [1,3)", re);
        ASSERT_MSG(im >= 2.0 && im < 5.0,
                   "RandomComplex range im=%g out of [2,5)", im);
    }
}

void test_randomcomplex_negative_range(void) {
    /* RandomComplex[{-2 - 3 I, -1 - 1 I}] -> re in [-2,-1), im in [-3,-1) */
    for (int i = 0; i < 100; i++) {
        double re, im;
        bool ok = eval_to_complex("RandomComplex[{-2 - 3 I, -1 - I}]", &re, &im);
        ASSERT_MSG(ok, "RandomComplex negative range did not return Complex");
        ASSERT_MSG(re >= -2.0 && re < -1.0,
                   "RandomComplex neg range re=%g", re);
        ASSERT_MSG(im >= -3.0 && im < -1.0,
                   "RandomComplex neg range im=%g", im);
    }
}

void test_randomcomplex_flat_list(void) {
    /* RandomComplex[{0, 1 + I}, 5] should give a list of 5 elements */
    char* s = eval_to_str("Length[RandomComplex[{0, 1 + I}, 5]]");
    ASSERT_STR_EQ(s, "5");
    free(s);

    /* RandomComplex[{0, 1 + I}, 0] should give an empty list */
    s = eval_to_str("RandomComplex[{0, 1 + I}, 0]");
    ASSERT_STR_EQ(s, "{}");
    free(s);
}

void test_randomcomplex_multidim(void) {
    /* RandomComplex[{0, 1+I}, {3, 4}] should give a 3x4 matrix */
    char* s = eval_to_str("Dimensions[RandomComplex[{0, 1 + I}, {3, 4}]]");
    ASSERT_STR_EQ(s, "{3, 4}");
    free(s);

    /* RandomComplex[{0, 1+I}, {2, 3, 4}] should give a 2x3x4 tensor */
    s = eval_to_str("Dimensions[RandomComplex[{0, 1 + I}, {2, 3, 4}]]");
    ASSERT_STR_EQ(s, "{2, 3, 4}");
    free(s);
}

void test_randomcomplex_seedrandom(void) {
    /* SeedRandom should make RandomComplex deterministic */
    assert_eval_eq("SeedRandom[42]", "Null", 0);
    double re1[5], im1[5];
    for (int i = 0; i < 5; i++) {
        bool ok = eval_to_complex("RandomComplex[]", &re1[i], &im1[i]);
        ASSERT(ok);
    }
    assert_eval_eq("SeedRandom[42]", "Null", 0);
    for (int i = 0; i < 5; i++) {
        double re, im;
        bool ok = eval_to_complex("RandomComplex[]", &re, &im);
        ASSERT(ok);
        ASSERT_MSG(re == re1[i] && im == im1[i],
                   "RandomComplex SeedRandom reproducibility failed at %d", i);
    }
}

void test_randomcomplex_symbolic(void) {
    assert_eval_eq("RandomComplex[x]", "RandomComplex[x]", 0);
    assert_eval_eq("RandomComplex[{a, b}]", "RandomComplex[{a, b}]", 0);
    assert_eval_eq("RandomComplex[{0, 1 + I}, x]", "RandomComplex[{0, 1 + I}, x]", 0);
}

void test_randomcomplex_attributes(void) {
    char* s = eval_to_str("MemberQ[Attributes[RandomComplex], Protected]");
    ASSERT_STR_EQ(s, "True");
    free(s);
}

void test_randomcomplex_info(void) {
    Expr* e = parse_expression("Information[RandomComplex]");
    Expr* r = evaluate(e);
    ASSERT(r->type == EXPR_STRING);
    expr_free(e);
    expr_free(r);
}

void test_randomcomplex_head(void) {
    /* RandomComplex[] should have head Complex */
    char* s = eval_to_str("Head[RandomComplex[]]");
    ASSERT_STR_EQ(s, "Complex");
    free(s);
}

void test_randomcomplex_coverage(void) {
    /* Both real and imaginary parts should spread across [0,1) */
    int re_low = 0, re_high = 0, im_low = 0, im_high = 0;
    for (int i = 0; i < 200; i++) {
        double re, im;
        bool ok = eval_to_complex("RandomComplex[]", &re, &im);
        ASSERT(ok);
        if (re < 0.5) re_low = 1; else re_high = 1;
        if (im < 0.5) im_low = 1; else im_high = 1;
    }
    ASSERT_MSG(re_low && re_high, "RandomComplex re not spread");
    ASSERT_MSG(im_low && im_high, "RandomComplex im not spread");
}

void test_randomcomplex_real_range(void) {
    /* RandomComplex[{1, 3}] with real corners -> Complex[re, 0.0] which
     * simplifies to just a real. Check that the result is a real in [1,3). */
    for (int i = 0; i < 50; i++) {
        double v = eval_to_real("RandomComplex[{1, 3}]");
        ASSERT_MSG(!isnan(v), "RandomComplex[{1,3}] did not return a number");
        ASSERT_MSG(v >= 1.0 && v < 3.0,
                   "RandomComplex[{1,3}] returned %g", v);
    }
}

void test_randomcomplex_single_dim_list(void) {
    char* s = eval_to_str("Dimensions[RandomComplex[{0, 1 + I}, {5}]]");
    ASSERT_STR_EQ(s, "{5}");
    free(s);
}

/* ==== WorkingPrecision tests ==== */

#ifdef USE_MPFR
#include "expr.h"

/* Evaluate an expression and return the result's ExprType, so callers can
 * verify that WorkingPrecision -> n produced an EXPR_MPFR atom. The
 * input expression must reduce to an atom (not a function). */
static int eval_to_type(const char* input) {
    Expr* e = parse_expression(input);
    Expr* r = evaluate(e);
    int t = (int)r->type;
    expr_free(e);
    expr_free(r);
    return t;
}

/* RandomReal[1, WorkingPrecision -> n] returns an MPFR atom whose
 * Precision is within 1 digit of n (Precision returns digits computed
 * from the MPFR bit count, which rounds up). */
void test_randomreal_workingprecision_basic(void) {
    /* Atom is MPFR-typed */
    ASSERT(eval_to_type("RandomReal[1, WorkingPrecision -> 30]") == (int)EXPR_MPFR);
    ASSERT(eval_to_type("RandomReal[WorkingPrecision -> 50]") == (int)EXPR_MPFR);

    /* Precision is at least the requested digit count. */
    char* s = eval_to_str("Precision[RandomReal[1, WorkingPrecision -> 30]] >= 30");
    ASSERT_STR_EQ(s, "True"); free(s);
    s = eval_to_str("Precision[RandomReal[1, WorkingPrecision -> 100]] >= 100");
    ASSERT_STR_EQ(s, "True"); free(s);

    /* Value remains in range. */
    for (int i = 0; i < 20; i++) {
        char* t = eval_to_str(
            "Module[{x = RandomReal[WorkingPrecision -> 40]},\n"
            "  And[x >= 0, x < 1]]");
        ASSERT_STR_EQ(t, "True"); free(t);
    }
}

void test_randomreal_workingprecision_range(void) {
    /* Bounds from a rational range should be honored exactly. */
    for (int i = 0; i < 20; i++) {
        char* t = eval_to_str(
            "Module[{x = RandomReal[{1/3, 2/3}, WorkingPrecision -> 40]},\n"
            "  And[x >= 1/3, x < 2/3]]");
        ASSERT_STR_EQ(t, "True"); free(t);
    }
}

void test_randomreal_workingprecision_list(void) {
    /* Returns a list of MPFR atoms at the requested precision. */
    char* s = eval_to_str("Length[RandomReal[1, 5, WorkingPrecision -> 40]]");
    ASSERT_STR_EQ(s, "5"); free(s);

    s = eval_to_str(
        "Module[{xs = RandomReal[1, 4, WorkingPrecision -> 50]},\n"
        "  And @@ ((# >= 50)& /@ (Precision /@ xs))]");
    ASSERT_STR_EQ(s, "True"); free(s);
}

void test_randomreal_workingprecision_array(void) {
    char* s = eval_to_str(
        "Dimensions[RandomReal[1, {2, 3}, WorkingPrecision -> 40]]");
    ASSERT_STR_EQ(s, "{2, 3}"); free(s);
}

void test_randomreal_machineprecision_option(void) {
    /* WorkingPrecision -> MachinePrecision (or a value at/below machine
     * precision) returns a plain EXPR_REAL. */
    ASSERT(eval_to_type("RandomReal[1, WorkingPrecision -> MachinePrecision]") == (int)EXPR_REAL);
    ASSERT(eval_to_type("RandomReal[1, WorkingPrecision -> 10]") == (int)EXPR_REAL);
}

void test_randomcomplex_workingprecision_basic(void) {
    /* Complex with MPFR real and imaginary components. */
    char* s = eval_to_str(
        "Module[{z = RandomComplex[1 + I, WorkingPrecision -> 40]},\n"
        "  {Head[z], Precision[Re[z]] >= 40, Precision[Im[z]] >= 40}]");
    ASSERT_STR_EQ(s, "{Complex, True, True}"); free(s);
}

void test_randomcomplex_workingprecision_list(void) {
    char* s = eval_to_str(
        "Length[RandomComplex[1 + I, 4, WorkingPrecision -> 30]]");
    ASSERT_STR_EQ(s, "4"); free(s);

    s = eval_to_str(
        "Dimensions[RandomComplex[1 + I, {2, 3}, WorkingPrecision -> 30]]");
    ASSERT_STR_EQ(s, "{2, 3}"); free(s);
}
#endif /* USE_MPFR */

/* ==== RandomChoice tests ==== */

void test_randomchoice_basic(void) {
    /* RandomChoice[{a, b, c}] should return one of a, b, c */
    for (int i = 0; i < 100; i++) {
        char* s = eval_to_str("RandomChoice[{a, b, c}]");
        ASSERT_MSG(strcmp(s, "a") == 0 || strcmp(s, "b") == 0 || strcmp(s, "c") == 0,
                   "RandomChoice[{a,b,c}] returned %s", s);
        free(s);
    }
}

void test_randomchoice_single_element(void) {
    /* RandomChoice[{42}] should always return 42 */
    for (int i = 0; i < 20; i++) {
        int64_t v = eval_to_int("RandomChoice[{42}]");
        ASSERT_MSG(v == 42, "RandomChoice[{42}] returned %lld", (long long)v);
    }
}

void test_randomchoice_flat_list(void) {
    /* RandomChoice[{a, b, c}, 5] should give a list of 5 elements */
    char* s = eval_to_str("Length[RandomChoice[{a, b, c}, 5]]");
    ASSERT_STR_EQ(s, "5");
    free(s);

    /* RandomChoice[{a, b}, 0] should give empty list */
    s = eval_to_str("RandomChoice[{a, b}, 0]");
    ASSERT_STR_EQ(s, "{}");
    free(s);
}

void test_randomchoice_multidim(void) {
    /* RandomChoice[{a, b, c}, {3, 4}] should give a 3x4 matrix */
    char* s = eval_to_str("Dimensions[RandomChoice[{a, b, c}, {3, 4}]]");
    ASSERT_STR_EQ(s, "{3, 4}");
    free(s);

    /* RandomChoice[{1, 2}, {2, 3, 4}] should give a 2x3x4 tensor */
    s = eval_to_str("Dimensions[RandomChoice[{1, 2}, {2, 3, 4}]]");
    ASSERT_STR_EQ(s, "{2, 3, 4}");
    free(s);
}

void test_randomchoice_coverage(void) {
    /* All elements should appear with enough trials */
    int seen_a = 0, seen_b = 0, seen_c = 0;
    for (int i = 0; i < 200; i++) {
        char* s = eval_to_str("RandomChoice[{a, b, c}]");
        if (strcmp(s, "a") == 0) seen_a = 1;
        if (strcmp(s, "b") == 0) seen_b = 1;
        if (strcmp(s, "c") == 0) seen_c = 1;
        free(s);
    }
    ASSERT_MSG(seen_a && seen_b && seen_c,
               "RandomChoice[{a,b,c}] didn't produce all elements");
}

void test_randomchoice_weighted(void) {
    /* RandomChoice[{1, 0, 0} -> {a, b, c}] should always return a */
    for (int i = 0; i < 50; i++) {
        char* s = eval_to_str("RandomChoice[{1, 0, 0} -> {a, b, c}]");
        ASSERT_STR_EQ(s, "a");
        free(s);
    }

    /* RandomChoice[{0, 0, 1} -> {a, b, c}] should always return c */
    for (int i = 0; i < 50; i++) {
        char* s = eval_to_str("RandomChoice[{0, 0, 1} -> {a, b, c}]");
        ASSERT_STR_EQ(s, "c");
        free(s);
    }
}

void test_randomchoice_weighted_list(void) {
    /* RandomChoice[{1, 0} -> {x, y}, 5] should give a list of 5 x's */
    char* s = eval_to_str("RandomChoice[{1, 0} -> {x, y}, 5]");
    ASSERT_STR_EQ(s, "{x, x, x, x, x}");
    free(s);
}

void test_randomchoice_weighted_multidim(void) {
    /* RandomChoice[{1, 0} -> {a, b}, {2, 3}] should give a 2x3 matrix of a's */
    char* s = eval_to_str("Dimensions[RandomChoice[{1, 0} -> {a, b}, {2, 3}]]");
    ASSERT_STR_EQ(s, "{2, 3}");
    free(s);

    s = eval_to_str("RandomChoice[{1, 0} -> {a, b}, {2, 3}]");
    ASSERT_STR_EQ(s, "{{a, a, a}, {a, a, a}}");
    free(s);
}

void test_randomchoice_seedrandom(void) {
    /* SeedRandom should make RandomChoice deterministic */
    assert_eval_eq("SeedRandom[42]", "Null", 0);
    char* seq1[10];
    for (int i = 0; i < 10; i++) {
        seq1[i] = eval_to_str("RandomChoice[{a, b, c, d, e}]");
    }
    assert_eval_eq("SeedRandom[42]", "Null", 0);
    for (int i = 0; i < 10; i++) {
        char* s = eval_to_str("RandomChoice[{a, b, c, d, e}]");
        ASSERT_MSG(strcmp(s, seq1[i]) == 0,
                   "RandomChoice SeedRandom failed at %d: %s != %s", i, s, seq1[i]);
        free(s);
    }
    for (int i = 0; i < 10; i++) free(seq1[i]);
}

void test_randomchoice_symbolic(void) {
    assert_eval_eq("RandomChoice[x]", "RandomChoice[x]", 0);
    assert_eval_eq("RandomChoice[{a, b}, x]", "RandomChoice[{a, b}, x]", 0);
}

void test_randomchoice_attributes(void) {
    char* s = eval_to_str("MemberQ[Attributes[RandomChoice], Protected]");
    ASSERT_STR_EQ(s, "True");
    free(s);
}

void test_randomchoice_info(void) {
    Expr* e = parse_expression("Information[RandomChoice]");
    Expr* r = evaluate(e);
    ASSERT(r->type == EXPR_STRING);
    expr_free(e);
    expr_free(r);
}

void test_randomchoice_integers(void) {
    /* RandomChoice with integer elements */
    for (int i = 0; i < 100; i++) {
        int64_t v = eval_to_int("RandomChoice[{1, 2, 3}]");
        ASSERT_MSG(v >= 1 && v <= 3,
                   "RandomChoice[{1,2,3}] returned %lld", (long long)v);
    }
}

void test_randomchoice_weighted_coverage(void) {
    /* With weights {1, 1}, both elements should appear */
    int seen_a = 0, seen_b = 0;
    for (int i = 0; i < 200; i++) {
        char* s = eval_to_str("RandomChoice[{1, 1} -> {a, b}]");
        if (strcmp(s, "a") == 0) seen_a = 1;
        if (strcmp(s, "b") == 0) seen_b = 1;
        free(s);
    }
    ASSERT_MSG(seen_a && seen_b,
               "RandomChoice weighted {1,1} didn't produce both elements");
}

void test_randomchoice_single_dim_list(void) {
    char* s = eval_to_str("Dimensions[RandomChoice[{a, b, c}, {5}]]");
    ASSERT_STR_EQ(s, "{5}");
    free(s);
}

/* ---- RandomSample[] ---- */

void test_randomsample_permutation(void) {
    /* RandomSample[{a, b, c}] gives a permutation: length 3, all elements present */
    char* s = eval_to_str("Length[RandomSample[{a, b, c}]]");
    ASSERT_STR_EQ(s, "3");
    free(s);

    /* Check that result is a permutation: Sort should give {a, b, c} */
    s = eval_to_str("Sort[RandomSample[{a, b, c}]]");
    ASSERT_STR_EQ(s, "{a, b, c}");
    free(s);
}

void test_randomsample_sample_n(void) {
    /* RandomSample[{a, b, c, d, e}, 3] gives 3 elements */
    char* s = eval_to_str("Length[RandomSample[{a, b, c, d, e}, 3]]");
    ASSERT_STR_EQ(s, "3");
    free(s);

    /* RandomSample[{a, b, c, d, e}, 1] gives 1 element */
    s = eval_to_str("Length[RandomSample[{a, b, c, d, e}, 1]]");
    ASSERT_STR_EQ(s, "1");
    free(s);
}

void test_randomsample_no_duplicates(void) {
    /* Sample of 5 from 5 elements: all must be unique (a permutation) */
    for (int trial = 0; trial < 20; trial++) {
        char* s = eval_to_str("Sort[RandomSample[{1, 2, 3, 4, 5}, 5]]");
        ASSERT_STR_EQ(s, "{1, 2, 3, 4, 5}");
        free(s);
    }

    /* Sample of 3 from 5: check uniqueness via DeleteDuplicates */
    for (int trial = 0; trial < 20; trial++) {
        char* s = eval_to_str("Module[{s = RandomSample[{a, b, c, d, e}, 3]}, Length[DeleteDuplicates[s]] == Length[s]]");
        ASSERT_STR_EQ(s, "True");
        free(s);
    }
}

void test_randomsample_single_element(void) {
    /* RandomSample[{42}] always returns {42} */
    for (int i = 0; i < 10; i++) {
        char* s = eval_to_str("RandomSample[{42}]");
        ASSERT_STR_EQ(s, "{42}");
        free(s);
    }

    /* RandomSample[{42}, 1] always returns {42} */
    for (int i = 0; i < 10; i++) {
        char* s = eval_to_str("RandomSample[{42}, 1]");
        ASSERT_STR_EQ(s, "{42}");
        free(s);
    }
}

void test_randomsample_sample_all(void) {
    /* Sampling all elements is a permutation */
    for (int trial = 0; trial < 20; trial++) {
        char* s = eval_to_str("Sort[RandomSample[{x, y, z}, 3]]");
        ASSERT_STR_EQ(s, "{x, y, z}");
        free(s);
    }
}

void test_randomsample_upto(void) {
    /* UpTo[3] from a 5-element list: gives 3 elements */
    char* s = eval_to_str("Length[RandomSample[{a, b, c, d, e}, UpTo[3]]]");
    ASSERT_STR_EQ(s, "3");
    free(s);

    /* UpTo[10] from a 5-element list: gives 5 elements (capped) */
    s = eval_to_str("Length[RandomSample[{a, b, c, d, e}, UpTo[10]]]");
    ASSERT_STR_EQ(s, "5");
    free(s);

    /* UpTo[0] gives empty list */
    s = eval_to_str("RandomSample[{a, b, c}, UpTo[0]]");
    ASSERT_STR_EQ(s, "{}");
    free(s);
}

void test_randomsample_upto_exceeds(void) {
    /* UpTo[100] from a 3-element list: should return all 3, sorted gives original */
    char* s = eval_to_str("Sort[RandomSample[{1, 2, 3}, UpTo[100]]]");
    ASSERT_STR_EQ(s, "{1, 2, 3}");
    free(s);
}

void test_randomsample_weighted(void) {
    /* Weighted sample: {1, 0, 0} -> {a, b, c} always picks a */
    for (int i = 0; i < 20; i++) {
        char* s = eval_to_str("RandomSample[{1, 0, 0} -> {a, b, c}, 1]");
        ASSERT_STR_EQ(s, "{a}");
        free(s);
    }

    /* Weighted sample of 2 from deterministic weights */
    char* s = eval_to_str("Sort[RandomSample[{1, 1, 0} -> {a, b, c}, 2]]");
    ASSERT_STR_EQ(s, "{a, b}");
    free(s);
}

void test_randomsample_weighted_deterministic(void) {
    /* {0, 0, 1} -> {x, y, z}: only z has weight, so sample of 1 is {z} */
    for (int i = 0; i < 10; i++) {
        char* s = eval_to_str("RandomSample[{0, 0, 1} -> {x, y, z}, 1]");
        ASSERT_STR_EQ(s, "{z}");
        free(s);
    }
}

void test_randomsample_weighted_upto(void) {
    /* Weighted sample with UpTo */
    char* s = eval_to_str("Length[RandomSample[{1, 2, 3} -> {a, b, c}, UpTo[2]]]");
    ASSERT_STR_EQ(s, "2");
    free(s);

    /* UpTo[10] caps to list length */
    s = eval_to_str("Length[RandomSample[{1, 2, 3} -> {a, b, c}, UpTo[10]]]");
    ASSERT_STR_EQ(s, "3");
    free(s);
}

void test_randomsample_weighted_permutation(void) {
    /* Weighted permutation: all elements present */
    char* s = eval_to_str("Sort[RandomSample[{1, 2, 3} -> {a, b, c}]]");
    ASSERT_STR_EQ(s, "{a, b, c}");
    free(s);
}

void test_randomsample_seedrandom(void) {
    /* SeedRandom produces reproducible results */
    char* s1 = eval_to_str("SeedRandom[42]; RandomSample[{a, b, c, d, e}, 3]");
    char* s2 = eval_to_str("SeedRandom[42]; RandomSample[{a, b, c, d, e}, 3]");
    ASSERT_MSG(strcmp(s1, s2) == 0,
               "SeedRandom[42] should give reproducible RandomSample, got %s vs %s", s1, s2);
    free(s1);
    free(s2);

    /* Different seeds give different results (probabilistically) */
    s1 = eval_to_str("SeedRandom[1]; RandomSample[Range[20], 10]");
    s2 = eval_to_str("SeedRandom[2]; RandomSample[Range[20], 10]");
    ASSERT_MSG(strcmp(s1, s2) != 0,
               "Different seeds should likely give different results");
    free(s1);
    free(s2);
}

void test_randomsample_coverage(void) {
    /* Over many trials, all elements should appear in position 1 at least once */
    int seen[5] = {0, 0, 0, 0, 0};
    for (int trial = 0; trial < 200; trial++) {
        char* s = eval_to_str("First[RandomSample[{1, 2, 3, 4, 5}]]");
        int64_t v = atoi(s);
        if (v >= 1 && v <= 5) seen[v - 1] = 1;
        free(s);
    }
    for (int i = 0; i < 5; i++) {
        ASSERT_MSG(seen[i], "Element %d never appeared first in RandomSample", i + 1);
    }
}

void test_randomsample_symbolic(void) {
    /* RandomSample[x] should remain unevaluated (x is not a list) */
    char* s = eval_to_str("RandomSample[x]");
    ASSERT_STR_EQ(s, "RandomSample[x]");
    free(s);

    /* RandomSample[{a, b}, x] should remain unevaluated (x is not integer or UpTo) */
    s = eval_to_str("RandomSample[{a, b}, x]");
    ASSERT_STR_EQ(s, "RandomSample[{a, b}, x]");
    free(s);
}

void test_randomsample_attributes(void) {
    char* s = eval_to_str("MemberQ[Attributes[RandomSample], Protected]");
    ASSERT_STR_EQ(s, "True");
    free(s);
}

void test_randomsample_info(void) {
    char* s = eval_to_str("Head[Information[RandomSample]]");
    ASSERT_STR_EQ(s, "String");
    free(s);
}

void test_randomsample_n_exceeds_length(void) {
    /* Requesting more elements than available without UpTo: stays unevaluated */
    char* s = eval_to_str("RandomSample[{a, b}, 5]");
    ASSERT_STR_EQ(s, "RandomSample[{a, b}, 5]");
    free(s);
}

void test_randomsample_empty_list_upto(void) {
    /* UpTo[5] on an empty-ish list: should return empty or stay unevaluated */
    /* Note: is_nonempty_list requires at least 1 element */
    char* s = eval_to_str("RandomSample[{a}, UpTo[5]]");
    ASSERT_STR_EQ(s, "{a}");
    free(s);
}

void test_randomsample_sample_zero(void) {
    /* Sample of 0 elements: gives empty list */
    char* s = eval_to_str("RandomSample[{a, b, c}, 0]");
    ASSERT_STR_EQ(s, "{}");
    free(s);
}

int main(void) {
    symtab_init();
    core_init();

    TEST(test_randominteger_no_args);
    TEST(test_randominteger_single_arg);
    TEST(test_randominteger_range);
    TEST(test_randominteger_flat_list);
    TEST(test_randominteger_multidim);
    TEST(test_seedrandom);
    TEST(test_seedrandom_no_args);
    TEST(test_seedrandom_returns_null);
    TEST(test_randominteger_symbolic);
    TEST(test_randominteger_bignum);
    TEST(test_randominteger_negative_range);
    TEST(test_randominteger_coverage);
    TEST(test_randominteger_attributes);
    TEST(test_randominteger_info);
    TEST(test_randominteger_single_dim_list);

    /* RandomReal tests */
    TEST(test_randomreal_no_args);
    TEST(test_randomreal_single_arg);
    TEST(test_randomreal_range);
    TEST(test_randomreal_rational_range);
    TEST(test_randomreal_flat_list);
    TEST(test_randomreal_multidim);
    TEST(test_randomreal_seedrandom);
    TEST(test_randomreal_symbolic);
    TEST(test_randomreal_attributes);
    TEST(test_randomreal_info);
    TEST(test_randomreal_negative_range);
    TEST(test_randomreal_coverage);
    TEST(test_randomreal_single_dim_list);
    TEST(test_randomreal_type);

    /* RandomComplex tests */
    TEST(test_randomcomplex_no_args);
    TEST(test_randomcomplex_single_complex_arg);
    TEST(test_randomcomplex_range);
    TEST(test_randomcomplex_negative_range);
    TEST(test_randomcomplex_flat_list);
    TEST(test_randomcomplex_multidim);
    TEST(test_randomcomplex_seedrandom);
    TEST(test_randomcomplex_symbolic);
    TEST(test_randomcomplex_attributes);
    TEST(test_randomcomplex_info);
    TEST(test_randomcomplex_head);
    TEST(test_randomcomplex_coverage);
    TEST(test_randomcomplex_real_range);
    TEST(test_randomcomplex_single_dim_list);

#ifdef USE_MPFR
    /* WorkingPrecision option (MPFR build only) */
    TEST(test_randomreal_workingprecision_basic);
    TEST(test_randomreal_workingprecision_range);
    TEST(test_randomreal_workingprecision_list);
    TEST(test_randomreal_workingprecision_array);
    TEST(test_randomreal_machineprecision_option);
    TEST(test_randomcomplex_workingprecision_basic);
    TEST(test_randomcomplex_workingprecision_list);
#endif

    /* RandomChoice tests */
    TEST(test_randomchoice_basic);
    TEST(test_randomchoice_single_element);
    TEST(test_randomchoice_flat_list);
    TEST(test_randomchoice_multidim);
    TEST(test_randomchoice_coverage);
    TEST(test_randomchoice_weighted);
    TEST(test_randomchoice_weighted_list);
    TEST(test_randomchoice_weighted_multidim);
    TEST(test_randomchoice_seedrandom);
    TEST(test_randomchoice_symbolic);
    TEST(test_randomchoice_attributes);
    TEST(test_randomchoice_info);
    TEST(test_randomchoice_integers);
    TEST(test_randomchoice_weighted_coverage);
    TEST(test_randomchoice_single_dim_list);

    /* RandomSample tests */
    TEST(test_randomsample_permutation);
    TEST(test_randomsample_sample_n);
    TEST(test_randomsample_no_duplicates);
    TEST(test_randomsample_single_element);
    TEST(test_randomsample_sample_all);
    TEST(test_randomsample_upto);
    TEST(test_randomsample_upto_exceeds);
    TEST(test_randomsample_weighted);
    TEST(test_randomsample_weighted_deterministic);
    TEST(test_randomsample_weighted_upto);
    TEST(test_randomsample_weighted_permutation);
    TEST(test_randomsample_seedrandom);
    TEST(test_randomsample_coverage);
    TEST(test_randomsample_symbolic);
    TEST(test_randomsample_attributes);
    TEST(test_randomsample_info);
    TEST(test_randomsample_n_exceeds_length);
    TEST(test_randomsample_empty_list_upto);
    TEST(test_randomsample_sample_zero);

    printf("All random tests passed!\n");
    return 0;
}
