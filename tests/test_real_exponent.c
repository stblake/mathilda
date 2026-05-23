/* Unit tests for RealExponent.
 *
 *   RealExponent[x]     -> Log[10, |x|]
 *   RealExponent[x, b]  -> Log[b,  |x|]
 *
 * Coverage:
 *   - Integer / BigInt (machine-int, 2^100, 10^50).
 *   - Rational (exact, |x| < 1, negative).
 *   - Machine Real (positive / negative / power-of-base / |x| < 1).
 *   - MPFR (N[Pi, 32], N[E, 32], MPFR base, mixed precisions).
 *   - Base argument (2, 3, 16, machine Real, MPFR, symbolic E).
 *   - Symbolic numeric (Pi, E, Pi^Pi) -> numericalized at machine precision.
 *   - Listable threading.
 *   - Identity:  RealExponent[b^p] == p * RealExponent[b].
 *   - Identity:  RealExponent[x y] == RealExponent[x] + RealExponent[y].
 *   - Zero (Integer 0, Rational 0/d, machine 0., MPFR 0``p)  ->  -Infinity.
 *   - Error paths:  0 args, 3+ args, Complex x, Complex b,
 *                   base 1, 0, -2, 0.5.
 *   - Symbolic non-numeric x left unevaluated.
 *   - Attributes / docstring / interned symbol pointer.
 *   - Memory-safety stress loop.
 */

#include "expr.h"
#include "eval.h"
#include "core.h"
#include "symtab.h"
#include "test_utils.h"
#include "parse.h"
#include "print.h"
#include "attr.h"
#include "sym_names.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

/* Capture stderr while `input` is parsed + evaluated.  Returns the
 * collected stderr text as a heap string (caller frees) and writes the
 * printed result into *out_result_str (also heap-allocated). */
static char* eval_capturing_stderr(const char* input, char** out_result_str) {
    const char* path = "/tmp/mathilda_real_exp_stderr.log";
    fflush(stderr);
    if (!freopen(path, "w+", stderr)) {
        if (out_result_str) *out_result_str = NULL;
        return NULL;
    }
    Expr* p = parse_expression(input);
    Expr* e = evaluate(p);
    if (out_result_str) *out_result_str = expr_to_string(e);
    expr_free(p);
    expr_free(e);

    fflush(stderr);
    freopen("/dev/tty", "w", stderr);

    FILE* f = fopen(path, "r");
    if (!f) return NULL;
    fseek(f, 0, SEEK_END);
    long n = ftell(f);
    if (n < 0) { fclose(f); return NULL; }
    fseek(f, 0, SEEK_SET);
    char* buf = malloc((size_t)n + 1);
    if (!buf) { fclose(f); return NULL; }
    size_t got = fread(buf, 1, (size_t)n, f);
    buf[got] = '\0';
    fclose(f);
    remove(path);
    return buf;
}

/* Evaluate `input` and ASSERT the resulting EXPR_REAL is within `tol` of
 * `expected`.  Frees both the parsed AST and the evaluated result. */
static void assert_eval_real_close(const char* input, double expected, double tol) {
    Expr* p = parse_expression(input);
    Expr* e = evaluate(p);
    expr_free(p);
    ASSERT_MSG(e->type == EXPR_REAL,
               "expected EXPR_REAL, got type %d for input: %s", e->type, input);
    ASSERT_MSG(fabs(e->data.real - expected) < tol,
               "expected %.15g (+- %g), got %.15g for input: %s",
               expected, tol, e->data.real, input);
    expr_free(e);
}

/* --- Integers / BigInts --------------------------------------------- */

static void test_int_small(void) {
    assert_eval_real_close("RealExponent[1]",  0.0,      1e-12);
    assert_eval_real_close("RealExponent[10]", 1.0,      1e-12);
    assert_eval_real_close("RealExponent[2]",  log10(2.0), 1e-12);
    assert_eval_real_close("RealExponent[100]",2.0,      1e-12);
}

static void test_int_negative(void) {
    /* Sign is discarded -- log_b(|x|). */
    assert_eval_real_close("RealExponent[-1000]", 3.0, 1e-12);
    assert_eval_real_close("RealExponent[-2]",    log10(2.0), 1e-12);
}

static void test_int_bigint_base10(void) {
    /* 10^50 in base 10: log_10(10^50) = 50. */
    assert_eval_real_close("RealExponent[10^50]", 50.0, 1e-9);
}

static void test_int_bigint_base2(void) {
    /* 2^100 in base 2: log_2(2^100) = 100. */
    assert_eval_real_close("RealExponent[2^100, 2]", 100.0, 1e-9);
}

/* --- Rationals ------------------------------------------------------ */

static void test_rat_small(void) {
    /* RealExponent[1/2] = log_10(1/2) = -log_10(2) ~= -0.30103. */
    assert_eval_real_close("RealExponent[1/2]", -log10(2.0), 1e-12);
}

static void test_rat_documented(void) {
    /* From Mathematica documentation. */
    assert_eval_real_close("RealExponent[987654321/123456789]",
                           log10(987654321.0 / 123456789.0), 1e-9);
}

static void test_rat_negative(void) {
    assert_eval_real_close("RealExponent[-3/2]", log10(1.5), 1e-12);
}

/* --- Machine reals -------------------------------------------------- */

static void test_real_documented(void) {
    /* From Mathematica documentation: y = RealExponent[123.456] -> 2.09151. */
    assert_eval_real_close("RealExponent[123.456]", log10(123.456), 1e-12);
}

static void test_real_base2_documented(void) {
    /* y = RealExponent[123.456, 2] -> 6.94785. */
    assert_eval_real_close("RealExponent[123.456, 2]", log10(123.456)/log10(2.0),
                           1e-12);
}

static void test_real_negative(void) {
    assert_eval_real_close("RealExponent[-123.456]", log10(123.456), 1e-12);
}

static void test_real_below_one(void) {
    assert_eval_real_close("RealExponent[0.001]", -3.0, 1e-12);
}

static void test_real_base_real(void) {
    /* RealExponent[10, 2.5] = log(10)/log(2.5). */
    assert_eval_real_close("RealExponent[10, 2.5]", log(10.0)/log(2.5), 1e-12);
}

/* --- Identities ----------------------------------------------------- */

static void test_identity_product(void) {
    /* RealExponent[x y] == RealExponent[x] + RealExponent[y]. */
    assert_eval_eq("RealExponent[3.4 * 5.6] == "
                   "RealExponent[3.4] + RealExponent[5.6]", "True", 0);
}

static void test_identity_power(void) {
    /* p * RealExponent[b] == RealExponent[b^p] for positive b, integer p. */
    assert_eval_eq("RealExponent[2.1^1000000000] == "
                   "1000000000 RealExponent[2.1]", "True", 0);
}

/* --- MPFR ----------------------------------------------------------- */

#ifdef USE_MPFR
static void test_mpfr_pi(void) {
    /* RealExponent[N[Pi, 32]] -> MPFR result at ~32 digits.
     * The expected value is log_10(Pi) = 0.4971498726941338... */
    Expr* p = parse_expression("RealExponent[N[Pi, 32]]");
    Expr* e = evaluate(p);
    expr_free(p);
    /* Must be an MPFR value, not a machine Real. */
    ASSERT(e->type == EXPR_MPFR);
    /* First 12 significant digits should match. */
    char* s = expr_to_string(e);
    ASSERT(s != NULL);
    ASSERT_MSG(strncmp(s, "0.4971498726", 12) == 0,
               "expected mantissa starting 0.4971498726, got: %s", s);
    free(s);
    expr_free(e);
}

static void test_mpfr_e_base2(void) {
    /* RealExponent[N[E, 32], 2] = log_2(E) ~= 1.4426950408889634... */
    Expr* p = parse_expression("RealExponent[N[E, 32], 2]");
    Expr* e = evaluate(p);
    expr_free(p);
    ASSERT(e->type == EXPR_MPFR);
    char* s = expr_to_string(e);
    ASSERT_MSG(strncmp(s, "1.4426950408", 12) == 0,
               "expected 1.4426950408..., got: %s", s);
    free(s);
    expr_free(e);
}

static void test_mpfr_machine_base(void) {
    /* MPFR x with machine-int base 7: result is MPFR at x's precision. */
    Expr* p = parse_expression("RealExponent[N[Pi, 40], 7]");
    Expr* e = evaluate(p);
    expr_free(p);
    ASSERT(e->type == EXPR_MPFR);
    /* log_7(Pi) ~= 0.5882747... */
    char* s = expr_to_string(e);
    ASSERT_MSG(strncmp(s, "0.5882747", 9) == 0,
               "expected 0.5882747..., got: %s", s);
    free(s);
    expr_free(e);
}

static void test_mpfr_identity(void) {
    /* MPFR analogue of the product identity:
     *   RealExponent[x*y] - (RealExponent[x] + RealExponent[y]) is zero
     *   modulo MPFR round-off at the working precision.
     *
     * We side-step Mathilda's not-yet-simplifying `Abs[<mpfr>]` by
     * computing the difference and exporting its machine-double value
     * via mpfr_get_d in the EXPR_MPFR leaf. */
    Expr* p = parse_expression(
        "Module[{x, y}, x = N[Pi, 40]; y = N[E, 40]; "
        "RealExponent[x*y] - (RealExponent[x] + RealExponent[y])]");
    Expr* e = evaluate(p);
    expr_free(p);
    double err;
    if (e->type == EXPR_REAL) err = e->data.real;
#ifdef USE_MPFR
    else if (e->type == EXPR_MPFR) err = mpfr_get_d(e->data.mpfr, MPFR_RNDN);
#endif
    else { ASSERT_MSG(false, "expected Real/MPFR, got type %d", e->type); err = 0; }
    ASSERT_MSG(fabs(err) < 1e-30,
               "product identity error %.6g exceeds 1e-30 tolerance", err);
    expr_free(e);
}
#endif

/* --- Symbolic numeric (numericalized) ------------------------------ */

static void test_symbolic_pi(void) {
    /* Bare Pi: numericalize at machine precision, return machine Real. */
    assert_eval_real_close("RealExponent[Pi]", log10(3.141592653589793), 1e-12);
}

static void test_symbolic_pi_power_pi(void) {
    /* Pi^Pi at machine precision: log_10(Pi^Pi) ~= 1.56184. */
    assert_eval_real_close("RealExponent[Pi^Pi]",
                           log10(pow(3.141592653589793, 3.141592653589793)),
                           1e-9);
}

static void test_symbolic_pi_base_e(void) {
    /* RealExponent[Pi, E] = ln(Pi) ~= 1.14473. */
    assert_eval_real_close("RealExponent[Pi, E]", log(3.141592653589793),
                           1e-12);
}

/* --- Listable ------------------------------------------------------- */

static void test_listable(void) {
    /* RealExponent[Range[5]] = {0., 0.301..., 0.477..., 0.602..., 0.699...}. */
    Expr* p = parse_expression("RealExponent[{1, 2, 3, 4, 5}]");
    Expr* e = evaluate(p);
    expr_free(p);
    ASSERT(e->type == EXPR_FUNCTION && e->data.function.arg_count == 5);
    double expected[] = { 0.0, log10(2.0), log10(3.0), log10(4.0), log10(5.0) };
    for (size_t i = 0; i < 5; i++) {
        Expr* v = e->data.function.args[i];
        ASSERT(v->type == EXPR_REAL);
        ASSERT_MSG(fabs(v->data.real - expected[i]) < 1e-12,
                   "RealExponent[{...}][%zu]: expected %.15g, got %.15g",
                   i, expected[i], v->data.real);
    }
    expr_free(e);
}

static void test_listable_base(void) {
    /* RealExponent[Pi, {2, 3, 5, 7, 10}] -- threads over the base. */
    Expr* p = parse_expression("RealExponent[Pi, {2, 3, 5, 7, 10}]");
    Expr* e = evaluate(p);
    expr_free(p);
    ASSERT(e->type == EXPR_FUNCTION && e->data.function.arg_count == 5);
    double pi = 3.141592653589793;
    double expected[] = {
        log(pi)/log(2.0),
        log(pi)/log(3.0),
        log(pi)/log(5.0),
        log(pi)/log(7.0),
        log10(pi),
    };
    for (size_t i = 0; i < 5; i++) {
        Expr* v = e->data.function.args[i];
        ASSERT(v->type == EXPR_REAL);
        ASSERT_MSG(fabs(v->data.real - expected[i]) < 1e-12,
                   "RealExponent[Pi, {...}][%zu]: expected %.15g, got %.15g",
                   i, expected[i], v->data.real);
    }
    expr_free(e);
}

/* --- Zero cases ----------------------------------------------------- */

static void test_zero_integer(void) {
    /* RealExponent[0] = -Infinity (per documentation). */
    assert_eval_eq("RealExponent[0]", "-Infinity", 0);
}

static void test_zero_machine_real(void) {
    /* Mathilda convention: Accuracy[0.] = Infinity, so RealExponent[0.] =
     * -Infinity.  This intentionally differs from Mathematica which
     * returns the finite ~-307.65. */
    assert_eval_eq("RealExponent[0.]", "-Infinity", 0);
}

static void test_zero_rational(void) {
    /* 0/5 normalizes to 0, but let's also check the explicit zero through
     * exact arithmetic. */
    assert_eval_eq("RealExponent[0/5]", "-Infinity", 0);
}

#ifdef USE_MPFR
static void test_zero_mpfr(void) {
    /* MPFR zero via N[0, 50] -- Mathilda Accuracy convention: -Infinity. */
    assert_eval_eq("RealExponent[N[0, 50]]", "-Infinity", 0);
}
#endif

/* --- Error paths ---------------------------------------------------- */

static void test_zero_args(void) {
    char* result = NULL;
    char* err = eval_capturing_stderr("RealExponent[]", &result);
    ASSERT(result != NULL);
    ASSERT_STR_EQ(result, "RealExponent[]");
    ASSERT(err != NULL);
    ASSERT_MSG(strstr(err, "RealExponent::argt") != NULL,
               "expected argt diagnostic, got: %s", err);
    ASSERT(strstr(err, "called with 0 arguments") != NULL);
    free(result);
    free(err);
}

static void test_four_args(void) {
    char* result = NULL;
    char* err = eval_capturing_stderr("RealExponent[x, 2, 3, 4]", &result);
    ASSERT(result != NULL);
    ASSERT(strstr(result, "RealExponent[x, 2, 3, 4]") != NULL);
    ASSERT(err != NULL);
    ASSERT(strstr(err, "RealExponent::argt") != NULL);
    ASSERT(strstr(err, "called with 4 arguments") != NULL);
    free(result);
    free(err);
}

static void test_base_one(void) {
    char* result = NULL;
    char* err = eval_capturing_stderr("RealExponent[5, 1]", &result);
    ASSERT(result != NULL);
    ASSERT_MSG(err && strstr(err, "RealExponent::ibase") != NULL,
               "expected ibase diagnostic, got: %s", err ? err : "(null)");
    free(result);
    free(err);
}

static void test_base_zero(void) {
    char* result = NULL;
    char* err = eval_capturing_stderr("RealExponent[5, 0]", &result);
    ASSERT(err && strstr(err, "RealExponent::ibase") != NULL);
    free(result);
    free(err);
}

static void test_base_negative(void) {
    char* result = NULL;
    char* err = eval_capturing_stderr("RealExponent[5, -2]", &result);
    ASSERT(err && strstr(err, "RealExponent::ibase") != NULL);
    free(result);
    free(err);
}

static void test_base_half(void) {
    /* Base = 1/2 < 1 should be rejected. */
    char* result = NULL;
    char* err = eval_capturing_stderr("RealExponent[5, 1/2]", &result);
    ASSERT(err && strstr(err, "RealExponent::ibase") != NULL);
    free(result);
    free(err);
}

static void test_base_real_le_one(void) {
    /* Base = 0.5 should be rejected. */
    char* result = NULL;
    char* err = eval_capturing_stderr("RealExponent[5, 0.5]", &result);
    ASSERT(err && strstr(err, "RealExponent::ibase") != NULL);
    free(result);
    free(err);
}

static void test_complex_x(void) {
    char* result = NULL;
    char* err = eval_capturing_stderr("RealExponent[2.34 + 9.1 I]", &result);
    /* Call is left unevaluated. */
    ASSERT(result && strstr(result, "RealExponent") != NULL);
    ASSERT_MSG(err && strstr(err, "RealExponent::realx") != NULL,
               "expected realx diagnostic, got: %s", err ? err : "(null)");
    free(result);
    free(err);
}

static void test_complex_base(void) {
    /* Complex base -> ::ibase (not ::realx). */
    char* result = NULL;
    char* err = eval_capturing_stderr("RealExponent[5, 1 + 2 I]", &result);
    ASSERT(err && strstr(err, "RealExponent::ibase") != NULL);
    free(result);
    free(err);
}

static void test_symbolic_x(void) {
    /* Plain unbound symbol -- left unevaluated, no error. */
    Expr* p = parse_expression("RealExponent[xyzfoo]");
    Expr* e = evaluate(p);
    char* s = expr_to_string(e);
    ASSERT(strstr(s, "RealExponent") != NULL);
    ASSERT(strstr(s, "xyzfoo") != NULL);
    free(s);
    expr_free(p);
    expr_free(e);
}

/* --- Attributes / docstring / interned symbol ----------------------- */

static void test_attributes(void) {
    SymbolDef* def = symtab_get_def("RealExponent");
    ASSERT(def != NULL);
    uint32_t a = get_attributes("RealExponent");
    ASSERT((a & ATTR_PROTECTED) != 0);
    ASSERT((a & ATTR_LISTABLE) != 0);
}

static void test_docstring(void) {
    SymbolDef* def = symtab_get_def("RealExponent");
    ASSERT(def != NULL);
    ASSERT(def->docstring != NULL);
    ASSERT(strstr(def->docstring, "RealExponent[x]") != NULL);
}

static void test_sym_pointer_interned(void) {
    ASSERT(SYM_RealExponent != NULL);
    ASSERT(strcmp(SYM_RealExponent, "RealExponent") == 0);
}

/* --- Memory-safety stress loop ------------------------------------- */

static void test_stress(void) {
    /* Mix every code path -- under valgrind this catches mpz / mpfr / Expr
     * leaks in any branch. */
    for (int k = 0; k < 30; k++) {
        assert_eval_real_close("RealExponent[123.456]",
                               log10(123.456), 1e-12);
        assert_eval_real_close("RealExponent[10^50]", 50.0, 1e-9);
        assert_eval_real_close("RealExponent[1/2]", -log10(2.0), 1e-12);
        assert_eval_real_close("RealExponent[Pi]",
                               log10(3.141592653589793), 1e-12);
        assert_eval_real_close("RealExponent[2^100, 2]", 100.0, 1e-9);
        assert_eval_eq("RealExponent[0]",  "-Infinity", 0);
        assert_eval_eq("RealExponent[0.]", "-Infinity", 0);
#ifdef USE_MPFR
        Expr* p = parse_expression("RealExponent[N[Pi, 30]]");
        Expr* e = evaluate(p);
        expr_free(p);
        expr_free(e);
#endif
    }
}

int main(void) {
    symtab_init();
    core_init();

    TEST(test_int_small);
    TEST(test_int_negative);
    TEST(test_int_bigint_base10);
    TEST(test_int_bigint_base2);

    TEST(test_rat_small);
    TEST(test_rat_documented);
    TEST(test_rat_negative);

    TEST(test_real_documented);
    TEST(test_real_base2_documented);
    TEST(test_real_negative);
    TEST(test_real_below_one);
    TEST(test_real_base_real);

    TEST(test_identity_product);
    TEST(test_identity_power);

#ifdef USE_MPFR
    TEST(test_mpfr_pi);
    TEST(test_mpfr_e_base2);
    TEST(test_mpfr_machine_base);
    TEST(test_mpfr_identity);
#endif

    TEST(test_symbolic_pi);
    TEST(test_symbolic_pi_power_pi);
    TEST(test_symbolic_pi_base_e);

    TEST(test_listable);
    TEST(test_listable_base);

    TEST(test_zero_integer);
    TEST(test_zero_machine_real);
    TEST(test_zero_rational);
#ifdef USE_MPFR
    TEST(test_zero_mpfr);
#endif

    TEST(test_zero_args);
    TEST(test_four_args);
    TEST(test_base_one);
    TEST(test_base_zero);
    TEST(test_base_negative);
    TEST(test_base_half);
    TEST(test_base_real_le_one);
    TEST(test_complex_x);
    TEST(test_complex_base);
    TEST(test_symbolic_x);

    TEST(test_attributes);
    TEST(test_docstring);
    TEST(test_sym_pointer_interned);

    TEST(test_stress);

    printf("All RealExponent tests passed!\n");
    return 0;
}
