/* Unit tests for MantissaExponent.
 *
 *   MantissaExponent[x]      -> {m, e} such that x = m * 10^e, 1/10 <= |m| < 1.
 *   MantissaExponent[x, b]   -> base-b mantissa / exponent.
 *
 * Coverage:
 *   - Integer (positive / negative / zero / one / power-of-base).
 *   - BigInt (2^100, 10^50).
 *   - Rational (terminating / recurring / |x| < 1 / |x| >= 1).
 *   - Machine real (positive / negative / very large / very small / zero).
 *   - MPFR real (N[Pi, 30], N[Pi, 50] in base 2).
 *   - Base argument: 2, 10, 16, 256.
 *   - Listable threading on x.
 *   - Mathematica identity: m * b^e == x (round-trip).
 *   - Error paths: 0 args, 3 args, base 0 / 1 / -2.
 *   - Symbolic x -- left unevaluated, no diagnostic.
 *   - Complex x -> ::realx diagnostic.
 *   - Attributes / docstring / interned symbol introspection.
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
 * printed result into *out_result_str (also heap-allocated).  */
static char* eval_capturing_stderr(const char* input, char** out_result_str) {
    const char* path = "/tmp/mathilda_mantissa_exp_stderr.log";
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

/* --- Integers -------------------------------------------------------- */

static void test_int_documented(void) {
    /* From the Mathematica documentation. */
    assert_eval_eq("MantissaExponent[123451]",
                   "{123451/1000000, 6}", 0);
}

static void test_int_one(void) {
    /* 1 = (1/10) * 10. */
    assert_eval_eq("MantissaExponent[1]", "{1/10, 1}", 0);
}

static void test_int_ten(void) {
    /* 10 = (1/10) * 10^2. */
    assert_eval_eq("MantissaExponent[10]", "{1/10, 2}", 0);
}

static void test_int_eight_base10(void) {
    /* 8 = (4/5) * 10. */
    assert_eval_eq("MantissaExponent[8]", "{4/5, 1}", 0);
}

static void test_int_negative(void) {
    /* Sign carries through to the mantissa. */
    assert_eval_eq("MantissaExponent[-123451]",
                   "{-123451/1000000, 6}", 0);
}

static void test_int_zero(void) {
    assert_eval_eq("MantissaExponent[0]", "{0, 0}", 0);
}

static void test_int_base2(void) {
    /* 1027 = (1027/2048) * 2^11; documented example. */
    assert_eval_eq("MantissaExponent[1027, 2]", "{1027/2048, 11}", 0);
}

static void test_int_base16(void) {
    /* 58127 = 0xE30F -> log_16(58127) ~= 3.97; e = 4, scale = 65536. */
    assert_eval_eq("MantissaExponent[58127, 16]",
                   "{58127/65536, 4}", 0);
}

static void test_int_base256(void) {
    /* 1000 = (125/32) * 256^? -- log_256(1000) ~= 1.246, e = 2, m = 1000/65536. */
    assert_eval_eq("MantissaExponent[1000, 256]",
                   "{125/8192, 2}", 0);
}

/* 2^100 in base 2: m = 1/2, e = 101 (since 2^100 = (1/2) * 2^101). */
static void test_int_bignum_base2(void) {
    assert_eval_eq("MantissaExponent[2^100, 2]", "{1/2, 101}", 0);
}

/* 10^50 in base 10: m = 1/10, e = 51. */
static void test_int_bignum_base10(void) {
    assert_eval_eq("MantissaExponent[10^50]", "{1/10, 51}", 0);
}

/* --- Rationals ------------------------------------------------------- */

static void test_rat_one_third(void) {
    /* 1/3 < 1, so e = 0 and mantissa = 1/3. */
    assert_eval_eq("MantissaExponent[1/3]", "{1/3, 0}", 0);
}

static void test_rat_three_halves(void) {
    /* 3/2 = (3/20) * 10. */
    assert_eval_eq("MantissaExponent[3/2]", "{3/20, 1}", 0);
}

static void test_rat_five_eighths_base2(void) {
    /* 5/8 in base 2: 5/8 = 0.101_2, so m = 5/8, e = 0. */
    assert_eval_eq("MantissaExponent[5/8, 2]", "{5/8, 0}", 0);
}

static void test_rat_negative(void) {
    /* -3/2 -> mantissa negative. */
    assert_eval_eq("MantissaExponent[-3/2]", "{-3/20, 1}", 0);
}

static void test_rat_tiny(void) {
    /* 1/1000 = (1/10) * 10^-2  -> {1/10, -2}. */
    assert_eval_eq("MantissaExponent[1/1000]", "{1/10, -2}", 0);
}

/* --- Machine reals --------------------------------------------------- */

static void test_real_documented(void) {
    /* 3.4 * 10^30 -> {0.34, 31}. */
    assert_eval_eq("MantissaExponent[3.4 10^30]", "{0.34, 31}", 0);
}

static void test_real_456(void) {
    /* 456.1414 -> {0.456141, 3}. The Mathilda printer shows ~6 sig figs. */
    Expr* p = parse_expression("MantissaExponent[456.1414]");
    Expr* e = evaluate(p);
    expr_free(p);
    ASSERT(e->type == EXPR_FUNCTION && e->data.function.arg_count == 2);
    Expr* m_e = e->data.function.args[0];
    Expr* ev  = e->data.function.args[1];
    ASSERT(m_e->type == EXPR_REAL);
    ASSERT(fabs(m_e->data.real - 0.4561414) < 1e-12);
    ASSERT(ev->type == EXPR_INTEGER && ev->data.integer == 3);
    expr_free(e);
}

static void test_real_negative(void) {
    /* Sign carries to mantissa. */
    Expr* p = parse_expression("MantissaExponent[-456.1414]");
    Expr* e = evaluate(p);
    expr_free(p);
    Expr* m_e = e->data.function.args[0];
    Expr* ev  = e->data.function.args[1];
    ASSERT(m_e->type == EXPR_REAL);
    ASSERT(fabs(m_e->data.real + 0.4561414) < 1e-12);
    ASSERT(ev->data.integer == 3);
    expr_free(e);
}

static void test_real_zero(void) {
    assert_eval_eq("MantissaExponent[0.]", "{0.0, 0}", 0);
}

static void test_real_below_one(void) {
    /* 0.5 = (0.5) * 10^0  -> {0.5, 0}. */
    Expr* p = parse_expression("MantissaExponent[0.5]");
    Expr* e = evaluate(p);
    expr_free(p);
    Expr* m_e = e->data.function.args[0];
    Expr* ev  = e->data.function.args[1];
    ASSERT(m_e->type == EXPR_REAL && fabs(m_e->data.real - 0.5) < 1e-15);
    ASSERT(ev->data.integer == 0);
    expr_free(e);
}

static void test_real_one_exactly(void) {
    /* 1.0 = (0.1) * 10^1  -> {0.1, 1}. */
    Expr* p = parse_expression("MantissaExponent[1.0]");
    Expr* e = evaluate(p);
    expr_free(p);
    Expr* m_e = e->data.function.args[0];
    Expr* ev  = e->data.function.args[1];
    ASSERT(m_e->type == EXPR_REAL && fabs(m_e->data.real - 0.1) < 1e-15);
    ASSERT(ev->data.integer == 1);
    expr_free(e);
}

static void test_real_base2(void) {
    /* 1.5 in base 2: 1.5 = (0.75) * 2^1 -> {0.75, 1}. */
    Expr* p = parse_expression("MantissaExponent[1.5, 2]");
    Expr* e = evaluate(p);
    expr_free(p);
    Expr* m_e = e->data.function.args[0];
    Expr* ev  = e->data.function.args[1];
    ASSERT(m_e->type == EXPR_REAL && fabs(m_e->data.real - 0.75) < 1e-15);
    ASSERT(ev->data.integer == 1);
    expr_free(e);
}

/* --- MPFR ------------------------------------------------------------ */

#ifdef USE_MPFR
static void test_mpfr_pi_30(void) {
    /* N[Pi, 30] -> MantissaExponent should give mantissa starting 0.31415... */
    Expr* p = parse_expression("MantissaExponent[N[Pi, 30]]");
    Expr* e = evaluate(p);
    expr_free(p);
    ASSERT(e->type == EXPR_FUNCTION && e->data.function.arg_count == 2);
    Expr* ev = e->data.function.args[1];
    ASSERT(ev->type == EXPR_INTEGER && ev->data.integer == 1);
    /* Verify the mantissa starts "0.3141592" via string rendering. */
    char* s = expr_to_string(e->data.function.args[0]);
    ASSERT(s != NULL);
    ASSERT_MSG(strncmp(s, "0.3141592", 9) == 0,
               "expected mantissa starting 0.3141592, got: %s", s);
    free(s);
    expr_free(e);
}

static void test_mpfr_pi_50_base2(void) {
    /* N[Pi, 50] in base 2: mantissa = Pi/4 ~= 0.7853981633974483... */
    Expr* p = parse_expression("MantissaExponent[N[Pi, 50], 2]");
    Expr* e = evaluate(p);
    expr_free(p);
    Expr* ev = e->data.function.args[1];
    ASSERT(ev->type == EXPR_INTEGER && ev->data.integer == 2);
    char* s = expr_to_string(e->data.function.args[0]);
    ASSERT_MSG(strncmp(s, "0.7853981", 9) == 0,
               "expected mantissa starting 0.7853981, got: %s", s);
    free(s);
    expr_free(e);
}

static void test_mpfr_identity(void) {
    /* The fundamental identity: m * b^e == x (within MPFR rounding). */
    Expr* p = parse_expression(
        "Module[{r, m, e}, r = MantissaExponent[N[Pi, 40], 7]; "
        "m = r[[1]]; e = r[[2]]; m * 7^e]");
    Expr* res = evaluate(p);
    expr_free(p);
    Expr* p2 = parse_expression("N[Pi, 40]");
    Expr* expected = evaluate(p2);
    expr_free(p2);
    /* The reconstructed value should be close to N[Pi, 40].  We check
     * via string prefix (first 20 decimal digits stay stable). */
    char* a = expr_to_string(res);
    char* b = expr_to_string(expected);
    ASSERT(a && b);
    ASSERT_MSG(strncmp(a, b, 20) == 0,
               "m * b^e mismatch:\n  recon: %s\n  exact: %s", a, b);
    free(a); free(b);
    expr_free(res);
    expr_free(expected);
}
#endif

/* --- Identity check for exact inputs -------------------------------- */

static void test_identity_exact(void) {
    /* m * b^e == x for an integer in an unusual base. */
    assert_eval_eq("Module[{r}, r = MantissaExponent[1027, 2]; "
                   "r[[1]] * 2^r[[2]]]",
                   "1027", 0);
    assert_eval_eq("Module[{r}, r = MantissaExponent[123451]; "
                   "r[[1]] * 10^r[[2]]]",
                   "123451", 0);
    assert_eval_eq("Module[{r}, r = MantissaExponent[-3/2]; "
                   "r[[1]] * 10^r[[2]]]",
                   "-3/2", 0);
}

/* --- Listable ------------------------------------------------------- */

static void test_listable(void) {
    assert_eval_eq(
        "MantissaExponent[{1.234, 5.678}]",
        "{{0.1234, 1}, {0.5678, 1}}", 0);
}

static void test_listable_ints(void) {
    assert_eval_eq(
        "MantissaExponent[{1027, 1027}, 2]",
        "{{1027/2048, 11}, {1027/2048, 11}}", 0);
}

/* --- Error paths --------------------------------------------------- */

static void test_zero_args(void) {
    char* result = NULL;
    char* err = eval_capturing_stderr("MantissaExponent[]", &result);
    ASSERT(result != NULL);
    ASSERT_STR_EQ(result, "MantissaExponent[]");
    ASSERT(err != NULL);
    ASSERT_MSG(strstr(err, "MantissaExponent::argt") != NULL,
               "expected argt diagnostic, got: %s", err);
    ASSERT(strstr(err, "called with 0 arguments") != NULL);
    free(result);
    free(err);
}

static void test_three_args(void) {
    char* result = NULL;
    char* err = eval_capturing_stderr("MantissaExponent[3, 4, 5]", &result);
    ASSERT(result != NULL);
    ASSERT(strstr(result, "MantissaExponent[3, 4, 5]") != NULL);
    ASSERT(err != NULL);
    ASSERT(strstr(err, "MantissaExponent::argt") != NULL);
    ASSERT(strstr(err, "called with 3 arguments") != NULL);
    free(result);
    free(err);
}

static void test_base_one(void) {
    char* result = NULL;
    char* err = eval_capturing_stderr("MantissaExponent[5, 1]", &result);
    ASSERT(result != NULL);
    ASSERT(err && strstr(err, "MantissaExponent::ibase") != NULL);
    free(result);
    free(err);
}

static void test_base_zero(void) {
    char* result = NULL;
    char* err = eval_capturing_stderr("MantissaExponent[5, 0]", &result);
    ASSERT(err && strstr(err, "MantissaExponent::ibase") != NULL);
    free(result);
    free(err);
}

static void test_base_negative(void) {
    char* result = NULL;
    char* err = eval_capturing_stderr("MantissaExponent[5, -2]", &result);
    ASSERT(err && strstr(err, "MantissaExponent::ibase") != NULL);
    free(result);
    free(err);
}

static void test_complex_x(void) {
    char* result = NULL;
    char* err = eval_capturing_stderr("MantissaExponent[2.34 + 9.1 I]", &result);
    /* Result is left unevaluated. */
    ASSERT(result && strstr(result, "MantissaExponent") != NULL);
    ASSERT_MSG(err && strstr(err, "MantissaExponent::realx") != NULL,
               "expected realx diagnostic, got: %s", err ? err : "(null)");
    free(result);
    free(err);
}

static void test_symbolic_x(void) {
    /* MantissaExponent[x] for a plain symbol: left unevaluated, no error. */
    Expr* p = parse_expression("MantissaExponent[xyzfoo]");
    Expr* e = evaluate(p);
    char* s = expr_to_string(e);
    ASSERT(strstr(s, "MantissaExponent") != NULL);
    ASSERT(strstr(s, "xyzfoo") != NULL);
    free(s);
    expr_free(p);
    expr_free(e);
}

static void test_symbolic_x_with_base(void) {
    /* MantissaExponent[x, 2] -- left unevaluated even with explicit base. */
    Expr* p = parse_expression("MantissaExponent[xyzfoo, 2]");
    Expr* e = evaluate(p);
    char* s = expr_to_string(e);
    ASSERT(strstr(s, "MantissaExponent") != NULL);
    free(s);
    expr_free(p);
    expr_free(e);
}

static void test_real_base_unevaluated(void) {
    /* Non-integer base: left unevaluated (no error). */
    Expr* p = parse_expression("MantissaExponent[5, 2.5]");
    Expr* e = evaluate(p);
    char* s = expr_to_string(e);
    ASSERT(strstr(s, "MantissaExponent[5, 2.5]") != NULL);
    free(s);
    expr_free(p);
    expr_free(e);
}

/* --- Attributes / docstring / interned symbol --------------------- */

static void test_attributes(void) {
    SymbolDef* def = symtab_get_def("MantissaExponent");
    ASSERT(def != NULL);
    uint32_t a = get_attributes("MantissaExponent");
    ASSERT((a & ATTR_PROTECTED) != 0);
    ASSERT((a & ATTR_LISTABLE) != 0);
}

static void test_docstring(void) {
    SymbolDef* def = symtab_get_def("MantissaExponent");
    ASSERT(def != NULL);
    ASSERT(def->docstring != NULL);
    ASSERT(strstr(def->docstring, "MantissaExponent[x]") != NULL);
}

static void test_sym_pointer_interned(void) {
    ASSERT(SYM_MantissaExponent != NULL);
    ASSERT(strcmp(SYM_MantissaExponent, "MantissaExponent") == 0);
}

/* --- Memory-safety stress loop ----------------------------------- */

static void test_stress(void) {
    /* Mix exact, machine real, MPFR, listable, error paths -- under
     * valgrind this catches mpz / mpfr / Expr leaks introduced anywhere
     * in the dispatch tree. */
    for (int k = 0; k < 30; k++) {
        assert_eval_eq("MantissaExponent[123451]",
                       "{123451/1000000, 6}", 0);
        assert_eval_eq("MantissaExponent[1027, 2]", "{1027/2048, 11}", 0);
        assert_eval_eq("MantissaExponent[2^100, 2]", "{1/2, 101}", 0);
        assert_eval_eq("MantissaExponent[3.4 10^30]", "{0.34, 31}", 0);
        assert_eval_eq("MantissaExponent[0]", "{0, 0}", 0);
        assert_eval_eq("MantissaExponent[0.]", "{0.0, 0}", 0);
        assert_eval_eq("MantissaExponent[-3/2]", "{-3/20, 1}", 0);
#ifdef USE_MPFR
        Expr* p = parse_expression("MantissaExponent[N[Pi, 30]]");
        Expr* e = evaluate(p);
        expr_free(p);
        expr_free(e);
#endif
    }
}

int main(void) {
    symtab_init();
    core_init();

    TEST(test_int_documented);
    TEST(test_int_one);
    TEST(test_int_ten);
    TEST(test_int_eight_base10);
    TEST(test_int_negative);
    TEST(test_int_zero);
    TEST(test_int_base2);
    TEST(test_int_base16);
    TEST(test_int_base256);
    TEST(test_int_bignum_base2);
    TEST(test_int_bignum_base10);

    TEST(test_rat_one_third);
    TEST(test_rat_three_halves);
    TEST(test_rat_five_eighths_base2);
    TEST(test_rat_negative);
    TEST(test_rat_tiny);

    TEST(test_real_documented);
    TEST(test_real_456);
    TEST(test_real_negative);
    TEST(test_real_zero);
    TEST(test_real_below_one);
    TEST(test_real_one_exactly);
    TEST(test_real_base2);

#ifdef USE_MPFR
    TEST(test_mpfr_pi_30);
    TEST(test_mpfr_pi_50_base2);
    TEST(test_mpfr_identity);
#endif

    TEST(test_identity_exact);

    TEST(test_listable);
    TEST(test_listable_ints);

    TEST(test_zero_args);
    TEST(test_three_args);
    TEST(test_base_one);
    TEST(test_base_zero);
    TEST(test_base_negative);
    TEST(test_complex_x);
    TEST(test_symbolic_x);
    TEST(test_symbolic_x_with_base);
    TEST(test_real_base_unevaluated);

    TEST(test_attributes);
    TEST(test_docstring);
    TEST(test_sym_pointer_interned);

    TEST(test_stress);

    printf("All MantissaExponent tests passed!\n");
    return 0;
}
