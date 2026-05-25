/* Unit tests for RealDigits.
 *
 *   RealDigits[x]                     default base 10, length from precision.
 *   RealDigits[x, b]                  base b.
 *   RealDigits[x, b, len]             exactly `len` digits.
 *   RealDigits[x, b, len, n]          `len` digits starting from b^n.
 *
 * Coverage:
 *   - Integer (terminating) decimal / binary / base-256.
 *   - Bignum integer (2^100).
 *   - Sign discarded for integers and reals.
 *   - Rational terminating (1/4, 1/8, 1/5, 1/10 in various bases).
 *   - Rational recurring decimal: 1/3, 19/7, 1/7, 22/7 -- recurring form.
 *   - Machine-precision reals 123.55555, 0.000012355555, 5.635, 1.234.
 *   - Indeterminate fill past precision.
 *   - 4-arg start offset (positive and negative).
 *   - MPFR via N[Pi, 50], N[Pi, 100]: length from precision.
 *   - Automatic length.
 *   - Symbolic numeric constant Pi gets numericalized when len given.
 *   - Zero special case: RealDigits[0] and RealDigits[0.].
 *   - Listable threading on x, b, both.
 *   - Error paths: 0 args, base 1, base 0, base -2, non-integer base,
 *     real base, negative len, real len, real position n.
 *   - Attribute / docstring / interned-symbol introspection.
 *   - Memory-safety stress loop (catches double-frees under valgrind).
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

/* Capture stderr while `input` is parsed + evaluated.  Returns the
 * collected stderr text as a heap string (caller frees) and writes the
 * printed result into *out_result_str (also heap-allocated).  Uses a
 * fixed temp file path; safe because tests run serially. */
static char* eval_capturing_stderr(const char* input, char** out_result_str) {
    const char* path = "/tmp/mathilda_realdigits_stderr.log";
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

static void test_int_basic(void) {
    /* Documented Mathematica example with the integer form. */
    assert_eval_eq("RealDigits[58127]", "{{5, 8, 1, 2, 7}, 5}", 0);
}

static void test_int_zero(void) {
    assert_eval_eq("RealDigits[0]", "{{0}, 0}", 0);
}

static void test_int_one(void) {
    assert_eval_eq("RealDigits[1]", "{{1}, 1}", 0);
}

static void test_int_ten(void) {
    assert_eval_eq("RealDigits[10]", "{{1, 0}, 2}", 0);
}

static void test_int_long_decimal(void) {
    assert_eval_eq("RealDigits[1234567890]",
                   "{{1, 2, 3, 4, 5, 6, 7, 8, 9, 0}, 10}", 0);
}

static void test_int_base2(void) {
    /* 13 in binary = 1101 -> {{1,1,0,1}, 4}. */
    assert_eval_eq("RealDigits[13, 2]", "{{1, 1, 0, 1}, 4}", 0);
}

static void test_int_base16(void) {
    /* 58127 = 0xE30F -> {{14, 3, 0, 15}, 4}. */
    assert_eval_eq("RealDigits[58127, 16]", "{{14, 3, 0, 15}, 4}", 0);
}

static void test_int_negative(void) {
    /* Sign discarded. */
    assert_eval_eq("RealDigits[-58127]", "{{5, 8, 1, 2, 7}, 5}", 0);
}

/* 2^100 in base 2 = one 1 then 100 zeros: 101 digits, exp = 101. */
static void test_int_bignum_base2(void) {
    Expr* p = parse_expression("RealDigits[2^100, 2]");
    Expr* e = evaluate(p);
    expr_free(p);
    ASSERT(e->type == EXPR_FUNCTION);
    ASSERT(e->data.function.arg_count == 2);
    Expr* digits = e->data.function.args[0];
    Expr* expv   = e->data.function.args[1];
    ASSERT(digits->type == EXPR_FUNCTION);
    ASSERT(digits->data.function.arg_count == 101);
    ASSERT(digits->data.function.args[0]->type == EXPR_INTEGER);
    ASSERT(digits->data.function.args[0]->data.integer == 1);
    for (size_t i = 1; i < 101; i++) {
        ASSERT(digits->data.function.args[i]->type == EXPR_INTEGER);
        ASSERT(digits->data.function.args[i]->data.integer == 0);
    }
    ASSERT(expv->type == EXPR_INTEGER && expv->data.integer == 101);
    expr_free(e);
}

/* --- Rationals (terminating) ---------------------------------------- */

static void test_rat_one_quarter(void) {
    /* 1/4 = 0.25 in base 10 -> {{2, 5}, 0}. */
    assert_eval_eq("RealDigits[1/4]", "{{2, 5}, 0}", 0);
}

static void test_rat_one_half(void) {
    /* 1/2 = 0.5 -> {{5}, 0}. */
    assert_eval_eq("RealDigits[1/2]", "{{5}, 0}", 0);
}

static void test_rat_one_tenth(void) {
    /* 1/10 = 0.1 -> {{1}, 0}. */
    assert_eval_eq("RealDigits[1/10]", "{{1}, 0}", 0);
}

static void test_rat_one_eighth_base2(void) {
    /* 1/8 = 0.001 in binary -> {{1}, -2}. */
    assert_eval_eq("RealDigits[1/8, 2]", "{{1}, -2}", 0);
}

static void test_rat_three_quarter(void) {
    /* 3/4 = 0.75 -> {{7, 5}, 0}. */
    assert_eval_eq("RealDigits[3/4]", "{{7, 5}, 0}", 0);
}

static void test_rat_five_halves(void) {
    /* 5/2 = 2.5 -> {{2, 5}, 1}. */
    assert_eval_eq("RealDigits[5/2]", "{{2, 5}, 1}", 0);
}

/* --- Rationals (recurring decimal) ---------------------------------- */

static void test_rat_one_third(void) {
    /* 1/3 = 0.(3) -> {{{3}}, 0}. */
    assert_eval_eq("RealDigits[1/3]", "{{{3}}, 0}", 0);
}

static void test_rat_one_seventh(void) {
    /* 1/7 = 0.(142857), period 6 -> {{{1,4,2,8,5,7}}, 0}. */
    assert_eval_eq("RealDigits[1/7]",
                   "{{{1, 4, 2, 8, 5, 7}}, 0}", 0);
}

static void test_rat_nineteen_sevenths(void) {
    /* 19/7 = 2.(714285) -> {{2, {7,1,4,2,8,5}}, 1}. */
    assert_eval_eq("RealDigits[19/7]",
                   "{{2, {7, 1, 4, 2, 8, 5}}, 1}", 0);
}

static void test_rat_one_sixth(void) {
    /* 1/6 = 0.1(6) -> {{1, {6}}, 0}. */
    assert_eval_eq("RealDigits[1/6]", "{{1, {6}}, 0}", 0);
}

static void test_rat_explicit_len(void) {
    /* Explicit len forces a flat list even for non-terminating rationals. */
    assert_eval_eq(
        "RealDigits[19/7, 10, 25]",
        "{{2, 7, 1, 4, 2, 8, 5, 7, 1, 4, 2, 8, 5, 7, 1, 4, 2, 8, 5, 7, 1, 4, 2, 8, 5}, 1}",
        0);
}

/* --- Machine reals -------------------------------------------------- */

static void test_real_123_55555(void) {
    /* Documented Mathematica example. */
    assert_eval_eq(
        "RealDigits[123.55555]",
        "{{1, 2, 3, 5, 5, 5, 5, 5, 0, 0, 0, 0, 0, 0, 0, 0}, 3}",
        0);
}

static void test_real_tiny(void) {
    /* 0.000012355555 -> exp = -4, same digits. */
    assert_eval_eq(
        "RealDigits[0.000012355555]",
        "{{1, 2, 3, 5, 5, 5, 5, 5, 0, 0, 0, 0, 0, 0, 0, 0}, -4}",
        0);
}

static void test_real_1234_base2(void) {
    /* Documented Mathematica example. */
    assert_eval_eq(
        "RealDigits[1.234, 2, 15]",
        "{{1, 0, 0, 1, 1, 1, 0, 1, 1, 1, 1, 1, 0, 0, 1}, 1}",
        0);
}

static void test_real_indeterminate_fill(void) {
    /* 5.635 has machine precision ~16 decimal digits; request 20 forces 4
     * Indeterminate slots at the LSB end. */
    assert_eval_eq(
        "RealDigits[5.635, 10, 20]",
        "{{5, 6, 3, 5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,"
        " Indeterminate, Indeterminate, Indeterminate, Indeterminate}, 1}",
        0);
}

static void test_real_indet_with_offset(void) {
    /* RealDigits[5.635, 10, 20, -15] -> {{0, Indeterminate * 19}, -14}. */
    Expr* p = parse_expression("RealDigits[5.635, 10, 20, -15]");
    Expr* e = evaluate(p);
    expr_free(p);
    ASSERT(e->type == EXPR_FUNCTION);
    Expr* digits = e->data.function.args[0];
    Expr* expv   = e->data.function.args[1];
    ASSERT(digits->data.function.arg_count == 20);
    ASSERT(digits->data.function.args[0]->type == EXPR_INTEGER);
    ASSERT(digits->data.function.args[0]->data.integer == 0);
    for (size_t i = 1; i < 20; i++) {
        Expr* d = digits->data.function.args[i];
        ASSERT(d->type == EXPR_SYMBOL && d->data.symbol == SYM_Indeterminate);
    }
    ASSERT(expv->type == EXPR_INTEGER && expv->data.integer == -14);
    expr_free(e);
}

static void test_real_negative_arg(void) {
    /* Sign discarded for inexact reals too. */
    assert_eval_eq(
        "RealDigits[-123.55555]",
        "{{1, 2, 3, 5, 5, 5, 5, 5, 0, 0, 0, 0, 0, 0, 0, 0}, 3}",
        0);
}

/* --- Pi numericalize, len ------------------------------------------ */

static void test_pi_25(void) {
    assert_eval_eq(
        "RealDigits[Pi, 10, 25]",
        "{{3, 1, 4, 1, 5, 9, 2, 6, 5, 3, 5, 8, 9, 7, 9, 3, 2, 3, 8, 4, 6, 2, 6, 4, 3}, 1}",
        0);
}

static void test_pi_20_neg5(void) {
    /* RealDigits[Pi, 10, 20, -5] -- 20 digits starting at 10^-5. */
    assert_eval_eq(
        "RealDigits[Pi, 10, 20, -5]",
        "{{9, 2, 6, 5, 3, 5, 8, 9, 7, 9, 3, 2, 3, 8, 4, 6, 2, 6, 4, 3}, -4}",
        0);
}

static void test_pi_20_pos5(void) {
    /* RealDigits[Pi, 10, 20, 5] -- 20 digits starting at 10^5. */
    assert_eval_eq(
        "RealDigits[Pi, 10, 20, 5]",
        "{{0, 0, 0, 0, 0, 3, 1, 4, 1, 5, 9, 2, 6, 5, 3, 5, 8, 9, 7, 9}, 6}",
        0);
}

/* --- MPFR ---------------------------------------------------------- */

static void test_mpfr_pi_50(void) {
    /* N[Pi, 50] gives a 50-digit MPFR; RealDigits should return ~50 digits. */
    Expr* p = parse_expression("RealDigits[N[Pi, 50]]");
    Expr* e = evaluate(p);
    expr_free(p);
    ASSERT(e->type == EXPR_FUNCTION);
    Expr* digits = e->data.function.args[0];
    Expr* expv   = e->data.function.args[1];
    /* The list should start with 3, 1, 4, 1, 5, ... */
    ASSERT(digits->data.function.arg_count >= 50);
    ASSERT(digits->data.function.args[0]->data.integer == 3);
    ASSERT(digits->data.function.args[1]->data.integer == 1);
    ASSERT(digits->data.function.args[2]->data.integer == 4);
    ASSERT(digits->data.function.args[3]->data.integer == 1);
    ASSERT(digits->data.function.args[4]->data.integer == 5);
    ASSERT(expv->type == EXPR_INTEGER && expv->data.integer == 1);
    expr_free(e);
}

static void test_mpfr_explicit_len(void) {
    /* RealDigits[N[Pi, 100], 10, 30] should give 30 digits. */
    Expr* p = parse_expression("RealDigits[N[Pi, 100], 10, 30]");
    Expr* e = evaluate(p);
    expr_free(p);
    Expr* digits = e->data.function.args[0];
    ASSERT(digits->data.function.arg_count == 30);
    /* Pi starts 3.14159265358979323846264338327950... */
    static const int expected[] = {3,1,4,1,5,9,2,6,5,3,5,8,9,7,9,3,
                                    2,3,8,4,6,2,6,4,3,3,8,3,2,7};
    for (int i = 0; i < 30; i++) {
        ASSERT(digits->data.function.args[i]->type == EXPR_INTEGER);
        ASSERT(digits->data.function.args[i]->data.integer == expected[i]);
    }
    expr_free(e);
}

/* --- Zero special cases -------------------------------------------- */

static void test_zero_exact(void) {
    assert_eval_eq("RealDigits[0]", "{{0}, 0}", 0);
}

static void test_zero_machine_real(void) {
    /* {{0}, -Floor[Accuracy[0.]]}: with Mathematica-compatible
     * Accuracy[0.] ≈ 323.607, the exponent is -323. */
    assert_eval_eq("RealDigits[0.]", "{{0}, -323}", 0);
}

/* --- Listable ------------------------------------------------------ */

static void test_listable_on_x(void) {
    assert_eval_eq(
        "RealDigits[{1.234, 5.678}]",
        "{{{1, 2, 3, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, 1},"
        " {{5, 6, 7, 8, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, 1}}",
        0);
}

static void test_listable_int(void) {
    assert_eval_eq("RealDigits[{12, 7}]",
                   "{{{1, 2}, 2}, {{7}, 1}}", 0);
}

/* --- Error paths --------------------------------------------------- */

static void test_zero_args(void) {
    char* result = NULL;
    char* err = eval_capturing_stderr("RealDigits[]", &result);
    ASSERT(result != NULL);
    ASSERT_STR_EQ(result, "RealDigits[]");
    ASSERT(err != NULL);
    ASSERT_MSG(strstr(err, "RealDigits::argb") != NULL,
               "expected argb diagnostic, got: %s", err);
    ASSERT(strstr(err, "called with 0 arguments") != NULL);
    free(result);
    free(err);
}

static void test_five_args(void) {
    char* result = NULL;
    char* err = eval_capturing_stderr("RealDigits[1, 2, 3, 4, 5]", &result);
    ASSERT(result != NULL);
    ASSERT(strstr(result, "RealDigits[1, 2, 3, 4, 5]") != NULL);
    ASSERT(err != NULL);
    ASSERT(strstr(err, "RealDigits::argb") != NULL);
    free(result);
    free(err);
}

static void test_base_one(void) {
    char* result = NULL;
    char* err = eval_capturing_stderr("RealDigits[5, 1]", &result);
    ASSERT(result != NULL);
    ASSERT(strstr(result, "RealDigits") != NULL);
    ASSERT(err != NULL);
    ASSERT_MSG(strstr(err, "RealDigits::ibase") != NULL,
               "expected ibase diagnostic, got: %s", err);
    free(result);
    free(err);
}

static void test_base_zero(void) {
    char* result = NULL;
    char* err = eval_capturing_stderr("RealDigits[5, 0]", &result);
    ASSERT(err && strstr(err, "RealDigits::ibase") != NULL);
    free(result);
    free(err);
}

static void test_base_negative(void) {
    char* result = NULL;
    char* err = eval_capturing_stderr("RealDigits[5, -2]", &result);
    ASSERT(err && strstr(err, "RealDigits::ibase") != NULL);
    free(result);
    free(err);
}

static void test_base_real(void) {
    /* Non-integer base -> ::ibase, not ::int. */
    char* result = NULL;
    char* err = eval_capturing_stderr("RealDigits[5, 2.5]", &result);
    ASSERT(err && strstr(err, "RealDigits::ibase") != NULL);
    free(result);
    free(err);
}

static void test_negative_len(void) {
    char* result = NULL;
    char* err = eval_capturing_stderr("RealDigits[5, 10, -1]", &result);
    ASSERT(err && strstr(err, "RealDigits::intnn") != NULL);
    free(result);
    free(err);
}

static void test_real_len(void) {
    char* result = NULL;
    char* err = eval_capturing_stderr("RealDigits[5, 10, 3.5]", &result);
    ASSERT(err && strstr(err, "RealDigits::intnn") != NULL);
    free(result);
    free(err);
}

static void test_real_pos_n(void) {
    char* result = NULL;
    char* err = eval_capturing_stderr("RealDigits[5, 10, 3, 1.5]", &result);
    ASSERT(err && strstr(err, "RealDigits::int") != NULL);
    free(result);
    free(err);
}

static void test_symbolic_x_no_len(void) {
    /* RealDigits[x] with symbolic x -- left unevaluated, no diagnostic. */
    Expr* p = parse_expression("RealDigits[xyzfoo]");
    Expr* e = evaluate(p);
    char* s = expr_to_string(e);
    ASSERT(strstr(s, "RealDigits") != NULL);
    free(s);
    expr_free(p);
    expr_free(e);
}

static void test_complex_x(void) {
    /* RealDigits of a complex number -> ::realx. */
    char* result = NULL;
    char* err = eval_capturing_stderr("RealDigits[1.0 + 2.0 I]", &result);
    ASSERT(err && strstr(err, "RealDigits::realx") != NULL);
    free(result);
    free(err);
}

/* --- Attributes / docstring / interned symbol --------------------- */

static void test_attributes(void) {
    SymbolDef* def = symtab_get_def("RealDigits");
    ASSERT(def != NULL);
    uint32_t a = get_attributes("RealDigits");
    ASSERT((a & ATTR_PROTECTED) != 0);
    ASSERT((a & ATTR_LISTABLE) != 0);
}

static void test_docstring(void) {
    SymbolDef* def = symtab_get_def("RealDigits");
    ASSERT(def != NULL);
    ASSERT(def->docstring != NULL);
    ASSERT(strstr(def->docstring, "RealDigits[x]") != NULL);
}

static void test_sym_pointer_interned(void) {
    ASSERT(SYM_RealDigits != NULL);
    ASSERT(strcmp(SYM_RealDigits, "RealDigits") == 0);
}

/* --- Memory-safety stress loop ----------------------------------- */

static void test_stress(void) {
    /* Mix exact, machine real, MPFR, listable, error -- all paths in
     * one loop so any leaked mpz/mpfr/Expr surfaces under valgrind. */
    for (int k = 0; k < 30; k++) {
        assert_eval_eq("RealDigits[58127]", "{{5, 8, 1, 2, 7}, 5}", 0);
        assert_eval_eq("RealDigits[19/7]", "{{2, {7, 1, 4, 2, 8, 5}}, 1}", 0);
        assert_eval_eq(
            "RealDigits[123.55555]",
            "{{1, 2, 3, 5, 5, 5, 5, 5, 0, 0, 0, 0, 0, 0, 0, 0}, 3}",
            0);
        assert_eval_eq("RealDigits[0]", "{{0}, 0}", 0);
        assert_eval_eq("RealDigits[0.]", "{{0}, -323}", 0);
        assert_eval_eq("RealDigits[1/8, 2]", "{{1}, -2}", 0);
        Expr* p = parse_expression("RealDigits[Pi, 10, 25]");
        Expr* e = evaluate(p);
        expr_free(p);
        expr_free(e);
    }
}

int main(void) {
    symtab_init();
    core_init();

    TEST(test_int_basic);
    TEST(test_int_zero);
    TEST(test_int_one);
    TEST(test_int_ten);
    TEST(test_int_long_decimal);
    TEST(test_int_base2);
    TEST(test_int_base16);
    TEST(test_int_negative);
    TEST(test_int_bignum_base2);

    TEST(test_rat_one_quarter);
    TEST(test_rat_one_half);
    TEST(test_rat_one_tenth);
    TEST(test_rat_one_eighth_base2);
    TEST(test_rat_three_quarter);
    TEST(test_rat_five_halves);

    TEST(test_rat_one_third);
    TEST(test_rat_one_seventh);
    TEST(test_rat_nineteen_sevenths);
    TEST(test_rat_one_sixth);
    TEST(test_rat_explicit_len);

    TEST(test_real_123_55555);
    TEST(test_real_tiny);
    TEST(test_real_1234_base2);
    TEST(test_real_indeterminate_fill);
    TEST(test_real_indet_with_offset);
    TEST(test_real_negative_arg);

    TEST(test_pi_25);
    TEST(test_pi_20_neg5);
    TEST(test_pi_20_pos5);

#ifdef USE_MPFR
    TEST(test_mpfr_pi_50);
    TEST(test_mpfr_explicit_len);
#endif

    TEST(test_zero_exact);
    TEST(test_zero_machine_real);

    TEST(test_listable_on_x);
    TEST(test_listable_int);

    TEST(test_zero_args);
    TEST(test_five_args);
    TEST(test_base_one);
    TEST(test_base_zero);
    TEST(test_base_negative);
    TEST(test_base_real);
    TEST(test_negative_len);
    TEST(test_real_len);
    TEST(test_real_pos_n);
    TEST(test_symbolic_x_no_len);
    TEST(test_complex_x);

    TEST(test_attributes);
    TEST(test_docstring);
    TEST(test_sym_pointer_interned);

    TEST(test_stress);

    printf("All RealDigits tests passed!\n");
    return 0;
}
