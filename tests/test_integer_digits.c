/* Unit tests for IntegerDigits.
 *
 *   IntegerDigits[n]            -- decimal digits, MSB-first.
 *   IntegerDigits[n, b]         -- base-b digits.
 *   IntegerDigits[n, b, len]    -- pad / truncate to length len.
 *
 * Coverage:
 *   - Decimal-digit decomposition: 58127, 1, 10, large multi-digit.
 *   - Binary, hex (base 16, includes digits > 9), base 8, large base 256.
 *   - Zero special-case (with and without padding).
 *   - Negative integers (sign discarded).
 *   - Padding: len > digit count, len == digit count, len == 0.
 *   - Truncation: len < digit count returns last len LSDs.
 *   - Bignum inputs: 10^30 and 2^100, including base-2 expansion of 2^100.
 *   - Listable threading on n, on b, on both, and with explicit len.
 *   - Error paths: base 1, base 0, base -2, negative len, real argument,
 *     symbolic argument, missing args.
 *   - Attribute, docstring, and interned-symbol introspection.
 *   - Repeated-evaluation stress loop to catch double-frees under valgrind.
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
    const char* path = "/tmp/mathilda_int_digits_stderr.log";
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

/* --- Decimal --------------------------------------------------------- */

static void test_decimal_basic(void) {
    assert_eval_eq("IntegerDigits[58127]", "{5, 8, 1, 2, 7}", 0);
}

static void test_decimal_one(void) {
    assert_eval_eq("IntegerDigits[1]", "{1}", 0);
}

static void test_decimal_ten(void) {
    assert_eval_eq("IntegerDigits[10]", "{1, 0}", 0);
}

static void test_decimal_long(void) {
    assert_eval_eq("IntegerDigits[1234567890]",
                   "{1, 2, 3, 4, 5, 6, 7, 8, 9, 0}", 0);
}

/* --- Other bases ----------------------------------------------------- */

static void test_base2(void) {
    /* 58127 in binary, 16 digits. */
    assert_eval_eq("IntegerDigits[58127, 2]",
                   "{1, 1, 1, 0, 0, 0, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1}", 0);
}

static void test_base16(void) {
    /* 58127_10 = E30F_16 -> digits {14, 3, 0, 15}. */
    assert_eval_eq("IntegerDigits[58127, 16]", "{14, 3, 0, 15}", 0);
}

static void test_base8(void) {
    /* 58127_10 = 161417_8 -> digits {1,6,1,4,1,7}. */
    assert_eval_eq("IntegerDigits[58127, 8]", "{1, 6, 1, 4, 1, 7}", 0);
}

static void test_base_large(void) {
    /* 58127 in base 256: 58127 = 227 * 256 + 15 = 0xE30F.
     * So the digits are {227, 15}. */
    assert_eval_eq("IntegerDigits[58127, 256]", "{227, 15}", 0);
}

static void test_small_in_base2(void) {
    assert_eval_eq("IntegerDigits[7, 2]", "{1, 1, 1}", 0);
}

/* --- Zero ----------------------------------------------------------- */

static void test_zero_decimal(void) {
    assert_eval_eq("IntegerDigits[0]", "{0}", 0);
}

static void test_zero_base2(void) {
    assert_eval_eq("IntegerDigits[0, 2]", "{0}", 0);
}

static void test_zero_padded(void) {
    /* Padding 0 to length 5 in any base gives 5 zeros. */
    assert_eval_eq("IntegerDigits[0, 2, 5]", "{0, 0, 0, 0, 0}", 0);
    assert_eval_eq("IntegerDigits[0, 10, 3]", "{0, 0, 0}", 0);
}

/* --- Sign discarded ------------------------------------------------- */

static void test_negative_int(void) {
    assert_eval_eq("IntegerDigits[-58127]", "{5, 8, 1, 2, 7}", 0);
}

static void test_negative_with_base(void) {
    assert_eval_eq("IntegerDigits[-7, 2]", "{1, 1, 1}", 0);
}

/* --- Padding & truncation ------------------------------------------ */

static void test_pad_when_longer(void) {
    /* IntegerDigits[Range[0,7], 2, 3] documented example. */
    assert_eval_eq(
        "IntegerDigits[Range[0, 7], 2, 3]",
        "{{0, 0, 0}, {0, 0, 1}, {0, 1, 0}, {0, 1, 1},"
        " {1, 0, 0}, {1, 0, 1}, {1, 1, 0}, {1, 1, 1}}",
        0);
}

static void test_pad_equal_length(void) {
    /* len exactly matches the natural digit count: no zeros prepended. */
    assert_eval_eq("IntegerDigits[58127, 10, 5]", "{5, 8, 1, 2, 7}", 0);
}

static void test_pad_to_zero(void) {
    /* len == 0 yields the empty list -- unusual but well-defined. */
    assert_eval_eq("IntegerDigits[58127, 10, 0]", "{}", 0);
}

static void test_truncate_last_n(void) {
    /* Documented example: keep the 4 least-significant digits of 6345354. */
    assert_eval_eq("IntegerDigits[6345354, 10, 4]", "{5, 3, 5, 4}", 0);
}

static void test_truncate_to_one(void) {
    /* The last digit of 58127 is 7. */
    assert_eval_eq("IntegerDigits[58127, 10, 1]", "{7}", 0);
}

static void test_truncate_base2(void) {
    /* 58127 in binary has 16 digits; keep the bottom 4 = 1111. */
    assert_eval_eq("IntegerDigits[58127, 2, 4]", "{1, 1, 1, 1}", 0);
}

static void test_pad_with_base16(void) {
    /* 58127_10 = E30F_16; padded to length 6 prepends two zeros. */
    assert_eval_eq("IntegerDigits[58127, 16, 6]", "{0, 0, 14, 3, 0, 15}", 0);
}

/* --- Bignum inputs -------------------------------------------------- */

static void test_bignum_decimal(void) {
    /* 10^30 has 31 digits in base 10: one '1' followed by thirty '0's. */
    assert_eval_eq(
        "IntegerDigits[10^30]",
        "{1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,"
        " 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}",
        0);
}

static void test_bignum_int64_overflow(void) {
    /* 2^63 = 9223372036854775808 is one past INT64_MAX, so the input
     * itself is a bigint.  Result is its 19 decimal digits. */
    assert_eval_eq(
        "IntegerDigits[2^63]",
        "{9, 2, 2, 3, 3, 7, 2, 0, 3, 6, 8, 5, 4, 7, 7, 5, 8, 0, 8}",
        0);
}

static void test_bignum_base2(void) {
    /* 2^100 in base 2 = one '1' followed by 100 zeros (101 digits). */
    Expr* p = parse_expression("IntegerDigits[2^100, 2]");
    Expr* e = evaluate(p);
    expr_free(p);
    /* Verify the structure: List with 101 args, first = 1, rest = 0. */
    ASSERT(e->type == EXPR_FUNCTION);
    ASSERT(e->data.function.head->type == EXPR_SYMBOL);
    ASSERT(e->data.function.head->data.symbol.name == SYM_List);
    ASSERT(e->data.function.arg_count == 101);
    ASSERT(e->data.function.args[0]->type == EXPR_INTEGER);
    ASSERT(e->data.function.args[0]->data.integer == 1);
    for (size_t i = 1; i < 101; i++) {
        ASSERT(e->data.function.args[i]->type == EXPR_INTEGER);
        ASSERT(e->data.function.args[i]->data.integer == 0);
    }
    expr_free(e);
}

static void test_bignum_negative(void) {
    /* Sign discarded on a bignum input. */
    assert_eval_eq(
        "IntegerDigits[-(10^15 + 1)]",
        "{1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1}",
        0);
}

static void test_bignum_truncate(void) {
    /* Keep only the 5 least-significant decimal digits of 10^30 + 12345. */
    assert_eval_eq("IntegerDigits[10^30 + 12345, 10, 5]",
                   "{1, 2, 3, 4, 5}", 0);
}

/* --- Listable threading --------------------------------------------- */

static void test_listable_on_n(void) {
    /* Documented example: IntegerDigits[{6, 7, 2}, 2]. */
    assert_eval_eq("IntegerDigits[{6, 7, 2}, 2]",
                   "{{1, 1, 0}, {1, 1, 1}, {1, 0}}", 0);
}

static void test_listable_on_base(void) {
    /* Documented example: IntegerDigits[7, {2, 3, 4}]. */
    assert_eval_eq("IntegerDigits[7, {2, 3, 4}]",
                   "{{1, 1, 1}, {2, 1}, {1, 3}}", 0);
}

static void test_listable_no_explicit_base(void) {
    /* 1-arg form threading. */
    assert_eval_eq("IntegerDigits[{0, 1, 12}]",
                   "{{0}, {1}, {1, 2}}", 0);
}

static void test_listable_with_len(void) {
    /* len is a scalar -- still threaded along n. */
    assert_eval_eq("IntegerDigits[{5, 13}, 2, 5]",
                   "{{0, 0, 1, 0, 1}, {0, 1, 1, 0, 1}}", 0);
}

/* --- Error paths ---------------------------------------------------- */

static void test_bad_base_one(void) {
    /* Base 1 prints IntegerDigits::ibase and leaves the call unevaluated. */
    const char* in = "IntegerDigits[5, 1]";
    Expr* p = parse_expression(in);
    Expr* e = evaluate(p);
    char* s = expr_to_string(e);
    ASSERT_MSG(strstr(s, "IntegerDigits") != NULL,
               "expected unevaluated form, got: %s", s);
    free(s);
    expr_free(p);
    expr_free(e);
}

static void test_bad_base_zero(void) {
    const char* in = "IntegerDigits[5, 0]";
    Expr* p = parse_expression(in);
    Expr* e = evaluate(p);
    char* s = expr_to_string(e);
    ASSERT_MSG(strstr(s, "IntegerDigits") != NULL,
               "expected unevaluated form, got: %s", s);
    free(s);
    expr_free(p);
    expr_free(e);
}

static void test_bad_base_negative(void) {
    const char* in = "IntegerDigits[5, -2]";
    Expr* p = parse_expression(in);
    Expr* e = evaluate(p);
    char* s = expr_to_string(e);
    ASSERT_MSG(strstr(s, "IntegerDigits") != NULL,
               "expected unevaluated form, got: %s", s);
    free(s);
    expr_free(p);
    expr_free(e);
}

static void test_bad_len_negative(void) {
    const char* in = "IntegerDigits[5, 10, -1]";
    Expr* p = parse_expression(in);
    Expr* e = evaluate(p);
    char* s = expr_to_string(e);
    ASSERT_MSG(strstr(s, "IntegerDigits") != NULL,
               "expected unevaluated form, got: %s", s);
    free(s);
    expr_free(p);
    expr_free(e);
}

static void test_real_n_left_unevaluated(void) {
    const char* in = "IntegerDigits[1.5]";
    Expr* p = parse_expression(in);
    Expr* e = evaluate(p);
    char* s = expr_to_string(e);
    ASSERT_MSG(strstr(s, "IntegerDigits") != NULL,
               "expected unevaluated form, got: %s", s);
    free(s);
    expr_free(p);
    expr_free(e);
}

static void test_symbolic_n_left_unevaluated(void) {
    const char* in = "IntegerDigits[x]";
    Expr* p = parse_expression(in);
    Expr* e = evaluate(p);
    char* s = expr_to_string(e);
    ASSERT_MSG(strstr(s, "IntegerDigits") != NULL,
               "expected unevaluated form, got: %s", s);
    free(s);
    expr_free(p);
    expr_free(e);
}

static void test_zero_args_left_unevaluated(void) {
    char* result = NULL;
    char* err = eval_capturing_stderr("IntegerDigits[]", &result);
    ASSERT(result != NULL);
    ASSERT_STR_EQ(result, "IntegerDigits[]");
    /* Mathematica diagnostic: IntegerDigits::argb with arg count 0. */
    ASSERT_MSG(err && strstr(err, "IntegerDigits::argb") != NULL,
               "expected argb diagnostic, got: %s", err ? err : "(null)");
    ASSERT(strstr(err, "called with 0 arguments") != NULL);
    ASSERT(strstr(err, "between 1 and 3 arguments are expected") != NULL);
    free(result);
    free(err);
}

static void test_four_args_left_unevaluated(void) {
    char* result = NULL;
    char* err = eval_capturing_stderr("IntegerDigits[5, 10, 3, 7]", &result);
    ASSERT(result != NULL);
    ASSERT(strstr(result, "IntegerDigits") != NULL);
    ASSERT_MSG(err && strstr(err, "IntegerDigits::argb") != NULL,
               "expected argb diagnostic, got: %s", err ? err : "(null)");
    ASSERT(strstr(err, "called with 4 arguments") != NULL);
    free(result);
    free(err);
}

static void test_five_args_argb(void) {
    /* Documented Mathematica example: 5-arg call. */
    char* result = NULL;
    char* err = eval_capturing_stderr("IntegerDigits[1, 2, 3, 4, 5]", &result);
    ASSERT(result != NULL);
    ASSERT(strstr(result, "IntegerDigits[1, 2, 3, 4, 5]") != NULL);
    ASSERT(err != NULL);
    ASSERT_MSG(strstr(err, "IntegerDigits::argb") != NULL,
               "expected argb diagnostic, got: %s", err);
    ASSERT(strstr(err, "called with 5 arguments") != NULL);
    free(result);
    free(err);
}

/* Singular "argument" (vs plural "arguments") for argc == 1.  We can't
 * actually invoke IntegerDigits with 1 argument as a failure (that's the
 * success case), so the singular form is exercised indirectly only.  This
 * test just guards the pluralisation logic in the 4+ case. */
static void test_argb_plural_form(void) {
    char* result = NULL;
    char* err = eval_capturing_stderr("IntegerDigits[1, 2, 3, 4]", &result);
    ASSERT(err != NULL);
    ASSERT(strstr(err, "called with 4 arguments;") != NULL);
    free(result);
    free(err);
}

static void test_real_n_int_diagnostic(void) {
    /* Mathematica's text: "Integer expected at position 1 in IntegerDigits[1.1234]." */
    char* result = NULL;
    char* err = eval_capturing_stderr("IntegerDigits[1.1234]", &result);
    ASSERT(result != NULL);
    ASSERT(strstr(result, "IntegerDigits[1.1234]") != NULL);
    ASSERT(err != NULL);
    ASSERT_MSG(strstr(err, "IntegerDigits::int") != NULL,
               "expected int diagnostic, got: %s", err);
    ASSERT(strstr(err, "Integer expected at position 1") != NULL);
    ASSERT(strstr(err, "IntegerDigits[1.1234]") != NULL);
    free(result);
    free(err);
}

static void test_complex_n_int_diagnostic(void) {
    /* Documented Mathematica example: complex first arg. */
    char* result = NULL;
    char* err = eval_capturing_stderr("IntegerDigits[1.1234 - 9 I]", &result);
    ASSERT(result != NULL);
    ASSERT(err != NULL);
    ASSERT_MSG(strstr(err, "IntegerDigits::int") != NULL,
               "expected int diagnostic, got: %s", err);
    ASSERT(strstr(err, "Integer expected at position 1") != NULL);
    free(result);
    free(err);
}

static void test_rational_n_int_diagnostic(void) {
    /* Exact rational also fails the integer test. */
    char* result = NULL;
    char* err = eval_capturing_stderr("IntegerDigits[3/2]", &result);
    ASSERT(result != NULL);
    ASSERT(err != NULL);
    ASSERT(strstr(err, "IntegerDigits::int") != NULL);
    ASSERT(strstr(err, "Integer expected at position 1") != NULL);
    free(result);
    free(err);
}

static void test_real_base_int_diagnostic(void) {
    char* result = NULL;
    char* err = eval_capturing_stderr("IntegerDigits[10, 2.5]", &result);
    ASSERT(result != NULL);
    ASSERT(err != NULL);
    ASSERT(strstr(err, "IntegerDigits::int") != NULL);
    ASSERT(strstr(err, "Integer expected at position 2") != NULL);
    free(result);
    free(err);
}

static void test_real_len_int_diagnostic(void) {
    char* result = NULL;
    char* err = eval_capturing_stderr("IntegerDigits[10, 2, 3.5]", &result);
    ASSERT(result != NULL);
    ASSERT(err != NULL);
    ASSERT(strstr(err, "IntegerDigits::int") != NULL);
    ASSERT(strstr(err, "Integer expected at position 3") != NULL);
    free(result);
    free(err);
}

/* --- Attribute / docstring / interned-symbol introspection ---------- */

static void test_attributes(void) {
    SymbolDef* def = symtab_get_def("IntegerDigits");
    ASSERT(def != NULL);
    /* attr.c table only carries the static base set; for a precise check
     * use the merged attribute resolver. */
    uint32_t a = get_attributes("IntegerDigits");
    ASSERT((a & ATTR_PROTECTED) != 0);
    ASSERT((a & ATTR_LISTABLE) != 0);
}

static void test_docstring_set(void) {
    SymbolDef* def = symtab_get_def("IntegerDigits");
    ASSERT(def != NULL);
    ASSERT(def->docstring != NULL);
    ASSERT(strstr(def->docstring, "decimal digits") != NULL);
    ASSERT(strstr(def->docstring, "base b") != NULL);
}

static void test_sym_pointer_interned(void) {
    ASSERT(SYM_IntegerDigits != NULL);
    ASSERT(strcmp(SYM_IntegerDigits, "IntegerDigits") == 0);
}

/* --- Memory-safety stress loop -------------------------------------- */

static void test_repeated_evaluation_does_not_corrupt(void) {
    /* Mix machine-int / bignum / threaded / truncated / padded / zero
     * paths so any mishandled mpz_t or args[] slot surfaces under valgrind. */
    for (int k = 0; k < 50; k++) {
        assert_eval_eq("IntegerDigits[58127]", "{5, 8, 1, 2, 7}", 0);
        assert_eval_eq("IntegerDigits[58127, 16]", "{14, 3, 0, 15}", 0);
        assert_eval_eq("IntegerDigits[0]", "{0}", 0);
        assert_eval_eq("IntegerDigits[0, 2, 5]", "{0, 0, 0, 0, 0}", 0);
        assert_eval_eq("IntegerDigits[6345354, 10, 4]", "{5, 3, 5, 4}", 0);
        assert_eval_eq("IntegerDigits[{6, 7, 2}, 2]",
                       "{{1, 1, 0}, {1, 1, 1}, {1, 0}}", 0);
        assert_eval_eq("IntegerDigits[2^63]",
                       "{9, 2, 2, 3, 3, 7, 2, 0, 3, 6, 8, 5, 4, 7, 7, 5, 8, 0, 8}",
                       0);
    }
}

int main(void) {
    symtab_init();
    core_init();

    TEST(test_decimal_basic);
    TEST(test_decimal_one);
    TEST(test_decimal_ten);
    TEST(test_decimal_long);

    TEST(test_base2);
    TEST(test_base16);
    TEST(test_base8);
    TEST(test_base_large);
    TEST(test_small_in_base2);

    TEST(test_zero_decimal);
    TEST(test_zero_base2);
    TEST(test_zero_padded);

    TEST(test_negative_int);
    TEST(test_negative_with_base);

    TEST(test_pad_when_longer);
    TEST(test_pad_equal_length);
    TEST(test_pad_to_zero);
    TEST(test_truncate_last_n);
    TEST(test_truncate_to_one);
    TEST(test_truncate_base2);
    TEST(test_pad_with_base16);

    TEST(test_bignum_decimal);
    TEST(test_bignum_int64_overflow);
    TEST(test_bignum_base2);
    TEST(test_bignum_negative);
    TEST(test_bignum_truncate);

    TEST(test_listable_on_n);
    TEST(test_listable_on_base);
    TEST(test_listable_no_explicit_base);
    TEST(test_listable_with_len);

    TEST(test_bad_base_one);
    TEST(test_bad_base_zero);
    TEST(test_bad_base_negative);
    TEST(test_bad_len_negative);
    TEST(test_real_n_left_unevaluated);
    TEST(test_symbolic_n_left_unevaluated);
    TEST(test_zero_args_left_unevaluated);
    TEST(test_four_args_left_unevaluated);
    TEST(test_five_args_argb);
    TEST(test_argb_plural_form);
    TEST(test_real_n_int_diagnostic);
    TEST(test_complex_n_int_diagnostic);
    TEST(test_rational_n_int_diagnostic);
    TEST(test_real_base_int_diagnostic);
    TEST(test_real_len_int_diagnostic);

    TEST(test_attributes);
    TEST(test_docstring_set);
    TEST(test_sym_pointer_interned);

    TEST(test_repeated_evaluation_does_not_corrupt);

    printf("All IntegerDigits tests passed!\n");
    return 0;
}
