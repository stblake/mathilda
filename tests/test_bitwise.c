/* Unit tests for the bit-level integer builtins (src/bitwise/).
 *
 *   BitLength[n] -- number of binary bits needed to represent the integer n.
 *
 * Coverage:
 *   - Small positives, powers of two (2^k -> k+1, 2^k-1 -> k), zero.
 *   - Negatives: BitLength[n] = BitLength[BitNot[n]], BitNot[n] = -n-1, so
 *     BitLength[-1]=0, BitLength[-2]=1, BitLength[-2^k]=k.
 *   - int64 boundary incl. INT64_MIN (-2^63 -> 63), where the int64 kernel's
 *     complement path is exact and IntegerLength's negate path would overflow.
 *   - Bignums: 2^100 -> 101, 2^100-1 -> 100, -2^100 -> 100, 100! -> 525.
 *   - Listable threading (flat + nested) and the packed / visible-NDArray fast
 *     paths (single-headed int64 result, packed == visible == List).
 *   - Compile[] round-trip at scalar and rank-1 array shapes.
 *   - Property tests: BitLength[n] == IntegerLength[n, 2] and
 *     == Floor[Log[2, n]] + 1 for n > 0.
 *   - Error paths: real / rational / complex arg (::int), wrong arity (::argx),
 *     symbolic arg (silent, unevaluated).
 *   - Attribute, docstring, interned-symbol introspection.
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

/* Capture stderr while `input` is parsed + evaluated.  Returns the collected
 * stderr text as a heap string (caller frees) and writes the printed result
 * into *out_result_str (also heap-allocated).  Fixed temp path; safe because
 * tests run serially. */
static char* eval_capturing_stderr(const char* input, char** out_result_str) {
    const char* path = "/tmp/mathilda_bitwise_stderr.log";
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

/* --- Small positives ------------------------------------------------- */

static void test_zero(void) {
    assert_eval_eq("BitLength[0]", "0", 0);
}

static void test_small_positives(void) {
    assert_eval_eq("BitLength[1]", "1", 0);
    assert_eval_eq("BitLength[2]", "2", 0);
    assert_eval_eq("BitLength[3]", "2", 0);
    assert_eval_eq("BitLength[4]", "3", 0);
    assert_eval_eq("BitLength[5]", "3", 0);
    assert_eval_eq("BitLength[7]", "3", 0);
    assert_eval_eq("BitLength[8]", "4", 0);
}

static void test_byte_boundaries(void) {
    assert_eval_eq("BitLength[255]", "8", 0);
    assert_eval_eq("BitLength[256]", "9", 0);
    assert_eval_eq("BitLength[1023]", "10", 0);
    assert_eval_eq("BitLength[1024]", "11", 0);
}

static void test_powers_of_two(void) {
    /* 2^k needs k+1 bits; 2^k - 1 needs k bits. */
    assert_eval_eq("BitLength[2^10]", "11", 0);
    assert_eval_eq("BitLength[2^10 - 1]", "10", 0);
    assert_eval_eq("BitLength[2^30]", "31", 0);
    assert_eval_eq("BitLength[2^30 - 1]", "30", 0);
}

/* --- Negatives (BitNot semantics) ------------------------------------ */

static void test_negatives(void) {
    /* BitLength[n] = BitLength[BitNot[n]], BitNot[n] = -n-1. */
    assert_eval_eq("BitLength[-1]", "0", 0);   /* BitNot[-1] = 0  */
    assert_eval_eq("BitLength[-2]", "1", 0);   /* BitNot[-2] = 1  */
    assert_eval_eq("BitLength[-3]", "2", 0);   /* BitNot[-3] = 2  */
    assert_eval_eq("BitLength[-4]", "3", 0);   /* BitNot[-4] = 3  */
    assert_eval_eq("BitLength[-8]", "3", 0);   /* BitNot[-8] = 7  */
    assert_eval_eq("BitLength[-256]", "8", 0); /* BitNot[-256] = 255 */
}

static void test_negative_powers_of_two(void) {
    /* BitLength[-2^k] = k, since BitNot[-2^k] = 2^k - 1. */
    assert_eval_eq("BitLength[-2^10]", "10", 0);
    assert_eval_eq("BitLength[-2^30]", "30", 0);
}

/* --- int64 boundary, including INT64_MIN ----------------------------- */

static void test_int64_boundary(void) {
    assert_eval_eq("BitLength[2^62]", "63", 0);
    assert_eval_eq("BitLength[2^63]", "64", 0);      /* just past int64 max */
    assert_eval_eq("BitLength[2^63 - 1]", "63", 0);  /* int64 max          */
}

static void test_int64_min(void) {
    /* -2^63 is INT64_MIN.  BitNot[-2^63] = 2^63 - 1, which has 63 bits.  The
     * int64 kernel complements the unsigned magnitude, so this is exact where
     * a negate-based kernel would overflow. */
    assert_eval_eq("BitLength[-9223372036854775808]", "63", 0);
    assert_eval_eq("BitLength[-2^63]", "63", 0);
}

/* --- Bignums --------------------------------------------------------- */

static void test_bignums(void) {
    assert_eval_eq("BitLength[2^100]", "101", 0);
    assert_eval_eq("BitLength[2^100 - 1]", "100", 0);
    assert_eval_eq("BitLength[-2^100]", "100", 0);
    /* Documented: IntegerLength[100!, 2] == 525, so BitLength[100!] == 525. */
    assert_eval_eq("BitLength[100!]", "525", 0);
}

/* --- Listable threading ---------------------------------------------- */

static void test_listable_flat(void) {
    assert_eval_eq("BitLength[{0, 1, 2, 7, 8, 255, 256}]",
                   "{0, 1, 2, 3, 4, 8, 9}", 0);
}

static void test_listable_negatives(void) {
    assert_eval_eq("BitLength[{-1, -2, -8, -256}]", "{0, 1, 3, 8}", 0);
}

static void test_listable_nested(void) {
    assert_eval_eq("BitLength[{{1, 10}, {100, 1000}}]",
                   "{{1, 4}, {7, 10}}", 0);
}

/* --- Packed / NDArray fast paths ------------------------------------- */

static void test_packed_range(void) {
    assert_eval_eq("BitLength[Range[8]]", "{1, 2, 2, 3, 3, 3, 3, 4}", 0);
}

static void test_visible_ndarray_int64(void) {
    /* A visible int64 NDArray must produce a single-headed int64 array whose
     * elements agree with the List path, INT64_MIN included. */
    assert_eval_eq(
        "BitLength[NDArray[{-9223372036854775808, -1, 0, 1, 255}, "
        "DataType -> \"int64\"]]",
        "NDArray[{63, 0, 0, 1, 8}]", 0);
}

static void test_packed_agrees_with_list(void) {
    /* The packed producer (Range) and the plain List must give the same
     * answer -- the surfaces are not allowed to diverge. */
    assert_eval_eq("BitLength[Range[16]] == BitLength[Table[n, {n, 1, 16}]]",
                   "True", 0);
}

/* --- Compile[] round-trip -------------------------------------------- */

static void test_compile_scalar(void) {
    assert_eval_eq("Compile[{{x, _Integer}}, BitLength[x]][13]", "4", 0);
    assert_eval_eq("Compile[{{x, _Integer}}, BitLength[x]][-8]", "3", 0);
    assert_eval_eq("Compile[{{x, _Integer}}, BitLength[x]][0]", "0", 0);
}

static void test_compile_array(void) {
    assert_eval_eq(
        "Compile[{{v, _Integer, 1}}, BitLength[v]][{1, 7, 8, 255}]",
        "{1, 3, 4, 8}", 0);
}

/* --- Property tests -------------------------------------------------- */

static void test_matches_integerlength_base2(void) {
    assert_eval_eq(
        "Table[BitLength[n] == IntegerLength[n, 2], {n, 1, 64}] === "
        "Table[True, {64}]",
        "True", 0);
}

static void test_floor_log_identity(void) {
    assert_eval_eq(
        "Table[BitLength[n] == Floor[Log[2, n]] + 1, {n, 1, 40}] === "
        "Table[True, {40}]",
        "True", 0);
}

/* --- Error paths ----------------------------------------------------- */

static void test_real_int_diagnostic(void) {
    char* result = NULL;
    char* err = eval_capturing_stderr("BitLength[3.5]", &result);
    ASSERT(result != NULL);
    ASSERT_STR_EQ(result, "BitLength[3.5]");
    ASSERT(err != NULL);
    ASSERT_MSG(strstr(err, "BitLength::int") != NULL,
               "expected int diagnostic, got: %s", err);
    ASSERT(strstr(err, "Integer expected at position 1") != NULL);
    free(result);
    free(err);
}

static void test_rational_int_diagnostic(void) {
    char* result = NULL;
    char* err = eval_capturing_stderr("BitLength[3/2]", &result);
    ASSERT(result != NULL);
    ASSERT(err != NULL);
    ASSERT(strstr(err, "BitLength::int") != NULL);
    ASSERT(strstr(err, "Integer expected at position 1") != NULL);
    free(result);
    free(err);
}

static void test_complex_int_diagnostic(void) {
    char* result = NULL;
    char* err = eval_capturing_stderr("BitLength[1 + 2 I]", &result);
    ASSERT(result != NULL);
    ASSERT(err != NULL);
    ASSERT_MSG(strstr(err, "BitLength::int") != NULL,
               "expected int diagnostic, got: %s", err);
    free(result);
    free(err);
}

static void test_zero_args_argx(void) {
    char* result = NULL;
    char* err = eval_capturing_stderr("BitLength[]", &result);
    ASSERT(result != NULL);
    ASSERT_STR_EQ(result, "BitLength[]");
    ASSERT(err != NULL);
    ASSERT_MSG(strstr(err, "BitLength::argx") != NULL,
               "expected argx diagnostic, got: %s", err);
    ASSERT(strstr(err, "called with 0 arguments") != NULL);
    ASSERT(strstr(err, "1 argument is expected") != NULL);
    free(result);
    free(err);
}

static void test_two_args_argx(void) {
    char* result = NULL;
    char* err = eval_capturing_stderr("BitLength[1, 2]", &result);
    ASSERT(result != NULL);
    ASSERT_STR_EQ(result, "BitLength[1, 2]");
    ASSERT(err != NULL);
    ASSERT_MSG(strstr(err, "BitLength::argx") != NULL,
               "expected argx diagnostic, got: %s", err);
    ASSERT(strstr(err, "called with 2 arguments") != NULL);
    free(result);
    free(err);
}

static void test_symbolic_left_unevaluated(void) {
    char* result = NULL;
    char* err = eval_capturing_stderr("BitLength[x]", &result);
    ASSERT(result != NULL);
    ASSERT(strstr(result, "BitLength[x]") != NULL);
    ASSERT_MSG(err == NULL || strstr(err, "BitLength::") == NULL,
               "expected no diagnostic for symbolic input, got: %s",
               err ? err : "(null)");
    free(result);
    free(err);
}

/* --- Introspection --------------------------------------------------- */

static void test_attributes(void) {
    SymbolDef* def = symtab_get_def("BitLength");
    ASSERT(def != NULL);
    uint32_t a = get_attributes("BitLength");
    ASSERT((a & ATTR_PROTECTED) != 0);
    ASSERT((a & ATTR_LISTABLE) != 0);
}

static void test_docstring_set(void) {
    SymbolDef* def = symtab_get_def("BitLength");
    ASSERT(def != NULL);
    ASSERT(def->docstring != NULL);
    ASSERT(strstr(def->docstring, "binary bits") != NULL);
}

static void test_sym_pointer_interned(void) {
    ASSERT(SYM_BitLength != NULL);
    ASSERT(strcmp(SYM_BitLength, "BitLength") == 0);
}

/* --- Memory-safety stress loop --------------------------------------- */

static void test_repeated_evaluation_does_not_corrupt(void) {
    /* Mix machine-int / bignum / negative / INT64_MIN / listable / packed /
     * error paths so any mishandled mpz_t or args[] slot surfaces under
     * valgrind. */
    for (int k = 0; k < 50; k++) {
        assert_eval_eq("BitLength[255]", "8", 0);
        assert_eval_eq("BitLength[-8]", "3", 0);
        assert_eval_eq("BitLength[0]", "0", 0);
        assert_eval_eq("BitLength[-2^63]", "63", 0);
        assert_eval_eq("BitLength[2^100]", "101", 0);
        assert_eval_eq("BitLength[{1, 10, 100, 1000}]", "{1, 4, 7, 10}", 0);
        assert_eval_eq("BitLength[Range[8]]", "{1, 2, 2, 3, 3, 3, 3, 4}", 0);
    }
}

int main(void) {
    symtab_init();
    core_init();

    TEST(test_zero);
    TEST(test_small_positives);
    TEST(test_byte_boundaries);
    TEST(test_powers_of_two);

    TEST(test_negatives);
    TEST(test_negative_powers_of_two);

    TEST(test_int64_boundary);
    TEST(test_int64_min);

    TEST(test_bignums);

    TEST(test_listable_flat);
    TEST(test_listable_negatives);
    TEST(test_listable_nested);

    TEST(test_packed_range);
    TEST(test_visible_ndarray_int64);
    TEST(test_packed_agrees_with_list);

    TEST(test_compile_scalar);
    TEST(test_compile_array);

    TEST(test_matches_integerlength_base2);
    TEST(test_floor_log_identity);

    TEST(test_real_int_diagnostic);
    TEST(test_rational_int_diagnostic);
    TEST(test_complex_int_diagnostic);
    TEST(test_zero_args_argx);
    TEST(test_two_args_argx);
    TEST(test_symbolic_left_unevaluated);

    TEST(test_attributes);
    TEST(test_docstring_set);
    TEST(test_sym_pointer_interned);

    TEST(test_repeated_evaluation_does_not_corrupt);

    printf("All BitLength tests passed!\n");
    return 0;
}
