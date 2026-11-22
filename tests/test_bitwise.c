/* Unit tests for the bit-level integer builtins (src/bitwise/).
 *
 *   BitLength[n] -- number of binary bits needed to represent the integer n.
 *
 * BitXor[n1, n2, ...] -- bitwise XOR (two's complement), Flat/Orderless/
 * OneIdentity/Listable: documented examples, bignums across the int64 boundary,
 * pair cancellation / integer folding / BitNot extraction, properties against
 * an IntegerDigits oracle, packed + NDArray (binary and n-ary), cross-surface
 * bit-equality incl. Compile, CompileDiagnostics, ::int + Quiet/Check.
 *
 * BitAnd[n1, n2, ...] -- bitwise AND (two's complement), same attributes:
 * documented examples, bignums, idempotence / integer folding / 0 absorbs /
 * -1 dropped, silent non-integers (no message), a
 * two's-complement oracle a + b == BitXor[a, b] + 2 BitAnd[a, b] across the
 * bignum boundary, packed + NDArray (binary and n-ary), cross-surface
 * bit-equality incl. Compile, CompileDiagnostics.
 *
 * BitOr[n1, n2, ...] -- bitwise OR, the dual of BitAnd: documented examples,
 * bignums (incl. a negative int64 shrinking a bignum result back to int64),
 * idempotence / integer folding / -1 absorbs / 0 dropped, silent non-integers,
 * oracles a + b == BitOr + BitAnd, BitOr == BitXor + BitAnd and De Morgan
 * across the bignum boundary, packed + NDArray (binary and n-ary),
 * cross-surface bit-equality incl. Compile, CompileDiagnostics.
 *
 * BitLength coverage:
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
#include "compile/autocompile.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>

/* Capture stderr while `input` is parsed + evaluated.  Returns the collected
 * stderr text as a heap string (caller frees) and writes the printed result
 * into *out_result_str (also heap-allocated).  Fixed temp path; safe because
 * tests run serially.  The original stderr fd is saved and restored with
 * dup/dup2 (as test_digit_sum.c does) rather than reopened as /dev/tty, which
 * does not exist in a headless run and left stderr closed for the next
 * capture. */
static char* eval_capturing_stderr(const char* input, char** out_result_str) {
    const char* path = "/tmp/mathilda_bitwise_stderr.log";
    fflush(stderr);
    int saved_stderr = dup(STDERR_FILENO);
    if (!freopen(path, "w+", stderr)) {
        if (saved_stderr != -1) close(saved_stderr);
        if (out_result_str) *out_result_str = NULL;
        return NULL;
    }

    Expr* p = parse_expression(input);
    Expr* e = evaluate(p);
    if (out_result_str) *out_result_str = expr_to_string(e);
    expr_free(p);
    expr_free(e);

    fflush(stderr);
    if (saved_stderr != -1) {
        dup2(saved_stderr, STDERR_FILENO);
        close(saved_stderr);
    }

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
    assert_eval_eq("BitLength[-4]", "2", 0);   /* BitNot[-4] = 3, BitLength[3] = 2 (WL: 2) */
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


/* ===================================================================== *
 *  BitXor[n1, n2, ...]                                                   *
 * ===================================================================== */

/* --- Documented examples --------------------------------------------- */

static void test_bitxor_basic(void) {
    assert_eval_eq("BitXor[61, 15]", "50", 0);
    assert_eval_eq("BitXor[3333, 5555, 7777, 9999]", "8664", 0);
    assert_eval_eq("BitXor[-2, -3]", "3", 0);
    assert_eval_eq("BitXor[-6, 3]", "-7", 0);
    assert_eval_eq("BitXor[]", "0", 0);
    assert_eval_eq("BitXor[5]", "5", 0);
    assert_eval_eq("BitXor[0, 0]", "0", 0);
    assert_eval_eq("BitXor[7, 7]", "0", 0);
    assert_eval_eq("BitXor[-1, -1]", "0", 0);
}

static void test_bitxor_range_threading(void) {
    assert_eval_eq("BitXor[Range[120], 5]",
        "{4, 7, 6, 1, 0, 3, 2, 13, 12, 15, 14, 9, 8, 11, 10, 21, 20, 23, 22, "
        "17, 16, 19, 18, 29, 28, 31, 30, 25, 24, 27, 26, 37, 36, 39, 38, 33, "
        "32, 35, 34, 45, 44, 47, 46, 41, 40, 43, 42, 53, 52, 55, 54, 49, 48, "
        "51, 50, 61, 60, 63, 62, 57, 56, 59, 58, 69, 68, 71, 70, 65, 64, 67, "
        "66, 77, 76, 79, 78, 73, 72, 75, 74, 85, 84, 87, 86, 81, 80, 83, 82, "
        "93, 92, 95, 94, 89, 88, 91, 90, 101, 100, 103, 102, 97, 96, 99, 98, "
        "109, 108, 111, 110, 105, 104, 107, 106, 117, 116, 119, 118, 113, "
        "112, 115, 114, 125}", 0);
}

static void test_bitxor_truth_table(void) {
    assert_eval_eq("Outer[BitXor, {1, 0}, {1, 0}]", "{{0, 1}, {1, 0}}", 0);
}

static void test_bitxor_empty_lists(void) {
    assert_eval_eq("BitXor[{}]", "{}", 0);
    assert_eval_eq("BitXor[{}, {}]", "{}", 0);
}

static void test_bitxor_listable_shapes(void) {
    assert_eval_eq("BitXor[{1, 2}, 3, {4, 5}]", "{6, 4}", 0);
    assert_eval_eq("BitXor[{{1, 2}, {3, 4}}, 1]", "{{0, 3}, {2, 5}}", 0);
    assert_eval_eq("BitXor[{1, 2}, {3, 4}]", "{2, 6}", 0);
}

/* --- Bignums --------------------------------------------------------- */

static void test_bitxor_bignums(void) {
    assert_eval_eq("BitXor[3^100, 5^100]",
        "7888609052210118054116793203790437104101594041461795136781437785275168", 0);
    /* Two's complement crosses the int64 boundary both ways. */
    assert_eval_eq("BitXor[2^63, -1]", "-9223372036854775809", 0);
    assert_eval_eq("BitXor[-2^63, 1]", "-9223372036854775807", 0);
    assert_eval_eq("BitXor[-2^63, -1]", "9223372036854775807", 0);
    assert_eval_eq("BitXor[2^63 - 1, -2^63]", "-1", 0);
    /* A bignum result that fits normalises back to a machine Integer. */
    assert_eval_eq("BitXor[2^70, 2^70 + 5]", "5", 0);
    assert_eval_eq("MachineIntegerQ[BitXor[2^70, 2^70 + 5]]", "True", 0);
    assert_eval_eq("BitXor[2^100, -2^100]", "-2535301200456458802993406410752", 0);
    assert_eval_eq("BitXor[-3^50, 7, -5^30]",
        "718817681789702120559719", 0);
    assert_eval_eq("BitXor[{2^64, 1}, 2^64]", "{0, 18446744073709551617}", 0);
}

static void test_bitxor_bignum_operands_untouched(void) {
    /* The bignum path XORs in place into a private result node: the operands
     * themselves -- bound, shared, or repeated -- must come out unchanged. */
    assert_eval_eq("p = 3^200; q = 5^150; r = BitXor[p, q, p, 7]; "
                   "{p == 3^200, q == 5^150, r == BitXor[q, 7]}",
                   "{True, True, True}", 0);
    assert_eval_eq("p = 2^100; {BitXor[p, p], p == 2^100}", "{0, True}", 0);
    assert_eval_eq("l = {2^90, 2^91}; BitXor[l, 2^90]; l == {2^90, 2^91}", "True", 0);
    /* machine constants folded into a bignum once, mixed signs */
    assert_eval_eq("BitXor[2^100, 1, 2, 4, -8]", "-1267650600228229401496703205377", 0);
}

/* --- Symbolic simplification ----------------------------------------- */

static void test_bitxor_pair_cancellation(void) {
    assert_eval_eq("BitXor[x, x]", "0", 0);
    assert_eval_eq("BitXor[x, x, x]", "x", 0);
    assert_eval_eq("BitXor[x, x, x, x]", "0", 0);
    assert_eval_eq("BitXor[x, x, x, x, x]", "x", 0);
    assert_eval_eq("BitXor[x, y, x, y]", "0", 0);
    assert_eval_eq("BitXor[x, y, x, y, x]", "x", 0);
    assert_eval_eq("BitXor[x, y, z, x]", "BitXor[y, z]", 0);
    assert_eval_eq("BitXor[f[a], f[a], 3]", "3", 0);
    assert_eval_eq("BitXor[a, 3, a]", "3", 0);
    /* A cancelled non-integer pair is gone, so it is never an error. */
    assert_eval_eq("BitXor[x, 2.5, 2.5]", "x", 0);
    assert_eval_eq("BitXor[2.5, 2.5]", "0", 0);
}

static void test_bitxor_integer_folding(void) {
    assert_eval_eq("BitXor[1, 2, x]", "BitXor[3, x]", 0);
    assert_eval_eq("BitXor[0, x]", "x", 0);
    assert_eval_eq("BitXor[x, 0, 0]", "x", 0);
    assert_eval_eq("BitXor[2, x, 2]", "x", 0);
    assert_eval_eq("BitXor[1, 2, x, x, y]", "BitXor[3, y]", 0);
    assert_eval_eq("BitXor[-5, x]", "BitXor[-5, x]", 0);
    assert_eval_eq("BitXor[x, 2^70, 2^70 + 1]", "BitXor[1, x]", 0);
}

static void test_bitxor_bitnot(void) {
    assert_eval_eq("BitXor[-1, x]", "BitNot[x]", 0);
    assert_eval_eq("BitXor[x, y, -1]", "BitNot[BitXor[x, y]]", 0);
    /* -1 is pulled out before the integers combine ... */
    assert_eval_eq("BitXor[1, -1, x]", "BitNot[BitXor[1, x]]", 0);
    assert_eval_eq("BitXor[3, -1, x]", "BitNot[BitXor[3, x]]", 0);
    assert_eval_eq("BitXor[-1, x, y, y, 3]", "BitNot[BitXor[3, x]]", 0);
    assert_eval_eq("BitXor[-1, 2^70, x]",
                   "BitNot[BitXor[1180591620717411303424, x]]", 0);
    /* ... and a combined constant of -1 is one too. */
    assert_eval_eq("BitXor[-2, 1, x]", "BitNot[x]", 0);
    assert_eval_eq("BitXor[-1, 0, x]", "BitNot[x]", 0);
    /* An even count cancels; everything else gone leaves the integer. */
    assert_eval_eq("BitXor[x, -1, -1]", "x", 0);
    assert_eval_eq("BitXor[-1, -1, -1, x]", "BitNot[x]", 0);
    assert_eval_eq("BitXor[-1, x, x]", "-1", 0);
    assert_eval_eq("BitXor[x, -1, y, y]", "BitNot[x]", 0);
    assert_eval_eq("BitXor[x, BitNot[x]]", "BitXor[x, BitNot[x]]", 0);
}

static void test_bitxor_orderless_flat(void) {
    assert_eval_eq("BitXor[y, x]", "BitXor[x, y]", 0);
    assert_eval_eq("BitXor[x]", "x", 0);
    assert_eval_eq("BitXor[BitXor[x, 1], BitXor[2, y]]", "BitXor[3, x, y]", 0);
    assert_eval_eq("BitXor[BitXor[x, y], x]", "y", 0);
}

/* --- Properties ------------------------------------------------------ */

static void test_bitxor_properties(void) {
    /* Self-inverse and associativity over a grid that includes negatives. */
    assert_eval_eq(
        "And @@ Flatten[Table[BitXor[BitXor[a, b], b] == a, "
        "{a, -20, 20}, {b, -20, 20}]]", "True", 0);
    assert_eval_eq(
        "And @@ Flatten[Table[BitXor[BitXor[a, b], c] == BitXor[a, BitXor[b, c]], "
        "{a, -6, 6}, {b, -6, 6}, {c, -6, 6}]]", "True", 0);
    /* Independent oracle: bitwise parity from the binary digits. */
    assert_eval_eq(
        "And @@ Flatten[Table[BitXor[a, b] == "
        "FromDigits[Mod[IntegerDigits[a, 2, 8] + IntegerDigits[b, 2, 8], 2], 2], "
        "{a, 0, 63}, {b, 0, 63}]]", "True", 0);
    /* Two's complement: BitXor[n, -1] == -n - 1 (BitNot). */
    assert_eval_eq(
        "And @@ Table[BitXor[n, -1] == -n - 1, {n, {-2^80, -7, 0, 9, 2^63, 3^60}}]",
        "True", 0);
}

/* --- Packed / NDArray surfaces --------------------------------------- */

static void test_bitxor_packed_and_ndarray(void) {
    assert_eval_eq("PackedArrayQ[BitXor[Range[10], 3]]", "True", 0);
    assert_eval_eq("BitXor[NDArray[{1, 2, 3}, DataType -> \"int64\"], 6]",
                   "NDArray[{7, 4, 5}]", 0);
    assert_eval_eq(
        "BitXor[NDArray[{-9223372036854775808, -1, 0}, DataType -> \"int64\"], "
        "NDArray[{1, -1, -1}, DataType -> \"int64\"]]",
        "NDArray[{-9223372036854775807, 0, -1}]", 0);
    /* n-ary: folded on the buffer, scalars combined first. */
    assert_eval_eq("BitXor[Range[10], 1, 2]",
                   "{2, 1, 0, 7, 6, 5, 4, 11, 10, 9}", 0);
    assert_eval_eq("PackedArrayQ[BitXor[Range[10], 1, 2]]", "True", 0);
    assert_eval_eq("BitXor[Range[5], Range[5], Reverse[Range[5]]]",
                   "{5, 4, 3, 2, 1}", 0);
    assert_eval_eq("BitXor[NDArray[{1, 2, 3}, DataType -> \"int64\"], 1, 2]",
                   "NDArray[{2, 1, 0}]", 0);
    assert_eval_eq("BitXor[Range[4], 1, 1]", "{1, 2, 3, 4}", 0);
    /* Off-buffer operands take the exact List path. */
    assert_eval_eq("BitXor[Range[3], 2^70]",
        "{1180591620717411303425, 1180591620717411303426, 1180591620717411303427}", 0);
    assert_eval_eq("BitXor[Range[3], x, x]", "{1, 2, 3}", 0);
    assert_eval_eq("BitXor[Range[3], -1]", "{-2, -3, -4}", 0);
}

static void test_bitxor_surfaces_agree(void) {
    assert_eval_eq(
        "SeedRandom[3]; a = Join[RandomInteger[{-2^62, 2^62}, 300], {-2^63, 2^63 - 1, 0, -1}]; "
        "b = Join[RandomInteger[{-2^62, 2^62}, 300], {-1, -2^63, 5, -2^63}]; "
        "ref = Table[BitXor[a[[i]], b[[i]]], {i, Length[a]}]; "
        "{BitXor[a, b] === ref, "
        " Normal[BitXor[NDArray[a, DataType -> \"int64\"], NDArray[b, DataType -> \"int64\"]]] === ref, "
        " Compile[{{v, _Integer, 1}, {w, _Integer, 1}}, BitXor[v, w]][a, b] === ref, "
        " Table[Compile[{{x, _Integer}, {y, _Integer}}, BitXor[x, y]][a[[i]], b[[i]]], {i, Length[a]}] === ref, "
        " BitXor[a, b, 99] === Table[BitXor[a[[i]], b[[i]], 99], {i, Length[a]}]}",
        "{True, True, True, True, True}", 0);
}

/* --- Compile[] ------------------------------------------------------- */

static void test_bitxor_compile(void) {
    assert_eval_eq("Compile[{{a, _Integer}, {b, _Integer}}, BitXor[a, b]][61, 15]", "50", 0);
    assert_eval_eq("Compile[{{a, _Integer}, {b, _Integer}}, BitXor[a, b]][-2, -3]", "3", 0);
    assert_eval_eq(
        "Compile[{{a, _Integer}, {b, _Integer}, {c, _Integer}, {d, _Integer}}, "
        "BitXor[a, b, c, d]][3333, 5555, 7777, 9999]", "8664", 0);
    assert_eval_eq("Compile[{{v, _Integer, 1}}, BitXor[v, 5]][Range[8]]",
                   "{4, 7, 6, 1, 0, 3, 2, 13}", 0);
    assert_eval_eq("Compile[{{v, _Integer, 1}}, BitXor[6, v, 4]][Range[4]]",
                   "{3, 0, 1, 6}", 0);
    assert_eval_eq(
        "Compile[{{v, _Integer, 1}, {w, _Integer, 1}}, BitXor[v, w]][{1, 2, 3}, {3, 2, 1}]",
        "{2, 0, 2}", 0);
}

static void test_bitxor_compile_diagnostics(void) {
    assert_eval_eq("\"Compiled\" /. CompileDiagnostics[{{a, _Integer}, {b, _Integer}}, BitXor[a, b]]", "True", 0);
    assert_eval_eq("\"Compiled\" /. CompileDiagnostics[{{a, _Integer}, {b, _Integer}, {c, _Integer}}, BitXor[a, b, c]]", "True", 0);
    assert_eval_eq("\"Compiled\" /. CompileDiagnostics[{{v, _Integer, 1}}, BitXor[v, 5]]", "True", 0);
    assert_eval_eq("\"Compiled\" /. CompileDiagnostics[{{v, _Integer, 1}, {w, _Integer, 1}}, BitXor[v, w]]", "True", 0);
    assert_eval_eq("\"Compiled\" /. CompileDiagnostics[{{v, _Integer, 1}}, BitXor[v, 6, 4]]", "True", 0);
}

/* --- Error paths ----------------------------------------------------- */

static void test_bitxor_real_int_diagnostic(void) {
    char* result = NULL;
    char* err = eval_capturing_stderr("BitXor[2.5, 3]", &result);
    ASSERT(result != NULL);
    ASSERT_STR_EQ(result, "BitXor[2.5, 3]");
    ASSERT(err != NULL);
    ASSERT_MSG(strstr(err, "BitXor::int") != NULL,
               "expected int diagnostic, got: %s", err);
    ASSERT(strstr(err, "Integer expected at position 1 in BitXor[2.5, 3]") != NULL);
    free(result);
    free(err);
}

static void test_bitxor_rational_complex_diagnostic(void) {
    const char* in[] = { "BitXor[1/2, x]", "BitXor[I, 3]", "BitXor[-1, x, 2.5]" };
    for (size_t i = 0; i < sizeof in / sizeof in[0]; i++) {
        char* result = NULL;
        char* err = eval_capturing_stderr(in[i], &result);
        ASSERT(result != NULL);
        ASSERT(strstr(result, "BitXor[") != NULL);    /* left unevaluated */
        ASSERT_MSG(err && strstr(err, "BitXor::int") != NULL,
                   "expected int diagnostic for %s, got: %s", in[i], err ? err : "(null)");
        free(result);
        free(err);
    }
}

static void test_bitxor_quiet_check(void) {
    assert_eval_eq("Check[BitXor[2.5, 3], failed]", "failed", 0);
    assert_eval_eq("Check[BitXor[x, x, 7], failed]", "7", 0);
    char* result = NULL;
    char* err = eval_capturing_stderr("Quiet[BitXor[2.5, 3]]", &result);
    ASSERT(result != NULL);
    ASSERT_STR_EQ(result, "BitXor[2.5, 3]");
    ASSERT_MSG(err == NULL || strstr(err, "BitXor::") == NULL,
               "Quiet leaked: %s", err ? err : "(null)");
    free(result);
    free(err);
}

static void test_bitxor_symbolic_silent(void) {
    char* result = NULL;
    char* err = eval_capturing_stderr("BitXor[x, y]", &result);
    ASSERT(result != NULL);
    ASSERT_STR_EQ(result, "BitXor[x, y]");
    ASSERT_MSG(err == NULL || strstr(err, "BitXor::") == NULL,
               "expected no diagnostic for symbolic input, got: %s", err ? err : "(null)");
    free(result);
    free(err);
}

/* --- Introspection --------------------------------------------------- */

static void test_bitxor_attributes(void) {
    assert_eval_eq("Attributes[BitXor]",
                   "{Flat, Listable, OneIdentity, Orderless, Protected}", 0);
}

static void test_bitxor_docstring_and_sym(void) {
    SymbolDef* def = symtab_get_def("BitXor");
    ASSERT(def != NULL);
    ASSERT(def->docstring != NULL);
    ASSERT(strstr(def->docstring, "bitwise XOR") != NULL);
    ASSERT(SYM_BitXor != NULL && strcmp(SYM_BitXor, "BitXor") == 0);
    ASSERT(SYM_BitNot != NULL && strcmp(SYM_BitNot, "BitNot") == 0);
}

static void test_bitxor_stress(void) {
    for (int k = 0; k < 50; k++) {
        assert_eval_eq("BitXor[61, 15]", "50", 0);
        assert_eval_eq("BitXor[3^100, 5^100, -7]", "-7888609052210118054116793203790437104101594041461795136781437785275175", 0);
        assert_eval_eq("BitXor[1, -1, x]", "BitNot[BitXor[1, x]]", 0);
        assert_eval_eq("BitXor[x, y, x, 2^70, 3]", "BitXor[1180591620717411303427, y]", 0);
        assert_eval_eq("BitXor[Range[6], 1, 2]", "{2, 1, 0, 7, 6, 5}", 0);
        assert_eval_eq("BitXor[Range[3], 2^70]",
            "{1180591620717411303425, 1180591620717411303426, 1180591620717411303427}", 0);
    }
}

/* ===================================================================== *
 *  BitAnd[n1, n2, ...]                                                   *
 * ===================================================================== */

/* --- Documented examples -------------------------------------------- */

static void test_bitand_documented_examples(void) {
    assert_eval_eq("BitAnd[61, 15]", "13", 0);
    assert_eval_eq("BitAnd[Range[70, 130], 17]",
        "{0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 16, 17, 16, 17, 16, 17, 16, 17, 16, 17, "
        "16, 17, 16, 17, 16, 17, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, "
        "16, 17, 16, 17, 16, 17, 16, 17, 16, 17, 16, 17, 16, 17, 16, 17, 0, 1, 0}", 0);
    assert_eval_eq("BitAnd[3^100, 5^100]",
        "503913279078601980753385379100266097026814443729", 0);
    assert_eval_eq("BitAnd[3333, 5555, 7777, 9999]", "1025", 0);
    assert_eval_eq("BitAnd[-2, -3]", "-4", 0);
    assert_eval_eq("BitAnd[x, y, y, x]", "BitAnd[x, y]", 0);
    assert_eval_eq("BitAnd[3, 5, x]", "BitAnd[1, x]", 0);
    assert_eval_eq("BitAnd[1, 2, x]", "0", 0);
    assert_eval_eq("Table[BitAnd[i, 2], {i, 10}]", "{0, 2, 2, 0, 0, 2, 2, 0, 0, 2}", 0);
    assert_eval_eq("Table[Sign[BitAnd[i, 2^2 + 2^3]], {i, 16}]",
                   "{0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0}", 0);
    assert_eval_eq("BitAnd[y, x]", "BitAnd[x, y]", 0);
    assert_eval_eq("BitAnd[]", "-1", 0);
    assert_eval_eq("BitAnd[{}]", "{}", 0);
    assert_eval_eq("BitAnd[{}, {}]", "{}", 0);
}

static void test_bitand_basic(void) {
    assert_eval_eq("BitAnd[-6, 3]", "2", 0);
    assert_eval_eq("BitAnd[0, 0]", "0", 0);
    assert_eval_eq("BitAnd[7, 7]", "7", 0);
    assert_eval_eq("BitAnd[5]", "5", 0);
    assert_eval_eq("BitAnd[-1, -1]", "-1", 0);
    assert_eval_eq("BitAnd[-7, -7]", "-7", 0);
    assert_eval_eq("BitAnd[-1]", "-1", 0);
    assert_eval_eq("Outer[BitAnd, {1, 0}, {1, 0}]", "{{1, 0}, {0, 0}}", 0);
}

static void test_bitand_listable_shapes(void) {
    assert_eval_eq("BitAnd[Range[16], 5]",
                   "{1, 0, 1, 4, 5, 4, 5, 0, 1, 0, 1, 4, 5, 4, 5, 0}", 0);
    assert_eval_eq("BitAnd[{1, 2}, 3, {4, 5}]", "{0, 0}", 0);
    assert_eval_eq("BitAnd[{{1, 2}, {3, 4}}, 1]", "{{1, 0}, {1, 0}}", 0);
    assert_eval_eq("BitAnd[{{1, 2}, {3, 4}}, 5]", "{{1, 0}, {1, 4}}", 0);
    assert_eval_eq("BitAnd[{1, 2}, {3, 4}]", "{1, 0}", 0);
    assert_eval_eq("BitAnd[{1, 2}, {3, 4}, {5, 6}]", "{1, 0}", 0);
    assert_eval_eq("BitAnd[0, {1, 2}]", "{0, 0}", 0);
    assert_eval_eq("BitAnd[{1, 2}, x]", "{BitAnd[1, x], BitAnd[2, x]}", 0);
    assert_eval_eq("BitAnd[x, {1, 2}, {3, 4}]", "{BitAnd[1, x], 0}", 0);
}

/* --- Bignums and the int64 boundary ---------------------------------- */

static void test_bitand_bignums(void) {
    assert_eval_eq("BitAnd[2^63, -1]", "9223372036854775808", 0);
    assert_eval_eq("BitAnd[-2^63, -1]", "-9223372036854775808", 0);
    assert_eval_eq("BitAnd[2^63 - 1, -2^63]", "0", 0);
    assert_eval_eq("BitAnd[-2^63, -2^63 + 1]", "-9223372036854775808", 0);
    assert_eval_eq("BitAnd[2^70, -1]", "1180591620717411303424", 0);
    assert_eval_eq("BitAnd[2^100, -2^100]", "1267650600228229401496703205376", 0);
    assert_eval_eq("BitAnd[-3^50, 7, -5^30]", "7", 0);
    assert_eval_eq("BitAnd[-3^50, 2^80 - 1]", "491027831922776585935927", 0);
    assert_eval_eq("BitAnd[-2^100 - 1, 2^101 - 1]",
                   "1267650600228229401496703205375", 0);
    assert_eval_eq("BitAnd[{2^64 + 3, 1}, 2^64 + 1]", "{18446744073709551617, 1}", 0);
    /* A bignum result that fits normalises back to a machine integer. */
    assert_eval_eq("BitAnd[-3^50, 7, -5^30] // MachineIntegerQ", "True", 0);
    assert_eval_eq("MachineIntegerQ[BitAnd[2^70 + 5, 7]]", "True", 0);
    /* Machine constants folded into a bignum once, mixed signs. */
    assert_eval_eq("BitAnd[2^100 - 1, 2, 4, -8]", "0", 0);
    assert_eval_eq("BitAnd[3^100, 5^100, -7]",
                   "503913279078601980753385379100266097026814443729", 0);
}

static void test_bitand_bignum_operands_untouched(void) {
    /* The bignum path ANDs in place into a private result node, set from the
     * first BigInt: the operands -- bound, shared or repeated -- are unchanged. */
    assert_eval_eq("p = 3^200; q = 5^150; r = BitAnd[p, q, p, -1]; "
                   "{p == 3^200, q == 5^150, r == BitAnd[3^200, 5^150]}",
                   "{True, True, True}", 0);
    assert_eval_eq("p = 2^100 + 7; {BitAnd[p, p], p == 2^100 + 7}",
                   "{1267650600228229401496703205383, True}", 0);
    assert_eval_eq("l = {2^90 + 1, 2^91 + 3}; BitAnd[l, 2^90 + 3]; l == {2^90 + 1, 2^91 + 3}",
                   "True", 0);
}

/* --- Symbolic simplification ----------------------------------------- */

static void test_bitand_idempotence(void) {
    assert_eval_eq("BitAnd[x, x]", "x", 0);
    assert_eval_eq("BitAnd[x, x, x]", "x", 0);
    assert_eval_eq("BitAnd[x, y, z, x]", "BitAnd[x, y, z]", 0);
    assert_eval_eq("BitAnd[f[u], f[u], 3]", "BitAnd[3, f[u]]", 0);
    assert_eval_eq("BitAnd[u, 3, u]", "BitAnd[3, u]", 0);
    assert_eval_eq("BitAnd[7, x, x, y]", "BitAnd[7, x, y]", 0);
    assert_eval_eq("BitAnd[2.5, 2.5]", "2.5", 0);
    assert_eval_eq("BitAnd[x, 2.5, 2.5]", "BitAnd[2.5, x]", 0);
    assert_eval_eq("BitAnd[I, I, x]", "BitAnd[I, x]", 0);
}

static void test_bitand_integer_folding(void) {
    assert_eval_eq("BitAnd[x, 2^70 + 3, 2^70 + 1]", "BitAnd[1180591620717411303425, x]", 0);
    assert_eval_eq("BitAnd[-5, x]", "BitAnd[-5, x]", 0);
    assert_eval_eq("BitAnd[-3^50, -5^40, x, x]",
                   "BitAnd[-9095627113487849460680488937, x]", 0);
    assert_eval_eq("BitAnd[2^70, 2^70 + 5, x]", "BitAnd[1180591620717411303424, x]", 0);
    assert_eval_eq("BitAnd[x, y, x, 2^70 + 3, 7]", "BitAnd[3, x, y]", 0);
    assert_eval_eq("BitAnd[3, 5, 2.5]", "BitAnd[1, 2.5]", 0);
    assert_eval_eq("BitAnd[3, 5, 2.5, x]", "BitAnd[1, 2.5, x]", 0);
}

static void test_bitand_zero_absorbs(void) {
    assert_eval_eq("BitAnd[0, x]", "0", 0);
    assert_eval_eq("BitAnd[1, 2, x, x, y]", "0", 0);
    assert_eval_eq("BitAnd[x, 2^70, 2^71]", "0", 0);
    assert_eval_eq("BitAnd[0, 2.5]", "0", 0);
    assert_eval_eq("BitAnd[0, 1/2]", "0", 0);
    assert_eval_eq("BitAnd[0, x, 1/2]", "0", 0);
    assert_eval_eq("BitAnd[0, x, y, 2.5, 1/2]", "0", 0);
    assert_eval_eq("BitAnd[1, 2, -1, x]", "0", 0);
}

static void test_bitand_minus_one_identity(void) {
    assert_eval_eq("BitAnd[-1, x]", "x", 0);
    assert_eval_eq("BitAnd[-1, -1, x]", "x", 0);
    assert_eval_eq("BitAnd[x, -1, -1]", "x", 0);
    assert_eval_eq("BitAnd[-1, x, y]", "BitAnd[x, y]", 0);
    assert_eval_eq("BitAnd[x, y, -1, y]", "BitAnd[x, y]", 0);
    assert_eval_eq("BitAnd[-1, 2.5]", "2.5", 0);
    assert_eval_eq("BitAnd[2^70, x, -1]", "BitAnd[1180591620717411303424, x]", 0);
    assert_eval_eq("BitAnd[x, -1, 2.5, 2.5]", "BitAnd[2.5, x]", 0);
}

static void test_bitand_no_extra_rules(void) {
    /* Neither complement nor absorption rules apply. */
    assert_eval_eq("BitAnd[x, BitNot[x]]", "BitAnd[x, BitNot[x]]", 0);
    assert_eval_eq("BitAnd[3, x, BitNot[x]]", "BitAnd[3, x, BitNot[x]]", 0);
    assert_eval_eq("BitAnd[BitNot[x], BitNot[y]]", "BitAnd[BitNot[x], BitNot[y]]", 0);
}

static void test_bitand_orderless_flat(void) {
    assert_eval_eq("BitAnd[x]", "x", 0);
    assert_eval_eq("BitAnd[BitAnd[x, 7], BitAnd[13, y]]", "BitAnd[5, x, y]", 0);
    assert_eval_eq("BitAnd[BitAnd[x, y], x]", "BitAnd[x, y]", 0);
}

/* --- Non-integers: left alone, silently ------------------------------ */

static void test_bitand_non_integers_silent(void) {
    /* (BitAnd[I, 3] is also silent; not pinned here because the canonical
     * order of I against an integer is a known core-comparator gap,
     * independent of BitAnd: it should sort I first.) */
    const char* in[]  = { "BitAnd[2.5, 3]", "BitAnd[1/2, x]", "BitAnd[I, x]",
                          "BitAnd[1.0, 3]", "BitAnd[x, y]" };
    const char* out[] = { "BitAnd[2.5, 3]", "BitAnd[1/2, x]", "BitAnd[I, x]",
                          "BitAnd[1.0, 3]", "BitAnd[x, y]" };
    for (size_t i = 0; i < sizeof in / sizeof in[0]; i++) {
        char* result = NULL;
        char* err = eval_capturing_stderr(in[i], &result);
        ASSERT(result != NULL);
        ASSERT_STR_EQ(result, out[i]);
        ASSERT_MSG(err == NULL || strstr(err, "BitAnd::") == NULL,
                   "unexpected diagnostic for %s: %s", in[i], err ? err : "(null)");
        free(result);
        free(err);
    }
    /* No message fires, so Check takes the value branch. */
    assert_eval_eq("Check[BitAnd[2.5, 3], failed]", "BitAnd[2.5, 3]", 0);
}

static void test_bitand_ragged_tdlen(void) {
    /* Ragged lists cannot thread: the evaluator's Listable step reports
     * Thread::tdlen and the call stays unevaluated (also through the n-ary
     * buffer fold, which hands mismatched shapes back to the List path). */
    assert_eval_eq("BitAnd[{1, 2}, {1, 2, 3}]", "BitAnd[{1, 2}, {1, 2, 3}]", 0);
    assert_eval_eq("BitAnd[Range[2], Range[3], 1]", "BitAnd[1, {1, 2}, {1, 2, 3}]", 0);
}

/* --- Properties ------------------------------------------------------ */

static void test_bitand_properties(void) {
    assert_eval_eq(
        "And @@ Flatten[Table[BitAnd[BitAnd[a, b], c] == BitAnd[a, BitAnd[b, c]], "
        "{a, -6, 6}, {b, -6, 6}, {c, -6, 6}]]", "True", 0);
    /* Independent oracle: bitwise product of the binary digits. */
    assert_eval_eq(
        "And @@ Flatten[Table[BitAnd[a, b] == "
        "FromDigits[IntegerDigits[a, 2, 8] * IntegerDigits[b, 2, 8], 2], "
        "{a, 0, 63}, {b, 0, 63}]]", "True", 0);
    /* a + b == BitXor[a, b] + 2 BitAnd[a, b] holds in two's complement for
     * every sign and size -- an oracle across the bignum boundary. */
    assert_eval_eq(
        "SeedRandom[7]; ns = Join[RandomInteger[{-2^200, 2^200}, 30], "
        "{-2^63, 2^63 - 1, 2^63, -2^63 - 1, 0, -1, 1, 3^100, -5^90}]; "
        "And @@ Flatten[Table[m + n == BitXor[m, n] + 2 BitAnd[m, n], {m, ns}, {n, ns}]]",
        "True", 0);
    /* BitAnd[n, -1] == n; BitAnd[n, 0] == 0. */
    assert_eval_eq(
        "And @@ Table[BitAnd[n, -1] == n && BitAnd[n, 0] == 0, "
        "{n, {-2^80, -7, 0, 9, 2^63, 3^60}}]", "True", 0);
}

/* --- Packed / NDArray surfaces --------------------------------------- */

static void test_bitand_packed_and_ndarray(void) {
    assert_eval_eq("PackedArrayQ[BitAnd[Range[10], 6]]", "True", 0);
    assert_eval_eq("BitAnd[NDArray[{1, 2, 3}, DataType -> \"int64\"], 6]",
                   "NDArray[{0, 2, 2}]", 0);
    assert_eval_eq(
        "BitAnd[NDArray[{-9223372036854775808, -1, 0}, DataType -> \"int64\"], "
        "NDArray[{-1, -9223372036854775808, -1}, DataType -> \"int64\"]]",
        "NDArray[{-9223372036854775808, -9223372036854775808, 0}]", 0);
    /* n-ary: folded on the buffer, scalars combined first. */
    assert_eval_eq("BitAnd[Range[10], 6, 3]", "{0, 2, 2, 0, 0, 2, 2, 0, 0, 2}", 0);
    assert_eval_eq("PackedArrayQ[BitAnd[Range[10], 6, 3]]", "True", 0);
    assert_eval_eq("BitAnd[Range[6], 7, 14]", "{0, 2, 2, 4, 4, 6}", 0);
    assert_eval_eq("BitAnd[Range[5], Range[5], Reverse[Range[5]]]", "{1, 0, 3, 0, 1}", 0);
    assert_eval_eq("BitAnd[NDArray[{7, 6, 5}, DataType -> \"int64\"], "
                   "NDArray[{3, 3, 3}, DataType -> \"int64\"], 2]",
                   "NDArray[{2, 2, 0}]", 0);
    assert_eval_eq("BitAnd[Range[4], -1, -1]", "{1, 2, 3, 4}", 0);
    assert_eval_eq("BitAnd[Range[4], 0]", "{0, 0, 0, 0}", 0);
    assert_eval_eq("BitAnd[Range[3], -2]", "{0, 2, 2}", 0);
    /* Off-buffer operands take the exact List path. */
    assert_eval_eq("BitAnd[Range[3], 2^70 + 2]", "{0, 2, 2}", 0);
    assert_eval_eq("BitAnd[Range[3], x, x]",
                   "{BitAnd[1, x], BitAnd[2, x], BitAnd[3, x]}", 0);
}

static void test_bitand_surfaces_agree(void) {
    assert_eval_eq(
        "SeedRandom[5]; ba = Join[RandomInteger[{-2^62, 2^62}, 300], {-2^63, 2^63 - 1, 0, -1}]; "
        "bb = Join[RandomInteger[{-2^62, 2^62}, 300], {-1, -2^63, 5, -2^63}]; "
        "ref = Table[BitAnd[ba[[i]], bb[[i]]], {i, Length[ba]}]; "
        "{BitAnd[ba, bb] === ref, "
        " Normal[BitAnd[NDArray[ba, DataType -> \"int64\"], NDArray[bb, DataType -> \"int64\"]]] === ref, "
        " Compile[{{v, _Integer, 1}, {w, _Integer, 1}}, BitAnd[v, w]][ba, bb] === ref, "
        " Table[Compile[{{s, _Integer}, {t, _Integer}}, BitAnd[s, t]][ba[[i]], bb[[i]]], {i, Length[ba]}] === ref, "
        " BitAnd[ba, bb, -99] === Table[BitAnd[ba[[i]], bb[[i]], -99], {i, Length[ba]}], "
        " Compile[{{v, _Integer, 1}, {w, _Integer, 1}}, BitAnd[v, w, -99]][ba, bb] === BitAnd[ba, bb, -99]}",
        "{True, True, True, True, True, True}", 0);
}

/* --- Compile[] ------------------------------------------------------- */

static void test_bitand_compile(void) {
    assert_eval_eq("Compile[{{a, _Integer}, {b, _Integer}}, BitAnd[a, b]][61, 15]", "13", 0);
    assert_eval_eq("Compile[{{a, _Integer}, {b, _Integer}}, BitAnd[a, b]][-2, -3]", "-4", 0);
    assert_eval_eq(
        "Compile[{{a, _Integer}, {b, _Integer}, {c, _Integer}, {d, _Integer}}, "
        "BitAnd[a, b, c, d]][3333, 5555, 7777, 9999]", "1025", 0);
    assert_eval_eq("Compile[{{v, _Integer, 1}}, BitAnd[v, 17]][Range[70, 80]]",
                   "{0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 16}", 0);
    assert_eval_eq("Compile[{{v, _Integer, 1}}, BitAnd[7, v, 14]][Range[6]]",
                   "{0, 2, 2, 4, 4, 6}", 0);
    assert_eval_eq(
        "Compile[{{v, _Integer, 1}, {w, _Integer, 1}}, BitAnd[v, w, 6]][{7, 15}, {3, 12}]",
        "{2, 4}", 0);
}

static void test_bitand_compile_diagnostics(void) {
    assert_eval_eq("\"Compiled\" /. CompileDiagnostics[{{a, _Integer}, {b, _Integer}}, BitAnd[a, b]]", "True", 0);
    assert_eval_eq("\"Compiled\" /. CompileDiagnostics[{{a, _Integer}, {b, _Integer}, {c, _Integer}}, BitAnd[a, b, c]]", "True", 0);
    assert_eval_eq("\"Compiled\" /. CompileDiagnostics[{{v, _Integer, 1}}, BitAnd[v, 255]]", "True", 0);
    assert_eval_eq("\"Compiled\" /. CompileDiagnostics[{{v, _Integer, 1}, {w, _Integer, 1}}, BitAnd[v, w]]", "True", 0);
    assert_eval_eq("\"Compiled\" /. CompileDiagnostics[{{v, _Integer, 1}}, BitAnd[v, 6, 4]]", "True", 0);
}

/* --- Introspection --------------------------------------------------- */

static void test_bitand_attributes(void) {
    assert_eval_eq("Attributes[BitAnd]",
                   "{Flat, Listable, OneIdentity, Orderless, Protected}", 0);
}

static void test_bitand_docstring_and_sym(void) {
    SymbolDef* def = symtab_get_def("BitAnd");
    ASSERT(def != NULL);
    ASSERT(def->docstring != NULL);
    ASSERT(strstr(def->docstring, "bitwise AND") != NULL);
    ASSERT(SYM_BitAnd != NULL && strcmp(SYM_BitAnd, "BitAnd") == 0);
}

static void test_bitand_stress(void) {
    for (int k = 0; k < 50; k++) {
        assert_eval_eq("BitAnd[61, 15]", "13", 0);
        assert_eval_eq("BitAnd[3^100, 5^100, -7]",
                       "503913279078601980753385379100266097026814443729", 0);
        assert_eval_eq("BitAnd[-1, x, x, y]", "BitAnd[x, y]", 0);
        assert_eval_eq("BitAnd[x, y, x, 2^70 + 3, 7]", "BitAnd[3, x, y]", 0);
        assert_eval_eq("BitAnd[Range[6], 7, 14]", "{0, 2, 2, 4, 4, 6}", 0);
        assert_eval_eq("BitAnd[Range[3], 2^70 + 2]", "{0, 2, 2}", 0);
    }
}

/* ===================================================================== *
 *  BitOr[n1, n2, ...]                                                    *
 * ===================================================================== */

/* --- Documented examples -------------------------------------------- */

static void test_bitor_documented_examples(void) {
    assert_eval_eq("BitOr[61, 15]", "63", 0);
    assert_eval_eq("BitOr[Range[-40, 40], 17]",
        "{-39, -39, -37, -37, -35, -35, -33, -33, -15, -15, -13, -13, -11, -11, "
        "-9, -9, -7, -7, -5, -5, -3, -3, -1, -1, -15, -15, -13, -13, -11, -11, "
        "-9, -9, -7, -7, -5, -5, -3, -3, -1, -1, 17, 17, 19, 19, 21, 21, 23, 23, "
        "25, 25, 27, 27, 29, 29, 31, 31, 17, 17, 19, 19, 21, 21, 23, 23, 25, 25, "
        "27, 27, 29, 29, 31, 31, 49, 49, 51, 51, 53, 53, 55, 55, 57}", 0);
    assert_eval_eq("BitOr[3^100, 5^100]",
        "7888609052210118054117297117069515706082347426840895402878464599718897", 0);
    assert_eval_eq("BitOr[3333, 5555, 7777, 9999]", "16383", 0);
    assert_eval_eq("BitOr[-2, -3]", "-1", 0);
    assert_eval_eq("BitOr[x, y, y, x]", "BitOr[x, y]", 0);
    assert_eval_eq("BitOr[1, 2, x]", "BitOr[3, x]", 0);
    assert_eval_eq("BitOr[-1, 2, x]", "-1", 0);
    assert_eval_eq("Outer[BitOr, {1, 0}, {1, 0}]", "{{1, 1}, {1, 0}}", 0);
    assert_eval_eq("BitOr[y, x]", "BitOr[x, y]", 0);
    assert_eval_eq("BitOr[]", "0", 0);
    assert_eval_eq("BitOr[{}]", "{}", 0);
    assert_eval_eq("BitOr[{}, {}]", "{}", 0);
}

static void test_bitor_basic(void) {
    assert_eval_eq("BitOr[-6, 3]", "-5", 0);
    assert_eval_eq("BitOr[0, 0]", "0", 0);
    assert_eval_eq("BitOr[7, 7]", "7", 0);
    assert_eval_eq("BitOr[5]", "5", 0);
    assert_eval_eq("BitOr[0]", "0", 0);
    assert_eval_eq("BitOr[-1, -1]", "-1", 0);
    assert_eval_eq("BitOr[-7, -7]", "-7", 0);
    assert_eval_eq("BitOr[1, 2, 4, 8]", "15", 0);
    assert_eval_eq("BitOr[-9223372036854775808, 9223372036854775807]", "-1", 0);
}

static void test_bitor_listable_shapes(void) {
    assert_eval_eq("BitOr[Range[16], 5]",
                   "{5, 7, 7, 5, 5, 7, 7, 13, 13, 15, 15, 13, 13, 15, 15, 21}", 0);
    assert_eval_eq("BitOr[{1, 2}, 8, {4, 5}]", "{13, 15}", 0);
    assert_eval_eq("BitOr[{{1, 2}, {3, 4}}, 1]", "{{1, 3}, {3, 5}}", 0);
    assert_eval_eq("BitOr[{1, 2}, {3, 4}]", "{3, 6}", 0);
    assert_eval_eq("BitOr[{1, 2}, {3, 4}, {8, 16}]", "{11, 22}", 0);
    assert_eval_eq("BitOr[-1, {1, 2}]", "{-1, -1}", 0);
    assert_eval_eq("BitOr[{1, 2}, x]", "{BitOr[1, x], BitOr[2, x]}", 0);
    assert_eval_eq("BitOr[x, {1, 0}, {-2, 0}]", "{-1, x}", 0);
}

/* --- Bignums and the int64 boundary ---------------------------------- */

static void test_bitor_bignums(void) {
    assert_eval_eq("BitOr[2^63, -1]", "-1", 0);
    assert_eval_eq("BitOr[2^63, 0]", "9223372036854775808", 0);
    assert_eval_eq("BitOr[-2^63, 2^63 - 1]", "-1", 0);
    assert_eval_eq("BitOr[2^70, 5]", "1180591620717411303429", 0);
    assert_eval_eq("BitOr[2^100, -2^100]", "-1267650600228229401496703205376", 0);
    assert_eval_eq("BitOr[-3^50, 7, -5^30]", "-5814238382973363081", 0);
    assert_eval_eq("BitOr[-3^50, 2^80 - 1]", "-1", 0);
    assert_eval_eq("BitOr[-2^100 - 1, 2^101 - 1]", "-1", 0);
    assert_eval_eq("BitOr[{2^64 + 3, 1}, 2^64 + 1]",
                   "{18446744073709551619, 18446744073709551617}", 0);
    assert_eval_eq("BitOr[-2^70, 1]", "-1180591620717411303423", 0);
    assert_eval_eq("BitOr[3^100, 5^100, 6]",
        "7888609052210118054117297117069515706082347426840895402878464599718903", 0);
    /* A negative machine integer ORed into a positive bignum sign-extends:
     * the result shrinks back into int64 and normalises to a machine integer. */
    assert_eval_eq("BitOr[-5, 2^70]", "-5", 0);
    assert_eval_eq("MachineIntegerQ[BitOr[-5, 2^70]]", "True", 0);
    assert_eval_eq("MachineIntegerQ[BitOr[-3^50, 7, -5^30]]", "True", 0);
    assert_eval_eq("MachineIntegerQ[BitOr[2^100, -2^100]]", "False", 0);
}

static void test_bitor_bignum_operands_untouched(void) {
    /* The bignum path ORs in place into a private result node, set from the
     * first BigInt: the operands -- bound, shared or repeated -- are unchanged. */
    assert_eval_eq("p = 3^200; q = 5^150; r = BitOr[p, q, p, 0]; "
                   "{p == 3^200, q == 5^150, r == BitOr[3^200, 5^150]}",
                   "{True, True, True}", 0);
    assert_eval_eq("p = 2^100 + 7; {BitOr[p, p], p == 2^100 + 7}",
                   "{1267650600228229401496703205383, True}", 0);
    assert_eval_eq("l = {2^90 + 1, 2^91 + 3}; BitOr[l, 2^90 + 3]; l == {2^90 + 1, 2^91 + 3}",
                   "True", 0);
    assert_eval_eq("p = -2^100; {BitOr[p, -7], p == -2^100}", "{-7, True}", 0);
}

/* --- Symbolic simplification ----------------------------------------- */

static void test_bitor_idempotence(void) {
    assert_eval_eq("BitOr[x, x]", "x", 0);
    assert_eval_eq("BitOr[x, x, x]", "x", 0);
    assert_eval_eq("BitOr[x, y, z, x]", "BitOr[x, y, z]", 0);
    assert_eval_eq("BitOr[f[u], f[u], 3]", "BitOr[3, f[u]]", 0);
    assert_eval_eq("BitOr[u, 3, u]", "BitOr[3, u]", 0);
    assert_eval_eq("BitOr[7, x, x, y]", "BitOr[7, x, y]", 0);
    assert_eval_eq("BitOr[2.5, 2.5]", "2.5", 0);
    assert_eval_eq("BitOr[x, 2.5, 2.5]", "BitOr[2.5, x]", 0);
    assert_eval_eq("BitOr[I, I, x]", "BitOr[I, x]", 0);
}

static void test_bitor_integer_folding(void) {
    assert_eval_eq("BitOr[x, 2^70 + 3, 2^70 + 1]", "BitOr[1180591620717411303427, x]", 0);
    assert_eval_eq("BitOr[-5, x]", "BitOr[-5, x]", 0);
    assert_eval_eq("BitOr[-3^50, -5^40, x, x]",
                   "BitOr[-37802229124771058671937, x]", 0);
    assert_eval_eq("BitOr[2^70, 2^70 + 5, x]", "BitOr[1180591620717411303429, x]", 0);
    assert_eval_eq("BitOr[x, y, x, 1, 2, 4]", "BitOr[7, x, y]", 0);
    assert_eval_eq("BitOr[1, 2, 2.5]", "BitOr[2.5, 3]", 0);
    assert_eval_eq("BitOr[1, 2, 2.5, x]", "BitOr[2.5, 3, x]", 0);
}

static void test_bitor_minus_one_absorbs(void) {
    assert_eval_eq("BitOr[-1, x]", "-1", 0);
    assert_eval_eq("BitOr[1, -2, x, x, y]", "-1", 0);
    assert_eval_eq("BitOr[x, -2^70, 2^70 - 1]", "-1", 0);
    assert_eval_eq("BitOr[-1, 2.5]", "-1", 0);
    assert_eval_eq("BitOr[-1, 1/2]", "-1", 0);
    assert_eval_eq("BitOr[-1, x, 1/2]", "-1", 0);
    assert_eval_eq("BitOr[-1, x, y, 2.5, 1/2]", "-1", 0);
    assert_eval_eq("BitOr[1, 2, 0, -1, x]", "-1", 0);
}

static void test_bitor_zero_identity(void) {
    assert_eval_eq("BitOr[0, x]", "x", 0);
    assert_eval_eq("BitOr[0, 0, x]", "x", 0);
    assert_eval_eq("BitOr[x, 0, 0]", "x", 0);
    assert_eval_eq("BitOr[0, x, y]", "BitOr[x, y]", 0);
    assert_eval_eq("BitOr[x, y, 0, y]", "BitOr[x, y]", 0);
    assert_eval_eq("BitOr[0, 2.5]", "2.5", 0);
    assert_eval_eq("BitOr[2^70, x, 0]", "BitOr[1180591620717411303424, x]", 0);
    assert_eval_eq("BitOr[x, 0, 2.5, 2.5]", "BitOr[2.5, x]", 0);
}

static void test_bitor_no_extra_rules(void) {
    /* Neither complement nor absorption rules apply. */
    assert_eval_eq("BitOr[x, BitNot[x]]", "BitOr[x, BitNot[x]]", 0);
    assert_eval_eq("BitOr[x, BitAnd[x, y]]", "BitOr[x, BitAnd[x, y]]", 0);
    assert_eval_eq("BitOr[BitNot[x], BitNot[y]]", "BitOr[BitNot[x], BitNot[y]]", 0);
}

static void test_bitor_orderless_flat(void) {
    assert_eval_eq("BitOr[x]", "x", 0);
    assert_eval_eq("BitOr[BitOr[x, 5], BitOr[10, y]]", "BitOr[15, x, y]", 0);
    assert_eval_eq("BitOr[BitOr[x, y], x]", "BitOr[x, y]", 0);
}

/* --- Non-integers: left alone, silently ------------------------------ */

static void test_bitor_non_integers_silent(void) {
    const char* in[]  = { "BitOr[2.5, 3]", "BitOr[1/2, x]", "BitOr[I, x]",
                          "BitOr[1.0, 3]", "BitOr[x, y]" };
    const char* out[] = { "BitOr[2.5, 3]", "BitOr[1/2, x]", "BitOr[I, x]",
                          "BitOr[1.0, 3]", "BitOr[x, y]" };
    for (size_t i = 0; i < sizeof in / sizeof in[0]; i++) {
        char* result = NULL;
        char* err = eval_capturing_stderr(in[i], &result);
        ASSERT(result != NULL);
        ASSERT_STR_EQ(result, out[i]);
        ASSERT_MSG(err == NULL || strstr(err, "BitOr::") == NULL,
                   "unexpected diagnostic for %s: %s", in[i], err ? err : "(null)");
        free(result);
        free(err);
    }
    /* No message fires, so Check takes the value branch. */
    assert_eval_eq("Check[BitOr[2.5, 3], failed]", "BitOr[2.5, 3]", 0);
}

static void test_bitor_ragged_tdlen(void) {
    assert_eval_eq("BitOr[{1, 2}, {1, 2, 3}]", "BitOr[{1, 2}, {1, 2, 3}]", 0);
    assert_eval_eq("BitOr[Range[2], Range[3], 4]", "BitOr[4, {1, 2}, {1, 2, 3}]", 0);
}

/* --- Properties ------------------------------------------------------ */

static void test_bitor_properties(void) {
    assert_eval_eq(
        "And @@ Flatten[Table[BitOr[BitOr[a, b], c] == BitOr[a, BitOr[b, c]], "
        "{a, -6, 6}, {b, -6, 6}, {c, -6, 6}]]", "True", 0);
    /* Independent oracle: bitwise max of the binary digits. */
    assert_eval_eq(
        "And @@ Flatten[Table[BitOr[a, b] == "
        "FromDigits[MapThread[Max, {IntegerDigits[a, 2, 8], IntegerDigits[b, 2, 8]}], 2], "
        "{a, 0, 63}, {b, 0, 63}]]", "True", 0);
    /* Two's-complement identities across every sign and size:
     *   a + b == BitOr[a, b] + BitAnd[a, b]
     *   BitOr[a, b] == BitXor[a, b] + BitAnd[a, b]  (disjoint bits)
     *   BitOr[a, b] == BitNot[BitAnd[BitNot[a], BitNot[b]]]  (De Morgan, -x-1) */
    assert_eval_eq(
        "SeedRandom[11]; ns = Join[RandomInteger[{-2^200, 2^200}, 30], "
        "{-2^63, 2^63 - 1, 2^63, -2^63 - 1, 0, -1, 1, 3^100, -5^90}]; "
        "And @@ Flatten[Table[m + n == BitOr[m, n] + BitAnd[m, n] && "
        "BitOr[m, n] == BitXor[m, n] + BitAnd[m, n] && "
        "BitOr[m, n] == -BitAnd[-m - 1, -n - 1] - 1, {m, ns}, {n, ns}]]",
        "True", 0);
    /* BitOr[n, 0] == n; BitOr[n, -1] == -1. */
    assert_eval_eq(
        "And @@ Table[BitOr[n, 0] == n && BitOr[n, -1] == -1, "
        "{n, {-2^80, -7, 0, 9, 2^63, 3^60}}]", "True", 0);
}

/* --- Packed / NDArray surfaces --------------------------------------- */

static void test_bitor_packed_and_ndarray(void) {
    assert_eval_eq("PackedArrayQ[BitOr[Range[10], 6]]", "True", 0);
    assert_eval_eq("BitOr[NDArray[{1, 2, 3}, DataType -> \"int64\"], 6]",
                   "NDArray[{7, 6, 7}]", 0);
    assert_eval_eq(
        "BitOr[NDArray[{-9223372036854775808, 0, 0}, DataType -> \"int64\"], "
        "NDArray[{9223372036854775807, -9223372036854775808, 0}, DataType -> \"int64\"]]",
        "NDArray[{-1, -9223372036854775808, 0}]", 0);
    /* n-ary: folded on the buffer, scalars combined first. */
    assert_eval_eq("BitOr[Range[6], 8, 16]", "{25, 26, 27, 28, 29, 30}", 0);
    assert_eval_eq("PackedArrayQ[BitOr[Range[6], 8, 16]]", "True", 0);
    assert_eval_eq("BitOr[Range[5], Range[5], Reverse[Range[5]]]", "{5, 6, 3, 6, 5}", 0);
    assert_eval_eq("BitOr[NDArray[{1, 2, 4}, DataType -> \"int64\"], "
                   "NDArray[{8, 8, 8}, DataType -> \"int64\"], 16]",
                   "NDArray[{25, 26, 28}]", 0);
    assert_eval_eq("BitOr[Range[4], 0, 0]", "{1, 2, 3, 4}", 0);
    assert_eval_eq("BitOr[Range[4], -1]", "{-1, -1, -1, -1}", 0);
    assert_eval_eq("BitOr[Range[3], -4]", "{-3, -2, -1}", 0);
    /* Off-buffer operands take the exact List path. */
    assert_eval_eq("BitOr[Range[3], 2^70]",
        "{1180591620717411303425, 1180591620717411303426, 1180591620717411303427}", 0);
    assert_eval_eq("BitOr[Range[3], x, x]",
                   "{BitOr[1, x], BitOr[2, x], BitOr[3, x]}", 0);
}

static void test_bitor_surfaces_agree(void) {
    assert_eval_eq(
        "SeedRandom[9]; oa = Join[RandomInteger[{-2^62, 2^62}, 300], {-2^63, 2^63 - 1, 0, -1}]; "
        "ob = Join[RandomInteger[{-2^62, 2^62}, 300], {0, 2^63 - 1, 5, -2^63}]; "
        "ref = Table[BitOr[oa[[i]], ob[[i]]], {i, Length[oa]}]; "
        "{BitOr[oa, ob] === ref, "
        " Normal[BitOr[NDArray[oa, DataType -> \"int64\"], NDArray[ob, DataType -> \"int64\"]]] === ref, "
        " Compile[{{v, _Integer, 1}, {w, _Integer, 1}}, BitOr[v, w]][oa, ob] === ref, "
        " Table[Compile[{{s, _Integer}, {t, _Integer}}, BitOr[s, t]][oa[[i]], ob[[i]]], {i, Length[oa]}] === ref, "
        " BitOr[oa, ob, 99] === Table[BitOr[oa[[i]], ob[[i]], 99], {i, Length[oa]}], "
        " Compile[{{v, _Integer, 1}, {w, _Integer, 1}}, BitOr[v, w, 99]][oa, ob] === BitOr[oa, ob, 99]}",
        "{True, True, True, True, True, True}", 0);
}

/* --- Compile[] ------------------------------------------------------- */

static void test_bitor_compile(void) {
    assert_eval_eq("Compile[{{a, _Integer}, {b, _Integer}}, BitOr[a, b]][61, 15]", "63", 0);
    assert_eval_eq("Compile[{{a, _Integer}, {b, _Integer}}, BitOr[a, b]][-2, -3]", "-1", 0);
    assert_eval_eq(
        "Compile[{{a, _Integer}, {b, _Integer}, {c, _Integer}, {d, _Integer}}, "
        "BitOr[a, b, c, d]][3333, 5555, 7777, 9999]", "16383", 0);
    assert_eval_eq("Compile[{{v, _Integer, 1}}, BitOr[v, 17]][Range[-4, 4]]",
                   "{-3, -3, -1, -1, 17, 17, 19, 19, 21}", 0);
    assert_eval_eq("Compile[{{v, _Integer, 1}}, BitOr[8, v, 16]][Range[6]]",
                   "{25, 26, 27, 28, 29, 30}", 0);
    assert_eval_eq(
        "Compile[{{v, _Integer, 1}, {w, _Integer, 1}}, BitOr[v, w, 64]][{1, 2}, {4, 8}]",
        "{69, 74}", 0);
}

static void test_bitor_compile_diagnostics(void) {
    assert_eval_eq("\"Compiled\" /. CompileDiagnostics[{{a, _Integer}, {b, _Integer}}, BitOr[a, b]]", "True", 0);
    assert_eval_eq("\"Compiled\" /. CompileDiagnostics[{{a, _Integer}, {b, _Integer}, {c, _Integer}}, BitOr[a, b, c]]", "True", 0);
    assert_eval_eq("\"Compiled\" /. CompileDiagnostics[{{v, _Integer, 1}}, BitOr[v, 255]]", "True", 0);
    assert_eval_eq("\"Compiled\" /. CompileDiagnostics[{{v, _Integer, 1}, {w, _Integer, 1}}, BitOr[v, w]]", "True", 0);
    assert_eval_eq("\"Compiled\" /. CompileDiagnostics[{{v, _Integer, 1}}, BitOr[v, 6, 4]]", "True", 0);
}

/* --- Introspection --------------------------------------------------- */

static void test_bitor_attributes(void) {
    assert_eval_eq("Attributes[BitOr]",
                   "{Flat, Listable, OneIdentity, Orderless, Protected}", 0);
}

static void test_bitor_docstring_and_sym(void) {
    SymbolDef* def = symtab_get_def("BitOr");
    ASSERT(def != NULL);
    ASSERT(def->docstring != NULL);
    ASSERT(strstr(def->docstring, "bitwise OR") != NULL);
    ASSERT(SYM_BitOr != NULL && strcmp(SYM_BitOr, "BitOr") == 0);
}

static void test_bitor_stress(void) {
    for (int k = 0; k < 50; k++) {
        assert_eval_eq("BitOr[61, 15]", "63", 0);
        assert_eval_eq("BitOr[3^100, 5^100, 6]",
            "7888609052210118054117297117069515706082347426840895402878464599718903", 0);
        assert_eval_eq("BitOr[0, x, x, y]", "BitOr[x, y]", 0);
        assert_eval_eq("BitOr[x, y, x, 1, 2, 4]", "BitOr[7, x, y]", 0);
        assert_eval_eq("BitOr[-5, 2^70]", "-5", 0);
        assert_eval_eq("BitOr[Range[6], 8, 16]", "{25, 26, 27, 28, 29, 30}", 0);
        assert_eval_eq("BitOr[Range[3], -2^70]",
            "{-1180591620717411303423, -1180591620717411303422, -1180591620717411303421}", 0);
    }
}

/* ===================================================================== *
 *  BitNot[n] = -1 - n (two's complement, unlimited digits)
 * ===================================================================== */

static void test_bitnot_documented_examples(void) {
    assert_eval_eq("BitNot[0]", "-1", 0);
    assert_eval_eq("BitNot[7]", "-8", 0);
    assert_eval_eq("BitNot[Range[-10, 10]]",
        "{9, 8, 7, 6, 5, 4, 3, 2, 1, 0, -1, -2, -3, -4, -5, -6, -7, -8, -9, -10, -11}", 0);
    assert_eval_eq("BitNot[BitNot[x]]", "x", 0);
    assert_eval_eq("BitNot[{}]", "{}", 0);

    char* result = NULL;
    char* err = eval_capturing_stderr("BitNot[]", &result);
    ASSERT(result != NULL);
    ASSERT_STR_EQ(result, "BitNot[]");
    ASSERT_MSG(err != NULL && strstr(err,
        "BitNot::argx: BitNot called with 0 arguments; 1 argument is expected.") != NULL,
        "expected argx diagnostic, got: %s", err ? err : "(null)");
    free(result);
    free(err);
}

static void test_bitnot_basic(void) {
    assert_eval_eq("BitNot[-1]", "0", 0);
    assert_eval_eq("BitNot[1]", "-2", 0);
    assert_eval_eq("BitNot[-8]", "7", 0);
    assert_eval_eq("BitNot[255]", "-256", 0);
    assert_eval_eq("BitNot[-256]", "255", 0);
}

static void test_bitnot_int64_extremes(void) {
    /* ~INT64_MIN = INT64_MAX and vice versa: total on int64, no overflow. */
    assert_eval_eq("BitNot[-9223372036854775808]", "9223372036854775807", 0);
    assert_eval_eq("BitNot[9223372036854775807]", "-9223372036854775808", 0);
    assert_eval_eq("BitNot[-2^63]", "9223372036854775807", 0);
    assert_eval_eq("Head[BitNot[2^63 - 1]]", "Integer", 0);
}

static void test_bitnot_bignums(void) {
    /* Crossing the int64 boundary in both directions. */
    assert_eval_eq("BitNot[2^63]", "-9223372036854775809", 0);
    assert_eval_eq("BitNot[-2^63 - 1]", "9223372036854775808", 0);
    assert_eval_eq("BitNot[2^70]", "-1180591620717411303425", 0);
    assert_eval_eq("BitNot[-2^70]", "1180591620717411303423", 0);
    assert_eval_eq("BitNot[100!] == -1 - 100!", "True", 0);
    assert_eval_eq("BitNot[-3^200] == 3^200 - 1", "True", 0);
    /* A bignum whose complement lands back in int64 range is demoted. */
    assert_eval_eq("BitNot[-9223372036854775809]", "9223372036854775808", 0);
    assert_eval_eq("BitNot[9223372036854775808]", "-9223372036854775809", 0);
    assert_eval_eq("BitNot[BitNot[2^63]] === 2^63", "True", 0);
}

static void test_bitnot_bignum_operands_untouched(void) {
    /* mpz_com writes a fresh node; a bound / shared operand is unchanged. */
    assert_eval_eq("p = 3^200; r = BitNot[p]; {p == 3^200, r == -1 - 3^200}",
                   "{True, True}", 0);
    assert_eval_eq("l = {2^90 + 1, -2^91}; BitNot[l]; l == {2^90 + 1, -2^91}",
                   "True", 0);
}

static void test_bitnot_involution(void) {
    assert_eval_eq("BitNot[BitNot[BitNot[x]]]", "BitNot[x]", 0);
    assert_eval_eq("BitNot[BitNot[BitNot[BitNot[x]]]]", "x", 0);
    assert_eval_eq("BitNot[BitNot[x + y]]", "x + y", 0);
    assert_eval_eq("BitNot[BitNot[{x, 2.5}]]", "{x, 2.5}", 0);
    /* The inner BitNot is shared, not mutated. */
    assert_eval_eq("b = BitNot[z]; {BitNot[b], b}", "{z, BitNot[z]}", 0);
    /* BitXor's -1 extraction now reduces against BitNot. */
    assert_eval_eq("BitXor[-1, BitNot[x]]", "x", 0);
    assert_eval_eq("BitNot[BitXor[-1, x]]", "x", 0);
}

static void test_bitnot_no_extra_rules(void) {
    assert_eval_eq("BitNot[-1 - x]", "BitNot[-1 - x]", 0);
    assert_eval_eq("BitAnd[x, BitNot[x]]", "BitAnd[x, BitNot[x]]", 0);
    assert_eval_eq("BitNot[x]", "BitNot[x]", 0);
}

static void test_bitnot_non_integers_silent(void) {
    const char* in[] = { "BitNot[2.5]", "BitNot[1.]", "BitNot[1/2]", "BitNot[I]",
                         "BitNot[Pi]", "BitNot[2 + 3 I]", "BitNot[x]" };
    const char* out[] = { "BitNot[2.5]", "BitNot[1.0]", "BitNot[1/2]", "BitNot[I]",
                          "BitNot[Pi]", "BitNot[2 + 3*I]", "BitNot[x]" };
    for (size_t i = 0; i < sizeof in / sizeof in[0]; i++) {
        char* result = NULL;
        char* err = eval_capturing_stderr(in[i], &result);
        ASSERT(result != NULL);
        ASSERT_STR_EQ(result, out[i]);
        ASSERT_MSG(err == NULL || strstr(err, "BitNot::") == NULL,
                   "unexpected diagnostic for %s: %s", in[i], err ? err : "(null)");
        free(result);
        free(err);
    }
    assert_eval_eq("Check[BitNot[2.5], failed]", "BitNot[2.5]", 0);
}

static void test_bitnot_argx(void) {
    char* result = NULL;
    char* err = eval_capturing_stderr("BitNot[x, y]", &result);
    ASSERT(result != NULL);
    ASSERT_STR_EQ(result, "BitNot[x, y]");
    ASSERT_MSG(err != NULL && strstr(err,
        "BitNot::argx: BitNot called with 2 arguments; 1 argument is expected.") != NULL,
        "expected argx diagnostic, got: %s", err ? err : "(null)");
    free(result);
    free(err);
    /* Routed through the message funnel: Quiet silences, Check sees it. */
    assert_eval_eq("Check[BitNot[1, 2], failed]", "failed", 0);
    err = eval_capturing_stderr("Quiet[BitNot[1, 2]]", &result);
    ASSERT(result != NULL);
    ASSERT_STR_EQ(result, "BitNot[1, 2]");
    ASSERT_MSG(err == NULL || strstr(err, "BitNot::") == NULL,
               "Quiet leaked: %s", err ? err : "(null)");
    free(result);
    free(err);
}

static void test_bitnot_listable_shapes(void) {
    assert_eval_eq("BitNot[{1, {2, x}}]", "{-2, {-3, BitNot[x]}}", 0);
    assert_eval_eq("BitNot[{{0, -1}, {2^70, 2.5}}]",
                   "{{-1, 0}, {-1180591620717411303425, BitNot[2.5]}}", 0);
    assert_eval_eq("BitNot[{{}, {}}]", "{{}, {}}", 0);
}

static void test_bitnot_properties(void) {
    /* BitNot[n] == -1 - n and involution, over every sign and size. */
    assert_eval_eq(
        "SeedRandom[13]; ns = Join[RandomInteger[{-2^200, 2^200}, 40], "
        "RandomInteger[{-2^62, 2^62}, 40], "
        "{-2^63, 2^63 - 1, 2^63, -2^63 - 1, 0, -1, 1, 3^100, -5^90}]; "
        "And @@ Table[BitNot[n] == -1 - n && BitNot[BitNot[n]] === n && "
        "BitXor[n, -1] == BitNot[n] && BitAnd[n, BitNot[n]] == 0 && "
        "BitOr[n, BitNot[n]] == -1, {n, ns}]", "True", 0);
    /* De Morgan, both ways. */
    assert_eval_eq(
        "SeedRandom[14]; ns = Join[RandomInteger[{-2^100, 2^100}, 15], {-2^63, 2^63 - 1, 0, -1}]; "
        "And @@ Flatten[Table[BitNot[BitAnd[m, n]] == BitOr[BitNot[m], BitNot[n]] && "
        "BitNot[BitOr[m, n]] == BitAnd[BitNot[m], BitNot[n]], {m, ns}, {n, ns}]]",
        "True", 0);
    /* Independent oracle: flip each of the low 8 binary digits. */
    assert_eval_eq(
        "And @@ Table[Mod[BitNot[a], 256] == FromDigits[1 - IntegerDigits[a, 2, 8], 2], "
        "{a, 0, 255}]", "True", 0);
}

static void test_bitnot_packed_and_ndarray(void) {
    assert_eval_eq("PackedArrayQ[BitNot[Range[10]]]", "True", 0);
    assert_eval_eq("BitNot[Range[5]]", "{-2, -3, -4, -5, -6}", 0);
    assert_eval_eq(
        "BitNot[NDArray[{-9223372036854775808, -1, 0, 9223372036854775807}, DataType -> \"int64\"]]",
        "NDArray[{9223372036854775807, 0, -1, -9223372036854775808}]", 0);
    assert_eval_eq("BitNot[NDArray[{{1, 2}, {3, 4}}, DataType -> \"int64\"]]",
                   "NDArray[{{-2, -3}, {-4, -5}}]", 0);
    assert_eval_eq("PackedArrayQ[BitNot[Partition[Range[12], 3]]]", "True", 0);
    /* A real buffer has no integer answer: elements stay unevaluated. */
    assert_eval_eq("BitNot[NDArray[{1.5, 2.}]]", "{BitNot[1.5], BitNot[2.0]}", 0);
    /* Off-buffer elements take the exact List path. */
    assert_eval_eq("BitNot[{1, 2^70, x}]", "{-2, -1180591620717411303425, BitNot[x]}", 0);
    assert_eval_eq("BitNot[BitNot[Range[6]]] === Range[6]", "True", 0);
}

static void test_bitnot_surfaces_agree(void) {
    assert_eval_eq(
        "SeedRandom[10]; oa = Join[RandomInteger[{-2^62, 2^62}, 500], "
        "{-2^63, 2^63 - 1, 0, -1}]; "
        "ref = Table[-1 - oa[[i]], {i, Length[oa]}]; "
        "{BitNot[oa] === ref, "
        " Normal[BitNot[NDArray[oa, DataType -> \"int64\"]]] === ref, "
        " Compile[{{v, _Integer, 1}}, BitNot[v]][oa] === ref, "
        " Table[Compile[{{s, _Integer}}, BitNot[s]][oa[[i]]], {i, Length[oa]}] === ref, "
        " Table[BitNot[oa[[i]]], {i, Length[oa]}] === ref}",
        "{True, True, True, True, True}", 0);
}

static void test_bitnot_compile(void) {
    assert_eval_eq("Compile[{{n, _Integer}}, BitNot[n]][5]", "-6", 0);
    assert_eval_eq("Compile[{{n, _Integer}}, BitNot[n]][-9223372036854775808]",
                   "9223372036854775807", 0);
    assert_eval_eq("Compile[{{n, _Integer}}, BitNot[BitNot[n] + 3]][10]", "7", 0);
    assert_eval_eq("Compile[{{v, _Integer, 1}}, BitNot[v]][Range[-3, 3]]",
                   "{2, 1, 0, -1, -2, -3, -4}", 0);
    assert_eval_eq("Compile[{{v, _Integer, 1}}, BitAnd[v, BitNot[6]]][Range[8]]",
                   "{1, 0, 1, 0, 1, 0, 1, 8}", 0);
    /* Not an integer: the interpreter leaves it unevaluated, so must Compile. */
    assert_eval_eq("Compile[{{x, _Real}}, BitNot[x]][2.5]", "BitNot[2.5]", 0);
}

static void test_bitnot_compile_diagnostics(void) {
    assert_eval_eq("\"Compiled\" /. CompileDiagnostics[{{n, _Integer}}, BitNot[n]]", "True", 0);
    assert_eval_eq("\"Compiled\" /. CompileDiagnostics[{{v, _Integer, 1}}, BitNot[v]]", "True", 0);
    assert_eval_eq("\"Compiled\" /. CompileDiagnostics[{{v, _Integer, 1}}, BitXor[v, BitNot[v]]]", "True", 0);
}

/* Auto-compilation (Table / Sum over a real iterator) lowers a body holding
 * BitNot through OP_BNOT_I; switching it off must not change the answer. */
static Expr* bitnot_eval_owned(const char* src) {
    Expr* p = parse_expression(src);
    Expr* e = evaluate(p);   /* evaluate borrows p */
    expr_free(p);
    return e;
}

static void test_bitnot_autocompile(void) {
    const char* tbl = "Table[x + BitNot[Floor[x]], {x, 0., 2000., 0.5}]";
    const char* sm  = "Sum[Sin[x] BitNot[Floor[x]], {x, 0., 500., 0.25}]";
    assert_eval_eq("Table[x + BitNot[Floor[x]], {x, 0., 2., 0.5}]",
                   "{-1.0, -0.5, -1.0, -0.5, -1.0}", 0);
    autocompile_set_enabled(false);
    Expr* t_off = bitnot_eval_owned(tbl);
    Expr* s_off = bitnot_eval_owned(sm);
    autocompile_set_enabled(true);
    Expr* t_on = bitnot_eval_owned(tbl);
    Expr* s_on = bitnot_eval_owned(sm);
    ASSERT(expr_eq(t_off, t_on));
    ASSERT(expr_eq(s_off, s_on));
    ASSERT(s_on->type == EXPR_REAL);
    expr_free(t_off); expr_free(s_off); expr_free(t_on); expr_free(s_on);
}

static void test_bitnot_attributes(void) {
    assert_eval_eq("Attributes[BitNot]", "{Listable, Protected}", 0);
}

static void test_bitnot_docstring_and_sym(void) {
    SymbolDef* def = symtab_get_def("BitNot");
    ASSERT(def != NULL);
    ASSERT(def->docstring != NULL);
    ASSERT(strstr(def->docstring, "bitwise NOT") != NULL);
    ASSERT(SYM_BitNot != NULL && strcmp(SYM_BitNot, "BitNot") == 0);
}

static void test_bitnot_stress(void) {
    for (int k = 0; k < 50; k++) {
        assert_eval_eq("BitNot[7]", "-8", 0);
        assert_eval_eq("BitNot[2^70]", "-1180591620717411303425", 0);
        assert_eval_eq("BitNot[-2^63 - 1]", "9223372036854775808", 0);
        assert_eval_eq("BitNot[BitNot[x]]", "x", 0);
        assert_eval_eq("BitNot[{1, {2, x}}]", "{-2, {-3, BitNot[x]}}", 0);
        assert_eval_eq("BitNot[Range[4]]", "{-2, -3, -4, -5}", 0);
    }
}

/* ---- BitSet ------------------------------------------------------------ */

static void test_bitset_documented_examples(void) {
    assert_eval_eq("IntegerDigits[42, 2]", "{1, 0, 1, 0, 1, 0}", 0);
    assert_eval_eq("BitSet[42, 0]", "43", 0);
    assert_eval_eq("IntegerDigits[BitSet[42, 0], 2]", "{1, 0, 1, 0, 1, 1}", 0);
    assert_eval_eq("BitSet[42, 4]", "58", 0);
    assert_eval_eq("IntegerDigits[BitSet[42, 4], 2]", "{1, 1, 1, 0, 1, 0}", 0);
    /* BitSet grows the integer to the size it needs. */
    assert_eval_eq("BitSet[42, 9]", "554", 0);
    assert_eval_eq("IntegerDigits[BitSet[42, 9], 2]",
                   "{1, 0, 0, 0, 1, 0, 1, 0, 1, 0}", 0);
}

static void test_bitset_basic(void) {
    assert_eval_eq("BitSet[0, 0]", "1", 0);
    assert_eval_eq("BitSet[42, 1]", "42", 0);   /* already set */
    assert_eval_eq("BitSet[42, 3]", "42", 0);
    assert_eval_eq("BitSet[0, 10]", "1024", 0);
    assert_eval_eq("BitSet[BitSet[0, 1], 2]", "6", 0);
    assert_eval_eq("Table[BitSet[0, k], {k, 0, 7}]",
                   "{1, 2, 4, 8, 16, 32, 64, 128}", 0);
}

static void test_bitset_negatives(void) {
    /* Two's complement, infinitely sign-extended. */
    assert_eval_eq("BitSet[-5, 1]", "-5", 0);
    assert_eval_eq("BitSet[-5, 2]", "-1", 0);
    assert_eval_eq("BitSet[-1, 0]", "-1", 0);
    assert_eval_eq("BitSet[-1, 100]", "-1", 0);
    assert_eval_eq("BitSet[-8, 0]", "-7", 0);
    assert_eval_eq("BitSet[-256, 7]", "-128", 0);
}

static void test_bitset_int64_boundary(void) {
    assert_eval_eq("BitSet[0, 62]", "4611686018427387904", 0);
    assert_eval_eq("Head[BitSet[0, 62]]", "Integer", 0);
    /* Bit 63 of a non-negative n is a bignum. */
    assert_eval_eq("BitSet[0, 63]", "9223372036854775808", 0);
    assert_eval_eq("BitSet[0, 64]", "18446744073709551616", 0);
    assert_eval_eq("BitSet[9223372036854775807, 63]", "18446744073709551615", 0);
    assert_eval_eq("BitSet[2^62, 62]", "4611686018427387904", 0);
    /* A negative int64 already has every bit from 63 up set. */
    assert_eval_eq("BitSet[-2^63, 63]", "-9223372036854775808", 0);
    assert_eval_eq("BitSet[-2^63, 64]", "-9223372036854775808", 0);
    assert_eval_eq("BitSet[-2^63, 62]", "-4611686018427387904", 0);
    assert_eval_eq("BitSet[-2^63, 0]", "-9223372036854775807", 0);
    assert_eval_eq("BitSet[-7, 1000]", "-7", 0);
}

static void test_bitset_bignums(void) {
    assert_eval_eq("BitSet[2^70, 0]", "1180591620717411303425", 0);
    assert_eval_eq("BitSet[2^70, 70]", "1180591620717411303424", 0);
    assert_eval_eq("BitSet[2^70, 71] == 3 2^70", "True", 0);
    assert_eval_eq("BitSet[-2^70, 0] == -2^70 + 1", "True", 0);
    assert_eval_eq("BitSet[-2^70, 70]", "-1180591620717411303424", 0);
    assert_eval_eq("BitSet[100!, 1000] == 100! + 2^1000", "True", 0);
    assert_eval_eq("BitLength[BitSet[1, 10^6]]", "1000001", 0);
    /* A negative bignum moved back into int64 range is demoted. */
    assert_eval_eq("BitSet[-2^64, 63]", "-9223372036854775808", 0);
    assert_eval_eq("Head[BitSet[-2^64, 63]]", "Integer", 0);
    assert_eval_eq("BitSet[-2^63 - 1, 63] == -1", "True", 0);
}

static void test_bitset_bignum_operands_untouched(void) {
    /* mpz_setbit writes a fresh node; a bound / shared operand is unchanged. */
    assert_eval_eq("p = 3^200; r = BitSet[p, 0]; {p == 3^200, r == 3^200}",
                   "{True, True}", 0);
    assert_eval_eq("p = 2^90; r = BitSet[p, 3]; {p == 2^90, r == 2^90 + 8}",
                   "{True, True}", 0);
    assert_eval_eq("l = {2^90, -2^91}; BitSet[l, 5]; l == {2^90, -2^91}", "True", 0);
}

static void test_bitset_unevaluated_silent(void) {
    /* No message for any of these, wrong argument counts included. */
    const char* in[] = { "BitSet[42, -1]", "BitSet[-10, -2]", "BitSet[2.5, 1]",
                         "BitSet[5, 1.]", "BitSet[x, 1]", "BitSet[5, x]",
                         "BitSet[1/2, 1]", "BitSet[5, 1/2]", "BitSet[I, 0]",
                         "BitSet[1, 2^70]", "BitSet[5]", "BitSet[5, 1, 2]",
                         "BitSet[]", "BitSet[BitSet[x, 1], 1]" };
    const char* out[] = { "BitSet[42, -1]", "BitSet[-10, -2]", "BitSet[2.5, 1]",
                          "BitSet[5, 1.0]", "BitSet[x, 1]", "BitSet[5, x]",
                          "BitSet[1/2, 1]", "BitSet[5, 1/2]", "BitSet[I, 0]",
                          "BitSet[1, 1180591620717411303424]", "BitSet[5]",
                          "BitSet[5, 1, 2]", "BitSet[]", "BitSet[BitSet[x, 1], 1]" };
    for (size_t i = 0; i < sizeof in / sizeof in[0]; i++) {
        char* result = NULL;
        char* err = eval_capturing_stderr(in[i], &result);
        ASSERT(result != NULL);
        ASSERT_STR_EQ(result, out[i]);
        ASSERT_MSG(err == NULL || strstr(err, "::") == NULL,
                   "unexpected diagnostic for %s: %s", in[i], err ? err : "(null)");
        free(result);
        free(err);
    }
    /* Check is the funnel-level oracle: no message of any kind fired. */
    for (size_t i = 0; i < sizeof in / sizeof in[0]; i++) {
        char src[128], want[128];
        snprintf(src, sizeof src, "Check[%s, failed]", in[i]);
        snprintf(want, sizeof want, "%s", out[i]);
        assert_eval_eq(src, want, 0);
    }
}

static void test_bitset_listable_shapes(void) {
    assert_eval_eq("BitSet[{1, 2, 3}, {0, 1, 2}]", "{1, 2, 7}", 0);
    assert_eval_eq("BitSet[{1, 2, 3}, 2]", "{5, 6, 7}", 0);
    assert_eval_eq("BitSet[3, {0, 1, 5}]", "{3, 3, 35}", 0);
    assert_eval_eq("BitSet[{1, {2, x}}, 3]", "{9, {10, BitSet[x, 3]}}", 0);
    assert_eval_eq("BitSet[{0, 2^70, -1, 2.5}, 64] === {2^64, 2^70 + 2^64, -1, BitSet[2.5, 64]}",
                   "True", 0);
    assert_eval_eq("BitSet[{}, 3]", "{}", 0);
    assert_eval_eq("BitSet[{{1, 2}, {3, 4}}, {1, 2}]", "{{3, 2}, {7, 4}}", 0);
}

static void test_bitset_ragged_tdlen(void) {
    /* Unequal lengths: Thread::tdlen from the evaluator's Listable threader,
     * the call left unevaluated (also from a packed operand). */
    assert_eval_eq("BitSet[{1, 2, 3}, {0, 1}]", "BitSet[{1, 2, 3}, {0, 1}]", 0);
    assert_eval_eq("BitSet[Range[3], Range[2]]", "BitSet[{1, 2, 3}, {1, 2}]", 0);
    /* Routed through the message funnel: printed on stderr, silenced by
     * Quiet, seen by Check. */
    char* result = NULL;
    char* err = eval_capturing_stderr("BitSet[{1, 2, 3}, {0, 1}]", &result);
    ASSERT(result != NULL);
    ASSERT_STR_EQ(result, "BitSet[{1, 2, 3}, {0, 1}]");
    ASSERT_MSG(err != NULL && strstr(err, "Thread::tdlen: Objects of unequal length in "
                                          "BitSet[{1, 2, 3}, {0, 1}] cannot be combined.") != NULL,
               "expected tdlen diagnostic, got: %s", err ? err : "(null)");
    free(result);
    free(err);
    err = eval_capturing_stderr("Quiet[BitSet[{1, 2, 3}, {0, 1}]]", &result);
    ASSERT(result != NULL);
    ASSERT_STR_EQ(result, "BitSet[{1, 2, 3}, {0, 1}]");
    ASSERT_MSG(err == NULL || strstr(err, "tdlen") == NULL,
               "Quiet leaked: %s", err ? err : "(null)");
    free(result);
    free(err);
    assert_eval_eq("Check[BitSet[{1, 2, 3}, {0, 1}], failed]", "failed", 0);
}

static void test_bitset_properties(void) {
    /* BitSet[n, k] == BitOr[n, 2^k]; idempotent; the bit reads back as set;
     * no other bit changes. Over every sign and size. */
    assert_eval_eq(
        "SeedRandom[21]; ns = Join[RandomInteger[{-2^200, 2^200}, 30], "
        "RandomInteger[{-2^62, 2^62}, 30], "
        "{-2^63, 2^63 - 1, 2^63, -2^63 - 1, 0, -1, 1, 3^100, -5^90}]; "
        "ks = {0, 1, 5, 31, 32, 61, 62, 63, 64, 65, 127, 128, 300}; "
        "And @@ Flatten[Table[BitSet[n, k] == BitOr[n, 2^k] && "
        "BitSet[BitSet[n, k], k] == BitSet[n, k] && "
        "BitAnd[BitSet[n, k], 2^k] == 2^k && "
        "BitXor[BitSet[n, k], n] == If[BitAnd[n, 2^k] == 0, 2^k, 0], "
        "{n, ns}, {k, ks}]]", "True", 0);
    /* Independent oracle: the binary digits of a non-negative n. */
    assert_eval_eq(
        "And @@ Flatten[Table[IntegerDigits[BitSet[a, k], 2, 12] == "
        "ReplacePart[IntegerDigits[a, 2, 12], 12 - k -> 1], {a, 0, 255}, {k, 0, 11}]]",
        "True", 0);
    /* Setting bits commutes. */
    assert_eval_eq(
        "And @@ Flatten[Table[BitSet[BitSet[n, j], k] == BitSet[BitSet[n, k], j], "
        "{n, {-9, 0, 77, 2^80}}, {j, {0, 3, 63, 90}}, {k, {1, 62, 64}}]]", "True", 0);
}

static void test_bitset_packed_and_ndarray(void) {
    assert_eval_eq("PackedArrayQ[BitSet[Range[10], 3]]", "True", 0);
    assert_eval_eq("BitSet[Range[10], 3]", "{9, 10, 11, 12, 13, 14, 15, 8, 9, 10}", 0);
    assert_eval_eq("PackedArrayQ[BitSet[5, Range[0, 9]]]", "True", 0);
    assert_eval_eq("BitSet[5, Range[0, 4]]", "{5, 7, 5, 13, 21}", 0);
    assert_eval_eq("PackedArrayQ[BitSet[Range[8], Range[0, 7]]]", "True", 0);
    assert_eval_eq("BitSet[Range[8], Range[0, 7]]", "{1, 2, 7, 12, 21, 38, 71, 136}", 0);
    assert_eval_eq("BitSet[NDArray[{1, 2, 3}, DataType -> \"int64\"], 4]",
                   "NDArray[{17, 18, 19}]", 0);
    assert_eval_eq("BitSet[5, NDArray[{0, 1, 2}, DataType -> \"int64\"]]",
                   "NDArray[{5, 7, 5}]", 0);
    assert_eval_eq("BitSet[NDArray[{1, 2, 3}, DataType -> \"int64\"], "
                   "NDArray[{0, 1, 62}, DataType -> \"int64\"]]",
                   "NDArray[{1, 2, 4611686018427387907}]", 0);
    assert_eval_eq("BitSet[NDArray[{{1, 2}, {3, 4}}, DataType -> \"int64\"], 3]",
                   "NDArray[{{9, 10}, {11, 12}}]", 0);
    assert_eval_eq("PackedArrayQ[BitSet[Partition[Range[12], 3], 6]]", "True", 0);
    /* A real buffer has no integer answer: elements stay unevaluated. */
    assert_eval_eq("BitSet[NDArray[{1.5, 2.}], 1]", "{BitSet[1.5, 1], BitSet[2.0, 1]}", 0);
}

static void test_bitset_buffer_declines(void) {
    /* Bit 63 of a non-negative buffer element is a bignum: the buffer is
     * abandoned and the List path answers exactly -- never a truncation. */
    assert_eval_eq("BitSet[Range[3], 63]",
                   "{9223372036854775809, 9223372036854775810, 9223372036854775811}", 0);
    assert_eval_eq("BitSet[NDArray[{0, 1}, DataType -> \"int64\"], 64]",
                   "{18446744073709551616, 18446744073709551617}", 0);
    /* All-negative buffer past bit 62: still int64, unchanged. */
    assert_eval_eq("BitSet[-Range[5], 63]", "{-1, -2, -3, -4, -5}", 0);
    assert_eval_eq("BitSet[-Range[5], 1000]", "{-1, -2, -3, -4, -5}", 0);
    assert_eval_eq("BitSet[-5, Range[60, 66]]", "{-5, -5, -5, -5, -5, -5, -5}", 0);
    /* Mixed signs past bit 62. */
    assert_eval_eq("BitSet[{-1, 0, 1}, 63]",
                   "{-1, 9223372036854775808, 9223372036854775809}", 0);
    /* A negative k element: that element stays unevaluated. */
    assert_eval_eq("BitSet[Range[3], {0, -1, 2}]", "{1, BitSet[2, -1], 7}", 0);
    assert_eval_eq("BitSet[1, {-1, 2}]", "{BitSet[1, -1], 5}", 0);
    assert_eval_eq("BitSet[NDArray[{1, 2}, DataType -> \"int64\"], -1]",
                   "{BitSet[1, -1], BitSet[2, -1]}", 0);
    assert_eval_eq("BitSet[{1, 2^70, x}, 1]",
                   "{3, 1180591620717411303426, BitSet[x, 1]}", 0);
}

static void test_bitset_surfaces_agree(void) {
    /* List == packed == NDArray == Compile (array and scalar), for every arm:
     * array n with scalar k, scalar n with array k, and two arrays. */
    assert_eval_eq(
        "SeedRandom[22]; na = Join[RandomInteger[{-2^62, 2^62}, 500], "
        "{-2^63, 2^63 - 1, 0, -1}]; ka = RandomInteger[{0, 62}, Length[na]]; "
        "ref1 = Table[BitOr[na[[i]], 2^37], {i, Length[na]}]; "
        "ref2 = Table[BitOr[-12345, 2^ka[[i]]], {i, Length[na]}]; "
        "ref3 = Table[BitOr[na[[i]], 2^ka[[i]]], {i, Length[na]}]; "
        "{BitSet[na, 37] === ref1, "
        " Normal[BitSet[NDArray[na, DataType -> \"int64\"], 37]] === ref1, "
        " Compile[{{v, _Integer, 1}}, BitSet[v, 37]][na] === ref1, "
        " BitSet[-12345, ka] === ref2, "
        " Normal[BitSet[-12345, NDArray[ka, DataType -> \"int64\"]]] === ref2, "
        " Compile[{{v, _Integer, 1}}, BitSet[-12345, v]][ka] === ref2, "
        " BitSet[na, ka] === ref3, "
        " Normal[BitSet[NDArray[na, DataType -> \"int64\"], NDArray[ka, DataType -> \"int64\"]]] === ref3, "
        " Compile[{{v, _Integer, 1}, {w, _Integer, 1}}, BitSet[v, w]][na, ka] === ref3, "
        " Table[Compile[{{s, _Integer}, {t, _Integer}}, BitSet[s, t]][na[[i]], ka[[i]]], {i, Length[na]}] === ref3, "
        " Table[BitSet[na[[i]], ka[[i]]], {i, Length[na]}] === ref3}",
        "{True, True, True, True, True, True, True, True, True, True, True}", 0);
}

static void test_bitset_compile(void) {
    assert_eval_eq("Compile[{{n, _Integer}, {k, _Integer}}, BitSet[n, k]][42, 9]", "554", 0);
    assert_eval_eq("Compile[{{n, _Integer}}, BitSet[n, 4]][42]", "58", 0);
    assert_eval_eq("Compile[{{n, _Integer}}, BitSet[n, 62]][0]", "4611686018427387904", 0);
    assert_eval_eq("Compile[{{n, _Integer}}, BitSet[n, 70]][-7]", "-7", 0);
    /* No int64 answer: the VM bails and the interpreter answers exactly. */
    assert_eval_eq("Compile[{{n, _Integer}}, BitSet[n, 63]][0]", "9223372036854775808", 0);
    assert_eval_eq("Compile[{{n, _Integer}, {k, _Integer}}, BitSet[n, k]][42, -1]",
                   "BitSet[42, -1]", 0);
    assert_eval_eq("Compile[{{v, _Integer, 1}}, BitSet[v, 2]][Range[5]]", "{5, 6, 7, 4, 5}", 0);
    assert_eval_eq("Compile[{{v, _Integer, 1}}, BitSet[v, 63]][Range[2]]",
                   "{9223372036854775809, 9223372036854775810}", 0);
    assert_eval_eq("Compile[{{n, _Integer}}, BitAnd[BitSet[n, 3], 15] + BitNot[n]][1]", "7", 0);
    /* Not an integer: the interpreter leaves it unevaluated, so must Compile. */
    assert_eval_eq("Compile[{{x, _Real}}, BitSet[x, 1]][2.5]", "BitSet[2.5, 1]", 0);
}

static void test_bitset_compile_diagnostics(void) {
    assert_eval_eq("\"Compiled\" /. CompileDiagnostics[{{n, _Integer}, {k, _Integer}}, BitSet[n, k]]", "True", 0);
    assert_eval_eq("\"Compiled\" /. CompileDiagnostics[{{v, _Integer, 1}}, BitSet[v, 3]]", "True", 0);
    assert_eval_eq("\"Compiled\" /. CompileDiagnostics[{{v, _Integer, 1}}, BitSet[7, v]]", "True", 0);
    assert_eval_eq("\"Compiled\" /. CompileDiagnostics[{{v, _Integer, 1}, {w, _Integer, 1}}, BitSet[v, w]]", "True", 0);
}

static void test_bitset_autocompile(void) {
    const char* tbl = "Table[x + BitSet[Floor[x], 3], {x, 0., 2000., 0.5}]";
    const char* sm  = "Sum[Sin[x] BitSet[Floor[x], 2], {x, 0., 500., 0.25}]";
    assert_eval_eq("Table[x + BitSet[Floor[x], 3], {x, 0., 2., 0.5}]",
                   "{8.0, 8.5, 10.0, 10.5, 12.0}", 0);
    autocompile_set_enabled(false);
    Expr* t_off = bitnot_eval_owned(tbl);
    Expr* s_off = bitnot_eval_owned(sm);
    autocompile_set_enabled(true);
    Expr* t_on = bitnot_eval_owned(tbl);
    Expr* s_on = bitnot_eval_owned(sm);
    ASSERT(expr_eq(t_off, t_on));
    ASSERT(expr_eq(s_off, s_on));
    ASSERT(s_on->type == EXPR_REAL);
    expr_free(t_off); expr_free(s_off); expr_free(t_on); expr_free(s_on);
}

static void test_bitset_attributes(void) {
    assert_eval_eq("Attributes[BitSet]", "{Listable, Protected}", 0);
}

static void test_bitset_docstring_and_sym(void) {
    SymbolDef* def = symtab_get_def("BitSet");
    ASSERT(def != NULL);
    ASSERT(def->docstring != NULL);
    ASSERT(strstr(def->docstring, "coefficient of 2^k") != NULL);
    ASSERT(SYM_BitSet != NULL && strcmp(SYM_BitSet, "BitSet") == 0);
}

static void test_bitset_stress(void) {
    for (int k = 0; k < 50; k++) {
        assert_eval_eq("BitSet[42, 9]", "554", 0);
        assert_eval_eq("BitSet[0, 64]", "18446744073709551616", 0);
        assert_eval_eq("BitSet[2^70, 0]", "1180591620717411303425", 0);
        assert_eval_eq("BitSet[-2^64, 63]", "-9223372036854775808", 0);
        assert_eval_eq("BitSet[x, 1]", "BitSet[x, 1]", 0);
        assert_eval_eq("BitSet[Range[4], 63] == Range[4] + 2^63", "True", 0);
        assert_eval_eq("BitSet[Range[4], 2]", "{5, 6, 7, 4}", 0);
    }
}

/* ---- BitClear ---------------------------------------------------------- */

static void test_bitclear_basic(void) {
    assert_eval_eq("BitClear[255, 0]", "254", 0);
    assert_eval_eq("BitClear[255, 7]", "127", 0);
    assert_eval_eq("BitClear[42, 0]", "42", 0);    /* already clear */
    assert_eval_eq("BitClear[42, 1]", "40", 0);
    assert_eval_eq("BitClear[0, 10]", "0", 0);
    assert_eval_eq("Table[BitClear[255, k], {k, 0, 7}]",
                   "{254, 253, 251, 247, 239, 223, 191, 127}", 0);
    /* Non-negative n: every high bit is already clear. */
    assert_eval_eq("BitClear[5, 100]", "5", 0);
    assert_eval_eq("BitClear[9223372036854775807, 62]", "4611686018427387903", 0);
    assert_eval_eq("BitClear[9223372036854775807, 63]", "9223372036854775807", 0);
}

static void test_bitclear_negatives_and_bignums(void) {
    assert_eval_eq("BitClear[-1, 0]", "-2", 0);
    assert_eval_eq("BitClear[-5, 0]", "-6", 0);
    assert_eval_eq("BitClear[-5, 1]", "-7", 0);
    assert_eval_eq("BitClear[-5, 2]", "-5", 0);    /* already clear */
    /* Clearing a high bit of a negative n subtracts 2^k: a bignum. */
    assert_eval_eq("BitClear[-1, 63]", "-9223372036854775809", 0);
    assert_eval_eq("BitClear[-1, 64]", "-18446744073709551617", 0);
    assert_eval_eq("BitClear[-1, 100] == -1 - 2^100", "True", 0);
    assert_eval_eq("BitClear[-2^63, 63]", "-18446744073709551616", 0);
    assert_eval_eq("BitClear[-2^64, 64]", "-36893488147419103232", 0);
    assert_eval_eq("BitClear[2^70 + 5, 70]", "5", 0);
    assert_eval_eq("Head[BitClear[2^70 + 5, 70]]", "Integer", 0);
    assert_eval_eq("BitClear[2^63, 63]", "0", 0);
    assert_eval_eq("BitClear[2^70, 71] == 2^70", "True", 0);
    assert_eval_eq("BitClear[2^70, 10^12] == 2^70", "True", 0);
    assert_eval_eq("BitClear[100! + 2^1000, 1000] == 100!", "True", 0);
    assert_eval_eq("p = 3^200; r = BitClear[p, 0]; {p == 3^200, r == 3^200 - 1}",
                   "{True, True}", 0);
}

static void test_bitclear_unevaluated_silent(void) {
    const char* in[] = { "BitClear[5, -1]", "BitClear[2.5, 1]", "BitClear[5, 1.]",
                         "BitClear[x, 1]", "BitClear[5, x]", "BitClear[1/2, 1]",
                         "BitClear[5, 2^70]", "BitClear[5]", "BitClear[5, 1, 2]" };
    const char* out[] = { "BitClear[5, -1]", "BitClear[2.5, 1]", "BitClear[5, 1.0]",
                          "BitClear[x, 1]", "BitClear[5, x]", "BitClear[1/2, 1]",
                          "BitClear[5, 1180591620717411303424]", "BitClear[5]",
                          "BitClear[5, 1, 2]" };
    for (size_t i = 0; i < sizeof in / sizeof in[0]; i++) {
        char src[128];
        snprintf(src, sizeof src, "Check[%s, failed]", in[i]);
        assert_eval_eq(src, out[i], 0);
    }
}

static void test_bitclear_properties(void) {
    /* BitClear[n, k] == BitAnd[n, BitNot[2^k]]; idempotent; the bit reads back
     * clear; BitSet undoes it when the bit was set. Over every sign and size. */
    assert_eval_eq(
        "SeedRandom[31]; ns = Join[RandomInteger[{-2^200, 2^200}, 30], "
        "RandomInteger[{-2^62, 2^62}, 30], "
        "{-2^63, 2^63 - 1, 2^63, -2^63 - 1, 0, -1, 1, 3^100, -5^90}]; "
        "ks = {0, 1, 5, 31, 32, 61, 62, 63, 64, 65, 127, 128, 300}; "
        "And @@ Flatten[Table[BitClear[n, k] == BitAnd[n, BitNot[2^k]] && "
        "BitClear[BitClear[n, k], k] == BitClear[n, k] && "
        "BitGet[BitClear[n, k], k] == 0 && "
        "BitSet[BitClear[n, k], k] == BitSet[n, k], "
        "{n, ns}, {k, ks}]]", "True", 0);
}

static void test_bitclear_packed_and_ndarray(void) {
    assert_eval_eq("PackedArrayQ[BitClear[Range[10], 0]]", "True", 0);
    assert_eval_eq("BitClear[Range[10], 0]", "{0, 2, 2, 4, 4, 6, 6, 8, 8, 10}", 0);
    assert_eval_eq("PackedArrayQ[BitClear[255, Range[0, 7]]]", "True", 0);
    assert_eval_eq("BitClear[255, Range[0, 3]]", "{254, 253, 251, 247}", 0);
    assert_eval_eq("BitClear[Range[8], Range[0, 7]]", "{0, 0, 3, 4, 5, 6, 7, 8}", 0);
    assert_eval_eq("BitClear[NDArray[{7, 15, 3}, DataType -> \"int64\"], 1]",
                   "NDArray[{5, 13, 1}]", 0);
    assert_eval_eq("BitClear[255, NDArray[{0, 1, 2}, DataType -> \"int64\"]]",
                   "NDArray[{254, 253, 251}]", 0);
    assert_eval_eq("BitClear[NDArray[{7, 7}, DataType -> \"int64\"], "
                   "NDArray[{0, 2}, DataType -> \"int64\"]]", "NDArray[{6, 3}]", 0);
    /* No int64 answer for a negative element past bit 62: the buffer is
     * abandoned and the List path answers exactly, never a truncation. */
    assert_eval_eq("BitClear[{-1, 0, 1}, 63]", "{-9223372036854775809, 0, 1}", 0);
    assert_eval_eq("BitClear[Range[3], 63]", "{1, 2, 3}", 0);
    assert_eval_eq("BitClear[Range[3], {0, -1, 2}]", "{0, BitClear[2, -1], 3}", 0);
}

static void test_bitclear_surfaces_agree(void) {
    assert_eval_eq(
        "SeedRandom[32]; na = Join[RandomInteger[{-2^62, 2^62}, 500], "
        "{-2^63, 2^63 - 1, 0, -1}]; ka = RandomInteger[{0, 62}, Length[na]]; "
        "ref1 = Table[BitAnd[na[[i]], BitNot[2^37]], {i, Length[na]}]; "
        "ref2 = Table[BitAnd[-12345, BitNot[2^ka[[i]]]], {i, Length[na]}]; "
        "ref3 = Table[BitAnd[na[[i]], BitNot[2^ka[[i]]]], {i, Length[na]}]; "
        "{BitClear[na, 37] === ref1, "
        " Normal[BitClear[NDArray[na, DataType -> \"int64\"], 37]] === ref1, "
        " Compile[{{v, _Integer, 1}}, BitClear[v, 37]][na] === ref1, "
        " BitClear[-12345, ka] === ref2, "
        " Normal[BitClear[-12345, NDArray[ka, DataType -> \"int64\"]]] === ref2, "
        " Compile[{{v, _Integer, 1}}, BitClear[-12345, v]][ka] === ref2, "
        " BitClear[na, ka] === ref3, "
        " BitClear[Developer`FromPackedArray[na], Developer`FromPackedArray[ka]] === ref3, "
        " Normal[BitClear[NDArray[na, DataType -> \"int64\"], NDArray[ka, DataType -> \"int64\"]]] === ref3, "
        " Compile[{{v, _Integer, 1}, {w, _Integer, 1}}, BitClear[v, w]][na, ka] === ref3, "
        " Table[Compile[{{s, _Integer}, {t, _Integer}}, BitClear[s, t]][na[[i]], ka[[i]]], {i, Length[na]}] === ref3, "
        " Table[BitClear[na[[i]], ka[[i]]], {i, Length[na]}] === ref3}",
        "{True, True, True, True, True, True, True, True, True, True, True, True}", 0);
}

static void test_bitclear_compile(void) {
    assert_eval_eq("Compile[{{n, _Integer}, {k, _Integer}}, BitClear[n, k]][255, 0]", "254", 0);
    assert_eval_eq("Compile[{{n, _Integer}}, BitClear[n, 70]][7]", "7", 0);
    /* No int64 answer: the VM bails and the interpreter answers exactly. */
    assert_eval_eq("Compile[{{n, _Integer}}, BitClear[n, 63]][-1]", "-9223372036854775809", 0);
    assert_eval_eq("Compile[{{n, _Integer}, {k, _Integer}}, BitClear[n, k]][42, -1]",
                   "BitClear[42, -1]", 0);
    assert_eval_eq("Compile[{{v, _Integer, 1}}, BitClear[v, 0]][Range[5]]", "{0, 2, 2, 4, 4}", 0);
    assert_eval_eq("\"Compiled\" /. CompileDiagnostics[{{x, _Integer}}, BitClear[x, 3]]", "True", 0);
    assert_eval_eq("\"Compiled\" /. CompileDiagnostics[{{n, _Integer}, {k, _Integer}}, BitClear[n, k]]", "True", 0);
    assert_eval_eq("\"Compiled\" /. CompileDiagnostics[{{v, _Integer, 1}}, BitClear[v, 3]]", "True", 0);
    assert_eval_eq("\"Compiled\" /. CompileDiagnostics[{{v, _Integer, 1}}, BitClear[7, v]]", "True", 0);
    assert_eval_eq("\"Compiled\" /. CompileDiagnostics[{{v, _Integer, 1}, {w, _Integer, 1}}, BitClear[v, w]]", "True", 0);
}

static void test_bitclear_attributes_docstring(void) {
    assert_eval_eq("Attributes[BitClear]", "{Listable, Protected}", 0);
    SymbolDef* def = symtab_get_def("BitClear");
    ASSERT(def != NULL);
    ASSERT(def->docstring != NULL);
    ASSERT(strstr(def->docstring, "coefficient of") != NULL);
    ASSERT(SYM_BitClear != NULL && strcmp(SYM_BitClear, "BitClear") == 0);
}

/* ---- BitGet ------------------------------------------------------------ */

static void test_bitget_documented_examples(void) {
    assert_eval_eq("IntegerDigits[42, 2]", "{1, 0, 1, 0, 1, 0}", 0);
    assert_eval_eq("BitGet[42, 0]", "0", 0);
    assert_eval_eq("BitGet[42, 3]", "1", 0);
    assert_eval_eq("BitGet[Range[30, 60], 3]",
                   "{1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, "
                   "0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1}", 0);
}

static void test_bitget_basic(void) {
    assert_eval_eq("BitGet[0, 0]", "0", 0);
    assert_eval_eq("BitGet[1, 0]", "1", 0);
    assert_eval_eq("BitGet[42, 1]", "1", 0);
    assert_eval_eq("BitGet[42, 5]", "1", 0);
    assert_eval_eq("BitGet[42, 6]", "0", 0);
    assert_eval_eq("BitGet[42, 1000]", "0", 0);
    assert_eval_eq("Table[BitGet[42, k], {k, 0, 7}]", "{0, 1, 0, 1, 0, 1, 0, 0}", 0);
    assert_eval_eq("Head[BitGet[42, 3]]", "Integer", 0);
}

static void test_bitget_negatives(void) {
    /* Two's complement, infinitely sign-extended. */
    assert_eval_eq("BitGet[-1, 0]", "1", 0);
    assert_eval_eq("BitGet[-1, 100]", "1", 0);
    assert_eval_eq("BitGet[-42, 1]", "1", 0);
    assert_eval_eq("BitGet[-42, 3]", "0", 0);
    assert_eval_eq("Table[BitGet[-42, k], {k, 0, 7}]", "{0, 1, 1, 0, 1, 0, 1, 1}", 0);
    assert_eval_eq("BitGet[-256, 7]", "0", 0);
    assert_eval_eq("BitGet[-256, 8]", "1", 0);
}

static void test_bitget_int64_boundary(void) {
    assert_eval_eq("BitGet[2^63 - 1, 62]", "1", 0);
    assert_eval_eq("BitGet[2^63 - 1, 63]", "0", 0);
    assert_eval_eq("BitGet[2^63 - 1, 64]", "0", 0);
    assert_eval_eq("BitGet[2^62, 62]", "1", 0);
    assert_eval_eq("BitGet[1, 63]", "0", 0);
    assert_eval_eq("BitGet[-2^63, 63]", "1", 0);
    assert_eval_eq("BitGet[-2^63, 62]", "0", 0);
    assert_eval_eq("BitGet[-2^63, 0]", "0", 0);
    assert_eval_eq("BitGet[-2^63, 64]", "1", 0);
    assert_eval_eq("BitGet[-7, 1000]", "1", 0);
    assert_eval_eq("BitGet[7, 2^62]", "0", 0);
    assert_eval_eq("BitGet[-7, 2^63 - 1]", "1", 0);
}

static void test_bitget_bignums(void) {
    assert_eval_eq("BitGet[2^100, 100]", "1", 0);
    assert_eval_eq("BitGet[2^100, 99]", "0", 0);
    assert_eval_eq("BitGet[2^100, 101]", "0", 0);
    assert_eval_eq("BitGet[-2^100, 100]", "1", 0);
    assert_eval_eq("BitGet[-2^100, 99]", "0", 0);
    assert_eval_eq("BitGet[-2^100, 101]", "1", 0);
    assert_eval_eq("BitGet[2^63, 63]", "1", 0);
    assert_eval_eq("BitGet[-2^63 - 1, 0]", "1", 0);
    assert_eval_eq("BitGet[-2^63 - 1, 63]", "0", 0);
    assert_eval_eq("BitGet[2^70 + 5, 2]", "1", 0);
    assert_eval_eq("BitGet[2^70, 2^62]", "0", 0);
    assert_eval_eq("BitGet[-3^200, 2^62]", "1", 0);
    assert_eval_eq("BitGet[3^200, 0]", "1", 0);
}

static void test_bitget_bignum_operands_untouched(void) {
    /* mpz_tstbit reads in place; a bound / shared operand is unchanged. */
    assert_eval_eq("p = 3^200; r = BitGet[p, 5]; p == 3^200", "True", 0);
    assert_eval_eq("p = -2^90; r = BitGet[p, 95]; {p == -2^90, r}", "{True, 1}", 0);
    assert_eval_eq("l = {2^90, -2^91}; BitGet[l, 5]; l == {2^90, -2^91}", "True", 0);
}

static void test_bitget_unevaluated_silent(void) {
    /* No message for any of these, wrong argument counts included. */
    const char* in[] = { "BitGet[42, -1]", "BitGet[-10, -2]", "BitGet[2.5, 1]",
                         "BitGet[5, 1.]", "BitGet[x, 1]", "BitGet[5, x]",
                         "BitGet[1/2, 1]", "BitGet[5, 1/2]", "BitGet[I, 0]",
                         "BitGet[1, 2^70]", "BitGet[-1, -2^70]", "BitGet[2^80, -1]",
                         "BitGet[5]", "BitGet[5, 1, 2]",
                         "BitGet[]", "BitGet[BitSet[x, 1], 1]" };
    const char* out[] = { "BitGet[42, -1]", "BitGet[-10, -2]", "BitGet[2.5, 1]",
                          "BitGet[5, 1.0]", "BitGet[x, 1]", "BitGet[5, x]",
                          "BitGet[1/2, 1]", "BitGet[5, 1/2]", "BitGet[I, 0]",
                          "BitGet[1, 1180591620717411303424]",
                          "BitGet[-1, -1180591620717411303424]",
                          "BitGet[1208925819614629174706176, -1]",
                          "BitGet[5]", "BitGet[5, 1, 2]", "BitGet[]",
                          "BitGet[BitSet[x, 1], 1]" };
    for (size_t i = 0; i < sizeof in / sizeof in[0]; i++) {
        char* result = NULL;
        char* err = eval_capturing_stderr(in[i], &result);
        ASSERT(result != NULL);
        ASSERT_STR_EQ(result, out[i]);
        ASSERT_MSG(err == NULL || strstr(err, "::") == NULL,
                   "unexpected diagnostic for %s: %s", in[i], err ? err : "(null)");
        free(result);
        free(err);
    }
    for (size_t i = 0; i < sizeof in / sizeof in[0]; i++) {
        char src[128];
        snprintf(src, sizeof src, "Check[%s, failed]", in[i]);
        assert_eval_eq(src, out[i], 0);
    }
}

static void test_bitget_listable_shapes(void) {
    assert_eval_eq("BitGet[{1, 2, 3}, {0, 1, 2}]", "{1, 1, 0}", 0);
    assert_eval_eq("BitGet[{1, 2, 3}, 1]", "{0, 1, 1}", 0);
    assert_eval_eq("BitGet[5, {0, 1, 2, 3}]", "{1, 0, 1, 0}", 0);
    assert_eval_eq("BitGet[{1, {2, x}}, 1]", "{0, {1, BitGet[x, 1]}}", 0);
    assert_eval_eq("BitGet[{0, 2^70, -1, 2.5}, 70]", "{0, 1, 1, BitGet[2.5, 70]}", 0);
    assert_eval_eq("BitGet[{}, 3]", "{}", 0);
    assert_eval_eq("BitGet[{{1, 2}, {3, 4}}, {0, 2}]", "{{1, 0}, {0, 1}}", 0);
}

static void test_bitget_ragged_tdlen(void) {
    assert_eval_eq("BitGet[Range[3], Range[2]]", "BitGet[{1, 2, 3}, {1, 2}]", 0);
    char* result = NULL;
    char* err = eval_capturing_stderr("BitGet[{1, 2, 3}, {0, 1}]", &result);
    ASSERT(result != NULL);
    ASSERT_STR_EQ(result, "BitGet[{1, 2, 3}, {0, 1}]");
    ASSERT_MSG(err != NULL && strstr(err, "Thread::tdlen: Objects of unequal length in "
                                          "BitGet[{1, 2, 3}, {0, 1}] cannot be combined.") != NULL,
               "expected tdlen diagnostic, got: %s", err ? err : "(null)");
    free(result);
    free(err);
    err = eval_capturing_stderr("Quiet[BitGet[{1, 2, 3}, {0, 1}]]", &result);
    ASSERT(result != NULL);
    ASSERT_STR_EQ(result, "BitGet[{1, 2, 3}, {0, 1}]");
    ASSERT_MSG(err == NULL || strstr(err, "tdlen") == NULL,
               "Quiet leaked: %s", err ? err : "(null)");
    free(result);
    free(err);
    assert_eval_eq("Check[BitGet[{1, 2, 3}, {0, 1}], failed]", "failed", 0);
}

static void test_bitget_properties(void) {
    /* Independent oracles over every sign and size: floor division, BitAnd
     * with 2^k, BitSet read-back, and reassembly of a non-negative n. */
    assert_eval_eq(
        "SeedRandom[31]; ns = Join[RandomInteger[{-2^200, 2^200}, 30], "
        "RandomInteger[{-2^62, 2^62}, 30], "
        "{-2^63, 2^63 - 1, 2^63, -2^63 - 1, 0, -1, 1, 3^100, -5^90}]; "
        "ks = {0, 1, 5, 31, 32, 61, 62, 63, 64, 65, 127, 128, 300, 1000}; "
        "And @@ Flatten[Table[BitGet[n, k] == Mod[Floor[n/2^k], 2] && "
        "BitGet[n, k] == If[BitAnd[n, 2^k] == 0, 0, 1] && "
        "BitGet[BitSet[n, k], k] == 1 && "
        "BitGet[BitNot[n], k] == 1 - BitGet[n, k], "
        "{n, ns}, {k, ks}]]", "True", 0);
    assert_eval_eq(
        "SeedRandom[32]; And @@ Table[Sum[BitGet[n, k] 2^k, {k, 0, BitLength[n]}] == n, "
        "{n, Join[RandomInteger[{0, 2^150}, 40], {0, 1, 2^63 - 1, 2^63, 2^64}]}]",
        "True", 0);
    assert_eval_eq(
        "And @@ Table[BitGet[a, Range[0, 11]] == Reverse[IntegerDigits[a, 2, 12]], "
        "{a, 0, 4095}]", "True", 0);
}

static void test_bitget_packed_and_ndarray(void) {
    assert_eval_eq("PackedArrayQ[BitGet[Range[10], 3]]", "True", 0);
    assert_eval_eq("BitGet[Range[10], 1]", "{0, 1, 1, 0, 0, 1, 1, 0, 0, 1}", 0);
    assert_eval_eq("PackedArrayQ[BitGet[5, Range[0, 9]]]", "True", 0);
    assert_eval_eq("BitGet[5, Range[0, 4]]", "{1, 0, 1, 0, 0}", 0);
    assert_eval_eq("PackedArrayQ[BitGet[Range[8], Range[0, 7]]]", "True", 0);
    assert_eval_eq("BitGet[2^Range[0, 7], Range[0, 7]]", "{1, 1, 1, 1, 1, 1, 1, 1}", 0);
    assert_eval_eq("BitGet[Range[8], Range[0, 7]]", "{1, 1, 0, 0, 0, 0, 0, 0}", 0);
    assert_eval_eq("BitGet[NDArray[{1, 2, 3}, DataType -> \"int64\"], 1]",
                   "NDArray[{0, 1, 1}]", 0);
    assert_eval_eq("BitGet[5, NDArray[{0, 1, 2}, DataType -> \"int64\"]]",
                   "NDArray[{1, 0, 1}]", 0);
    assert_eval_eq("BitGet[NDArray[{1, 2, -3}, DataType -> \"int64\"], "
                   "NDArray[{0, 1, 62}, DataType -> \"int64\"]]",
                   "NDArray[{1, 1, 1}]", 0);
    assert_eval_eq("BitGet[NDArray[{{1, 2}, {3, 4}}, DataType -> \"int64\"], 0]",
                   "NDArray[{{1, 0}, {1, 0}}]", 0);
    assert_eval_eq("PackedArrayQ[BitGet[Partition[Range[12], 3], 2]]", "True", 0);
    /* Past bit 62 the buffer still answers: the sign bit. */
    assert_eval_eq("BitGet[{-1, 0, 5, -5}, 63]", "{1, 0, 0, 1}", 0);
    assert_eval_eq("BitGet[{-1, 0, 5, -5}, 10^6]", "{1, 0, 0, 1}", 0);
    assert_eval_eq("BitGet[-5, Range[60, 66]]", "{1, 1, 1, 1, 1, 1, 1}", 0);
    assert_eval_eq("BitGet[{-1, 3}, {100, 2^62}]", "{1, 0}", 0);
    assert_eval_eq("PackedArrayQ[BitGet[-Range[5], 1000]]", "True", 0);
    /* A real buffer has no integer answer: elements stay unevaluated. */
    assert_eval_eq("BitGet[NDArray[{1.5, 2.}], 1]", "{BitGet[1.5, 1], BitGet[2.0, 1]}", 0);
}

static void test_bitget_buffer_declines(void) {
    /* A negative k: those elements stay unevaluated, the others answer. */
    assert_eval_eq("BitGet[Range[3], {0, -1, 2}]", "{1, BitGet[2, -1], 0}", 0);
    assert_eval_eq("BitGet[1, {-1, 0}]", "{BitGet[1, -1], 1}", 0);
    assert_eval_eq("BitGet[NDArray[{1, 2}, DataType -> \"int64\"], -1]",
                   "{BitGet[1, -1], BitGet[2, -1]}", 0);
    assert_eval_eq("BitGet[Range[2], -1]", "{BitGet[1, -1], BitGet[2, -1]}", 0);
    assert_eval_eq("BitGet[{1, 2^70, x}, 70]", "{0, 1, BitGet[x, 70]}", 0);
    assert_eval_eq("BitGet[{-2^63, 2^63 - 1}, {-2^63, 0}]",
                   "{BitGet[-9223372036854775808, -9223372036854775808], 1}", 0);
}

static void test_bitget_surfaces_agree(void) {
    /* List == packed == NDArray == Compile (array and scalar), for every arm. */
    assert_eval_eq(
        "SeedRandom[33]; na = Join[RandomInteger[{-2^62, 2^62}, 500], "
        "{-2^63, 2^63 - 1, 0, -1}]; ka = RandomInteger[{0, 100}, Length[na]]; "
        "ref1 = Table[Mod[Floor[na[[i]]/2^37], 2], {i, Length[na]}]; "
        "ref2 = Table[Mod[Floor[-12345/2^ka[[i]]], 2], {i, Length[na]}]; "
        "ref3 = Table[Mod[Floor[na[[i]]/2^ka[[i]]], 2], {i, Length[na]}]; "
        "{BitGet[na, 37] === ref1, "
        " Normal[BitGet[NDArray[na, DataType -> \"int64\"], 37]] === ref1, "
        " Compile[{{v, _Integer, 1}}, BitGet[v, 37]][na] === ref1, "
        " BitGet[-12345, ka] === ref2, "
        " Normal[BitGet[-12345, NDArray[ka, DataType -> \"int64\"]]] === ref2, "
        " Compile[{{v, _Integer, 1}}, BitGet[-12345, v]][ka] === ref2, "
        " BitGet[na, ka] === ref3, "
        " Normal[BitGet[NDArray[na, DataType -> \"int64\"], NDArray[ka, DataType -> \"int64\"]]] === ref3, "
        " Compile[{{v, _Integer, 1}, {w, _Integer, 1}}, BitGet[v, w]][na, ka] === ref3, "
        " Table[Compile[{{s, _Integer}, {t, _Integer}}, BitGet[s, t]][na[[i]], ka[[i]]], {i, Length[na]}] === ref3, "
        " Table[BitGet[na[[i]], ka[[i]]], {i, Length[na]}] === ref3, "
        " Table[BitGet[FromPackedArray[na][[i]], ka[[i]]], {i, Length[na]}] === ref3}",
        "{True, True, True, True, True, True, True, True, True, True, True, True}", 0);
}

static void test_bitget_compile(void) {
    assert_eval_eq("Compile[{{n, _Integer}, {k, _Integer}}, BitGet[n, k]][42, 3]", "1", 0);
    assert_eval_eq("Compile[{{n, _Integer}}, BitGet[n, 0]][42]", "0", 0);
    assert_eval_eq("Compile[{{n, _Integer}}, BitGet[n, 63]][-1]", "1", 0);
    assert_eval_eq("Compile[{{n, _Integer}}, BitGet[n, 1000]][7]", "0", 0);
    /* k < 0: the VM bails and the interpreter leaves it unevaluated. */
    assert_eval_eq("Compile[{{n, _Integer}, {k, _Integer}}, BitGet[n, k]][42, -1]",
                   "BitGet[42, -1]", 0);
    assert_eval_eq("Compile[{{v, _Integer, 1}}, BitGet[v, 1]][Range[5]]", "{0, 1, 1, 0, 0}", 0);
    assert_eval_eq("Compile[{{n, _Integer}}, Sum[BitGet[n, k] 2^k, {k, 0, 10}]][1234]", "1234", 0);
    assert_eval_eq("Compile[{{x, _Real}}, BitGet[x, 1]][2.5]", "BitGet[2.5, 1]", 0);
}

static void test_bitget_compile_diagnostics(void) {
    assert_eval_eq("\"Compiled\" /. CompileDiagnostics[{{n, _Integer}, {k, _Integer}}, BitGet[n, k]]", "True", 0);
    assert_eval_eq("\"Compiled\" /. CompileDiagnostics[{{v, _Integer, 1}}, BitGet[v, 3]]", "True", 0);
    assert_eval_eq("\"Compiled\" /. CompileDiagnostics[{{v, _Integer, 1}}, BitGet[7, v]]", "True", 0);
    assert_eval_eq("\"Compiled\" /. CompileDiagnostics[{{v, _Integer, 1}, {w, _Integer, 1}}, BitGet[v, w]]", "True", 0);
}

static void test_bitget_autocompile(void) {
    const char* tbl = "Table[x + BitGet[Floor[x], 3], {x, 0., 2000., 0.5}]";
    const char* sm  = "Sum[Sin[x] BitGet[Floor[x], 2], {x, 0., 500., 0.25}]";
    assert_eval_eq("Table[x + BitGet[Floor[x], 1], {x, 0., 2., 0.5}]",
                   "{0.0, 0.5, 1.0, 1.5, 3.0}", 0);
    autocompile_set_enabled(false);
    Expr* t_off = bitnot_eval_owned(tbl);
    Expr* s_off = bitnot_eval_owned(sm);
    autocompile_set_enabled(true);
    Expr* t_on = bitnot_eval_owned(tbl);
    Expr* s_on = bitnot_eval_owned(sm);
    ASSERT(expr_eq(t_off, t_on));
    ASSERT(expr_eq(s_off, s_on));
    ASSERT(s_on->type == EXPR_REAL);
    expr_free(t_off); expr_free(s_off); expr_free(t_on); expr_free(s_on);
}

static void test_bitget_attributes(void) {
    assert_eval_eq("Attributes[BitGet]", "{Listable, Protected}", 0);
}

static void test_bitget_docstring_and_sym(void) {
    SymbolDef* def = symtab_get_def("BitGet");
    ASSERT(def != NULL);
    ASSERT(def->docstring != NULL);
    ASSERT(strstr(def->docstring, "coefficient of 2^k") != NULL);
    ASSERT(SYM_BitGet != NULL && strcmp(SYM_BitGet, "BitGet") == 0);
}

static void test_bitget_stress(void) {
    for (int k = 0; k < 50; k++) {
        assert_eval_eq("BitGet[42, 3]", "1", 0);
        assert_eval_eq("BitGet[-1, 100]", "1", 0);
        assert_eval_eq("BitGet[2^100, 100]", "1", 0);
        assert_eval_eq("BitGet[-2^100, 99]", "0", 0);
        assert_eval_eq("BitGet[x, 1]", "BitGet[x, 1]", 0);
        assert_eval_eq("BitGet[Range[4], 1]", "{0, 1, 1, 0}", 0);
        assert_eval_eq("BitGet[-Range[4], 100]", "{1, 1, 1, 1}", 0);
    }
}

/* ---- BitShiftLeft / BitShiftRight --------------------------------------- */

static void test_bitshift_documented_examples(void) {
    assert_eval_eq("BitShiftLeft[32, 3]", "256", 0);
    assert_eval_eq("BitShiftLeft[2^200 - 1, 7]",
                   "205688069665150755269371147819668813122841983204197482918576000", 0);
    assert_eval_eq("BitShiftLeft[{1, 2, 4, 16}, 2]", "{4, 8, 16, 64}", 0);
    assert_eval_eq("BitShiftRight[32, 3]", "4", 0);
    assert_eval_eq("BitShiftRight[32, 7]", "0", 0);
    assert_eval_eq("BitShiftRight[2^200 - 1, 7]",
                   "12554203470773361527671578846415332832204710888928069025791", 0);
    assert_eval_eq("BitShiftRight[{1, 2, 4, 16}, 2]", "{0, 0, 1, 4}", 0);
}

static void test_bitshift_basic(void) {
    assert_eval_eq("BitShiftLeft[5]", "10", 0);
    assert_eval_eq("BitShiftRight[5]", "2", 0);
    assert_eval_eq("BitShiftLeft[0, 1000]", "0", 0);
    assert_eval_eq("BitShiftLeft[0, 10^9]", "0", 0);
    assert_eval_eq("BitShiftLeft[0, -5]", "0", 0);
    assert_eval_eq("BitShiftLeft[7, 0]", "7", 0);
    assert_eval_eq("BitShiftRight[7, 0]", "7", 0);
    assert_eval_eq("Table[BitShiftLeft[1, k], {k, 0, 8}]", "{1, 2, 4, 8, 16, 32, 64, 128, 256}", 0);
    assert_eval_eq("Table[BitShiftRight[255, k], {k, 0, 9}]",
                   "{255, 127, 63, 31, 15, 7, 3, 1, 0, 0}", 0);
    /* A negative k shifts the other way. */
    assert_eval_eq("BitShiftLeft[5, -1]", "2", 0);
    assert_eval_eq("BitShiftRight[5, -2]", "20", 0);
    assert_eval_eq("BitShiftRight[1, -62]", "4611686018427387904", 0);
    assert_eval_eq("Head[BitShiftLeft[3, 2]]", "Integer", 0);
}

static void test_bitshift_negatives_truncate(void) {
    /* A right shift truncates toward zero (IntegerPart), unlike the
     * arithmetic shift's floor. */
    assert_eval_eq("BitShiftRight[-5, 1]", "-2", 0);
    assert_eval_eq("BitShiftRight[-7, 2]", "-1", 0);
    assert_eval_eq("BitShiftRight[-1, 1]", "0", 0);
    assert_eval_eq("BitShiftRight[-1, 100]", "0", 0);
    assert_eval_eq("BitShiftLeft[-5, -1]", "-2", 0);
    assert_eval_eq("BitShiftRight[-8, 3]", "-1", 0);
    assert_eval_eq("BitShiftRight[-9, 3]", "-1", 0);
    assert_eval_eq("BitShiftLeft[-3, 2]", "-12", 0);
    assert_eval_eq("BitShiftRight[Range[-10, 10], 2]",
                   "{-2, -2, -2, -1, -1, -1, -1, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2}", 0);
}

static void test_bitshift_int64_boundary(void) {
    assert_eval_eq("BitShiftLeft[1, 62]", "4611686018427387904", 0);
    assert_eval_eq("BitShiftLeft[1, 63]", "9223372036854775808", 0);
    assert_eval_eq("BitShiftLeft[-1, 63]", "-9223372036854775808", 0);
    assert_eval_eq("Head[BitShiftLeft[-1, 63]]", "Integer", 0);
    assert_eval_eq("BitShiftLeft[-1, 64]", "-18446744073709551616", 0);
    assert_eval_eq("BitShiftLeft[2^62, 1]", "9223372036854775808", 0);
    assert_eval_eq("BitShiftLeft[-2^62, 1]", "-9223372036854775808", 0);
    assert_eval_eq("BitShiftLeft[2^63 - 1, 1]", "18446744073709551614", 0);
    assert_eval_eq("BitShiftRight[2^63 - 1, -1]", "18446744073709551614", 0);
    assert_eval_eq("BitShiftRight[-2^63, 1]", "-4611686018427387904", 0);
    assert_eval_eq("BitShiftRight[-2^63, 63]", "-1", 0);
    assert_eval_eq("BitShiftRight[-2^63, 64]", "0", 0);
    assert_eval_eq("BitShiftRight[2^63 - 1, 62]", "1", 0);
    assert_eval_eq("BitShiftRight[2^63 - 1, 63]", "0", 0);
    assert_eval_eq("BitShiftLeft[-2^63, 0]", "-9223372036854775808", 0);
    assert_eval_eq("BitShiftLeft[5, -2^63]", "0", 0);
    assert_eval_eq("BitShiftRight[5, 2^63 - 1]", "0", 0);
    /* Zero shifts to 0 by any count, even -2^63, whose negation is no int64. */
    assert_eval_eq("BitShiftRight[0, -2^63]", "0", 0);
}

static void test_bitshift_bignums(void) {
    assert_eval_eq("BitShiftLeft[2^100, -100]", "1", 0);
    assert_eval_eq("Head[BitShiftLeft[2^100, -100]]", "Integer", 0);
    assert_eval_eq("BitShiftLeft[-(2^100) + 1, -99]", "-1", 0);
    assert_eval_eq("BitShiftRight[2^64, 1]", "9223372036854775808", 0);
    assert_eval_eq("BitShiftRight[2^64 + 5, -3]", "147573952589676412968", 0);
    assert_eval_eq("BitShiftRight[-2^64, 1]", "-9223372036854775808", 0);
    assert_eval_eq("BitShiftRight[-2^64 - 1, 1]", "-9223372036854775808", 0);
    assert_eval_eq("BitShiftRight[-2^64 - 1, 2]", "-4611686018427387904", 0);
    assert_eval_eq("BitShiftRight[3^200, 1000]", "0", 0);
    assert_eval_eq("BitShiftRight[-3^200, 1000]", "0", 0);
    assert_eval_eq("BitShiftRight[-3^200, 2^62]", "0", 0);
    assert_eval_eq("BitShiftLeft[3, 70] == 3 2^70", "True", 0);
    assert_eval_eq("BitShiftLeft[-3^50, 300] == -3^50 2^300", "True", 0);
    assert_eval_eq("BitShiftRight[-3^50, 7] == IntegerPart[-3^50 / 2^7]", "True", 0);
    assert_eval_eq("IntegerLength[BitShiftLeft[1, 100000], 2]", "100001", 0);
}

static void test_bitshift_operands_untouched(void) {
    assert_eval_eq("p = 3^200; r = BitShiftLeft[p, 5]; p == 3^200", "True", 0);
    assert_eval_eq("p = -2^90; r = BitShiftRight[p, 3]; {p == -2^90, r == -2^87}",
                   "{True, True}", 0);
    assert_eval_eq("l = {2^90, -2^91}; BitShiftRight[l, 5]; l == {2^90, -2^91}", "True", 0);
    assert_eval_eq("q = 7; BitShiftLeft[q, 70]; q", "7", 0);
}

static void test_bitshift_non_integers(void) {
    /* n * 2^k for any other number; a complex n shifts each part, so integer
     * parts still truncate. */
    assert_eval_eq("BitShiftLeft[2.5, 1]", "5.0", 0);
    assert_eval_eq("BitShiftLeft[2.5]", "5.0", 0);
    assert_eval_eq("BitShiftRight[2.5]", "1.25", 0);
    assert_eval_eq("BitShiftRight[7.5, 1]", "3.75", 0);
    assert_eval_eq("BitShiftRight[-7.5, 1]", "-3.75", 0);
    assert_eval_eq("BitShiftLeft[2.5, -1]", "1.25", 0);
    assert_eval_eq("BitShiftRight[1.5, -3]", "12.0", 0);
    assert_eval_eq("BitShiftLeft[1.5, 1000] == 1.5 * 2.^1000", "True", 0);
    assert_eval_eq("BitShiftLeft[0.0, 5000]", "0.0", 0);
    assert_eval_eq("BitShiftLeft[1/2, 1]", "1", 0);
    assert_eval_eq("BitShiftLeft[7/3, 1]", "14/3", 0);
    assert_eval_eq("BitShiftRight[1/2, 1]", "1/4", 0);
    assert_eval_eq("BitShiftRight[7/2, 1]", "7/4", 0);
    assert_eval_eq("BitShiftRight[-7/2, 1]", "-7/4", 0);
    assert_eval_eq("BitShiftLeft[1/3, -2]", "1/12", 0);
    assert_eval_eq("BitShiftLeft[3 + I, 1]", "6 + 2*I", 0);
    assert_eval_eq("BitShiftRight[3 + 5 I, 1]", "1 + 2*I", 0);
    assert_eval_eq("BitShiftRight[-3 - 5 I, 1]", "-1 - 2*I", 0);
    assert_eval_eq("BitShiftRight[Complex[1/2, 3], 1]", "1/4 + I", 0);
    assert_eval_eq("BitShiftLeft[Complex[1/2, 3], 1]", "1 + 6*I", 0);
    assert_eval_eq("BitShiftRight[1 + I, 1]", "0", 0);
    assert_eval_eq("BitShiftRight[3.5 + 5 I, 1]", "1.75 + 2.5*I", 0);
    /* Arbitrary precision: the scale is exact and keeps the precision. */
    assert_eval_eq("Precision[BitShiftLeft[N[Pi, 30], 2]] == Precision[N[Pi, 30]]", "True", 0);
    assert_eval_eq("BitShiftLeft[N[Pi, 30], 2] == 4 N[Pi, 30]", "True", 0);
    assert_eval_eq("BitShiftRight[N[1, 25], 3] == N[1/8, 25]", "True", 0);
    /* n * 2^k leaves these as they are. */
    assert_eval_eq("BitShiftLeft[Infinity, 1]", "Infinity", 0);
    assert_eval_eq("BitShiftRight[-Infinity, 3]", "-Infinity", 0);
    assert_eval_eq("BitShiftLeft[ComplexInfinity]", "ComplexInfinity", 0);
    assert_eval_eq("BitShiftRight[Indeterminate, 1]", "Indeterminate", 0);
}

/* Evaluate `in` capturing stderr; assert the printed result is `out` and that
 * stderr carries `diag` (or no `::` at all when diag is NULL). */
static void bitshift_expect_diag(const char* in, const char* out, const char* diag) {
    char* result = NULL;
    char* err = eval_capturing_stderr(in, &result);
    ASSERT(result != NULL);
    ASSERT_STR_EQ(result, out);
    if (diag)
        ASSERT_MSG(err != NULL && strstr(err, diag) != NULL,
                   "expected '%s' for %s, got: %s", diag, in, err ? err : "(null)");
    else
        ASSERT_MSG(err == NULL || strstr(err, "::") == NULL,
                   "unexpected diagnostic for %s: %s", in, err ? err : "(null)");
    free(result);
    free(err);
}

static void test_bitshift_messages(void) {
    bitshift_expect_diag("BitShiftLeft[5, 1.5]", "BitShiftLeft[5, 1.5]",
        "BitShiftLeft::intm: Machine-sized integer expected at position 2 in "
        "BitShiftLeft[5, 1.5].");
    bitshift_expect_diag("BitShiftRight[5, 2.]", "BitShiftRight[5, 2.0]",
        "BitShiftRight::intm: Machine-sized integer expected at position 2 in");
    bitshift_expect_diag("BitShiftLeft[1.5, 2.]", "BitShiftLeft[1.5, 2.0]", "BitShiftLeft::intm");
    bitshift_expect_diag("BitShiftLeft[5, 1/2]", "BitShiftLeft[5, 1/2]", "BitShiftLeft::intm");
    bitshift_expect_diag("BitShiftLeft[x, 1.5]", "BitShiftLeft[x, 1.5]", "BitShiftLeft::intm");
    bitshift_expect_diag("BitShiftLeft[0, 2^70]", "BitShiftLeft[0, 1180591620717411303424]",
                         "BitShiftLeft::intm");
    bitshift_expect_diag("BitShiftRight[5, -2^70]",
                         "BitShiftRight[5, -1180591620717411303424]", "BitShiftRight::intm");
    bitshift_expect_diag("BitShiftLeft[5, 1, 2]", "BitShiftLeft[5, 1, 2]",
        "BitShiftLeft::argt: BitShiftLeft called with 3 arguments; 1 or 2 arguments "
        "are expected.");
    bitshift_expect_diag("BitShiftRight[]", "BitShiftRight[]",
        "BitShiftRight::argt: BitShiftRight called with 0 arguments; 1 or 2 arguments "
        "are expected.");
    /* A shift too large to represent overflows instead of aborting. */
    bitshift_expect_diag("BitShiftLeft[5, 2^63 - 1]", "Overflow[]",
                         "General::ovfl: Overflow occurred in computation.");
    bitshift_expect_diag("BitShiftRight[5, -2^63]", "Overflow[]", "General::ovfl");
    bitshift_expect_diag("BitShiftLeft[1/2, 2^63 - 1]", "Overflow[]", "General::ovfl");
    /* A non-integer shifted right that far underflows; a machine real is 0. */
    bitshift_expect_diag("BitShiftRight[3/2, 2^63 - 1]", "Underflow[]",
                         "General::unfl: Underflow occurred in computation.");
    bitshift_expect_diag("BitShiftRight[N[3, 30], 2^63 - 1]", "Underflow[]", "General::unfl");
    bitshift_expect_diag("BitShiftRight[3., 2^63 - 1]", "0.0", NULL);
    bitshift_expect_diag("BitShiftLeft[0., 2^63 - 1]", "0.0", NULL);
    bitshift_expect_diag("BitShiftLeft[0, 2^63 - 1]", "0", NULL);
    /* Symbolic arguments wait, silently. */
    const char* silent[] = { "BitShiftLeft[x, 2]", "BitShiftLeft[5, x]", "BitShiftRight[x]",
                             "BitShiftLeft[Pi, 1]", "BitShiftLeft[True, 1]",
                             "BitShiftRight[f[1], 2]", "BitShiftLeft[x, y]" };
    for (size_t i = 0; i < sizeof silent / sizeof silent[0]; i++)
        bitshift_expect_diag(silent[i], silent[i], NULL);
    /* Every diagnostic honours Quiet and Check. */
    const char* noisy[] = { "BitShiftLeft[5, 1.5]", "BitShiftLeft[5, 1, 2]",
                            "BitShiftRight[]", "BitShiftLeft[5, 2^63 - 1]",
                            "BitShiftRight[3/2, 2^63 - 1]",
                            "BitShiftLeft[{1, 2}, {1, 2, 3}]" };
    for (size_t i = 0; i < sizeof noisy / sizeof noisy[0]; i++) {
        char src[128], *result = NULL, *err;
        snprintf(src, sizeof src, "Quiet[%s]", noisy[i]);
        err = eval_capturing_stderr(src, &result);
        ASSERT_MSG(err == NULL || strstr(err, "::") == NULL,
                   "Quiet leaked for %s: %s", noisy[i], err ? err : "(null)");
        free(result);
        free(err);
        snprintf(src, sizeof src, "Check[%s, failed]", noisy[i]);
        assert_eval_eq(src, "failed", 0);
    }
    assert_eval_eq("Check[BitShiftLeft[x, 2], failed]", "BitShiftLeft[x, 2]", 0);
}

static void test_bitshift_listable_shapes(void) {
    assert_eval_eq("BitShiftLeft[{1, 2}, {1, 2}]", "{2, 8}", 0);
    assert_eval_eq("BitShiftLeft[5, {0, 1, -1, 70}]", "{5, 10, 2, 5902958103587056517120}", 0);
    assert_eval_eq("BitShiftLeft[{1, x, 2.5}, 3]", "{8, BitShiftLeft[x, 3], 20.0}", 0);
    assert_eval_eq("BitShiftRight[{{1, 2}, {3, 4}}, {0, 1}]", "{{1, 2}, {1, 2}}", 0);
    assert_eval_eq("BitShiftLeft[{}, 3]", "{}", 0);
    assert_eval_eq("BitShiftLeft[{1, {2, 3}}]", "{2, {4, 6}}", 0);
    assert_eval_eq("BitShiftLeft[{1.5, 2}, 1]", "{3.0, 4}", 0);
    bitshift_expect_diag("BitShiftLeft[{1, 2}, {1, 2, 3}]", "BitShiftLeft[{1, 2}, {1, 2, 3}]",
        "Thread::tdlen: Objects of unequal length in BitShiftLeft[{1, 2}, {1, 2, 3}] "
        "cannot be combined.");
}

static void test_bitshift_properties(void) {
    /* Independent oracles over every sign and size: n 2^k, IntegerPart, the
     * two heads as each other's negation, and the round trip. */
    assert_eval_eq(
        "SeedRandom[41]; ns = Join[RandomInteger[{-2^200, 2^200}, 30], "
        "RandomInteger[{-2^62, 2^62}, 30], "
        "{-2^63, 2^63 - 1, 2^63, -2^63 - 1, 0, -1, 1, 3^100, -5^90}]; "
        "ks = {0, 1, 5, 31, 32, 61, 62, 63, 64, 65, 127, 128, 300, 1000}; "
        "And @@ Flatten[Table[BitShiftLeft[n, k] == n 2^k && "
        "BitShiftRight[n, k] == IntegerPart[n / 2^k] && "
        "BitShiftLeft[n, -k] == BitShiftRight[n, k] && "
        "BitShiftRight[n, -k] == BitShiftLeft[n, k] && "
        "BitShiftRight[BitShiftLeft[n, k], k] == n, "
        "{n, ns}, {k, ks}]]", "True", 0);
    /* Truncation toward zero is odd symmetry; floor would not be. */
    assert_eval_eq(
        "SeedRandom[42]; And @@ Flatten[Table[BitShiftRight[-n, k] == -BitShiftRight[n, k], "
        "{n, RandomInteger[{0, 2^100}, 50]}, {k, {0, 1, 3, 63, 64, 99}}]]",
        "True", 0);
}

static void test_bitshift_packed_and_ndarray(void) {
    assert_eval_eq("PackedArrayQ[BitShiftLeft[Range[10], 2]]", "True", 0);
    assert_eval_eq("PackedArrayQ[BitShiftRight[Range[-10, 10], 2]]", "True", 0);
    assert_eval_eq("PackedArrayQ[BitShiftLeft[3, Range[0, 9]]]", "True", 0);
    assert_eval_eq("PackedArrayQ[BitShiftRight[Range[8], Range[0, 7]]]", "True", 0);
    assert_eval_eq("PackedArrayQ[BitShiftLeft[Range[8]]]", "True", 0);
    assert_eval_eq("PackedArrayQ[BitShiftRight[-Range[8]]]", "True", 0);
    assert_eval_eq("PackedArrayQ[BitShiftLeft[Partition[Range[12], 3], -1]]", "True", 0);
    assert_eval_eq("BitShiftLeft[3, Range[-2, 4]]", "{0, 1, 3, 6, 12, 24, 48}", 0);
    assert_eval_eq("BitShiftRight[-Range[8]]", "{0, -1, -1, -2, -2, -3, -3, -4}", 0);
    assert_eval_eq("BitShiftLeft[NDArray[{1, 2, 3}, DataType -> \"int64\"], 2]",
                   "NDArray[{4, 8, 12}]", 0);
    assert_eval_eq("BitShiftRight[NDArray[{-5, 7}, DataType -> \"int64\"]]", "NDArray[{-2, 3}]", 0);
    assert_eval_eq("BitShiftLeft[NDArray[{-5, 7}, DataType -> \"int64\"]]", "NDArray[{-10, 14}]", 0);
    assert_eval_eq("BitShiftRight[5, NDArray[{0, 1, -1}, DataType -> \"int64\"]]",
                   "NDArray[{5, 2, 10}]", 0);
    assert_eval_eq("BitShiftLeft[NDArray[{1, -2, 3}, DataType -> \"int64\"], "
                   "NDArray[{0, 1, -1}, DataType -> \"int64\"]]", "NDArray[{1, -4, 1}]", 0);
    assert_eval_eq("BitShiftRight[NDArray[{{-9, 9}, {17, -17}}, DataType -> \"int64\"], 3]",
                   "NDArray[{{-1, 1}, {2, -2}}]", 0);
    /* Every right shift is an int64, however far. */
    assert_eval_eq("BitShiftRight[{-1, 0, 5, -5, -2^63}, 10^6]", "{0, 0, 0, 0, 0}", 0);
    assert_eval_eq("PackedArrayQ[BitShiftRight[-Range[5], 10^6]]", "True", 0);
    assert_eval_eq("BitShiftLeft[{0, 0}, 10^6]", "{0, 0}", 0);
    /* A real buffer: the List path answers each element as n 2^k. */
    assert_eval_eq("BitShiftLeft[NDArray[{1.5, 2.}], 1]", "{3.0, 4.0}", 0);
    assert_eval_eq("BitShiftLeft[{1.5, 2.5}, 1]", "{3.0, 5.0}", 0);
}

static void test_bitshift_buffer_declines(void) {
    /* A left shift that leaves int64 abandons the buffer -- never wraps. */
    assert_eval_eq("BitShiftLeft[{1, 2^62}, 2]", "{4, 18446744073709551616}", 0);
    assert_eval_eq("BitShiftLeft[NDArray[{1, 2^62}, DataType -> \"int64\"], 2]",
                   "{4, 18446744073709551616}", 0);
    assert_eval_eq("BitShiftLeft[{1, -1, 0}, 63]", "{9223372036854775808, -9223372036854775808, 0}", 0);
    assert_eval_eq("BitShiftLeft[{-1, 0}, 63]", "{-9223372036854775808, 0}", 0);
    /* Past 62 the constant arm declines, and the per-element loop still
     * answers a buffer of -1s and 0s in int64. */
    assert_eval_eq("PackedArrayQ[BitShiftLeft[Join[0 Range[5] - 1, 0 Range[5]], 63]]", "True", 0);
    assert_eval_eq("BitShiftRight[{1, 3}, -63]", "{9223372036854775808, 27670116110564327424}", 0);
    assert_eval_eq("BitShiftLeft[3, {1, 62, 2}]", "{6, 13835058055282163712, 12}", 0);
    assert_eval_eq("BitShiftLeft[{1, 2}, {1, 63}]", "{2, 18446744073709551616}", 0);
    assert_eval_eq("BitShiftRight[{1, 2}, {-1, -63}]", "{2, 18446744073709551616}", 0);
    assert_eval_eq("BitShiftLeft[{2^62, 1}, {1, -1}]", "{9223372036854775808, 0}", 0);
    assert_eval_eq("BitShiftLeft[5, {-70, 70}]", "{0, 5902958103587056517120}", 0);
    /* A count of -2^63 (BitShiftRight) is a left shift by 2^63. */
    assert_eval_eq("BitShiftRight[{0, 0}, -2^63]", "{0, 0}", 0);
    assert_eval_eq("Quiet[BitShiftRight[{0, 1}, -2^63]]", "{0, Overflow[]}", 0);
}

static void test_bitshift_surfaces_agree(void) {
    /* List == packed == NDArray == Compile (array and scalar), for every arm:
     * scalar k, scalar n, two buffers, one argument -- with counts of both
     * signs, so the fast and the general per-element passes both run. */
    assert_eval_eq(
        "SeedRandom[43]; na = Join[RandomInteger[{-2^40, 2^40}, 500], "
        "{-2^63, 2^63 - 1, 0, -1, 1}]; ka = RandomInteger[{-70, 20}, Length[na]]; "
        "kp = RandomInteger[{0, 20}, Length[na]]; "
        "nd[x_] := NDArray[x, DataType -> \"int64\"]; "
        "refL[a_, b_] := Table[If[b[[i]] >= 0, a[[i]] 2^b[[i]], IntegerPart[a[[i]] / 2^-b[[i]]]], "
        "  {i, Length[a]}]; "
        "refR[a_, b_] := refL[a, -b]; "
        "sm = Select[na, Abs[#] < 2^40 &]; ks = RandomInteger[{-70, 20}, Length[sm]]; "
        "kq = Take[kp, Length[sm]]; "
        "Flatten[{"
        " BitShiftLeft[sm, ks] === refL[sm, ks], "
        " Normal[BitShiftLeft[nd[sm], nd[ks]]] === refL[sm, ks], "
        " Compile[{{v, _Integer, 1}, {w, _Integer, 1}}, BitShiftLeft[v, w]][sm, ks] === refL[sm, ks], "
        " BitShiftRight[sm, ks] === refR[sm, ks], "
        " Normal[BitShiftRight[nd[sm], nd[ks]]] === refR[sm, ks], "
        " Compile[{{v, _Integer, 1}, {w, _Integer, 1}}, BitShiftRight[v, w]][sm, ks] === refR[sm, ks], "
        " BitShiftRight[na, kp] === refR[na, kp], "
        " Normal[BitShiftRight[nd[na], nd[kp]]] === refR[na, kp], "
        " BitShiftLeft[sm, kq] === refL[sm, kq], "
        " BitShiftLeft[-12345, ka] === refL[Table[-12345, Length[ka]], ka], "
        " Normal[BitShiftLeft[-12345, nd[ka]]] === refL[Table[-12345, Length[ka]], ka], "
        " BitShiftRight[-12345, ka] === refR[Table[-12345, Length[ka]], ka], "
        " Normal[BitShiftRight[-12345, nd[ka]]] === refR[Table[-12345, Length[ka]], ka], "
        " Compile[{{v, _Integer, 1}}, BitShiftRight[-12345, v]][ka] === refR[Table[-12345, Length[ka]], ka], "
        " Table[BitShiftRight[na, s] === refR[na, Table[s, Length[na]]], {s, {0, 1, 7, 62, 63, 64, -1, -20}}], "
        " Table[Normal[BitShiftRight[nd[na], s]] === refR[na, Table[s, Length[na]]], {s, {0, 1, 7, 62, 63, 64}}], "
        " Table[BitShiftLeft[sm, s] === refL[sm, Table[s, Length[sm]]], {s, {0, 1, 7, 22, 23, -1, -63, -64}}], "
        " Table[Compile[{{v, _Integer, 1}}, BitShiftLeft[v, 7]][sm] === refL[sm, Table[7, Length[sm]]], {1}], "
        " BitShiftLeft[sm] === refL[sm, Table[1, Length[sm]]], "
        " BitShiftRight[na] === refR[na, Table[1, Length[na]]], "
        " Normal[BitShiftRight[nd[na]]] === refR[na, Table[1, Length[na]]], "
        " Compile[{{v, _Integer, 1}}, BitShiftRight[v]][na] === refR[na, Table[1, Length[na]]], "
        " Table[Compile[{{s, _Integer}, {t, _Integer}}, BitShiftLeft[s, t]][sm[[i]], ks[[i]]], {i, Length[sm]}] === refL[sm, ks], "
        " Table[BitShiftRight[FromPackedArray[sm][[i]], ks[[i]]], {i, Length[sm]}] === refR[sm, ks]"
        "}] // Union",
        "{True}", 0);
}

static void test_bitshift_compile(void) {
    assert_eval_eq("Compile[{{n, _Integer}, {k, _Integer}}, BitShiftLeft[n, k]][5, 3]", "40", 0);
    assert_eval_eq("Compile[{{n, _Integer}, {k, _Integer}}, BitShiftRight[n, k]][-5, 3]", "0", 0);
    assert_eval_eq("Compile[{{n, _Integer}, {k, _Integer}}, BitShiftRight[n, k]][-7, 2]", "-1", 0);
    assert_eval_eq("Compile[{{n, _Integer}, {k, _Integer}}, BitShiftRight[n, k]][5, -2]", "20", 0);
    assert_eval_eq("Compile[{{n, _Integer}, {k, _Integer}}, BitShiftLeft[n, k]][5, -2]", "1", 0);
    assert_eval_eq("Compile[{{n, _Integer}, {k, _Integer}}, BitShiftRight[n, k]][-1, 64]", "0", 0);
    assert_eval_eq("Compile[{{n, _Integer}}, BitShiftLeft[n]][-3]", "-6", 0);
    assert_eval_eq("Compile[{{n, _Integer}}, BitShiftRight[n]][-3]", "-1", 0);
    /* No int64 answer: the VM bails and the interpreter answers exactly. */
    assert_eval_eq("Compile[{{n, _Integer}, {k, _Integer}}, BitShiftLeft[n, k]][1, 64]",
                   "18446744073709551616", 0);
    assert_eval_eq("Compile[{{n, _Integer}}, BitShiftLeft[n, 70]][3]", "3541774862152233910272", 0);
    assert_eval_eq("Compile[{{v, _Integer, 1}}, BitShiftLeft[v, 2]][Range[5]]", "{4, 8, 12, 16, 20}", 0);
    assert_eval_eq("Compile[{{n, _Integer}}, Sum[BitShiftRight[n, k], {k, 0, 10}]][1000]", "1994", 0);
    assert_eval_eq("Compile[{{x, _Real}}, BitShiftLeft[x, 1]][2.5]", "5.0", 0);
}

static void test_bitshift_compile_diagnostics(void) {
    const char* in[] = {
        "CompileDiagnostics[{{n, _Integer}, {k, _Integer}}, BitShiftLeft[n, k]]",
        "CompileDiagnostics[{{n, _Integer}, {k, _Integer}}, BitShiftRight[n, k]]",
        "CompileDiagnostics[{{n, _Integer}}, BitShiftLeft[n]]",
        "CompileDiagnostics[{{n, _Integer}}, BitShiftRight[n]]",
        "CompileDiagnostics[{{v, _Integer, 1}}, BitShiftLeft[v, 3]]",
        "CompileDiagnostics[{{v, _Integer, 1}}, BitShiftRight[v, 3]]",
        "CompileDiagnostics[{{v, _Integer, 1}}, BitShiftLeft[7, v]]",
        "CompileDiagnostics[{{v, _Integer, 1}}, BitShiftRight[v]]",
        "CompileDiagnostics[{{v, _Integer, 1}, {w, _Integer, 1}}, BitShiftLeft[v, w]]",
        "CompileDiagnostics[{{v, _Integer, 1}, {w, _Integer, 1}}, BitShiftRight[v, w]]" };
    for (size_t i = 0; i < sizeof in / sizeof in[0]; i++) {
        char src[256];
        snprintf(src, sizeof src, "\"Compiled\" /. %s", in[i]);
        assert_eval_eq(src, "True", 0);
    }
}

static void test_bitshift_autocompile(void) {
    const char* tbl = "Table[x + BitShiftRight[Floor[x], 2], {x, 0., 2000., 0.5}]";
    const char* sm  = "Sum[Sin[x] BitShiftLeft[Floor[x], -1], {x, 0., 500., 0.25}]";
    assert_eval_eq("Table[x + BitShiftLeft[Floor[x], 1], {x, 0., 2., 0.5}]",
                   "{0.0, 0.5, 3.0, 3.5, 6.0}", 0);
    autocompile_set_enabled(false);
    Expr* t_off = bitnot_eval_owned(tbl);
    Expr* s_off = bitnot_eval_owned(sm);
    autocompile_set_enabled(true);
    Expr* t_on = bitnot_eval_owned(tbl);
    Expr* s_on = bitnot_eval_owned(sm);
    ASSERT(expr_eq(t_off, t_on));
    ASSERT(expr_eq(s_off, s_on));
    ASSERT(s_on->type == EXPR_REAL);
    expr_free(t_off); expr_free(s_off); expr_free(t_on); expr_free(s_on);
}

static void test_bitshift_attributes_docstrings(void) {
    assert_eval_eq("Attributes[BitShiftLeft]", "{Listable, Protected}", 0);
    assert_eval_eq("Attributes[BitShiftRight]", "{Listable, Protected}", 0);
    SymbolDef* def = symtab_get_def("BitShiftLeft");
    ASSERT(def != NULL && def->docstring != NULL);
    ASSERT(strstr(def->docstring, "shifts the binary bits in the integer n to the left") != NULL);
    def = symtab_get_def("BitShiftRight");
    ASSERT(def != NULL && def->docstring != NULL);
    ASSERT(strstr(def->docstring, "IntegerPart[n / 2^k]") != NULL);
    ASSERT(SYM_BitShiftLeft != NULL && strcmp(SYM_BitShiftLeft, "BitShiftLeft") == 0);
    ASSERT(SYM_BitShiftRight != NULL && strcmp(SYM_BitShiftRight, "BitShiftRight") == 0);
}

static void test_bitshift_stress(void) {
    for (int k = 0; k < 50; k++) {
        assert_eval_eq("BitShiftLeft[32, 3]", "256", 0);
        assert_eval_eq("BitShiftRight[-5, 1]", "-2", 0);
        assert_eval_eq("BitShiftLeft[2^100, -100]", "1", 0);
        assert_eval_eq("BitShiftLeft[3, 70] == 3 2^70", "True", 0);
        assert_eval_eq("BitShiftRight[3 + 5 I, 1]", "1 + 2*I", 0);
        assert_eval_eq("BitShiftLeft[x, 1]", "BitShiftLeft[x, 1]", 0);
        assert_eval_eq("BitShiftLeft[Range[4], {1, -1, 2, 70}]",
                       "{2, 1, 12, 4722366482869645213696}", 0);
        assert_eval_eq("BitShiftRight[-Range[4]]", "{0, -1, -1, -2}", 0);
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

    TEST(test_bitxor_basic);
    TEST(test_bitxor_range_threading);
    TEST(test_bitxor_truth_table);
    TEST(test_bitxor_empty_lists);
    TEST(test_bitxor_listable_shapes);
    TEST(test_bitxor_bignums);
    TEST(test_bitxor_bignum_operands_untouched);
    TEST(test_bitxor_pair_cancellation);
    TEST(test_bitxor_integer_folding);
    TEST(test_bitxor_bitnot);
    TEST(test_bitxor_orderless_flat);
    TEST(test_bitxor_properties);
    TEST(test_bitxor_packed_and_ndarray);
    TEST(test_bitxor_surfaces_agree);
    TEST(test_bitxor_compile);
    TEST(test_bitxor_compile_diagnostics);
    TEST(test_bitxor_real_int_diagnostic);
    TEST(test_bitxor_rational_complex_diagnostic);
    TEST(test_bitxor_quiet_check);
    TEST(test_bitxor_symbolic_silent);
    TEST(test_bitxor_attributes);
    TEST(test_bitxor_docstring_and_sym);
    TEST(test_bitxor_stress);

    TEST(test_bitand_documented_examples);
    TEST(test_bitand_basic);
    TEST(test_bitand_listable_shapes);
    TEST(test_bitand_bignums);
    TEST(test_bitand_bignum_operands_untouched);
    TEST(test_bitand_idempotence);
    TEST(test_bitand_integer_folding);
    TEST(test_bitand_zero_absorbs);
    TEST(test_bitand_minus_one_identity);
    TEST(test_bitand_no_extra_rules);
    TEST(test_bitand_orderless_flat);
    TEST(test_bitand_non_integers_silent);
    TEST(test_bitand_ragged_tdlen);
    TEST(test_bitand_properties);
    TEST(test_bitand_packed_and_ndarray);
    TEST(test_bitand_surfaces_agree);
    TEST(test_bitand_compile);
    TEST(test_bitand_compile_diagnostics);
    TEST(test_bitand_attributes);
    TEST(test_bitand_docstring_and_sym);
    TEST(test_bitand_stress);

    TEST(test_bitor_documented_examples);
    TEST(test_bitor_basic);
    TEST(test_bitor_listable_shapes);
    TEST(test_bitor_bignums);
    TEST(test_bitor_bignum_operands_untouched);
    TEST(test_bitor_idempotence);
    TEST(test_bitor_integer_folding);
    TEST(test_bitor_minus_one_absorbs);
    TEST(test_bitor_zero_identity);
    TEST(test_bitor_no_extra_rules);
    TEST(test_bitor_orderless_flat);
    TEST(test_bitor_non_integers_silent);
    TEST(test_bitor_ragged_tdlen);
    TEST(test_bitor_properties);
    TEST(test_bitor_packed_and_ndarray);
    TEST(test_bitor_surfaces_agree);
    TEST(test_bitor_compile);
    TEST(test_bitor_compile_diagnostics);
    TEST(test_bitor_attributes);
    TEST(test_bitor_docstring_and_sym);
    TEST(test_bitor_stress);

    TEST(test_bitnot_documented_examples);
    TEST(test_bitnot_basic);
    TEST(test_bitnot_int64_extremes);
    TEST(test_bitnot_bignums);
    TEST(test_bitnot_bignum_operands_untouched);
    TEST(test_bitnot_involution);
    TEST(test_bitnot_no_extra_rules);
    TEST(test_bitnot_non_integers_silent);
    TEST(test_bitnot_argx);
    TEST(test_bitnot_listable_shapes);
    TEST(test_bitnot_properties);
    TEST(test_bitnot_packed_and_ndarray);
    TEST(test_bitnot_surfaces_agree);
    TEST(test_bitnot_compile);
    TEST(test_bitnot_compile_diagnostics);
    TEST(test_bitnot_autocompile);
    TEST(test_bitnot_attributes);
    TEST(test_bitnot_docstring_and_sym);
    TEST(test_bitnot_stress);

    TEST(test_bitset_documented_examples);
    TEST(test_bitset_basic);
    TEST(test_bitset_negatives);
    TEST(test_bitset_int64_boundary);
    TEST(test_bitset_bignums);
    TEST(test_bitset_bignum_operands_untouched);
    TEST(test_bitset_unevaluated_silent);
    TEST(test_bitset_listable_shapes);
    TEST(test_bitset_ragged_tdlen);
    TEST(test_bitset_properties);
    TEST(test_bitset_packed_and_ndarray);
    TEST(test_bitset_buffer_declines);
    TEST(test_bitset_surfaces_agree);
    TEST(test_bitset_compile);
    TEST(test_bitset_compile_diagnostics);
    TEST(test_bitset_autocompile);
    TEST(test_bitset_attributes);
    TEST(test_bitset_docstring_and_sym);
    TEST(test_bitset_stress);
    TEST(test_bitclear_basic);
    TEST(test_bitclear_negatives_and_bignums);
    TEST(test_bitclear_unevaluated_silent);
    TEST(test_bitclear_properties);
    TEST(test_bitclear_packed_and_ndarray);
    TEST(test_bitclear_surfaces_agree);
    TEST(test_bitclear_compile);
    TEST(test_bitclear_attributes_docstring);
    TEST(test_bitget_documented_examples);
    TEST(test_bitget_basic);
    TEST(test_bitget_negatives);
    TEST(test_bitget_int64_boundary);
    TEST(test_bitget_bignums);
    TEST(test_bitget_bignum_operands_untouched);
    TEST(test_bitget_unevaluated_silent);
    TEST(test_bitget_listable_shapes);
    TEST(test_bitget_ragged_tdlen);
    TEST(test_bitget_properties);
    TEST(test_bitget_packed_and_ndarray);
    TEST(test_bitget_buffer_declines);
    TEST(test_bitget_surfaces_agree);
    TEST(test_bitget_compile);
    TEST(test_bitget_compile_diagnostics);
    TEST(test_bitget_autocompile);
    TEST(test_bitget_attributes);
    TEST(test_bitget_docstring_and_sym);
    TEST(test_bitget_stress);

    TEST(test_bitshift_documented_examples);
    TEST(test_bitshift_basic);
    TEST(test_bitshift_negatives_truncate);
    TEST(test_bitshift_int64_boundary);
    TEST(test_bitshift_bignums);
    TEST(test_bitshift_operands_untouched);
    TEST(test_bitshift_non_integers);
    TEST(test_bitshift_messages);
    TEST(test_bitshift_listable_shapes);
    TEST(test_bitshift_properties);
    TEST(test_bitshift_packed_and_ndarray);
    TEST(test_bitshift_buffer_declines);
    TEST(test_bitshift_surfaces_agree);
    TEST(test_bitshift_compile);
    TEST(test_bitshift_compile_diagnostics);
    TEST(test_bitshift_autocompile);
    TEST(test_bitshift_attributes_docstrings);
    TEST(test_bitshift_stress);

    printf("All BitLength / BitXor / BitAnd / BitOr / BitNot / BitSet / BitGet / BitShiftLeft / BitShiftRight tests passed!\n");
    return 0;
}
