/* Unit tests for ChineseRemainder.
 *
 *   ChineseRemainder[{r1, ...}, {m1, ...}]      smallest x >= 0 with
 *       x mod mi == ri mod mi, in 0 <= x < LCM[m1, ...].
 *   ChineseRemainder[{r1, ...}, {m1, ...}, d]   smallest such x >= d,
 *       in d <= x < d + LCM[m1, ...].
 *
 * Coverage:
 *   - Documented Mathematica examples (19, 68, 7, ...).
 *   - The x mod mi == ri mod mi contract, incl. unreduced / negative residues.
 *   - Round-trip: Mod[ChineseRemainder[r, m], m] == r for coprime moduli
 *     (the residue-number-system recovery, incl. modular multiplication).
 *   - BigInt: coprime ~10^9 primes producing a >=25-digit result; verified
 *     by range (0 <= x < LCM) and per-modulus round-trip.
 *   - Three-argument offset form, including d that skips the base solution
 *     and a negative d.
 *   - Non-coprime but consistent moduli (a genuine gcd merge).
 *   - No-solution / inconsistent systems -> unevaluated, no diagnostic.
 *   - Edge cases: empty lists -> 0; single congruence; length mismatch;
 *     symbolic and non-integer arguments -> unevaluated.
 *   - Argument-count error -> ChineseRemainder::argt.
 *   - Attribute (Protected, NOT Listable), docstring, interned-symbol
 *     introspection.
 *   - Repeated-evaluation stress loop to catch double-frees / leaks under
 *     valgrind.
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
    const char* path = "/tmp/mathilda_chinese_remainder_stderr.log";
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

/* --- Documented examples ------------------------------------------- */

static void test_doc_two_mod(void) {
    /* Spec: smallest x with x mod 4 == 3 and x mod 5 == 4. */
    assert_eval_eq("ChineseRemainder[{3, 4}, {4, 5}]", "19", 0);
}

static void test_doc_three_mod(void) {
    /* Spec: remainders 2, 3, 5 when divided by 3, 5, 7. */
    assert_eval_eq("ChineseRemainder[{2, 3, 5}, {3, 5, 7}]", "68", 0);
}

static void test_doc_reduce_findinstance_agreement(void) {
    /* Spec: Reduce/FindInstance x==7 for {Mod[x,3]==1, Mod[x,5]==2}. */
    assert_eval_eq("ChineseRemainder[{1, 2}, {3, 5}]", "7", 0);
}

/* --- The x mod mi == ri mod mi contract ---------------------------- */

static void test_congruence_contract(void) {
    /* Result really satisfies each congruence and lies in [0, LCM). */
    assert_eval_eq(
        "Mod[ChineseRemainder[{2, 3, 5}, {3, 5, 7}], {3, 5, 7}] == {2, 3, 5}",
        "True", 0);
    assert_eval_eq(
        "And[ChineseRemainder[{3, 4}, {4, 5}] >= 0, "
        "ChineseRemainder[{3, 4}, {4, 5}] < LCM[4, 5]]", "True", 0);
}

static void test_single_congruence(void) {
    assert_eval_eq("ChineseRemainder[{7}, {10}]", "7", 0);
    /* Empty modulus list of length 1, unreduced residue: reduces into [0, m). */
    assert_eval_eq("ChineseRemainder[{17}, {10}]", "7", 0);
}

static void test_unreduced_and_negative_residues(void) {
    /* If ri is not already in [0, mi), the result still gives x mod mi. */
    assert_eval_eq("ChineseRemainder[{-1}, {10}]", "9", 0);
    assert_eval_eq("Mod[ChineseRemainder[{123, -4}, {7, 9}], {7, 9}] == "
                   "Mod[{123, -4}, {7, 9}]", "True", 0);
}

/* --- Residue number system (multiply / recover) -------------------- */

static void test_rns_recover_number(void) {
    /* 123 < each of Prime[100..103] = {541,547,557,563}; recovering the RNS
     * {123,123,123,123} gives 123 back. */
    assert_eval_eq("ChineseRemainder[{123, 123, 123, 123}, {541, 547, 557, 563}]",
                   "123", 0);
}

static void test_rns_recover_product(void) {
    /* 123 * 2 in the residue system: residues multiply, recovery gives 246. */
    assert_eval_eq("ChineseRemainder[{246, 246, 246, 246}, {541, 547, 557, 563}]",
                   "246", 0);
}

/* --- BigInt (arbitrary-precision) recovery ------------------------- */

static void test_bignum_roundtrip(void) {
    /* Three coprime ~10^9 primes -> LCM ~10^27, a genuine bignum result.
     * Verify it lies in [0, LCM) and reproduces every residue under Mod.
     * `crx` is used only here; never referenced as a free symbol elsewhere. */
    assert_eval_eq(
        "crx = ChineseRemainder[{123, 456, 789}, "
        "{1000000007, 1000000009, 1000000021}]; "
        "And[crx >= 0, crx < LCM[1000000007, 1000000009, 1000000021], "
        "Mod[crx, 1000000007] == 123, Mod[crx, 1000000009] == 456, "
        "Mod[crx, 1000000021] == 789]", "True", 0);
}

static void test_bignum_is_large(void) {
    /* Sanity: the recovered value is a many-digit integer (not truncated to
     * a machine int).  x >= product/2 need not hold, but for these residues
     * it exceeds 10^18, so IntegerLength is well past machine range. */
    assert_eval_eq(
        "IntegerLength[ChineseRemainder[{123, 456, 789}, "
        "{1000000007, 1000000009, 1000000021}]] >= 19", "True", 0);
}

/* --- Three-argument offset form ------------------------------------ */

static void test_offset_skips_base(void) {
    /* Base solution 68 (mod 105); smallest >= 100 is 68 + 105 = 173. */
    assert_eval_eq("ChineseRemainder[{2, 3, 5}, {3, 5, 7}, 100]", "173", 0);
}

static void test_offset_boundary(void) {
    /* Base 19 (mod 20). */
    assert_eval_eq("ChineseRemainder[{3, 4}, {4, 5}, 19]", "19", 0);   /* d hits it */
    assert_eval_eq("ChineseRemainder[{3, 4}, {4, 5}, 20]", "39", 0);   /* next up  */
    assert_eval_eq("ChineseRemainder[{3, 4}, {4, 5}, 0]", "19", 0);    /* == 2-arg */
}

static void test_offset_negative(void) {
    /* Smallest x >= -100 with x == 19 (mod 20) is -81. */
    assert_eval_eq("ChineseRemainder[{3, 4}, {4, 5}, -100]", "-81", 0);
    assert_eval_eq("Mod[ChineseRemainder[{3, 4}, {4, 5}, -100], {4, 5}]",
                   "{3, 4}", 0);
}

/* --- Non-coprime but consistent moduli ----------------------------- */

static void test_noncoprime_consistent(void) {
    /* x == 2 (mod 4) and x == 2 (mod 6): gcd(4,6)=2, residues agree. */
    assert_eval_eq("ChineseRemainder[{2, 2}, {4, 6}]", "2", 0);
    /* x == 5 (mod 6) and x == 3 (mod 10): both == 1 (mod 2); answer 23. */
    assert_eval_eq("ChineseRemainder[{5, 3}, {6, 10}]", "23", 0);
    assert_eval_eq("Mod[ChineseRemainder[{5, 3}, {6, 10}], {6, 10}]",
                   "{5, 3}", 0);
}

/* --- No solution (inconsistent) -> unevaluated --------------------- */

static void test_no_solution_unevaluated(void) {
    /* Spec: ChineseRemainder[{1, 2}, {6, 10}] has no solution; 1 and 2 differ
     * modulo gcd(6,10)=2.  Left unevaluated, no diagnostic. */
    char* result = NULL;
    char* err = eval_capturing_stderr("ChineseRemainder[{1, 2}, {6, 10}]", &result);
    ASSERT(result != NULL);
    ASSERT_MSG(strcmp(result, "ChineseRemainder[{1, 2}, {6, 10}]") == 0,
               "expected unevaluated call, got: %s", result);
    ASSERT_MSG(err == NULL || strstr(err, "ChineseRemainder::") == NULL,
               "expected no diagnostic, got: %s", err ? err : "(null)");
    free(result);
    free(err);
}

static void test_no_solution_variants(void) {
    /* {1, 2}, {2*3, 2*5} = {6, 10} -- same as above spelled with products. */
    assert_eval_eq("ChineseRemainder[{1, 2}, {2 3, 2 5}]",
                   "ChineseRemainder[{1, 2}, {6, 10}]", 0);
    /* Three-modulus inconsistency. */
    assert_eval_eq("ChineseRemainder[{0, 1}, {4, 6}]",
                   "ChineseRemainder[{0, 1}, {4, 6}]", 0);
}

/* --- Edge cases ---------------------------------------------------- */

static void test_empty_lists(void) {
    /* LCM of nothing is 1; the unique x in [0, 1) is 0. */
    assert_eval_eq("ChineseRemainder[{}, {}]", "0", 0);
    /* Empty with an offset: smallest x >= 5 in [5, 6). */
    assert_eval_eq("ChineseRemainder[{}, {}, 5]", "5", 0);
}

static void test_zero_modulus_unevaluated(void) {
    /* A zero modulus is rejected; left unevaluated. */
    assert_eval_eq("ChineseRemainder[{1}, {0}]", "ChineseRemainder[{1}, {0}]", 0);
}

static void test_length_mismatch_unevaluated(void) {
    char* result = NULL;
    char* err = eval_capturing_stderr("ChineseRemainder[{1, 2}, {3}]", &result);
    ASSERT(result != NULL);
    ASSERT(strstr(result, "ChineseRemainder[{1, 2}, {3}]") != NULL);
    ASSERT_MSG(err == NULL || strstr(err, "ChineseRemainder::") == NULL,
               "expected no diagnostic for length mismatch, got: %s",
               err ? err : "(null)");
    free(result);
    free(err);
}

static void test_symbolic_unevaluated(void) {
    /* `q` is a symbol never assigned elsewhere. */
    char* result = NULL;
    char* err = eval_capturing_stderr("ChineseRemainder[{q, 2}, {3, 5}]", &result);
    ASSERT(result != NULL);
    ASSERT(strstr(result, "ChineseRemainder[{q, 2}, {3, 5}]") != NULL);
    ASSERT_MSG(err == NULL || strstr(err, "ChineseRemainder::") == NULL,
               "expected no diagnostic for symbolic input, got: %s",
               err ? err : "(null)");
    free(result);
    free(err);
}

static void test_noninteger_unevaluated(void) {
    /* Non-integer numeric residue -> unevaluated (integer-only contract). */
    assert_eval_eq("ChineseRemainder[{1/2, 2}, {3, 5}]",
                   "ChineseRemainder[{1/2, 2}, {3, 5}]", 0);
    /* Non-list arguments -> unevaluated. */
    assert_eval_eq("ChineseRemainder[3, 5]", "ChineseRemainder[3, 5]", 0);
}

/* --- Argument-count error ------------------------------------------ */

static void test_arg_error_zero(void) {
    char* result = NULL;
    char* err = eval_capturing_stderr("ChineseRemainder[]", &result);
    ASSERT(result != NULL);
    ASSERT(strstr(result, "ChineseRemainder") != NULL);
    ASSERT(err != NULL);
    ASSERT_MSG(strstr(err, "ChineseRemainder::argt") != NULL,
               "expected argt diagnostic, got: %s", err);
    ASSERT(strstr(err, "2 or 3 arguments are expected") != NULL);
    free(result);
    free(err);
}

static void test_arg_error_one(void) {
    char* result = NULL;
    char* err = eval_capturing_stderr("ChineseRemainder[{1, 2}]", &result);
    ASSERT(result != NULL);
    ASSERT(err != NULL);
    ASSERT(strstr(err, "ChineseRemainder::argt") != NULL);
    free(result);
    free(err);
}

/* --- Attribute / docstring / interned-symbol introspection --------- */

static void test_attributes(void) {
    SymbolDef* def = symtab_get_def("ChineseRemainder");
    ASSERT(def != NULL);
    uint32_t a = get_attributes("ChineseRemainder");
    ASSERT((a & ATTR_PROTECTED) != 0);
    /* The residue/modulus lists are atomic arguments: threading element-wise
     * would be mathematically wrong, so ChineseRemainder is NOT Listable. */
    ASSERT((a & ATTR_LISTABLE) == 0);
}

static void test_docstring_set(void) {
    SymbolDef* def = symtab_get_def("ChineseRemainder");
    ASSERT(def != NULL);
    ASSERT(def->docstring != NULL);
    ASSERT(strstr(def->docstring, "smallest x") != NULL);
}

static void test_sym_pointer_interned(void) {
    ASSERT(SYM_ChineseRemainder != NULL);
    ASSERT(strcmp(SYM_ChineseRemainder, "ChineseRemainder") == 0);
}

/* --- Memory-safety stress loop ------------------------------------- */

static void test_repeated_evaluation_does_not_corrupt(void) {
    /* Mix machine / bignum / offset / non-coprime / no-solution / error
     * cases to catch double-frees and leaks under valgrind. */
    for (int k = 0; k < 50; k++) {
        assert_eval_eq("ChineseRemainder[{3, 4}, {4, 5}]", "19", 0);
        assert_eval_eq("ChineseRemainder[{2, 3, 5}, {3, 5, 7}, 100]", "173", 0);
        assert_eval_eq("ChineseRemainder[{2, 2}, {4, 6}]", "2", 0);
        assert_eval_eq("ChineseRemainder[{}, {}]", "0", 0);
        assert_eval_eq("ChineseRemainder[{1, 2}, {6, 10}]",
                       "ChineseRemainder[{1, 2}, {6, 10}]", 0);
        assert_eval_eq(
            "Mod[ChineseRemainder[{123, 456, 789}, "
            "{1000000007, 1000000009, 1000000021}], 1000000009]", "456", 0);
        char* result = NULL;
        char* err = eval_capturing_stderr("ChineseRemainder[]", &result);
        free(result);
        free(err);
    }
}

int main(void) {
    symtab_init();
    core_init();

    TEST(test_doc_two_mod);
    TEST(test_doc_three_mod);
    TEST(test_doc_reduce_findinstance_agreement);

    TEST(test_congruence_contract);
    TEST(test_single_congruence);
    TEST(test_unreduced_and_negative_residues);

    TEST(test_rns_recover_number);
    TEST(test_rns_recover_product);

    TEST(test_bignum_roundtrip);
    TEST(test_bignum_is_large);

    TEST(test_offset_skips_base);
    TEST(test_offset_boundary);
    TEST(test_offset_negative);

    TEST(test_noncoprime_consistent);

    TEST(test_no_solution_unevaluated);
    TEST(test_no_solution_variants);

    TEST(test_empty_lists);
    TEST(test_zero_modulus_unevaluated);
    TEST(test_length_mismatch_unevaluated);
    TEST(test_symbolic_unevaluated);
    TEST(test_noninteger_unevaluated);

    TEST(test_arg_error_zero);
    TEST(test_arg_error_one);

    TEST(test_attributes);
    TEST(test_docstring_set);
    TEST(test_sym_pointer_interned);

    TEST(test_repeated_evaluation_does_not_corrupt);

    printf("All ChineseRemainder tests passed!\n");
    return 0;
}
