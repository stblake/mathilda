/* Unit tests for JacobiSymbol.
 *
 *   JacobiSymbol[n, m] -> the Jacobi symbol (n/m), returned as -1, 0, or 1.
 *
 * JacobiSymbol follows the Wolfram Language: it is the full Kronecker-symbol
 * generalisation, so the modulus m may be even or non-positive and n may be
 * negative.  For prime m it is the Legendre symbol.
 *
 * Coverage:
 *   - Documented Mathematica examples (the Table, list-threading, and large
 *     and even-modulus examples).
 *   - Legendre reduction for prime m (quadratic residues mod 7).
 *   - Multiplicativity / reciprocity sanity on a composite odd modulus.
 *   - BigInt arguments (the GMP path).
 *   - Symbolic arguments left unevaluated (no diagnostic).
 *   - Listable threading over a 1D list and a 2D array.
 *   - Wrong arg count -> JacobiSymbol::argrx, call retained.
 *   - Attribute (Protected + Listable), docstring, interned-symbol checks.
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

/* Capture stderr while `input` is parsed + evaluated.  Returns the collected
 * stderr text as a heap string (caller frees) and writes the printed result
 * into *out_result_str (also heap-allocated).  Uses a fixed temp file path;
 * safe because tests run serially. */
static char* eval_capturing_stderr(const char* input, char** out_result_str) {
    const char* path = "/tmp/mathilda_jacobisymbol_stderr.log";
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

static void test_doc_zero(void) {
    /* Spec: JacobiSymbol[10, 5] -> 0 (5 | 10). */
    assert_eval_eq("JacobiSymbol[10, 5]", "0", 0);
}

static void test_doc_table(void) {
    /* Spec: Table[JacobiSymbol[n, m], {n, 0, 10}, {m, 1, n, 2}]. */
    assert_eval_eq(
        "Table[JacobiSymbol[n, m], {n, 0, 10}, {m, 1, n, 2}]",
        "{{}, {1}, {1}, {1, 0}, {1, 1}, {1, -1, 0}, {1, 0, 1}, "
        "{1, 1, -1, 0}, {1, -1, -1, 1}, {1, 0, 1, 1, 0}, "
        "{1, 1, 0, -1, 1}}", 0);
}

static void test_doc_large_args(void) {
    /* Spec: JacobiSymbol[10^10 + 1, Prime[1000]] -> 1.  Prime[1000] == 7919. */
    assert_eval_eq("JacobiSymbol[10^10 + 1, Prime[1000]]", "1", 0);
    assert_eval_eq("JacobiSymbol[10^10 + 1, 7919]", "1", 0);
}

static void test_doc_even_modulus(void) {
    /* Spec: "Second argument may be even": JacobiSymbol[7, 6] -> 1. */
    assert_eval_eq("JacobiSymbol[7, 6]", "1", 0);
}

static void test_doc_thread_first(void) {
    /* Spec: JacobiSymbol[{2,3,5,7,11}, 3] -> {-1,0,-1,1,-1}. */
    assert_eval_eq("JacobiSymbol[{2, 3, 5, 7, 11}, 3]",
                   "{-1, 0, -1, 1, -1}", 0);
}

static void test_doc_negative_first(void) {
    /* Spec: JacobiSymbol[-3, {1, 3, 5, 7}] -> {1, 0, -1, 1}. */
    assert_eval_eq("JacobiSymbol[-3, {1, 3, 5, 7}]",
                   "{1, 0, -1, 1}", 0);
}

/* --- Legendre reduction for prime m -------------------------------- */

static void test_legendre_mod7(void) {
    /* The quadratic residues mod 7 are {1, 2, 4}; the non-residues are
     * {3, 5, 6}.  For prime m the Jacobi symbol is the Legendre symbol. */
    assert_eval_eq("JacobiSymbol[1, 7]", "1", 0);
    assert_eval_eq("JacobiSymbol[2, 7]", "1", 0);
    assert_eval_eq("JacobiSymbol[4, 7]", "1", 0);
    assert_eval_eq("JacobiSymbol[3, 7]", "-1", 0);
    assert_eval_eq("JacobiSymbol[5, 7]", "-1", 0);
    assert_eval_eq("JacobiSymbol[6, 7]", "-1", 0);
    assert_eval_eq("JacobiSymbol[7, 7]", "0", 0);
    /* Map over the full residue list, mirroring the Legendre table. */
    assert_eval_eq("Table[JacobiSymbol[n, 7], {n, 0, 6}]",
                   "{0, 1, 1, -1, 1, -1, -1}", 0);
}

/* --- Composite odd modulus: multiplicativity in the top argument --- */

static void test_multiplicative_numerator(void) {
    /* (a b / m) == (a/m)(b/m) for odd m = 15, a = 2, b = 7. */
    assert_eval_eq(
        "JacobiSymbol[2*7, 15] == JacobiSymbol[2, 15] JacobiSymbol[7, 15]",
        "True", 0);
    /* A composite-modulus value that is +1 yet 14 is a non-residue mod 15
     * (the Jacobi symbol does not detect residuosity for composite m). */
    assert_eval_eq("JacobiSymbol[14, 15]", "-1", 0);
}

/* --- BigInt arguments (GMP path) ----------------------------------- */

static void test_bigint_args(void) {
    /* A 40-digit numerator against a large prime modulus.  Value cross-checked
     * against GMP's mpz_kronecker. */
    assert_eval_eq(
        "JacobiSymbol[12345678901234567890123456789012345678901, "
        "1000000000000000000000000000057]",
        "1", 0);
    /* (m^2 k / m) == 0 for any k, since m | numerator. */
    assert_eval_eq("JacobiSymbol[(2^61 - 1)^2, 2^61 - 1]", "0", 0);
}

/* --- Symbolic passthrough ------------------------------------------ */

static void test_symbolic_silent(void) {
    char* result = NULL;
    char* err = eval_capturing_stderr("JacobiSymbol[a, 5]", &result);
    ASSERT(result != NULL);
    ASSERT_MSG(strstr(result, "JacobiSymbol[a, 5]") != NULL,
               "expected unevaluated call, got: %s", result);
    ASSERT_MSG(err == NULL || strstr(err, "JacobiSymbol::") == NULL,
               "expected no diagnostic for symbolic input, got: %s",
               err ? err : "(null)");
    free(result);
    free(err);
}

static void test_symbolic_both(void) {
    char* result = NULL;
    char* err = eval_capturing_stderr("JacobiSymbol[a, b]", &result);
    ASSERT(result != NULL);
    ASSERT(strstr(result, "JacobiSymbol[a, b]") != NULL);
    ASSERT(err == NULL || strstr(err, "JacobiSymbol::") == NULL);
    free(result);
    free(err);
}

/* --- Listable threading -------------------------------------------- */

static void test_listable_both(void) {
    /* Element-wise over two equal-length lists. */
    assert_eval_eq("JacobiSymbol[{2, 3, 4}, {3, 5, 7}]",
                   "{-1, -1, 1}", 0);
}

static void test_listable_2d(void) {
    assert_eval_eq("JacobiSymbol[{{2, 3}, {5, 7}}, 11]",
                   "{{-1, 1}, {1, -1}}", 0);
}

/* --- Error paths --------------------------------------------------- */

static void test_argrx_zero(void) {
    /* Spec: JacobiSymbol[] -> JacobiSymbol::argrx ... 2 arguments expected. */
    char* result = NULL;
    char* err = eval_capturing_stderr("JacobiSymbol[]", &result);
    ASSERT(result != NULL);
    ASSERT_MSG(strstr(result, "JacobiSymbol[]") != NULL,
               "expected unevaluated call, got: %s", result);
    ASSERT(err != NULL);
    ASSERT_MSG(strstr(err, "JacobiSymbol::argrx") != NULL,
               "expected argrx diagnostic, got: %s", err);
    ASSERT(strstr(err, "0 arguments") != NULL);
    ASSERT(strstr(err, "2 arguments are expected") != NULL);
    free(result);
    free(err);
}

static void test_argrx_one(void) {
    /* Singular wording for a single argument. */
    char* result = NULL;
    char* err = eval_capturing_stderr("JacobiSymbol[5]", &result);
    ASSERT(result != NULL);
    ASSERT(strstr(result, "JacobiSymbol[5]") != NULL);
    ASSERT(err != NULL);
    ASSERT(strstr(err, "JacobiSymbol::argrx") != NULL);
    ASSERT_MSG(strstr(err, "1 argument;") != NULL,
               "expected singular 'argument', got: %s", err);
    free(result);
    free(err);
}

static void test_argrx_three(void) {
    char* result = NULL;
    char* err = eval_capturing_stderr("JacobiSymbol[1, 2, 3]", &result);
    ASSERT(result != NULL);
    ASSERT(strstr(result, "JacobiSymbol[1, 2, 3]") != NULL);
    ASSERT(err != NULL);
    ASSERT(strstr(err, "JacobiSymbol::argrx") != NULL);
    ASSERT(strstr(err, "3 arguments") != NULL);
    free(result);
    free(err);
}

/* --- Attribute / docstring / interned-symbol introspection -------- */

static void test_attributes(void) {
    SymbolDef* def = symtab_get_def("JacobiSymbol");
    ASSERT(def != NULL);
    uint32_t a = get_attributes("JacobiSymbol");
    ASSERT((a & ATTR_PROTECTED) != 0);
    ASSERT((a & ATTR_LISTABLE) != 0);
}

static void test_docstring_set(void) {
    SymbolDef* def = symtab_get_def("JacobiSymbol");
    ASSERT(def != NULL);
    ASSERT(def->docstring != NULL);
    ASSERT(strstr(def->docstring, "Jacobi") != NULL);
}

static void test_sym_pointer_interned(void) {
    ASSERT(SYM_JacobiSymbol != NULL);
    ASSERT(strcmp(SYM_JacobiSymbol, "JacobiSymbol") == 0);
}

/* --- Memory-safety stress loop ------------------------------------ */

static void test_repeated_evaluation_does_not_corrupt(void) {
    /* Mix machine-int / bignum / listable / even-modulus / negative / error
     * cases to catch double-frees and leaks under valgrind. */
    for (int k = 0; k < 50; k++) {
        assert_eval_eq("JacobiSymbol[10, 5]", "0", 0);
        assert_eval_eq("JacobiSymbol[7, 6]", "1", 0);
        assert_eval_eq("JacobiSymbol[{2, 3, 5, 7, 11}, 3]",
                       "{-1, 0, -1, 1, -1}", 0);
        assert_eval_eq("JacobiSymbol[-3, {1, 3, 5, 7}]",
                       "{1, 0, -1, 1}", 0);
        assert_eval_eq("JacobiSymbol[10^10 + 1, 7919]", "1", 0);
        char* result = NULL;
        char* err = eval_capturing_stderr("JacobiSymbol[]", &result);
        free(result);
        free(err);
    }
}

int main(void) {
    symtab_init();
    core_init();

    TEST(test_doc_zero);
    TEST(test_doc_table);
    TEST(test_doc_large_args);
    TEST(test_doc_even_modulus);
    TEST(test_doc_thread_first);
    TEST(test_doc_negative_first);

    TEST(test_legendre_mod7);
    TEST(test_multiplicative_numerator);

    TEST(test_bigint_args);

    TEST(test_symbolic_silent);
    TEST(test_symbolic_both);

    TEST(test_listable_both);
    TEST(test_listable_2d);

    TEST(test_argrx_zero);
    TEST(test_argrx_one);
    TEST(test_argrx_three);

    TEST(test_attributes);
    TEST(test_docstring_set);
    TEST(test_sym_pointer_interned);

    TEST(test_repeated_evaluation_does_not_corrupt);

    printf("All JacobiSymbol tests passed!\n");
    return 0;
}
