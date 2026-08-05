/* Unit tests for the Compile[] B1 read-only Association parameter bag.
 *
 * Coverage matrix (each head, ASSOC_ARG and ASSOC_CONST paths):
 *   - Lookup[a,k] / Lookup[a,k,default], KeyExistsQ/KeyMemberQ/KeyFreeQ,
 *     Length, Values — over a declared `{p, _Association[, _valtype]}` argument
 *     bag and over a compile-time-constant association (literal or folded global).
 *
 * Each test asserts:
 *   1. "it lowers":  CompileDiagnostics[spec, body] contains Compiled -> True.
 *   2. "compiled == interpreted":  the compiled value === the interpreter's OWN
 *      result on the SAME association input (never a hand-computed constant).
 *      For a scalar that is a direct `===`; for the packed Values result it is
 *      Normal[...] === the interpreter's Values (Normal removes packed-vs-plain).
 *   3. Cliff:   a COMPOSED body compiles whole (Lookup[p,"a"]+x, Total[Values[p]]).
 *   4. Decline: a missing key with no default declines at runtime and the
 *      interpreter answers (Missing[...]), never a crash or a wrong number.
 *   5. B2: a runtime-varying integer/real key lowers to ASSOC_LOOKUP_DYN (the
 *      Lookup[a,i]-in-a-loop pattern); a key past the machine-scalar boundary
 *      (an array) still bails the whole body.
 *   6. Purity/hoist: a loop-invariant constant-key Lookup appears ONCE in the
 *      compiled listing (CSE + LICM), proving the op is pure.
 *   7. Const-fold: Lookup of a constant association folds to a constant — no
 *      ASSOC_LOOKUP in the listing — which is the auto-compilation story.
 *   8. A repeated-evaluation loop to surface double-free / leak regressions.
 */

#include "expr.h"
#include "eval.h"
#include "core.h"
#include "symtab.h"
#include "test_utils.h"
#include "parse.h"
#include "print.h"
#include <string.h>
#include <stdlib.h>

/* Evaluate `input`; assert its printed result CONTAINS `substr`. */
static void assert_eval_contains(const char* input, const char* substr) {
    Expr* p = parse_expression(input);
    Expr* e = evaluate(p);
    char* s = expr_to_string(e);
    ASSERT_MSG(strstr(s, substr) != NULL,
               "expected result of %s to contain \"%s\", got: %s", input, substr, s);
    free(s); expr_free(p); expr_free(e);
}
/* Evaluate `input`; assert its printed result does NOT contain `substr`. */
static void assert_eval_lacks(const char* input, const char* substr) {
    Expr* p = parse_expression(input);
    Expr* e = evaluate(p);
    char* s = expr_to_string(e);
    ASSERT_MSG(strstr(s, substr) == NULL,
               "expected result of %s to NOT contain \"%s\", got: %s", input, substr, s);
    free(s); expr_free(p); expr_free(e);
}
static void assert_true(const char* input) { assert_eval_eq(input, "True", 0); }
static void assert_lowers(const char* diag)     { assert_eval_contains(diag, "Compiled\" -> True"); }
static void assert_not_lowers(const char* diag) { assert_eval_contains(diag, "Compiled\" -> False"); }

/* ------------------------------------------------------------------ *
 *  ASSOC_ARG — an association passed as a compiled argument bag       *
 * ------------------------------------------------------------------ */

static void test_lookup_arg(void) {
    assert_lowers("CompileDiagnostics[{{p, _Association}, {x, _Real}}, Lookup[p, \"a\"] + x]");
    /* compiled value === interpreter's own value on the same association */
    assert_true("Compile[{{p, _Association}, {x, _Real}}, Lookup[p, \"a\"] + x]"
                "[<|\"a\" -> 2., \"b\" -> 5.|>, 3.] === (Lookup[<|\"a\" -> 2., \"b\" -> 5.|>, \"a\"] + 3.)");
    /* integer keys, several probes on the SAME bag (CSE-safe: distinct keys) */
    assert_true("Compile[{{p, _Association}}, Lookup[p, 1] - Lookup[p, 3]]"
                "[<|1 -> 10., 2 -> 20., 3 -> 4.|>] === 6.");
}

static void test_lookup_default(void) {
    /* present key: default ignored; absent key: default returned (no decline) */
    assert_true("Compile[{{p, _Association}}, Lookup[p, \"k\", -1.]]"
                "[<|\"k\" -> 9.|>] === 9.");
    assert_true("Compile[{{p, _Association}}, Lookup[p, \"k\", -1.]]"
                "[<|\"x\" -> 1.|>] === -1.");
}

static void test_membership(void) {
    assert_lowers("CompileDiagnostics[{{p, _Association}}, KeyExistsQ[p, \"k\"]]");
    assert_true("Compile[{{p, _Association}}, KeyExistsQ[p, \"k\"]][<|\"k\" -> 1.|>] === True");
    assert_true("Compile[{{p, _Association}}, KeyExistsQ[p, \"k\"]][<|\"z\" -> 1.|>] === False");
    assert_true("Compile[{{p, _Association}}, KeyFreeQ[p, \"k\"]][<|\"z\" -> 1.|>] === True");
    /* KeyMemberQ is the same op */
    assert_true("Compile[{{p, _Association}}, KeyMemberQ[p, 2]][<|1 -> 1., 2 -> 2.|>] === True");
}

static void test_length_arg(void) {
    assert_lowers("CompileDiagnostics[{{p, _Association}}, Length[p]]");
    assert_true("Compile[{{p, _Association}}, Length[p]][<|\"a\" -> 1., \"b\" -> 2., \"c\" -> 3.|>] === 3");
    /* Length[array] must still lower via the array path — the shared head is not hijacked. */
    assert_lowers("CompileDiagnostics[{{v, _Real, 1}}, Length[v]]");
    assert_true("Compile[{{v, _Real, 1}}, Length[v]][{1., 2., 3., 4.}] === 4");
}

static void test_value_types(void) {
    /* integer bag */
    assert_true("Compile[{{p, _Association, _Integer}}, Length[p] + Lookup[p, \"n\"]]"
                "[<|\"n\" -> 10, \"m\" -> 20|>] === 12");
    /* complex bag */
    assert_true("Compile[{{p, _Association, _Complex}}, Re[Lookup[p, \"z\"]]]"
                "[<|\"z\" -> 3. + 4. I|>] === 3.");
}

static void test_values_arg(void) {
    assert_lowers("CompileDiagnostics[{{p, _Association, _Real}}, Total[Values[p]]]");
    assert_true("Compile[{{p, _Association, _Real}}, Total[Values[p]]]"
                "[<|\"a\" -> 1., \"b\" -> 2., \"c\" -> 3.|>] === 6.");
    /* the packed result Normal-equals the interpreter's own Values on the same bag */
    assert_true("Normal[Compile[{{p, _Association, _Real}}, Values[p]][<|\"a\" -> 1., \"b\" -> 2.|>]]"
                " === Values[<|\"a\" -> 1., \"b\" -> 2.|>]");
}

/* ------------------------------------------------------------------ *
 *  Cliff — a composed body compiles whole (a bail costs the WHOLE body)*
 * ------------------------------------------------------------------ */

static void test_cliff(void) {
    assert_lowers("CompileDiagnostics[{{p, _Association}, {x, _Real}},"
                  " If[KeyExistsQ[p, \"k\"], Lookup[p, \"k\"], Lookup[p, \"d\", 0.]] * x]");
    assert_true("Compile[{{p, _Association}, {x, _Real}},"
                " If[KeyExistsQ[p, \"k\"], Lookup[p, \"k\"], Lookup[p, \"d\", 0.]] * x]"
                "[<|\"k\" -> 4.|>, 2.] === 8.");
    /* Sum over Values, all inside one compiled body. */
    assert_true("Compile[{{p, _Association, _Real}}, Total[Values[p]]^2]"
                "[<|\"a\" -> 1., \"b\" -> 3.|>] === 16.");
}

/* ------------------------------------------------------------------ *
 *  Decline — missing key, no default: runtime falls back, no crash    *
 * ------------------------------------------------------------------ */

static void test_decline_missing_key(void) {
    /* It still LOWERS (the shape is compilable); at runtime the absent key with
     * no default declines and the interpreter answers with Missing[...]. */
    assert_lowers("CompileDiagnostics[{{p, _Association}}, Lookup[p, \"zzz\"] + 1.]");
    assert_eval_contains("Compile[{{p, _Association}}, Lookup[p, \"zzz\"] + 1.][<|\"a\" -> 2.|>]",
                         "Missing");
}

/* ------------------------------------------------------------------ *
 *  B2 — a runtime-varying integer / real key compiles (ASSOC_LOOKUP_DYN)*
 * ------------------------------------------------------------------ */

static void test_runtime_key(void) {
    /* The Lookup[a, i]-in-a-loop pattern that B1 bailed on now lowers. */
    assert_lowers("CompileDiagnostics[{{p, _Association, _Real}, {k, _Integer}}, Lookup[p, k, -1.]]");
    /* parity vs interpreter, present and absent (default). */
    assert_true("Compile[{{p, _Association, _Real}, {k, _Integer}}, Lookup[p, k, -1.]]"
                "[<|1 -> 10., 2 -> 20.|>, 2] === Lookup[<|1 -> 10., 2 -> 20.|>, 2, -1.]");
    assert_true("Compile[{{p, _Association, _Real}, {k, _Integer}}, Lookup[p, k, -1.]]"
                "[<|1 -> 10., 2 -> 20.|>, 9] === -1.");
    /* the workhorse: sum Lookup[a, i] over a loop of runtime keys. */
    assert_true("Compile[{{p, _Association, _Real}, {n, _Integer}},"
                " Module[{s = 0.}, Do[s = s + Lookup[p, i, 0.], {i, 1, n}]; s]]"
                "[<|1 -> 10., 2 -> 20., 3 -> 30.|>, 5] === 60.");
    /* integer value type, runtime key. */
    assert_true("Compile[{{p, _Association, _Integer}, {k, _Integer}}, Lookup[p, k, 0]]"
                "[<|10 -> 100, 20 -> 200|>, 20] === 200");
    /* the op is a real runtime read — a CompilePrint shows the DYN lookup, and it
     * is NOT const-folded away (unlike a literal key). */
    assert_eval_contains("CompilePrint[{{p, _Association, _Real}, {k, _Integer}}, Lookup[p, k]]", "Lookup[");
    /* Still bails past the machine-scalar boundary: an array-typed key is not a
     * key the register model can carry, so the whole body falls back. */
    assert_not_lowers("CompileDiagnostics[{{p, _Association}, {v, _Integer, 1}}, Lookup[p, v]]");
}

/* ------------------------------------------------------------------ *
 *  Purity / hoist — a loop-invariant lookup appears once              *
 * ------------------------------------------------------------------ */

static void test_lookup_is_hoisted(void) {
    /* The listing must contain exactly one Lookup[ — CSE + LICM lift the pure,
     * loop-invariant constant-key probe out of the Sum body. */
    Expr* p = parse_expression(
        "StringCount[ToString[CompilePrint[{{p, _Association}, {n, _Integer}},"
        " Sum[Lookup[p, \"a\"], {i, 1, n}]]], \"Lookup[\"]");
    Expr* e = evaluate(p);
    char* s = expr_to_string(e);
    ASSERT_MSG(strcmp(s, "1") == 0, "loop-invariant Lookup was not hoisted (count=%s, expected 1)", s);
    free(s); expr_free(p); expr_free(e);
}

/* ------------------------------------------------------------------ *
 *  Const-fold — a constant association folds (the auto-compile story) *
 * ------------------------------------------------------------------ */

static void test_const_fold_literal(void) {
    /* Lookup of a literal association folds to a constant: no ASSOC_LOOKUP op. */
    assert_lowers("CompileDiagnostics[{{x, _Real}}, Lookup[<|\"a\" -> 2.|>, \"a\"] + x]");
    assert_eval_lacks("CompilePrint[{{x, _Real}}, Lookup[<|\"a\" -> 2.|>, \"a\"] + x]", "Lookup[");
    assert_true("Compile[{{x, _Real}}, Lookup[<|\"a\" -> 2.|>, \"a\"] + x][3.] === 5.");
    /* KeyExistsQ / Length of a constant fold too. */
    assert_true("Compile[{{x, _Real}}, If[KeyExistsQ[<|\"a\" -> 1.|>, \"a\"], x, -x]][2.] === 2.");
    assert_true("Compile[{{x, _Real}}, Length[<|\"a\" -> 1., \"b\" -> 2.|>] + x][1.] === 3.");
}

/* Auto-compilation end-to-end: a global association captured in a Table body. */
static void test_autocompile_global(void) {
    expr_free(eval_and_free(parse_expression("gassoc = <|\"k\" -> 3.|>")));
    /* Table folds the global lookup and matches the interpreter. */
    assert_true("Table[Lookup[gassoc, \"k\"] * i, {i, 1, 4}] === {3., 6., 9., 12.}");
    expr_free(eval_and_free(parse_expression("Clear[gassoc]")));
}

/* ------------------------------------------------------------------ *
 *  B3 — KeyDrop / KeyTake produce OWNED association values            *
 * ------------------------------------------------------------------ */

static void test_keydrop_keytake(void) {
    /* lowers, and the compiled function RETURNS an association */
    assert_lowers("CompileDiagnostics[{{p, _Association, _Real}}, KeyDrop[p, \"b\"]]");
    assert_lowers("CompileDiagnostics[{{p, _Association, _Real}}, KeyTake[p, {\"a\", \"c\"}]]");
    assert_true("AssociationQ[Compile[{{p, _Association, _Real}}, KeyDrop[p, \"b\"]][<|\"a\" -> 1., \"b\" -> 2.|>]]");

    /* parity vs interpreter — single key, key-list, order preserved */
    assert_true("Compile[{{p, _Association, _Real}}, KeyDrop[p, \"b\"]][<|\"a\" -> 1., \"b\" -> 2., \"c\" -> 3.|>]"
                " === KeyDrop[<|\"a\" -> 1., \"b\" -> 2., \"c\" -> 3.|>, \"b\"]");
    assert_true("Compile[{{p, _Association, _Real}}, KeyTake[p, {\"a\", \"c\"}]]"
                "[<|\"a\" -> 1., \"b\" -> 2., \"c\" -> 3., \"d\" -> 4.|>]"
                " === KeyTake[<|\"a\" -> 1., \"b\" -> 2., \"c\" -> 3., \"d\" -> 4.|>, {\"a\", \"c\"}]");

    /* integer keys / integer value type */
    assert_true("Compile[{{p, _Association, _Integer}}, KeyDrop[p, 2]][<|1 -> 100, 2 -> 200, 3 -> 300|>]"
                " === <|1 -> 100, 3 -> 300|>");

    /* edge cases: dropping an absent key is a no-op; dropping every key empties */
    assert_true("Compile[{{p, _Association, _Real}}, KeyDrop[p, \"zzz\"]][<|\"a\" -> 1.|>] === <|\"a\" -> 1.|>");
    assert_true("Compile[{{p, _Association, _Real}}, KeyDrop[p, {\"a\", \"b\"}]][<|\"a\" -> 1., \"b\" -> 2.|>] === <||>");
    assert_true("Compile[{{p, _Association, _Real}}, KeyTake[p, {\"zzz\"}]][<|\"a\" -> 1.|>] === <||>");

    /* the returned association is a genuine, usable association (index intact) */
    assert_true("Lookup[Compile[{{p, _Association, _Real}}, KeyDrop[p, \"b\"]][<|\"a\" -> 5., \"b\" -> 9.|>], \"a\"] === 5.");
    assert_true("Length[Compile[{{p, _Association, _Real}}, KeyTake[p, {\"a\", \"c\"}]][<|\"a\" -> 1., \"b\" -> 2., \"c\" -> 3.|>]] === 2");

    /* constant source (a literal association) — the transform runs from the spec */
    assert_true("Compile[{{x, _Real}}, KeyTake[<|1 -> 10., 2 -> 20., 3 -> 30.|>, {1, 3}]][0.]"
                " === <|1 -> 10., 3 -> 30.|>");
}

static void test_keysel_boundaries(void) {
    /* a runtime-varying key set is not in B3 — the whole body falls back */
    assert_not_lowers("CompileDiagnostics[{{p, _Association, _Real}, {k, _Integer}}, KeyDrop[p, k]]");
    /* KeyDrop on a non-association operand is not an association op -> bails */
    assert_not_lowers("CompileDiagnostics[{{v, _Real, 1}}, KeyDrop[v, 1]]");
}

/* ------------------------------------------------------------------ *
 *  B3 composition — a produced association is a first-class operand   *
 * ------------------------------------------------------------------ */

static void test_composition(void) {
    /* a scalar reader over a produced association */
    assert_lowers("CompileDiagnostics[{{p, _Association, _Real}}, Lookup[KeyDrop[p, \"b\"], \"a\"]]");
    assert_true("Compile[{{p, _Association, _Real}}, Lookup[KeyDrop[p, \"b\"], \"a\"]]"
                "[<|\"a\" -> 5., \"b\" -> 9., \"c\" -> 3.|>] === 5.");
    /* nested transforms — chained filtering, returns an association */
    assert_lowers("CompileDiagnostics[{{p, _Association, _Real}}, KeyTake[KeyDrop[p, \"a\"], {\"b\"}]]");
    assert_true("Compile[{{p, _Association, _Real}}, KeyTake[KeyDrop[p, \"a\"], {\"b\"}]]"
                "[<|\"a\" -> 1., \"b\" -> 2., \"c\" -> 3.|>]"
                " === KeyTake[KeyDrop[<|\"a\" -> 1., \"b\" -> 2., \"c\" -> 3.|>, \"a\"], {\"b\"}]");
    /* Values of a produced association, aggregated */
    assert_true("Compile[{{p, _Association, _Real}}, Total[Values[KeyDrop[p, \"b\"]]]]"
                "[<|\"a\" -> 1., \"b\" -> 99., \"c\" -> 3.|>] === 4.");
    /* Length / membership over a produced association */
    assert_true("Compile[{{p, _Association, _Real}}, Length[KeyDrop[KeyTake[p, {\"a\", \"b\", \"c\"}], \"b\"]]]"
                "[<|\"a\" -> 1., \"b\" -> 2., \"c\" -> 3., \"d\" -> 4.|>] === 2");
    assert_true("Compile[{{p, _Association, _Real}}, KeyExistsQ[KeyDrop[p, \"a\"], \"a\"]]"
                "[<|\"a\" -> 1., \"b\" -> 2.|>] === False");
    /* a runtime (B2) key over a produced association */
    assert_true("Compile[{{p, _Association, _Real}, {i, _Integer}}, Lookup[KeyDrop[p, 2], i, -1.]]"
                "[<|1 -> 10., 2 -> 20., 3 -> 30.|>, 3] === 30.");
    assert_true("Compile[{{p, _Association, _Real}, {i, _Integer}}, Lookup[KeyDrop[p, 2], i, -1.]]"
                "[<|1 -> 10., 2 -> 20., 3 -> 30.|>, 2] === -1.");
}

/* ------------------------------------------------------------------ *
 *  B3 — Counts[machine array] -> an association of element -> count    *
 * ------------------------------------------------------------------ */

static void test_counts(void) {
    assert_lowers("CompileDiagnostics[{{v, _Integer, 1}}, Counts[v]]");
    /* parity vs interpreter, integer and real element arrays */
    assert_true("Compile[{{v, _Integer, 1}}, Counts[v]][{1, 2, 1, 3, 1, 2}]"
                " === Counts[{1, 2, 1, 3, 1, 2}]");
    assert_true("Compile[{{v, _Real, 1}}, Counts[v]][{1., 2., 2., 3.}] === Counts[{1., 2., 2., 3.}]");

    /* composition: look a count up (values are integers, default 0 for absent) */
    assert_true("Compile[{{v, _Integer, 1}, {k, _Integer}}, Lookup[Counts[v], k, 0]][{1, 2, 1, 3, 1}, 1] === 3");
    assert_true("Compile[{{v, _Integer, 1}, {k, _Integer}}, Lookup[Counts[v], k, 0]][{1, 2, 1}, 9] === 0");

    /* Length[Counts] = number of distinct elements; Total[Values[Counts]] = Length[v] */
    assert_true("Compile[{{v, _Real, 1}}, Length[Counts[v]]][{1., 2., 2., 3., 3., 3.}] === 3");
    assert_true("Compile[{{v, _Integer, 1}}, Total[Values[Counts[v]]]][{5, 5, 7, 9, 9}] === 5");

    /* a rank-2 (non-vector) operand is outside the subset -> the body falls back */
    assert_not_lowers("CompileDiagnostics[{{m, _Integer, 2}}, Counts[m]]");
}

/* ------------------------------------------------------------------ *
 *  Leak / double-free surface — build, apply, free in a loop          *
 * ------------------------------------------------------------------ */

static void test_repeated_eval_no_leak(void) {
    for (int i = 0; i < 50; i++) {
        Expr* r = eval_and_free(parse_expression(
            "Compile[{{p, _Association}, {x, _Real}}, Lookup[p, \"a\"] + Lookup[p, \"b\", 0.] + x]"
            "[<|\"a\" -> 1., \"b\" -> 2.|>, 3.]"));
        char* s = expr_to_string(r);
        ASSERT_STR_EQ(s, "6.0");
        free(s); expr_free(r);
    }
    /* Values path in the loop too (owned packed result freed each call). */
    for (int i = 0; i < 50; i++) {
        Expr* r = eval_and_free(parse_expression(
            "Compile[{{p, _Association, _Real}}, Total[Values[p]]][<|\"a\" -> 1., \"b\" -> 2., \"c\" -> 3.|>]"));
        char* s = expr_to_string(r);
        ASSERT_STR_EQ(s, "6.0");
        free(s); expr_free(r);
    }
    /* B3: an owned association RESULT built and freed every iteration. */
    for (int i = 0; i < 50; i++) {
        Expr* r = eval_and_free(parse_expression(
            "Compile[{{p, _Association, _Real}}, KeyDrop[p, \"b\"]][<|\"a\" -> 1., \"b\" -> 2., \"c\" -> 3.|>]"));
        char* s = expr_to_string(r);
        ASSERT_STR_EQ(s, "<|\"a\" -> 1.0, \"c\" -> 3.0|>");
        free(s); expr_free(r);
    }
    /* B3 composition: a produced association is consumed AND freed each iteration
     * (a leak or double-free in the free-source discipline shows up here). */
    for (int i = 0; i < 50; i++) {
        Expr* r = eval_and_free(parse_expression(
            "Compile[{{p, _Association, _Real}}, Total[Values[KeyTake[KeyDrop[p, \"a\"], {\"b\", \"c\"}]]]]"
            "[<|\"a\" -> 1., \"b\" -> 2., \"c\" -> 3.|>]"));
        char* s = expr_to_string(r);
        ASSERT_STR_EQ(s, "5.0");
        free(s); expr_free(r);
    }
    /* B3 Counts: array -> association built and freed each iteration. */
    for (int i = 0; i < 50; i++) {
        Expr* r = eval_and_free(parse_expression(
            "Compile[{{v, _Integer, 1}, {k, _Integer}}, Lookup[Counts[v], k, 0]][{1, 2, 1, 3, 1}, 1]"));
        char* s = expr_to_string(r);
        ASSERT_STR_EQ(s, "3");
        free(s); expr_free(r);
    }
}

int main(void) {
    symtab_init();
    core_init();

    TEST(test_lookup_arg);
    TEST(test_lookup_default);
    TEST(test_membership);
    TEST(test_length_arg);
    TEST(test_value_types);
    TEST(test_values_arg);
    TEST(test_cliff);
    TEST(test_decline_missing_key);
    TEST(test_runtime_key);
    TEST(test_lookup_is_hoisted);
    TEST(test_const_fold_literal);
    TEST(test_autocompile_global);
    TEST(test_keydrop_keytake);
    TEST(test_keysel_boundaries);
    TEST(test_composition);
    TEST(test_counts);
    TEST(test_repeated_eval_no_leak);

    printf("All Compile[] Association (B1/B2/B3) tests passed!\n");
    return 0;
}
