/* test_condition_downvalue.c
 *
 * Pins the behaviour of /;-guarded and ?-guarded DownValues — i.e.
 * rules whose LHS is wrapped in Condition[...] or PatternTest[...]
 * after parsing.  These are critical for the CRC integral table
 * (src/internal/CRCMathTablesIntegrals.m), where virtually every
 * non-trivial formula is `IntegrateTable[lhs, x_] /; FreeQ[...]`.
 *
 * Prior to 2026-05-15 the §3.5 dispatch filter in symtab.c only
 * stripped HoldPattern/Verbatim; a top-level Condition wrapper
 * caused dispatch_arity and first_arg_head_canon to be computed
 * from the WRONG node (the Condition's args instead of the inner
 * function pattern's args).  The filter then skipped every such
 * rule before match() was even invoked, even though the matcher
 * itself handles Condition correctly post-dispatch.  Result:
 *   f[a_, b_] /; True := rhs   —  never fired.
 *   f[a_] /; FreeQ[a, x] := rhs  — never fired.
 *
 * These tests:
 *   - Pin every Condition / PatternTest shape we care about
 *     (single-arg, multi-arg, compound first slot, Optional
 *     slots, nested /;, false-condition negative cases).
 *   - Guard against re-introducing the dispatch-skip behaviour.
 *   - Sanity-check the existing dispatch wins (literal-head
 *     filtering) that must keep working after the fix.
 *   - Smoke-test that CRC formulas now actually fire end-to-end
 *     once the .m table is loaded.
 */

#include "test_utils.h"
#include "core.h"
#include "expr.h"
#include "eval.h"
#include "parse.h"
#include "symtab.h"
#include "match.h"

#include <stdio.h>
#include <string.h>

/* ----- helpers ---------------------------------------------------- */

/* Evaluate `src` and check whether the result's outermost head (as a
 * symbol name) equals `head_name`.  Used as a cheap "did the rule fire
 * vs. stay unevaluated" probe — if the rule fired we get whatever the
 * RHS evaluates to, otherwise we get the original head back. */
static bool result_head_is(const char* src, const char* head_name) {
    Expr* p = parse_expression(src);
    Expr* r = evaluate(p);
    expr_free(p);
    if (!r) return false;
    bool ok = false;
    if (r->type == EXPR_FUNCTION && r->data.function.head
        && r->data.function.head->type == EXPR_SYMBOL
        && strcmp(r->data.function.head->data.symbol, head_name) == 0) {
        ok = true;
    } else if (r->type == EXPR_SYMBOL
        && strcmp(r->data.symbol, head_name) == 0) {
        ok = true;
    }
    expr_free(r);
    return ok;
}

static bool rule_did_fire(const char* src, const char* original_head) {
    /* "fired" = result's outermost head is NOT the original input head,
     * OR the original input was already a symbol/atom (in which case
     * the RHS replaces it).  We treat "result still has original head"
     * as the failure-to-fire signal. */
    return !result_head_is(src, original_head);
}

/* Quiet helper that just runs `src` (e.g. a rule definition) and
 * discards the value, since we don't care what `SetDelayed` returns. */
static void eval_and_drop(const char* src) {
    Expr* p = parse_expression(src);
    if (!p) {
        fprintf(stderr, "eval_and_drop: failed to parse: %s\n", src);
        exit(1);
    }
    Expr* r = evaluate(p);
    expr_free(p);
    if (r) expr_free(r);
}

/* ----- 1. Single-arg /; rules -------------------------------------- */

static void test_cond_single_arg_true(void) {
    eval_and_drop("ClearAll[CT1]");
    eval_and_drop("CT1[a_] /; True := \"hit\"");
    assert_eval_eq("CT1[5]", "\"hit\"", 0);
    assert_eval_eq("CT1[x]", "\"hit\"", 0);
    assert_eval_eq("CT1[f[1,2,3]]", "\"hit\"", 0);
}

static void test_cond_single_arg_false(void) {
    eval_and_drop("ClearAll[CT2]");
    eval_and_drop("CT2[a_] /; False := \"hit\"");
    /* Condition is False, rule must NOT fire — must remain unevaluated. */
    ASSERT(result_head_is("CT2[5]", "CT2"));
    ASSERT(result_head_is("CT2[x]", "CT2"));
}

static void test_cond_single_arg_freeq(void) {
    eval_and_drop("ClearAll[CT3]");
    eval_and_drop("CT3[a_] /; FreeQ[a, x] := {a, \"free\"}");
    assert_eval_eq("CT3[5]", "{5, \"free\"}", 0);
    assert_eval_eq("CT3[y]", "{y, \"free\"}", 0);
    assert_eval_eq("CT3[Sin[y]]", "{Sin[y], \"free\"}", 0);
    /* x appears inside — Condition is False, rule does NOT fire. */
    ASSERT(result_head_is("CT3[x]", "CT3"));
    ASSERT(result_head_is("CT3[Sin[x]]", "CT3"));
}

/* ----- 2. Multi-arg /; rules --------------------------------------- */

static void test_cond_two_arg_true(void) {
    eval_and_drop("ClearAll[CT4]");
    eval_and_drop("CT4[a_, b_] /; True := {a, b, \"hit\"}");
    assert_eval_eq("CT4[1, 2]", "{1, 2, \"hit\"}", 0);
    assert_eval_eq("CT4[x, y]", "{x, y, \"hit\"}", 0);
    assert_eval_eq("CT4[Sin[z], w]", "{Sin[z], w, \"hit\"}", 0);
}

static void test_cond_two_arg_freeq(void) {
    eval_and_drop("ClearAll[CT5]");
    /* CRC-table-style: FreeQ over the second argument. */
    eval_and_drop("CT5[a_, x_] /; FreeQ[a, x] := {a, x, \"hit\"}");
    /* 5 does not contain y. */
    assert_eval_eq("CT5[5, y]", "{5, y, \"hit\"}", 0);
    /* Sin[w] does not contain y. */
    assert_eval_eq("CT5[Sin[w], y]", "{Sin[w], y, \"hit\"}", 0);
    /* y DOES contain y — Condition False, rule must not fire. */
    ASSERT(result_head_is("CT5[y, y]", "CT5"));
    ASSERT(result_head_is("CT5[Sin[y], y]", "CT5"));
}

static void test_cond_three_arg(void) {
    eval_and_drop("ClearAll[CT6]");
    eval_and_drop("CT6[a_, b_, c_] /; FreeQ[{a, b}, c] := {a, b, c, \"hit\"}");
    assert_eval_eq("CT6[1, 2, x]", "{1, 2, x, \"hit\"}", 0);
    assert_eval_eq("CT6[Sin[u], Cos[v], w]", "{Sin[u], Cos[v], w, \"hit\"}", 0);
    /* c (=x) appears in a (=x) — Condition False. */
    ASSERT(result_head_is("CT6[x, 2, x]", "CT6"));
}

/* ----- 3. Compound first slot -------------------------------------- */

static void test_cond_compound_first_arg_true(void) {
    eval_and_drop("ClearAll[CT7]");
    eval_and_drop("CT7[Sin[x_], y_] /; True := {x, y, \"hit\"}");
    assert_eval_eq("CT7[Sin[3], 5]", "{3, 5, \"hit\"}", 0);
    assert_eval_eq("CT7[Sin[a], b]", "{a, b, \"hit\"}", 0);
    /* First arg is Cos, not Sin — must not fire. */
    ASSERT(result_head_is("CT7[Cos[3], 5]", "CT7"));
    /* First arg is a bare integer — must not fire. */
    ASSERT(result_head_is("CT7[42, 5]", "CT7"));
}

static void test_cond_compound_first_arg_freeq(void) {
    eval_and_drop("ClearAll[CT8]");
    /* Sin's inner is allowed to contain x (it's bound here); the
     * Condition just guards on a separate symbol. */
    eval_and_drop("CT8[Sin[a_ x_], x_] /; FreeQ[a, x] := -Cos[a x]/a");
    /* Sin[2 y]: a=2, x=y. FreeQ[2, y] is True. */
    assert_eval_eq("CT8[Sin[2 y], y]", "-1/2 Cos[2 y]", 0);
    /* Sin[u y]: a=u, x=y. FreeQ[u, y] is True. */
    assert_eval_eq("CT8[Sin[u y], y]", "-Cos[u y]/u", 0);
    /* Sin[y y] = Sin[y^2]: doesn't fit Times[a_, x_] cleanly — must
     * fall through. (sanity check: NOT firing here is fine.) */
}

/* ----- 4. Optional-bearing slot under Condition -------------------- */

static void test_cond_with_optional_slot(void) {
    eval_and_drop("ClearAll[CT9]");
    /* CRC-table-style: `a_. f_` is `Optional[a_] · f_`.  Without
     * Condition this already worked; the bug was that adding a
     * Condition broke it. */
    eval_and_drop("CT9[a_. f_, x_] /; FreeQ[a, x] := {a, f, x, \"hit\"}");
    /* Explicit factor: 2 y → a=2, f=y. */
    assert_eval_eq("CT9[2 y, y]", "{2, y, y, \"hit\"}", 0);
    assert_eval_eq("CT9[3 Sin[y], y]", "{3, Sin[y], y, \"hit\"}", 0);
    /* Symbol factor: a y → a=a, f=y. FreeQ[a, y] is True. */
    assert_eval_eq("CT9[a y, y]", "{a, y, y, \"hit\"}", 0);
}

static void test_cond_with_optional_inside_compound(void) {
    eval_and_drop("ClearAll[CT10]");
    /* Sin[a_. x_] — the canonical CRC pattern (Formula 446-ish). */
    eval_and_drop("CT10[Sin[a_. x_], x_] /; FreeQ[a, x] := -Cos[a x]/a");
    assert_eval_eq("CT10[Sin[2 y], y]", "-1/2 Cos[2 y]", 0);
    assert_eval_eq("CT10[Sin[3 z], z]", "-1/3 Cos[3 z]", 0);
}

/* ----- 5. Nested Condition ---------------------------------------- */

static void test_cond_nested_both_true(void) {
    eval_and_drop("ClearAll[CT11]");
    eval_and_drop("CT11[a_, b_] /; True /; True := {a, b, \"hit\"}");
    assert_eval_eq("CT11[1, 2]", "{1, 2, \"hit\"}", 0);
}

static void test_cond_nested_one_false(void) {
    eval_and_drop("ClearAll[CT12]");
    eval_and_drop("CT12[a_, b_] /; True /; False := {a, b, \"hit\"}");
    ASSERT(result_head_is("CT12[1, 2]", "CT12"));
}

/* ----- 6. PatternTest at the LHS root ----------------------------- */

static void test_pattern_test_inside_slot(void) {
    /* PatternTest INSIDE a slot (a_?NumberQ) already worked before the
     * fix because pattern_arg_head_canon's recursion handles it.  We
     * keep the test as a regression guard. */
    eval_and_drop("ClearAll[CT13]");
    eval_and_drop("CT13[a_?NumberQ, b_] := {a, b, \"hit\"}");
    assert_eval_eq("CT13[42, y]", "{42, y, \"hit\"}", 0);
    /* x is symbolic — PatternTest fails, rule doesn't fire. */
    ASSERT(result_head_is("CT13[x, y]", "CT13"));
}

/* ----- 7. Dispatch invariants (must still hold) ------------------- */

static void test_dispatch_first_arg_head_still_filters(void) {
    /* The dispatch filter must continue to skip rules whose literal
     * first-arg head doesn't match the input.  This is the perf win
     * we don't want to lose. */
    eval_and_drop("ClearAll[CT14]");
    eval_and_drop("CT14[Sin[x_], y_] /; True := {x, y, \"sin-hit\"}");
    eval_and_drop("CT14[Cos[x_], y_] /; True := {x, y, \"cos-hit\"}");
    assert_eval_eq("CT14[Sin[3], 5]", "{3, 5, \"sin-hit\"}", 0);
    assert_eval_eq("CT14[Cos[3], 5]", "{3, 5, \"cos-hit\"}", 0);
    /* Tan should match neither — bubbles back. */
    ASSERT(result_head_is("CT14[Tan[3], 5]", "CT14"));
}

static void test_dispatch_arity_still_filters(void) {
    eval_and_drop("ClearAll[CT15]");
    eval_and_drop("CT15[a_, b_, c_] /; True := \"three-hit\"");
    eval_and_drop("CT15[a_, b_] /; True := \"two-hit\"");
    assert_eval_eq("CT15[1, 2, 3]", "\"three-hit\"", 0);
    assert_eval_eq("CT15[1, 2]", "\"two-hit\"", 0);
    /* One-arg call: neither rule fits; must bubble back. */
    ASSERT(result_head_is("CT15[1]", "CT15"));
}

/* ----- 8. Cross-check: Replace path was always fine --------------- */

static void test_replace_with_condition_still_works(void) {
    /* Sanity: Replace bypasses the dispatch filter entirely and
     * always invoked the matcher directly.  These already worked
     * pre-fix; we keep them as a regression guard. */
    assert_eval_eq("Replace[5, a_ /; True -> \"hit\"]", "\"hit\"", 0);
    assert_eval_eq("ReplaceAll[Sin[x], Sin[a_] /; FreeQ[a, q] -> a]", "x", 0);
    assert_eval_eq("{1, 2, 3} /. a_Integer /; a > 1 -> \"big\"",
                   "{1, \"big\", \"big\"}", 0);
}

/* ----- 9. CRC table integration smoke ----------------------------- */

static void test_crc_table_formula_fires(void) {
    /* End-to-end: lazy-load the CRC table via Method->"CRCTable" and
     * verify a /;-guarded formula now produces an answer.  This is
     * the user-visible reason the dispatch fix matters.  We use the
     * `Sin[a_. x_]` family (Formula 446 in the .m file). */
    Expr* p = parse_expression(
        "Integrate[Sin[2 x], x, Method -> \"CRCTable\"]");
    Expr* r = evaluate(p);
    expr_free(p);
    ASSERT(r != NULL);
    /* "fired" = the outer head is NOT Integrate any more (the table
     * looked it up successfully). */
    bool still_integrate = (r->type == EXPR_FUNCTION
        && r->data.function.head
        && r->data.function.head->type == EXPR_SYMBOL
        && strcmp(r->data.function.head->data.symbol, "Integrate") == 0);
    if (still_integrate) {
        char* s = expr_to_string(r);
        fprintf(stderr, "CRC formula did not fire: result = %s\n", s);
        free(s);
    }
    ASSERT(!still_integrate);
    expr_free(r);

    /* Polynomial x^3 should also fire via Formula 19 (or similar
     * monomial reducer): `IntegrateTable[x_^n_., x_] /; FreeQ[n, x]
     * && n =!= -1 := x^(n+1)/(n+1)`. */
    ASSERT(!result_head_is(
        "Integrate[x^3, x, Method -> \"CRCTable\"]", "Integrate"));

    /* a x (linear) — Formula 17ish: `IntegrateTable[a_ f_, x_] /;
     * FreeQ[a, x] := a IntegrateTable[f, x]`. */
    ASSERT(!result_head_is(
        "Integrate[a x, x, Method -> \"CRCTable\"]", "Integrate"));
}

/* ----- 9b. §3.4 LHS pattern canonicalization ----------------------
 *
 * Pin the behaviour that DownValue LHSs are canonicalized at insertion
 * time so the §3.5 dispatch filter and the matcher see the same shape
 * that the evaluator produces for runtime inputs.  Specifically:
 *
 *   - `1/(x_ Sqrt[a_ + b_. x_ + c_. x_^2])` parses to
 *     `Power[Times[x_, Sqrt[...]], -1]`, but every input of that
 *     mathematical form is evaluated to `Times[Power[x,-1],
 *     Power[Sqrt[...],-1]]`.  Without canonicalization the dispatch
 *     filter cached `Power` and rejected `Times` before the matcher
 *     even ran (~80 CRC integral formulas silently failed).
 *   - `Power[Power[base, e1], e2]` with numeric exponents must collapse
 *     to `Power[base, e1*e2]` (matches the evaluator's nested-Power
 *     normalization).
 */
static void test_canon_power_times_lhs_dispatches(void) {
    /* CRC formula 246 shape — the LHS that previously failed to
     * dispatch.  Use a synthetic head `CTcn1` so we don't depend on
     * the CRC table being loaded. */
    eval_and_drop("ClearAll[CTcn1]");
    eval_and_drop(
        "CTcn1[1/(x_ Sqrt[a_ + b_. x_ + c_. x_^2]), x_] /; "
        "FreeQ[{a,b,c}, x] && a > 0 := {a, b, c, \"FIRED\"}");
    /* The rule MUST fire after canonicalization.  Pre-fix this returned
     * the original CTcn1[...] unchanged. */
    ASSERT(rule_did_fire(
        "CTcn1[1/(x Sqrt[2 + 3 x + 5 x^2]), x]", "CTcn1"));
}

static void test_canon_power_of_power_lhs_dispatches(void) {
    /* `Sqrt[Sqrt[x_]]` parses to Power[Power[x_, 1/2], 1/2]; the
     * evaluator collapses this to Power[x_, 1/4] for any concrete
     * input.  Without canonicalization the stored pattern would keep
     * the nested form and never match. */
    eval_and_drop("ClearAll[CTcn2]");
    eval_and_drop(
        "CTcn2[Sqrt[Sqrt[x_]]] := {x, \"FOLDED\"}");
    ASSERT(rule_did_fire("CTcn2[Sqrt[Sqrt[5]]]", "CTcn2"));
    /* And an input written directly as a quarter root must hit the
     * same canonical form. */
    ASSERT(rule_did_fire("CTcn2[5^(1/4)]", "CTcn2"));
}

static void test_canon_does_not_alter_pattern_semantics(void) {
    /* Pattern wrappers must NOT be folded.  HoldPattern explicitly
     * opts out of canonicalization, so a pattern wrapped in
     * HoldPattern keeps its parser shape and still matches inputs of
     * the same parser shape (this is exactly the role HoldPattern
     * exists to play). */
    eval_and_drop("ClearAll[CTcn3]");
    eval_and_drop(
        "CTcn3[HoldPattern[1/(x_ y_)]] := {x, y, \"HELD\"}");
    /* The held LHS literally requires the input to arrive as
     * Power[Times[x,y],-1].  Inputs that auto-distribute won't match;
     * this is fine — HoldPattern explicitly preserves literal shape. */
    ASSERT(!rule_did_fire("CTcn3[1/(p q)]", "CTcn3"));

    /* /;-guards still evaluate at match time with bindings, not at
     * rule-insertion time.  This rule's guard `x > 0` would error if
     * canonicalization eagerly evaluated it (x would be unbound), but
     * because we leave Condition guards alone it works correctly and
     * fires for positive inputs only. */
    eval_and_drop("ClearAll[CTcn4]");
    eval_and_drop("CTcn4[x_] /; x > 0 := \"POS\"");
    ASSERT(rule_did_fire("CTcn4[5]", "CTcn4"));
    ASSERT(!rule_did_fire("CTcn4[-5]", "CTcn4"));
}

/* ----- 10. Negative cases for the fix ----------------------------- */

static void test_cond_false_with_compound_first_arg(void) {
    /* The fix must NOT cause Condition-False rules to spuriously fire. */
    eval_and_drop("ClearAll[CT16]");
    eval_and_drop("CT16[Sin[x_], y_] /; False := \"hit\"");
    ASSERT(result_head_is("CT16[Sin[3], 5]", "CT16"));
}

static void test_cond_freeq_false_blocks_compound(void) {
    eval_and_drop("ClearAll[CT17]");
    eval_and_drop("CT17[Sin[a_. x_], x_] /; FreeQ[a, x] := -Cos[a x]/a");
    /* a binds to x in Sin[x x] — FreeQ[a, x] is False here.  Actually
     * Sin[x x] auto-evaluates Times[x, x] to x^2; let's use a clearer
     * negative test. */
    /* Use Sin[x^2] — doesn't match Sin[a_. x_] structurally (because
     * the inner is x^2, not a_ x_). */
    ASSERT(result_head_is(
        "CT17[Sin[x^2], x]", "CT17"));
}

/* ----- driver ----------------------------------------------------- */

int main(void) {
    symtab_init();
    core_init();

    TEST(test_cond_single_arg_true);
    TEST(test_cond_single_arg_false);
    TEST(test_cond_single_arg_freeq);
    TEST(test_cond_two_arg_true);
    TEST(test_cond_two_arg_freeq);
    TEST(test_cond_three_arg);
    TEST(test_cond_compound_first_arg_true);
    TEST(test_cond_compound_first_arg_freeq);
    TEST(test_cond_with_optional_slot);
    TEST(test_cond_with_optional_inside_compound);
    TEST(test_cond_nested_both_true);
    TEST(test_cond_nested_one_false);
    TEST(test_pattern_test_inside_slot);
    TEST(test_dispatch_first_arg_head_still_filters);
    TEST(test_dispatch_arity_still_filters);
    TEST(test_replace_with_condition_still_works);
    TEST(test_canon_power_times_lhs_dispatches);
    TEST(test_canon_power_of_power_lhs_dispatches);
    TEST(test_canon_does_not_alter_pattern_semantics);
    TEST(test_cond_false_with_compound_first_arg);
    TEST(test_cond_freeq_false_blocks_compound);
    TEST(test_crc_table_formula_fires);

    printf("All Condition-DownValue dispatch tests passed!\n");
    return 0;
}
