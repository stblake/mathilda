/* test_rule_dispatch.c -- M4 (§3.5 + §3.6) tests.
 *
 * §3.5 invariants:
 *   - Each Rule carries a precomputed dispatch_arity and
 *     first_arg_head_canon. Variable-arity patterns (top-level
 *     BlankSequence / Optional / Repeated / OptionsPattern) get
 *     dispatch_arity = -1 and first_arg_head_canon = NULL so they
 *     always pass the filter.
 *   - apply_down_values uses these fields to skip rules whose top-level
 *     shape cannot match the input before calling the matcher. Skipped
 *     rules must never affect evaluation -- they are unreachable by
 *     definition.
 *
 * §3.6 invariants:
 *   - Rules in down_values are stored in descending specificity order.
 *     Specificity is a function of the pattern alone (literals beat
 *     blanks beat sequence blanks, with extra credit for Condition /
 *     PatternTest constraints).
 *   - Insertion order is the tie-breaker for equal specificity.
 *   - Re-defining the same pattern replaces the RHS in place; the
 *     rule's position in the list is unchanged.
 *
 * Functional behaviour:
 *   - Specific patterns fire before general patterns of the same head,
 *     even when the general one was inserted first.
 *   - Many DownValues on a single head still dispatch correctly under
 *     the filtered scan.
 */

#include "core.h"
#include "expr.h"
#include "symtab.h"
#include "eval.h"
#include "parse.h"
#include "print.h"
#include "sym_intern.h"
#include "test_utils.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdbool.h>

/* Helper: count rules in a list. */
static size_t rule_list_length(Rule* r) {
    size_t n = 0;
    while (r) { n++; r = r->next; }
    return n;
}

/* Helper: nth rule in list (0-indexed); returns NULL if out of range. */
static Rule* rule_at(Rule* r, size_t idx) {
    while (r && idx > 0) { r = r->next; idx--; }
    return r;
}

/* Helper: parse a pattern -> replacement pair and call symtab_add_down_value. */
static void add_down(const char* head, const char* pattern_src, const char* repl_src) {
    Expr* p = parse_expression(pattern_src);
    Expr* r = parse_expression(repl_src);
    symtab_add_down_value(head, p, r);
    expr_free(p);
    expr_free(r);
}

/* Convenience: evaluate a string and return its printed form. Caller frees. */
static char* eval_str(const char* src) {
    Expr* parsed = parse_expression(src);
    Expr* result = evaluate(parsed);
    expr_free(parsed);
    char* s = expr_to_string(result);
    expr_free(result);
    return s;
}

/* ------------------------------------------------------------------ */
/* 1. Dispatch metadata is computed correctly.                         */
/* ------------------------------------------------------------------ */

static void test_dispatch_meta_fixed_arity_no_head(void) {
    symtab_clear_symbol("f1");
    /* f1[x_, y_] := x + y -- arity 2, first arg wildcard head */
    add_down("f1", "f1[x_, y_]", "x + y");

    Rule* r = symtab_get_down_values("f1");
    ASSERT(r != NULL);
    ASSERT(r->dispatch_arity == 2);
    ASSERT(r->first_arg_head_canon == NULL);
}

static void test_dispatch_meta_typed_blank(void) {
    symtab_clear_symbol("f2");
    /* f2[n_Integer] := n^2 -- arity 1, first arg head Integer */
    add_down("f2", "f2[n_Integer]", "n^2");

    Rule* r = symtab_get_down_values("f2");
    ASSERT(r != NULL);
    ASSERT(r->dispatch_arity == 1);
    ASSERT(r->first_arg_head_canon == intern_symbol("Integer"));
}

static void test_dispatch_meta_function_head(void) {
    symtab_clear_symbol("f3");
    /* f3[Sin[x_], y_] := y -- arity 2, first arg head Sin */
    add_down("f3", "f3[Sin[x_], y_]", "y");

    Rule* r = symtab_get_down_values("f3");
    ASSERT(r != NULL);
    ASSERT(r->dispatch_arity == 2);
    ASSERT(r->first_arg_head_canon == intern_symbol("Sin"));
}

static void test_dispatch_meta_sequence_pattern(void) {
    symtab_clear_symbol("f4");
    /* f4[x__] := x -- variable arity, wildcard head */
    add_down("f4", "f4[x__]", "x");

    Rule* r = symtab_get_down_values("f4");
    ASSERT(r != NULL);
    ASSERT(r->dispatch_arity == -1);
    ASSERT(r->first_arg_head_canon == NULL);
}

static void test_dispatch_meta_optional(void) {
    symtab_clear_symbol("f5");
    /* f5[x_, y_:0] := x + y -- variable arity (Optional makes it 1 or 2) */
    add_down("f5", "f5[x_, y_:0]", "x + y");

    Rule* r = symtab_get_down_values("f5");
    ASSERT(r != NULL);
    ASSERT(r->dispatch_arity == -1);
}

static void test_dispatch_meta_nested_optional_head(void) {
    /* Regression: a function-call slot whose head is fixed-arity but whose
     * argument is Optional must NOT specialise on that head. An omitted
     * Optional collapses the call to its base form, so a concrete value
     * matching the slot need not be headed by the call's head:
     *   Power[b_, n_.] matches a bare `b` as b^1   (head Symbol, not Power)
     *   Times[a_., x_] matches `x`                 (head Symbol, not Times)
     * If first_arg_head_canon were "Power"/"Times", the §3.5 filter would
     * skip the rule for those collapsed inputs. */
    symtab_clear_symbol("f5b");
    /* f5b[Power[b_, n_.]] := {b, n} -- arity 1, first slot must be wildcard. */
    add_down("f5b", "f5b[Power[b_, n_.]]", "{b, n}");
    Rule* r = symtab_get_down_values("f5b");
    ASSERT(r != NULL);
    ASSERT(r->dispatch_arity == 1);
    ASSERT(r->first_arg_head_canon == NULL);

    symtab_clear_symbol("f5c");
    /* f5c[Times[a_., x_]] := {a, x} -- same, Optional in a Times slot. */
    add_down("f5c", "f5c[Times[a_., x_]]", "{a, x}");
    r = symtab_get_down_values("f5c");
    ASSERT(r != NULL);
    ASSERT(r->dispatch_arity == 1);
    ASSERT(r->first_arg_head_canon == NULL);

    /* Control: NO Optional -> head must still be specialised. */
    symtab_clear_symbol("f5d");
    add_down("f5d", "f5d[Power[b_, n_]]", "{b, n}");
    r = symtab_get_down_values("f5d");
    ASSERT(r != NULL);
    ASSERT(r->dispatch_arity == 1);
    ASSERT(r->first_arg_head_canon == intern_symbol("Power"));
}

static void test_functional_nested_optional_collapses(void) {
    /* The rule must actually FIRE on a collapsed-Optional input that does
     * not share the pattern's call head -- this is what the filter fix
     * unblocks end-to-end. */
    symtab_clear_symbol("fp");
    add_down("fp", "fp[Power[b_, n_.]]", "{b, n}");

    char* s = eval_str("fp[k]");          /* k as k^1 */
    ASSERT_STR_EQ(s, "{k, 1}");
    free(s);

    s = eval_str("fp[w^3]");              /* explicit Power still matches */
    ASSERT_STR_EQ(s, "{w, 3}");
    free(s);

    symtab_clear_symbol("ft");
    add_down("ft", "ft[Times[a_., x_]]", "{a, x}");

    s = eval_str("ft[z]");                /* z as 1*z */
    ASSERT_STR_EQ(s, "{1, z}");
    free(s);
}

static void test_dispatch_meta_literal_int(void) {
    symtab_clear_symbol("f6");
    /* f6[0] := 0 -- arity 1, first arg head Integer (literal atom) */
    add_down("f6", "f6[0]", "0");

    Rule* r = symtab_get_down_values("f6");
    ASSERT(r != NULL);
    ASSERT(r->dispatch_arity == 1);
    ASSERT(r->first_arg_head_canon == intern_symbol("Integer"));
}

static void test_dispatch_meta_holdpattern_transparency(void) {
    symtab_clear_symbol("f7");
    /* HoldPattern[f7[Plus[a_, b_]]] := a*b -- arity 1, first arg head Plus */
    add_down("f7", "HoldPattern[f7[Plus[a_, b_]]]", "a*b");

    Rule* r = symtab_get_down_values("f7");
    ASSERT(r != NULL);
    ASSERT(r->dispatch_arity == 1);
    ASSERT(r->first_arg_head_canon == intern_symbol("Plus"));
}

/* ------------------------------------------------------------------ */
/* 2. Specificity ordering -- §3.6.                                    */
/* ------------------------------------------------------------------ */

static void test_specificity_literal_beats_blank(void) {
    symtab_clear_symbol("g1");
    /* Insert the GENERAL rule first; specificity sort must move the
     * specific one ahead of it. */
    add_down("g1", "g1[x_]", "g1[x_]_general");
    add_down("g1", "g1[0]", "g1[0]_specific");

    Rule* head = symtab_get_down_values("g1");
    ASSERT(rule_list_length(head) == 2);
    /* The specific (literal-arg) rule must come first. */
    Rule* first = head;
    ASSERT(first->specificity > first->next->specificity);
    /* Sanity: pattern of first rule is the literal one. */
    Expr* zero = parse_expression("g1[0]");
    ASSERT(expr_eq(first->pattern, zero));
    expr_free(zero);
}

static void test_specificity_typed_beats_blank(void) {
    symtab_clear_symbol("g2");
    add_down("g2", "g2[x_]", "general");
    add_down("g2", "g2[n_Integer]", "typed");

    Rule* head = symtab_get_down_values("g2");
    ASSERT(rule_list_length(head) == 2);
    /* Typed Blank[Integer] is more specific than plain Blank[]. */
    ASSERT(head->specificity > head->next->specificity);
}

static void test_specificity_blank_beats_blanksequence(void) {
    symtab_clear_symbol("g3");
    add_down("g3", "g3[x__]", "seq");
    add_down("g3", "g3[x_]", "single");

    Rule* head = symtab_get_down_values("g3");
    ASSERT(rule_list_length(head) == 2);
    /* Single blank is more specific than a sequence blank. */
    ASSERT(head->specificity > head->next->specificity);
}

static void test_specificity_condition_adds_score(void) {
    symtab_clear_symbol("g4");
    /* Plain x_ vs x_ /; cond: a condition does not raise the score (the
     * reference does not promote a guarded rule past a more specific one),
     * but between two rules that differ only in a guard the guarded one goes
     * first. */
    add_down("g4", "g4[x_]", "plain");
    add_down("g4", "g4[x_ /; x > 0]", "with_cond");

    Rule* head = symtab_get_down_values("g4");
    ASSERT(rule_list_length(head) == 2);
    ASSERT(head->specificity == head->next->specificity);
    /* the guarded rule is the first one: the plain g4[x_] comes second */
    Expr* plain = parse_expression("g4[x_]");
    ASSERT(!expr_eq(head->pattern, plain));
    ASSERT(expr_eq(head->next->pattern, plain));
    expr_free(plain);
}

static void test_specificity_insertion_order_tie_break(void) {
    symtab_clear_symbol("g5");
    /* Three rules with identical specificity (one blank, one literal each)
     * that are NOT the same up to renaming, so none replaces another.
     * (g5[a_], g5[b_], g5[c_] would be one rule: see the alpha test below.) */
    add_down("g5", "g5[a_, 1]", "first");
    add_down("g5", "g5[b_, 2]", "second");
    add_down("g5", "g5[c_, 3]", "third");

    Rule* head = symtab_get_down_values("g5");
    ASSERT(rule_list_length(head) == 3);
    /* All three have the same specificity, so insertion order should hold. */
    ASSERT(head->specificity == head->next->specificity);
    ASSERT(head->next->specificity == head->next->next->specificity);

    /* Verify the patterns are in insertion order. */
    Expr* p0 = parse_expression("g5[a_, 1]");
    Expr* p1 = parse_expression("g5[b_, 2]");
    Expr* p2 = parse_expression("g5[c_, 3]");
    ASSERT(expr_eq(rule_at(head, 0)->pattern, p0));
    ASSERT(expr_eq(rule_at(head, 1)->pattern, p1));
    ASSERT(expr_eq(rule_at(head, 2)->pattern, p2));
    expr_free(p0); expr_free(p1); expr_free(p2);
}

/* The list is ordered by pattern generality, as in the reference: a new rule
 * goes ahead of a stored one only when the stored pattern matches everything
 * the new one does; incomparable patterns keep definition order whatever their
 * scores; pattern-free definitions lead; of two patterns that cover each other
 * the newer goes first.  Every expected order was read off the reference. */
static void expect_order(const char* defs, const char* want) {
    char* got = eval_str(defs);
    if (strcmp(got, want) != 0) {
        printf("  %s\n    want %s\n    got  %s\n", defs, want, got);
        ASSERT(0);
    }
    free(got);
}

static void test_generality_order(void) {
    /* FeynCalc DiracTrick: neither covers the other, so the guarded,
     * earlier rule stays first (it had been moved behind on a score). */
    expect_order("o1[h[b___G, di_G, c__G]] := 1 /; c1 /; c2; "
                 "o1[h[b___, G[5], c:G[_[__]].., d___]] := 2; Last /@ DownValues[o1]",
                 "{1, 2}");
    expect_order("o2[x_] := 1; o2[x_Integer] := 2; o2[1] := 3; o2[x_String] := 4; "
                 "o2[x_?EvenQ] := 5; Last /@ DownValues[o2]", "{3, 2, 4, 5, 1}");
    expect_order("o3[x__] := 1; o3[x_, y_] := 2; o3[x_, y__] := 3; o3[1, y_] := 4; "
                 "Last /@ DownValues[o3]", "{4, 2, 3, 1}");
    expect_order("o4[x_Integer] := 1; o4[x_?EvenQ] := 2; Last /@ DownValues[o4]", "{1, 2}");
    expect_order("o5[x_Integer] := 1; o5[x_Integer?EvenQ] := 2; Last /@ DownValues[o5]", "{2, 1}");
    expect_order("o6[h[x___]] := 1; o6[h[x___G]] := 2; Last /@ DownValues[o6]", "{2, 1}");
    expect_order("o7[h[x___G]] := 1; o7[h[G[1], x___]] := 2; Last /@ DownValues[o7]", "{1, 2}");
    expect_order("o8[h[x__, y_]] := 1; o8[h[x_, y__]] := 2; Last /@ DownValues[o8]", "{1, 2}");
    expect_order("o9[x_?NumericQ] := 1; o9[1] := 2; Last /@ DownValues[o9]", "{2, 1}");
    expect_order("o10[h[x__]] := 1; o10[h[]] := 2; Last /@ DownValues[o10]", "{2, 1}");
    expect_order("o11[x_] := 1; o11[x_^n_.] := 2; Last /@ DownValues[o11]", "{2, 1}");
    expect_order("o12[x_^n_.] := 1; o12[x_] := 2; Last /@ DownValues[o12]", "{2, 1}");
    expect_order("o13[x:(_)..] := 1; o13[x__] := 2; Last /@ DownValues[o13]", "{1, 2}");
    expect_order("o14[x_ /; x > 0] := 1; o14[x_ /; x < 0] := 2; Last /@ DownValues[o14]", "{1, 2}");
    expect_order("o15[x_, OptionsPattern[]] := 1; o15[x_, a -> 1] := 2; "
                 "Last /@ DownValues[o15]", "{2, 1}");
}

/* Pattern-free definitions are found through the literal index: memoisation
 * is linear, and redefinition, removal, Block, DownValues[f] = ... and
 * Flat/Orderless heads behave as with the scan (values from the reference). */
static void test_literal_index(void) {
    char* t = eval_str("ClearAll[mf]; mf[n_] := mf[n] = n; "
                       "First[AbsoluteTiming[Do[mf[i], {i, 20000}]]] < 2");
    ASSERT_STR_EQ(t, "True"); free(t);
    t = eval_str("Length[DownValues[mf]]"); ASSERT_STR_EQ(t, "20001"); free(t);
    t = eval_str("ClearAll[lf]; lf[x_] := g[x]; lf[1] = 10; lf[2] = 20; lf[1] =.; "
                 "lf[2] = 21; {lf[1], lf[2], lf[3], Length[DownValues[lf]]}");
    ASSERT_STR_EQ(t, "{g[1], 21, g[3], 2}"); free(t);
    /* Block's value is evaluated again after lf's rules are restored, so the
     * lf[2] left unevaluated inside becomes 21 (checked on 13.2 and 14; the
     * earlier lf[2] was a stale fixed-point stamp surviving the restore). */
    t = eval_str("Block[{lf}, lf[1] = 99; {lf[1], lf[2]}]");
    ASSERT_STR_EQ(t, "{99, 21}"); free(t);
    t = eval_str("{lf[1], lf[2]}"); ASSERT_STR_EQ(t, "{g[1], 21}"); free(t);
    t = eval_str("DownValues[lf] = {HoldPattern[lf[5]] :> 55, HoldPattern[lf[y_]] :> h[y]}; "
                 "lf[7] = 77; {lf[5], lf[7], lf[2]}");
    ASSERT_STR_EQ(t, "{55, 77, h[2]}"); free(t);
    t = eval_str("ClearAll[lk]; lk[1] := 1 /; False; lk[1] := 2; lk[x_] := 3; {lk[1], lk[2]}");
    ASSERT_STR_EQ(t, "{2, 3}"); free(t);
    t = eval_str("ClearAll[lo]; SetAttributes[lo, Orderless]; lo[a, b] = 1; lo[b, a]");
    ASSERT_STR_EQ(t, "1"); free(t);
    t = eval_str("ClearAll[lp]; lp[{1, 2}] = 3; lp[1/2] = 6; {lp[{1, 2}], lp[1/2], lp[0.5]}");
    ASSERT_STR_EQ(t, "{3, 6, lp[0.5]}"); free(t);
}

static void test_specificity_overwrite_keeps_position(void) {
    symtab_clear_symbol("g6");
    add_down("g6", "g6[a_, 1]", "v1");
    add_down("g6", "g6[b_, 2]", "v2");
    add_down("g6", "g6[a_, 1]", "v1_new"); /* overwrite first rule's RHS */

    Rule* head = symtab_get_down_values("g6");
    ASSERT(rule_list_length(head) == 2);
    /* Order preserved: g6[a_, 1] still at index 0 (was already there). */
    Expr* pa = parse_expression("g6[a_, 1]");
    ASSERT(expr_eq(rule_at(head, 0)->pattern, pa));
    expr_free(pa);
}

/* A left-hand side that is the same up to renaming of its pattern variables
 * replaces the rule -- pattern and right-hand side, since the RHS refers to the
 * new names -- in place; a different condition keeps both rules. */
static void test_alpha_equivalent_redefinition_replaces(void) {
    symtab_clear_symbol("g7");
    add_down("g7", "g7[a_, z]", "v1");
    add_down("g7", "g7[b_, 9]", "other");
    add_down("g7", "g7[c_, z]", "v2");          /* alpha-equivalent to the first */
    Rule* head = symtab_get_down_values("g7");
    ASSERT(rule_list_length(head) == 2);
    Expr* pc = parse_expression("g7[c_, z]");
    ASSERT(expr_eq(rule_at(head, 0)->pattern, pc));   /* new names, same slot */
    expr_free(pc);

    char* s = eval_str("g7[1, z]");
    ASSERT_STR_EQ(s, "v2");
    free(s);
    s = eval_str("g7[x_] := x^2; g7[y_] := y^3; g7[2]");
    ASSERT_STR_EQ(s, "8");
    free(s);
    s = eval_str("g7[u_] := pos /; u > 0; g7[w_] := neg /; w < 0; {g7[1], g7[-1]}");
    ASSERT_STR_EQ(s, "{pos, neg}");
    free(s);
}

/* ------------------------------------------------------------------ */
/* 3. Functional dispatch -- specific rules win over general ones.     */
/* ------------------------------------------------------------------ */

static void test_functional_specific_wins(void) {
    symtab_clear_symbol("h1");
    /* General first, specific second -- specific must still fire. */
    add_down("h1", "h1[x_]", "general");
    add_down("h1", "h1[0]", "zero");

    char* s = eval_str("h1[0]");
    ASSERT_STR_EQ(s, "zero");
    free(s);

    s = eval_str("h1[5]");
    ASSERT_STR_EQ(s, "general");
    free(s);
}

static void test_functional_arity_dispatch(void) {
    symtab_clear_symbol("h2");
    add_down("h2", "h2[x_]", "one_arg");
    add_down("h2", "h2[x_, y_]", "two_args");
    add_down("h2", "h2[x_, y_, z_]", "three_args");

    char* s = eval_str("h2[a]");
    ASSERT_STR_EQ(s, "one_arg");
    free(s);

    s = eval_str("h2[a, b]");
    ASSERT_STR_EQ(s, "two_args");
    free(s);

    s = eval_str("h2[a, b, c]");
    ASSERT_STR_EQ(s, "three_args");
    free(s);
}

static void test_functional_head_dispatch(void) {
    symtab_clear_symbol("h3");
    add_down("h3", "h3[Sin[x_]]", "got_sin");
    add_down("h3", "h3[Cos[x_]]", "got_cos");
    add_down("h3", "h3[x_]", "fallback");

    char* s = eval_str("h3[Sin[a]]");
    ASSERT_STR_EQ(s, "got_sin");
    free(s);

    s = eval_str("h3[Cos[a]]");
    ASSERT_STR_EQ(s, "got_cos");
    free(s);

    s = eval_str("h3[Tan[a]]");
    ASSERT_STR_EQ(s, "fallback");
    free(s);

    /* Atomic input falls through both head-specialized rules. */
    s = eval_str("h3[42]");
    ASSERT_STR_EQ(s, "fallback");
    free(s);
}

static void test_functional_sequence_pattern_still_fires(void) {
    symtab_clear_symbol("h4");
    /* Variable-arity rule must NOT be skipped by the dispatch filter. */
    add_down("h4", "h4[x__]", "got_seq");

    char* s = eval_str("h4[1, 2, 3, 4, 5]");
    ASSERT_STR_EQ(s, "got_seq");
    free(s);

    s = eval_str("h4[a]");
    ASSERT_STR_EQ(s, "got_seq");
    free(s);
}

static void test_functional_typed_blank_dispatch(void) {
    symtab_clear_symbol("h5");
    add_down("h5", "h5[n_Integer]", "is_int");
    add_down("h5", "h5[r_Real]", "is_real");
    add_down("h5", "h5[s_Symbol]", "is_sym");
    add_down("h5", "h5[x_]", "any");

    char* s = eval_str("h5[5]");
    ASSERT_STR_EQ(s, "is_int");
    free(s);

    s = eval_str("h5[3.14]");
    ASSERT_STR_EQ(s, "is_real");
    free(s);

    s = eval_str("h5[xyz]");
    ASSERT_STR_EQ(s, "is_sym");
    free(s);

    /* Function input falls through to the wildcard. */
    s = eval_str("h5[Plus[a, b]]");
    ASSERT_STR_EQ(s, "any");
    free(s);
}

/* ------------------------------------------------------------------ */
/* 4. Wildcard rules are never filtered out.                           */
/* ------------------------------------------------------------------ */

static void test_filter_skips_wildcard_arg_rules(void) {
    symtab_clear_symbol("h6");
    /* General rule first, several head-specialised rules after.
     * Specific rules should fire when applicable; the general rule is
     * the fallback for everything else. */
    add_down("h6", "h6[x_]", "fallback");
    add_down("h6", "h6[Sin[x_]]", "specific_sin");
    add_down("h6", "h6[Cos[x_]]", "specific_cos");

    char* s = eval_str("h6[Sin[q]]");
    ASSERT_STR_EQ(s, "specific_sin");
    free(s);

    s = eval_str("h6[Cos[q]]");
    ASSERT_STR_EQ(s, "specific_cos");
    free(s);

    s = eval_str("h6[Tan[q]]");
    ASSERT_STR_EQ(s, "fallback");
    free(s);
}

/* ------------------------------------------------------------------ */
/* 5. Stress: many DownValues on one symbol.                           */
/* ------------------------------------------------------------------ */

static void test_heavy_downvalue_dispatch(void) {
    symtab_clear_symbol("hv");
    /* Add 100 specific rules: hv[k] -> Integer[k+1000] for k in 0..99.
     * Then a wildcard fallback. */
    for (int k = 0; k < 100; k++) {
        char pat[32], repl[32];
        snprintf(pat, sizeof(pat), "hv[%d]", k);
        snprintf(repl, sizeof(repl), "%d", k + 1000);
        add_down("hv", pat, repl);
    }
    add_down("hv", "hv[x_]", "general");

    /* Each numeric input must hit its specific rule, NOT the general one. */
    for (int k = 0; k < 100; k += 7) {
        char inp[32], exp[32];
        snprintf(inp, sizeof(inp), "hv[%d]", k);
        snprintf(exp, sizeof(exp), "%d", k + 1000);
        char* s = eval_str(inp);
        ASSERT_STR_EQ(s, exp);
        free(s);
    }
    /* Symbolic input falls through to general. */
    char* s = eval_str("hv[abc]");
    ASSERT_STR_EQ(s, "general");
    free(s);
}

/* ------------------------------------------------------------------ */
/* 6. Cross-arity rules don't fire for wrong-arity inputs.             */
/* ------------------------------------------------------------------ */

static void test_arity_mismatch_does_not_fire(void) {
    symtab_clear_symbol("h7");
    add_down("h7", "h7[x_]", "one");
    /* Calling h7[a, b] should leave it unevaluated -- no rule matches. */
    char* s = eval_str("h7[a, b]");
    ASSERT_STR_EQ(s, "h7[a, b]");
    free(s);

    /* And h7[a] should evaluate. */
    s = eval_str("h7[a]");
    ASSERT_STR_EQ(s, "one");
    free(s);
}

/* ------------------------------------------------------------------ */
/* 7. Filter is sound for all atomic input head classes.               */
/* ------------------------------------------------------------------ */

static void test_filter_sound_for_all_atom_kinds(void) {
    symtab_clear_symbol("h8");
    /* Only a Sin-headed first arg should fire. Other inputs (atoms of
     * various kinds, other functions) should leave the call
     * unevaluated. */
    add_down("h8", "h8[Sin[x_]]", "got_sin");

    /* Function input matching head: fires. */
    char* s = eval_str("h8[Sin[1]]");
    ASSERT_STR_EQ(s, "got_sin");
    free(s);

    /* Integer input: filtered out, no rule fires. */
    s = eval_str("h8[5]");
    ASSERT_STR_EQ(s, "h8[5]");
    free(s);

    /* Real input. */
    s = eval_str("h8[3.14]");
    ASSERT_STR_EQ(s, "h8[3.14]");
    free(s);

    /* Symbol input. */
    s = eval_str("h8[xyz]");
    ASSERT_STR_EQ(s, "h8[xyz]");
    free(s);

    /* Wrong-headed function input. */
    s = eval_str("h8[Cos[a]]");
    ASSERT_STR_EQ(s, "h8[Cos[a]]");
    free(s);
}

/* ------------------------------------------------------------------ */
/* 8. Specific-rule preference under the order user actually expects.  */
/* ------------------------------------------------------------------ */

static void test_classical_recursion_dispatch(void) {
    /* Classical fib definition: base case must beat recursive case
     * regardless of insertion order. */
    symtab_clear_symbol("myfib");
    /* Insert recursive case first. */
    add_down("myfib", "myfib[n_]", "myfib[n - 1] + myfib[n - 2]");
    add_down("myfib", "myfib[0]", "0");
    add_down("myfib", "myfib[1]", "1");

    char* s = eval_str("myfib[5]");
    ASSERT_STR_EQ(s, "5");
    free(s);

    s = eval_str("myfib[10]");
    ASSERT_STR_EQ(s, "55");
    free(s);
}

int main(void) {
    symtab_init();
    core_init();

    TEST(test_dispatch_meta_fixed_arity_no_head);
    TEST(test_dispatch_meta_typed_blank);
    TEST(test_dispatch_meta_function_head);
    TEST(test_dispatch_meta_sequence_pattern);
    TEST(test_dispatch_meta_optional);
    TEST(test_dispatch_meta_nested_optional_head);
    TEST(test_dispatch_meta_literal_int);
    TEST(test_dispatch_meta_holdpattern_transparency);

    TEST(test_specificity_literal_beats_blank);
    TEST(test_specificity_typed_beats_blank);
    TEST(test_specificity_blank_beats_blanksequence);
    TEST(test_specificity_condition_adds_score);
    TEST(test_specificity_insertion_order_tie_break);
    TEST(test_specificity_overwrite_keeps_position);
    TEST(test_generality_order);
    TEST(test_literal_index);
    TEST(test_alpha_equivalent_redefinition_replaces);

    TEST(test_functional_specific_wins);
    TEST(test_functional_arity_dispatch);
    TEST(test_functional_head_dispatch);
    TEST(test_functional_sequence_pattern_still_fires);
    TEST(test_functional_typed_blank_dispatch);
    TEST(test_functional_nested_optional_collapses);

    TEST(test_filter_skips_wildcard_arg_rules);
    TEST(test_heavy_downvalue_dispatch);
    TEST(test_arity_mismatch_does_not_fire);
    TEST(test_filter_sound_for_all_atom_kinds);
    TEST(test_classical_recursion_dispatch);

    printf("All test_rule_dispatch tests passed!\n");
    symtab_clear();
    intern_clear();
    return 0;
}
