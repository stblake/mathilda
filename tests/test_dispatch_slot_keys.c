/* test_dispatch_slot_keys.c
 *
 * Pins three FeynCalc-driven speed changes whose risk is a silently skipped
 * rule or a changed match set.  Every expected value below was computed by
 * the reference kernel.
 *
 *  1. Per-slot DownValue dispatch keys (symtab.c, Rule.slot_keys): a rule is
 *     skipped before the matcher when an argument's head (or call arity) can
 *     not fit its slot -- by position, or for an Orderless head by some
 *     assignment.  FeynCalc's Pair (Orderless, 36 rules) tried every rule on
 *     every call because only the first argument was keyed, and never for an
 *     Orderless head.  The tests check that rules still fire through every
 *     shape the key reads: Orderless permutations, pattern-named head
 *     alternatives, call arities, Flat argument heads, Verbatim, Optional,
 *     a trailing sequence.
 *
 *  2. A call pattern on another head never takes a GROUP of a Flat
 *     operand sequence, OneIdentity or not (match.c): MatchQ[2 x y, y_^n_. x]
 *     is False.  Offering group spans made n_?NumberQ o1_. + o2_:0 exponential
 *     in the length of a sum (FeynCalc's DotSimplify).
 *
 *  3. FreeQ with an atom / _h / Alternatives form walks the tree comparing
 *     pointers instead of calling the matcher at every node (funcprog.c).
 */

#include "test_utils.h"
#include "core.h"
#include "expr.h"
#include "eval.h"
#include "parse.h"
#include "symtab.h"

#include <stdio.h>
#include <string.h>
#include <time.h>

static void run(const char* src) {
    Expr* p = parse_expression(src);
    Expr* r = evaluate(p);
    expr_free(p);
    expr_free(r);
}

/* ---- 1. slot keys ------------------------------------------------- */

static void test_orderless_assignment(void) {
    run("SetAttributes[po, Orderless]; po[a_g, b_h] := {a, b}");
    assert_eval_eq("po[h[1], g[2]]", "{g[2], h[1]}", 0);
    assert_eval_eq("po[g[2], k[1]]", "po[g[2], k[1]]", 0);
}

static void test_flat_argument_head_arity(void) {
    /* The slot's arity is a lower bound for a Flat argument head. */
    run("SetAttributes[fl, Flat]; rr[fl[x_, y_]] := {x, y}");
    assert_eval_eq("Length[rr[fl[1, 2, 3]]]", "2", 0);
    assert_eval_eq("rr[fl[1]]", "rr[fl[1]]", 0);
}

static void test_head_alternatives_and_arity(void) {
    run("q[(a : g | h)[x_]] := {a, x}; q2[g[x_, y_]] := {x, y}");
    assert_eval_eq("q[h[1]]", "{h, 1}", 0);
    assert_eval_eq("q[k[1]]", "q[k[1]]", 0);
    assert_eval_eq("q[g[1, 2]]", "q[g[1, 2]]", 0);
    assert_eval_eq("q2[g[1]]", "q2[g[1]]", 0);
    assert_eval_eq("q2[g[1, 2]]", "{1, 2}", 0);
}

static void test_verbatim_optional_sequence(void) {
    run("v[Verbatim[_Integer]] := 1; o[x_, y_ : 0] := {x, y}; m[g[x_], r___] := {x, r}");
    assert_eval_eq("v[_Integer]", "1", 0);
    assert_eval_eq("v[5]", "v[5]", 0);
    assert_eval_eq("o[1]", "{1, 0}", 0);
    assert_eval_eq("m[g[1], 2, 3]", "{1, 2, 3}", 0);
    assert_eval_eq("m[k[1], 2]", "m[k[1], 2]", 0);
}

static void test_orderless_times_slot(void) {
    /* Pair[n_Integer x_, y_] fires through either argument of an Orderless
     * head that was Orderless when the rule was made. */
    run("SetAttributes[pq, Orderless]; pq[n_Integer x_, y_] := {n, x, y}");
    assert_eval_eq("pq[s, 2 t]", "{2, t, s}", 0);
    assert_eval_eq("pq[2 t, s]", "{2, t, s}", 0);
    assert_eval_eq("pq[s, t]", "pq[s, t]", 0);
}

/* ---- 2. no group spans for a call on another head ---------------- */

static void test_no_group_for_foreign_call(void) {
    assert_eval_eq("MatchQ[2 x y, y_^n_. x]", "False", 0);
    assert_eval_eq("ReplaceList[a + b + c, n_ o1_. + o2_ :> {n, o1, o2}]",
                   "{{a, 1, b + c}, {b, 1, a + c}, {c, 1, a + b}}", 0);
    assert_eval_eq("ReplaceList[a b c d, y_^n_. z_ :> {y, n, z}]",
                   "{{a, 1, b c d}, {b, 1, a c d}, {c, 1, a b d}, {d, 1, a b c}}", 0);
    assert_eval_eq("ReplaceList[Hold[a + b + c + d], Hold[k_. y_ + z_] :> {k, y, z}]",
                   "{{1, a, b + c + d}, {1, b, a + c + d}, {1, c, a + b + d}, {1, d, a + b + c}}", 0);
    /* The whole subject is still reachable through the defaults (f[f[...]]
     * with every Optional defaulted), the one group the reference allows. */
    assert_eval_eq("Sort[ReplaceList[2 x, x_. y_^n_. :> {x, y, n}]]",
                   "{{1, 2 x, 1}, {2, x, 1}, {x, 2, 1}}", 0);
    assert_eval_eq("ReplaceList[a + b + c, n_ o1_. + o2_:0 :> {n, o1, o2}]",
                   "{{a, 1, b + c}, {b, 1, a + c}, {c, 1, a + b}, {a + b + c, 1, 0}}", 0);
}

static void test_numeric_coefficient_pattern_is_linear(void) {
    /* 40 terms: was ~3^n (20 s at 14 terms); now microseconds. */
    clock_t t0 = clock();
    assert_eval_eq("MatchQ[Sum[c[i] f[i], {i, 40}] + 3 g, n_?NumberQ o1_. + o2_:0]", "True", 0);
    assert_eval_eq("MatchQ[Sum[c[i] f[i], {i, 40}], n_?NumberQ o1_. + o2_:0]", "False", 0);
    assert_eval_eq("MatchQ[Sum[c[i] f[i], {i, 40}], n_?NumberQ o1_ + o2_]", "False", 0);
    double secs = (double)(clock() - t0) / CLOCKS_PER_SEC;
    ASSERT_MSG(secs < 2.0, "numeric-coefficient Plus pattern took %.2f s", secs);
}

/* ---- UpValue dispatch: outer head and slot keys ------------------- */

static void test_upvalue_outer_head(void) {
    /* An UpValue pattern's own head must be the enclosing call's; a pattern
     * construct at the top (u | v)[...] or f_[...] keys nothing. */
    run("g /: Set[g[a_], c_] := (gs[a] = c); g[1] = 5;"
        "h /: (u | v)[h[x_]] := {x};"
        "SetAttributes[oo, Orderless]; m /: oo[m[x_], n_Integer] := x + n;"
        "q /: f_[q[x_]] := {f, x};"
        "z /: z[x_] + z[y_] := z[x + y]");
    assert_eval_eq("gs[1]", "5", 0);
    assert_eval_eq("{u[h[2]], v[h[3]], w[h[4]]}", "{{2}, {3}, w[h[4]]}", 0);
    assert_eval_eq("oo[7, m[1]]", "8", 0);
    assert_eval_eq("r[q[1]]", "{r, 1}", 0);
    assert_eval_eq("z[1] + z[2] + a", "a + z[3]", 0);
    assert_eval_eq("Times[2, g[3]]", "2 g[3]", 0);
}

/* ---- OptionValue injection shares untouched subtrees --------------- */

static void test_optionvalue_injection(void) {
    /* Only the path to each OptionValue is rebuilt (options_builtin.c); the
     * bound arguments are shared, so a big argument keeps its fixed-point
     * stamps.  Resolution must still reach every depth, the head, and a guard. */
    run("Options[ov] = {a -> 1, b -> 2};"
        "ov[x_, OptionsPattern[]] := {x, OptionValue[a], g[h[OptionValue[b]]], OptionValue[ov, b]};"
        "ow[x_, OptionsPattern[ov]] := OptionValue[a][x];"
        "oc[x_, OptionsPattern[ov]] := x /; OptionValue[a] === 3");
    assert_eval_eq("ov[k[1, 2], a -> 5]", "{k[1, 2], 5, g[h[2]], 2}", 0);
    assert_eval_eq("ov[k[1, 2]]", "{k[1, 2], 1, g[h[2]], 2}", 0);
    assert_eval_eq("ow[q, a -> s]", "s[q]", 0);
    assert_eval_eq("{oc[1, a -> 3], oc[1]}", "{1, oc[1]}", 0);
    clock_t t0 = clock();
    run("bigx = Sum[c[i] f[i], {i, 20000}]; Do[ov[bigx], {20}]");
    double secs = (double)(clock() - t0) / CLOCKS_PER_SEC;
    ASSERT_MSG(secs < 3.0, "OptionsPattern call on a 20000-term argument took %.2f s", secs);
}

/* ---- 3. FreeQ fast path ------------------------------------------ */

static void test_freeq_simple_forms(void) {
    assert_eval_eq("FreeQ[f[g[x]], x]", "False", 0);
    assert_eval_eq("FreeQ[f[1/2], Rational]", "False", 0);
    assert_eval_eq("FreeQ[f[1/2], 2]", "True", 0);
    assert_eval_eq("FreeQ[<|a -> 1|>, a]", "True", 0);
    assert_eval_eq("FreeQ[f[x], f, Heads -> False]", "True", 0);
    assert_eval_eq("FreeQ[f[x], f]", "False", 0);
    assert_eval_eq("FreeQ[{1, 2}, _List]", "False", 0);
    assert_eval_eq("FreeQ[a + b, a | c]", "False", 0);
    assert_eval_eq("FreeQ[a + b, c | d]", "True", 0);
    assert_eval_eq("FreeQ[a + b, _Plus]", "False", 0);
    assert_eval_eq("FreeQ[a + b, _Times]", "True", 0);
    assert_eval_eq("FreeQ[f[x], _f]", "False", 0);
    assert_eval_eq("FreeQ[f[I], Complex]", "False", 0);
    assert_eval_eq("FreeQ[f[3], 3]", "False", 0);
    assert_eval_eq("FreeQ[f[\"s\"], \"s\"]", "False", 0);
    assert_eval_eq("FreeQ[f[x_], x]", "False", 0);
    assert_eval_eq("FreeQ[f[2^100], 2^100]", "False", 0);
    assert_eval_eq("FreeQ[f[2.5], 2.5]", "False", 0);
    assert_eval_eq("FreeQ[Hold[1 + 1], 2]", "True", 0);
    assert_eval_eq("FreeQ[f[x, g[y]], y, {1}]", "True", 0);
    assert_eval_eq("FreeQ[f[x, g[y]], y, {2}]", "False", 0);
}

int main(void) {
    symtab_init();
    core_init();

    TEST(test_orderless_assignment);
    TEST(test_flat_argument_head_arity);
    TEST(test_head_alternatives_and_arity);
    TEST(test_verbatim_optional_sequence);
    TEST(test_orderless_times_slot);
    TEST(test_no_group_for_foreign_call);
    TEST(test_numeric_coefficient_pattern_is_linear);
    TEST(test_upvalue_outer_head);
    TEST(test_optionvalue_injection);
    TEST(test_freeq_simple_forms);

    printf("All dispatch slot-key / Flat-span / FreeQ tests passed.\n");
    return 0;
}
