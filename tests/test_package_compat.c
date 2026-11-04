/*
 * Legacy-package compatibility (issue #125): the kernel features the
 * DifferentialForms package (Zizza/Jentschura, 2007) relies on, reduced to
 * plain-kernel repros so they are guarded without the package itself (which
 * is not vendored; tools/check_differentialforms.py runs the full corpus).
 *
 *   - Literal, the obsolete System` synonym of HoldPattern;
 *   - the positional ReplacePart[expr, new, pos] / [expr, new, pos, npos];
 *   - Dt of a symbol with the Constant attribute;
 *   - Attributes[{s1, s2, ...}];
 *   - named operator heads (\[Wedge] -> Wedge) being System` symbols, so a
 *     package's rules for them are reached from the infix form;
 *   - Context["name"] resolved through $ContextPath;
 *   - display: SequenceForm, Format[..., InputForm | TeXForm], HoldForm in
 *     InputForm, TeX strings.
 */

#include "test_utils.h"
#include "expr.h"
#include "eval.h"
#include "parse.h"
#include "core.h"
#include "symtab.h"
#include <stdio.h>
#include <string.h>

/* Evaluate and discard (a definition / setup statement). */
static void run(const char* src) {
    Expr* e = parse_expression(src);
    ASSERT(e != NULL);
    Expr* r = evaluate(e);
    expr_free(e);
    if (r) expr_free(r);
}

/* ---- Literal ---- */

static void test_literal_is_holdpattern_for_matching(void) {
    /* The package's ExteriorDerivative: Dt[f] /. Literal[Dt[x_Symbol]] -> ... */
    assert_eval_eq("Dt[x^2 y] /. Literal[Dt[v_Symbol]] -> dd[v]", "2 x y dd[x] + x^2 dd[y]", 0);
    assert_eval_eq("f[2] /. Literal[f[a_]] -> a", "2", 0);
    assert_eval_eq("MatchQ[f[1], Literal[f[_]]]", "True", 0);
    assert_eval_eq("Cases[{f[1], g[2]}, Literal[f[_]]]", "{f[1]}", 0);
    assert_eval_eq("Hold[f[1]] /. Literal[f[a_]] :> a", "Hold[1]", 0);
    /* Held: Literal[1 + 1] stays as written. */
    assert_eval_eq("FullForm[Literal[1 + 1]]", "Literal[Plus[1, 1]]", 0);
}

static void test_literal_assignment_and_attributes(void) {
    run("Literal[litk[a_]] := a^2");
    assert_eval_eq("litk[3]", "9", 0);
    assert_eval_eq("Attributes[Literal]", "{HoldAll, Protected}", 0);
    assert_eval_eq("Context[Literal]", "\"System`\"", 0);
    /* A distinct symbol: ReleaseHold leaves it alone. */
    assert_eval_eq("ReleaseHold[Literal[1 + 1]]", "Literal[1 + 1]", 0);
}

/* ---- positional ReplacePart ---- */

static void test_replacepart_positional(void) {
    assert_eval_eq("ReplacePart[{a, b, c}, x, {1}]", "{x, b, c}", 0);
    assert_eval_eq("ReplacePart[{a, b, c}, x, 2]", "{a, x, c}", 0);
    assert_eval_eq("ReplacePart[{a, b, c}, x, -1]", "{a, b, x}", 0);
    assert_eval_eq("ReplacePart[{a, b, c}, x, {{1}, {3}}]", "{x, b, x}", 0);
    assert_eval_eq("ReplacePart[{{a, b}, {c, d}}, x, {2, 1}]", "{{a, b}, {x, d}}", 0);
    assert_eval_eq("ReplacePart[{{a, b}, {c, d}}, x, {{1, 2}, {2, 1}}]", "{{a, x}, {x, d}}", 0);
    assert_eval_eq("ReplacePart[f[a, b, c], x, 0]", "x[a, b, c]", 0);
    /* {} names no position; {{}} names the whole expression. */
    assert_eval_eq("ReplacePart[{a, b, c}, x, {}]", "{a, b, c}", 0);
    assert_eval_eq("ReplacePart[{a, b, c}, x, {{}}]", "x", 0);
    /* Parts of new at npos. */
    assert_eval_eq("ReplacePart[{a, b, c}, {u, v}, {{1}, {2}}, {{2}, {1}}]", "{v, u, c}", 0);
    /* The DifferentialForms Indices idiom. */
    assert_eval_eq("ReplacePart[1 + RotateRight[{i, j, k}], 1, {1}]", "{1, 1 + i, 1 + j}", 0);
}

static void test_replacepart_positional_errors(void) {
    test_mute_messages();
    /* A missing part is an error (::partw) and the call stays unevaluated;
     * a malformed spec is ::psl. */
    assert_eval_eq("ReplacePart[{a, b, c}, x, {5}]", "ReplacePart[{a, b, c}, x, {5}]", 0);
    assert_eval_eq("ReplacePart[{a, b, c}, x, y]", "ReplacePart[{a, b, c}, x, y]", 0);
    assert_eval_eq("Check[ReplacePart[{a, b, c}, x, {5}], failed]", "failed", 0);
}

/* ---- Dt and the Constant attribute ---- */

static void test_dt_constant_attribute(void) {
    run("SetAttributes[kc, Constant]");
    assert_eval_eq("Dt[kc]", "0", 0);
    assert_eval_eq("Dt[kc x]", "kc Dt[x]", 0);
    assert_eval_eq("Dt[kc x, x]", "kc", 0);
    assert_eval_eq("Dt[kc^2 + y]", "Dt[y]", 0);
    assert_eval_eq("Dt[kc, x]", "0", 0);
    /* Without the attribute a symbol is a variable. */
    assert_eval_eq("Dt[kn x]", "x Dt[kn] + kn Dt[x]", 0);
}

/* ---- Attributes over a list ---- */

static void test_attributes_list(void) {
    run("SetAttributes[kc2, Constant]");
    assert_eval_eq("Attributes[{kc2, plainsym}]", "{{Constant}, {}}", 0);
    assert_eval_eq("FreeQ[Attributes[{kc2, plainsym}], Constant]", "False", 0);
    assert_eval_eq("Attributes[{}]", "{}", 0);
    test_mute_messages();
    assert_eval_eq("Attributes[{1}]", "{Attributes[1]}", 0);
}

/* ---- operator heads are System`, Context["name"] ---- */

static void test_operator_heads_are_system(void) {
    run("BeginPackage[\"WedgePkg`\"]");
    run("Wedge::usage = \"w\"");
    run("CircleTimes::usage = \"c\"");
    run("Begin[\"`Private`\"]");
    run("Wedge[a_, b_] := wedged[a, b]");
    run("End[]");
    run("EndPackage[]");
    /* The package's Wedge IS System`Wedge, so the infix form reaches it. */
    assert_eval_eq("Context[Wedge]", "\"System`\"", 0);
    assert_eval_eq("Context[CircleTimes]", "\"System`\"", 0);
    assert_eval_eq("p \\[Wedge] q", "WedgePkg`Private`wedged[p, q]", 0);
}

static void test_context_of_string_name(void) {
    run("BeginPackage[\"StrPkg`\"]");
    run("strPkgFn::usage = \"u\"");
    run("Begin[\"`Private`\"]");
    run("strPkgFn[x_] := x");
    run("End[]");
    run("EndPackage[]");
    assert_eval_eq("Context[\"strPkgFn\"]", "\"StrPkg`\"", 0);
    assert_eval_eq("Context[\"StrPkg`strPkgFn\"]", "\"StrPkg`\"", 0);
    assert_eval_eq("Context[\"Plus\"]", "\"System`\"", 0);
}

/* ---- display ---- */

static void test_sequenceform_and_format_values(void) {
    run("Format[fq[a_]] := SequenceForm[\"(\", a, \") \", \"dq\"]");
    run("Format[fq[a_], InputForm] := HoldForm[a qq]");
    run("Format[fq[a_], TeXForm] := \"TEX\"");
    assert_eval_eq("ToString[fq[1]]", "\"(1) dq\"", 0);
    assert_eval_eq("ToString[{fq[2], fq[3] + 1}]", "\"{(2) dq, 1 + (3) dq}\"", 0);
    assert_eval_eq("ToString[fq[1], InputForm]", "\"HoldForm[1*qq]\"", 0);
    assert_eval_eq("ToString[TeXForm[fq[1]]]", "\"\\\\text{TEX}\"", 0);
    assert_eval_eq("ToString[SequenceForm[\"a\", 1, \" b\"]]", "\"a1 b\"", 0);
    assert_eval_eq("ToString[SequenceForm[\"a\", 1], InputForm]", "\"\\\"a\\\"1\"", 0);
}

static void test_holdform_inputform_and_tex_strings(void) {
    assert_eval_eq("ToString[HoldForm[1 + 1], InputForm]", "\"HoldForm[1 + 1]\"", 0);
    assert_eval_eq("ToString[HoldForm[1*x], InputForm]", "\"HoldForm[1*x]\"", 0);
    assert_eval_eq("ToString[{HoldForm[2*3]}, InputForm]", "\"{HoldForm[2*3]}\"", 0);
    assert_eval_eq("ToString[HoldForm[1*x]]", "\"1 x\"", 0);
    assert_eval_eq("ToString[TeXForm[\"abc\"]]", "\"\\\\text{abc}\"", 0);
}

int main(void) {
    symtab_init();
    core_init();

    TEST(test_literal_is_holdpattern_for_matching);
    TEST(test_literal_assignment_and_attributes);
    TEST(test_replacepart_positional);
    TEST(test_replacepart_positional_errors);
    TEST(test_dt_constant_attribute);
    TEST(test_attributes_list);
    TEST(test_operator_heads_are_system);
    TEST(test_context_of_string_name);
    TEST(test_sequenceform_and_format_values);
    TEST(test_holdform_inputform_and_tex_strings);

    printf("All package-compat tests passed.\n");
    return 0;
}
