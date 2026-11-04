/*
 * Tests for the Context system: Context[], BeginPackage[], Begin[], End[],
 * EndPackage[], and the $Context / $ContextPath system variables.
 */

#include "test_utils.h"
#include "expr.h"
#include "eval.h"
#include "parse.h"
#include "core.h"
#include "symtab.h"
#include "context.h"
#include <stdio.h>
#include <string.h>

/* ---- Direct API (non-parser) tests ---- */

static void test_initial_state(void) {
    /* core_init() loads the FullSimplify package (BeginPackage["FullSimplify`"]
     * ... EndPackage[]) and then the SeparablePolynomialDecomposition package.
     * Each EndPackage prepends its package context to $ContextPath, so the
     * exported symbols (FullSimplify, RegisterTransforms,
     * SeparablePolynomialDecomposition) resolve under their short names. The
     * later load comes first. The current context returns to Global` and
     * System` stays last. */
    ASSERT_STR_EQ(context_current(), "Global`");
    ASSERT(context_path_count() == 5);
    ASSERT_STR_EQ(context_path_at(0), "SeparablePolynomialDecomposition`");
    ASSERT_STR_EQ(context_path_at(1), "FullSimplify`");
    ASSERT_STR_EQ(context_path_at(2), "Global`");
    ASSERT_STR_EQ(context_path_at(3), "Mathilda`");   /* Mathilda-only builtins */
    ASSERT_STR_EQ(context_path_at(4), "System`");
}

static void test_resolve_bare_builtin(void) {
    /* A builtin like "Plus" is stored bare; resolve should keep it bare. */
    char* r = context_resolve_name("Plus");
    ASSERT_STR_EQ(r, "Plus");
    free(r);
}

static void test_resolve_system_prefix_stripped(void) {
    /* "System`Plus" should canonicalize to bare "Plus". */
    char* r = context_resolve_name("System`Plus");
    ASSERT_STR_EQ(r, "Plus");
    free(r);
}

static void test_resolve_dollar_unqualified(void) {
    /* In Global`, a new $Foo stays bare like any other new name. (Inside a
     * package it is qualified -- see test_package_dollar_symbol.) */
    char* r = context_resolve_name("$Foo");
    ASSERT_STR_EQ(r, "$Foo");
    free(r);
}

static void test_resolve_absolute(void) {
    char* r = context_resolve_name("MyPkg`f");
    ASSERT_STR_EQ(r, "MyPkg`f");
    free(r);
}

static void test_begin_end_roundtrip(void) {
    ASSERT(context_begin("MyCtx`"));
    ASSERT_STR_EQ(context_current(), "MyCtx`");
    char* closed = NULL;
    ASSERT(context_end(&closed));
    ASSERT_STR_EQ(closed, "MyCtx`");
    free(closed);
    ASSERT_STR_EQ(context_current(), "Global`");
}

static void test_begin_relative(void) {
    ASSERT(context_begin("Outer`"));
    ASSERT_STR_EQ(context_current(), "Outer`");
    ASSERT(context_begin("`Inner`"));        /* relative: Outer`Inner` */
    ASSERT_STR_EQ(context_current(), "Outer`Inner`");
    char* closed = NULL;
    ASSERT(context_end(&closed)); free(closed);
    ASSERT_STR_EQ(context_current(), "Outer`");
    ASSERT(context_end(&closed)); free(closed);
    ASSERT_STR_EQ(context_current(), "Global`");
}

static void test_begin_package_path(void) {
    ASSERT(context_begin_package("Pkg`", NULL, 0));
    ASSERT_STR_EQ(context_current(), "Pkg`");
    ASSERT(context_path_count() == 2);
    ASSERT_STR_EQ(context_path_at(0), "Pkg`");
    ASSERT_STR_EQ(context_path_at(1), "System`");
    char* closed = NULL;
    ASSERT(context_end_package(&closed));
    free(closed);
    ASSERT_STR_EQ(context_current(), "Global`");
    /* Pkg` should have been prepended to the restored path. */
    ASSERT_STR_EQ(context_path_at(0), "Pkg`");
    /* Clean up so subsequent tests start fresh(ish): pop the prepended entry
     * by running another BeginPackage/EndPackage cycle — not ideal. Instead,
     * we accept the residual Pkg` and continue; other tests are tolerant. */
}

static void test_begin_package_with_needs(void) {
    const char* needs[] = { "Needed1`", "Needed2`" };
    ASSERT(context_begin_package("NPkg`", needs, 2));
    ASSERT_STR_EQ(context_current(), "NPkg`");
    /* Expected path: {NPkg`, System`, Needed1`, Needed2`}. */
    ASSERT(context_path_count() == 4);
    ASSERT_STR_EQ(context_path_at(0), "NPkg`");
    ASSERT_STR_EQ(context_path_at(1), "System`");
    ASSERT_STR_EQ(context_path_at(2), "Needed1`");
    ASSERT_STR_EQ(context_path_at(3), "Needed2`");
    char* closed = NULL;
    ASSERT(context_end_package(&closed));
    free(closed);
}

/* ---- End-to-end (parser + evaluator) tests ---- */

static void test_context_builtin_no_args(void) {
    assert_eval_eq("Context[]", "\"Global`\"", 0);
}

static void test_context_builtin_of_builtin(void) {
    assert_eval_eq("Context[Plus]", "\"System`\"", 0);
    assert_eval_eq("Context[Sin]", "\"System`\"", 0);
}

static void test_dollar_context_and_path(void) {
    assert_eval_eq("$Context", "\"Global`\"", 0);
    /* The FullSimplify and SeparablePolynomialDecomposition packages, loaded
     * by core_init, leave their contexts on $ContextPath (see
     * test_initial_state). */
    assert_eval_eq("$ContextPath",
                   "{\"SeparablePolynomialDecomposition`\", \"FullSimplify`\", \"Global`\", \"Mathilda`\", \"System`\"}", 0);
}

static void test_begin_changes_context(void) {
    assert_eval_eq("Begin[\"Cell$$001`\"]", "\"Cell$$001`\"", 0);
    assert_eval_eq("$Context", "\"Cell$$001`\"", 0);
    assert_eval_eq("End[]", "\"Cell$$001`\"", 0);
    assert_eval_eq("$Context", "\"Global`\"", 0);
}

static void test_package_symbol_qualification(void) {
    /* Inside BeginPackage, a new symbol is qualified into the package. */
    assert_eval_eq("BeginPackage[\"MyPkgA`\"]", "\"MyPkgA`\"", 0);
    /* After BeginPackage the path is {MyPkgA`, System`} but not Global`. */
    assert_eval_eq("$ContextPath", "{\"MyPkgA`\", \"System`\"}", 0);
    /* Define a package function. */
    assert_eval_eq("f[x_] := x^2 + 1", "Null", 0);
    assert_eval_eq("EndPackage[]", "Null", 0);
    /* Inside the package, "f" resolved to MyPkgA`f. After EndPackage, MyPkgA`
     * is on the path, so the short name "f" still finds MyPkgA`f. */
    assert_eval_eq("f[3]", "10", 0);
    /* And we can always call it by its fully qualified name. */
    assert_eval_eq("MyPkgA`f[5]", "26", 0);
}

static void test_nested_begin_private(void) {
    /* Exercises Begin["`Private`"] within a package. Note: symbol references
     * are resolved at parse time, so a public stub (priv) needs to exist in
     * the public context before code that uses it is parsed; otherwise the
     * later Private-context definition creates a distinct symbol. We define
     * priv in the public context first, then re-enter Private. */
    assert_eval_eq("BeginPackage[\"UtilPkg`\"]", "\"UtilPkg`\"", 0);
    /* Seed priv in the public context so pub sees it. */
    assert_eval_eq("priv[x_] := 10 x", "Null", 0);
    assert_eval_eq("pub[x_] := priv[x] + 1", "Null", 0);
    /* Now drop into the Private context: verify the context string updates. */
    assert_eval_eq("Begin[\"`Private`\"]", "\"UtilPkg`Private`\"", 0);
    assert_eval_eq("$Context", "\"UtilPkg`Private`\"", 0);
    /* A helper defined here gets the Private context. */
    assert_eval_eq("helper[x_] := x + x", "Null", 0);
    assert_eval_eq("End[]", "\"UtilPkg`Private`\"", 0);
    assert_eval_eq("EndPackage[]", "Null", 0);
    /* Short-name lookup finds the public symbols. */
    assert_eval_eq("pub[2]", "21", 0);
    /* The helper in Private is reachable only by its full name. */
    assert_eval_eq("UtilPkg`Private`helper[4]", "8", 0);
}

static void test_system_prefix_stripped_in_parser(void) {
    /* Typing "System`Plus" should resolve to the bare builtin. */
    assert_eval_eq("System`Plus[2, 3]", "5", 0);
}

static void test_context_of_user_symbol(void) {
    assert_eval_eq("myvar = 5", "5", 0);
    assert_eval_eq("Context[myvar]", "\"Global`\"", 0);
    assert_eval_eq("Context[\"myvar\"]", "\"Global`\"", 0);
}

static void test_package_dollar_symbol(void) {
    /* Issue #104: a leading '$' is not special -- $T2 inside a package is
     * Test104`$T2, exactly like T1. */
    assert_eval_eq("BeginPackage[\"Test104`\"]", "\"Test104`\"", 0);
    assert_eval_eq("T1 = 137", "137", 0);
    assert_eval_eq("$T2 = 42", "42", 0);
    assert_eval_eq("Context[T1]", "\"Test104`\"", 0);
    assert_eval_eq("Context[$T2]", "\"Test104`\"", 0);
    /* System $-variables still resolve to System inside the package. */
    assert_eval_eq("Context[$Failed]", "\"System`\"", 0);
    assert_eval_eq("Context[$Assumptions]", "\"System`\"", 0);
    assert_eval_eq("$RecursionLimit", "1024", 0);
    assert_eval_eq("EndPackage[]", "Null", 0);
    /* After EndPackage the short name finds the package symbol. */
    assert_eval_eq("$T2", "42", 0);
    assert_eval_eq("Test104`$T2", "42", 0);
    assert_eval_eq("Context[$T2]", "\"Test104`\"", 0);
    /* Begin[] qualifies $-names the same way. */
    assert_eval_eq("Begin[\"Foo104`\"]", "\"Foo104`\"", 0);
    assert_eval_eq("$y104 = 1", "1", 0);
    assert_eval_eq("Context[$y104]", "\"Foo104`\"", 0);
    assert_eval_eq("End[]", "\"Foo104`\"", 0);
}

static void test_context_of_system_symbols(void) {
    /* Non-function System symbols (constants, $-variables, structural heads)
     * report System`, not Global`. Scratch names the kernel interns (x, t)
     * remain Global`. */
    /* Context is HoldFirst; call it directly (Context /@ {...} would
     * evaluate $Version to its string first). */
    assert_eval_eq("{Context[$Version], Context[$Failed], Context[Pi], Context[Automatic], Context[List]}",
                   "{\"System`\", \"System`\", \"System`\", \"System`\", \"System`\"}", 0);
    assert_eval_eq("{Context[x], Context[t]}", "{\"Global`\", \"Global`\"}", 0);
    assert_eval_eq("$g104 = 1; Context[$g104]", "\"Global`\"", 0);
}

/* ---- inert front-end / box names are System`, even inside a package ---- */
static void test_inert_frontend_names_are_system(void) {
    /* A package that mentions these must reference the System` symbol rather
     * than mint a private copy in its own context. */
    /* Separate statements: symbols are resolved when a statement is parsed. */
    assert_eval_eq("BeginPackage[\"ZZfe`\"]; Begin[\"`Private`\"];", "Null", 0);
    assert_eval_eq("zzfe = {Context[Bold], Context[RowBox], Context[$FrontEndSession], "
                   "Context[SubscriptBox], Context[TableForm], Context[Integral], Context[zzfeLocal]};",
                   "Null", 0);
    assert_eval_eq("End[]; EndPackage[];", "Null", 0);
    assert_eval_eq("ZZfe`Private`zzfe",
                   "{\"System`\", \"System`\", \"System`\", \"System`\", \"System`\", "
                   "\"System`\", \"ZZfe`Private`\"}", 0);
    /* They stay inert. */
    assert_eval_eq("{Bold, RowBox[{\"a\", \"b\"}], Darker[x]}", "{Bold, RowBox[{\"a\", \"b\"}], Darker[x]}", 0);
}

/* ---- $ContextPath / $Context are assignable (they ARE the state) ---- */
static void test_contextpath_assignment(void) {
    size_t n0 = context_path_count();
    assert_eval_eq("CpFoo`bar = 1; Names[\"bar\"]", "{}", 0);
    assert_eval_eq("$ContextPath = Append[$ContextPath, \"CpFoo`\"]; "
                   "{Names[\"bar\"], ToString[Hold[CpFoo`bar]]}",
                   "{{\"bar\"}, \"Hold[bar]\"}", 0);
    ASSERT(context_path_count() == n0 + 1);
    ASSERT_STR_EQ(context_path_at(n0), "CpFoo`");
    /* A bare name now resolves through the new path entry. */
    assert_eval_eq("bar", "1", 0);
    /* Invalid values: message, old value kept, Set still returns the rhs. */
    assert_eval_eq("Quiet[$ContextPath = 5]", "5", 0);
    assert_eval_eq("Check[$ContextPath = {\"A`\", 3}, \"msg\"]", "\"msg\"", 0);
    assert_eval_eq("Last[$ContextPath]", "\"CpFoo`\"", 0);
    /* AppendTo / DeleteCases round trip. */
    assert_eval_eq("AppendTo[$ContextPath, \"CpZed`\"]; Last[$ContextPath]", "\"CpZed`\"", 0);
    assert_eval_eq("$ContextPath = DeleteCases[$ContextPath, \"CpZed`\" | \"CpFoo`\"]; "
                   "{Names[\"bar\"], ToString[Hold[CpFoo`bar]]}",
                   "{{}, \"Hold[CpFoo`bar]\"}", 0);
    ASSERT(context_path_count() == n0);
    /* Block scopes it, and restores it on exit. */
    assert_eval_eq("Block[{$ContextPath = {\"System`\", \"CpFoo`\"}}, Names[\"bar\"]]",
                   "{\"bar\"}", 0);
    assert_eval_eq("Names[\"bar\"]", "{}", 0);
    assert_eval_eq("Attributes[$ContextPath]", "{}", 0);
}

static void test_context_assignment(void) {
    assert_eval_eq("$Context = \"CpQ`\"", "\"CpQ`\"", 0);
    ASSERT_STR_EQ(context_current(), "CpQ`");
    /* New symbols are created in the new context. */
    assert_eval_eq("cpNew = 1; Context[cpNew]", "\"CpQ`\"", 0);
    assert_eval_eq("$Context = \"Global`\"; Names[\"CpQ`*\"]", "{\"CpQ`cpNew\"}", 0);
    assert_eval_eq("Quiet[$Context = \"nobacktick\"]; $Context", "\"Global`\"", 0);
    assert_eval_eq("Check[$Context = 3, \"msg\"]", "\"msg\"", 0);
    assert_eval_eq("Attributes[$Context]", "{}", 0);
    /* Begin/End still nest on top of a directly assigned value. */
    assert_eval_eq("Begin[\"CpB`\"]; c1 = $Context; End[]; {c1, $Context}",
                   "{\"CpB`\", \"Global`\"}", 0);
}

/* ---- L7: the Mathilda` context, and package resolution ----
 * Mathilda-only kernel names live in Mathilda` (on the default $ContextPath);
 * genuine System names stay System`.  Inside BeginPackage -- whose path is
 * {ctx`, System`} -- neither Mathilda-only builtins nor Global` symbols are
 * visible, so a package gets its own symbols of those names. */
static void test_mathilda_context(void) {
    assert_eval_eq("{Context[DataType], Context[NDArray], Context[CompileDiagnostics], Context[$SimplifyDebug]}",
                   "{\"Mathilda`\", \"Mathilda`\", \"Mathilda`\", \"Mathilda`\"}", 0);
    assert_eval_eq("{Context[Plus], Context[List], Context[Sin], Context[$Version]}",
                   "{\"System`\", \"System`\", \"System`\", \"System`\"}", 0);
    assert_eval_eq("MemberQ[$ContextPath, \"Mathilda`\"] && MemberQ[$Packages, \"Mathilda`\"]", "True", 0);
    /* Qualified spellings reach the same builtin. */
    assert_eval_eq("{Mathilda`DataType === DataType, System`Plus === Plus}", "{True, True}", 0);
    assert_eval_eq("Head[Mathilda`ToNDArray[{1., 2.}]] === Head[ToNDArray[{1., 2.}]]", "True", 0);
    /* Needs["Mathilda`"] is satisfied already: no load, no message. */
    assert_eval_eq("Check[Needs[\"Mathilda`\"], \"msg\"]", "Null", 0);
}

static void test_package_does_not_see_global_or_mathilda(void) {
    assert_eval_eq("l7foo = 1;", "Null", 0);
    assert_eval_eq("BeginPackage[\"L7P`\"];", "Null", 0);
    /* A valued Global` symbol, an unvalued Global` symbol the kernel itself made
     * (x), and a Mathilda-only builtin: all private to the package now. */
    assert_eval_eq("l7ctx = {Context[l7foo], Context[x], Context[DataType], Context[Plus]};", "Null", 0);
    /* Defining the package's own DataType: no Set::wrsym, no Protected clash. */
    assert_eval_eq("Check[DataType[q_] := {\"pkg\", q}; DataType::usage = \"pkg\"; \"ok\", \"msg\"]",
                   "\"ok\"", 0);
    assert_eval_eq("EndPackage[];", "Null", 0);
    assert_eval_eq("L7P`l7ctx", "{\"L7P`\", \"L7P`\", \"L7P`\", \"System`\"}", 0);
    /* After EndPackage the package's symbol comes first; the kernel's is still
     * reachable qualified, and the Global` one is untouched. */
    assert_eval_eq("{DataType[1], Context[DataType], l7foo}", "{{\"pkg\", 1}, \"L7P`\", 1}", 0);
    assert_eval_eq("Head[Mathilda`NDArray[{1, 2}]]", "NDArray", 0);
}

static void test_package_can_request_mathilda(void) {
    /* A package that wants the kernel's own builtins asks for Mathilda`. */
    assert_eval_eq("BeginPackage[\"L7Q`\", {\"Mathilda`\"}];", "Null", 0);
    assert_eval_eq("l7q = {Context[ToNDArray], Context[CompileDiagnostics]};", "Null", 0);
    assert_eval_eq("EndPackage[];", "Null", 0);
    assert_eval_eq("L7Q`l7q", "{\"Mathilda`\", \"Mathilda`\"}", 0);
    /* Begin without BeginPackage keeps the path: Global` and Mathilda` stay visible. */
    assert_eval_eq("l7g = 5; Begin[\"L7R`\"];", "Null", 0);
    assert_eval_eq("l7r = {Context[l7g], Context[DataType]};", "Null", 0);
    assert_eval_eq("End[]", "\"L7R`\"", 0);
    assert_eval_eq("L7R`l7r", "{\"Global`\", \"L7P`\"}", 0);
    /* The kernel's own FullSimplify package declared its Mathilda` dependency:
     * no private copies of the kernel builtins it uses. */
    assert_eval_eq("Select[Names[\"FullSimplify`Private`*\"], "
                   "MemberQ[{\"SimplifyCount\", \"LoadModule\"}, Last[StringSplit[#, \"`\"]]] &]", "{}", 0);
    /* Nor the lazily-loaded ParallelMixed package (issue #111): without the
     * {"Mathilda`"} request its LoadModule["mixed/logrewrite.m"] became an
     * inert ParallelMixed`Private`LoadModule and Integrate printed
     * "logrewrite.m not found". */
    assert_eval_eq("LoadModule[\"mixed/ParallelMixed.m\"]", "True", 0);
    assert_eval_eq("Select[Names[\"ParallelMixed`Private`*\"], "
                   "MemberQ[Names[\"Mathilda`*\"], Last[StringSplit[#, \"`\"]]] &]", "{}", 0);
}

/* Global`x and x are the same symbol, whichever spelling came first. */
static void test_global_qualified_is_bare(void) {
    assert_eval_eq("gq2 = 2; {Global`gq2, Names[\"Global`gq2\"]}", "{2, {\"gq2\"}}", 0);
    assert_eval_eq("Global`gq1 = 1; {gq1, Names[\"Global`gq1\"]}", "{1, {\"gq1\"}}", 0);
    assert_eval_eq("Global`gq3 := 3; Global`gq3 =.; gq3", "gq3", 0);
    /* A System name is not aliased by its Global` spelling. */
    assert_eval_eq("Global`Sin === Sin", "False", 0);
}

/* Developer`ToPackedArray / FromPackedArray / PackedArrayQ are the bare
 * builtins: the Developer` spelling reaches the same symbol. */
static void test_developer_packing_aliases(void) {
    assert_eval_eq("Developer`PackedArrayQ[Developer`ToPackedArray[{1, 2, 3}]]", "True", 0);
    assert_eval_eq("Developer`PackedArrayQ[Developer`FromPackedArray[ToPackedArray[{1, 2}]]]", "False", 0);
    assert_eval_eq("Developer`ToPackedArray === ToPackedArray", "True", 0);
    assert_eval_eq("Context[Developer`PackedArrayQ] === Context[PackedArrayQ]", "True", 0);
    /* Only the listed names alias: any other Developer` name is its own symbol. */
    assert_eval_eq("Context[Developer`NotAnAliasName]", "\"Developer`\"", 0);
}

int main(void) {
    symtab_init();
    core_init();

    TEST(test_initial_state);
    TEST(test_resolve_bare_builtin);
    TEST(test_resolve_system_prefix_stripped);
    TEST(test_resolve_dollar_unqualified);
    TEST(test_resolve_absolute);
    TEST(test_begin_end_roundtrip);
    TEST(test_begin_relative);
    /* test_dollar_context_and_path expects $ContextPath to be exactly
     * {"Global`", "System`"}. It must run before test_begin_package_path
     * and test_begin_package_with_needs, which permanently prepend their
     * package context to $ContextPath via EndPackage (see the comment in
     * test_begin_package_path about "other tests are tolerant"). */
    TEST(test_context_builtin_no_args);
    TEST(test_context_builtin_of_builtin);
    TEST(test_dollar_context_and_path);
    TEST(test_begin_package_path);
    TEST(test_begin_package_with_needs);
    TEST(test_begin_changes_context);
    TEST(test_package_symbol_qualification);
    TEST(test_nested_begin_private);
    TEST(test_system_prefix_stripped_in_parser);
    TEST(test_context_of_user_symbol);
    TEST(test_package_dollar_symbol);
    TEST(test_context_of_system_symbols);
    TEST(test_inert_frontend_names_are_system);
    TEST(test_mathilda_context);
    TEST(test_package_does_not_see_global_or_mathilda);
    TEST(test_package_can_request_mathilda);
    TEST(test_global_qualified_is_bare);
    TEST(test_contextpath_assignment);
    TEST(test_developer_packing_aliases);
    TEST(test_context_assignment);

    printf("All context tests passed!\n");
    return 0;
}
