/*
 * test_needs.c -- Needs, $Packages, $ContextAliases, $Path, and the
 * context-name / $Path search in Get.
 *
 * Every test writes small packages into a private scratch directory
 * (/tmp/Mathilda_needs_<pid>/) and points $Path at it, so the suite is
 * independent of the CTest working directory.  Packages bump a Global`
 * counter when they are read, which is how "loaded exactly once" is checked.
 * Messages are observed through Check[..., "msg"], which also proves they
 * route through the Quiet/Check funnel.
 */

#include "test_utils.h"
#include "expr.h"
#include "eval.h"
#include "parse.h"
#include "core.h"
#include "symtab.h"
#include "context.h"
#include "readwrite.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <sys/stat.h>

static char g_dir[256];

/* Write `body` to <g_dir>/<rel>, creating up to two levels of parent dirs. */
static void write_file(const char* rel, const char* body) {
    char path[512];
    snprintf(path, sizeof(path), "%s/%s", g_dir, rel);
    for (char* p = path + strlen(g_dir) + 1; *p; p++) {
        if (*p == '/') {
            *p = '\0';
            mkdir(path, 0755);
            *p = '/';
        }
    }
    FILE* fp = fopen(path, "w");
    assert(fp != NULL);
    fputs(body, fp);
    fclose(fp);
}

/* A package source for context `ctx` exporting `fn` (fn[x] = x^2) that
 * increments the Global` <counter> every time it is read.  The counter is
 * written bare, before BeginPackage, where $Context is still Global`.
 *
 * Note: a statement is parsed whole before any of it runs, so a package
 * symbol (or an alias) is only usable from the NEXT statement after the
 * Get/Needs that creates it -- exactly as in WL. */
static void write_package(const char* rel, const char* ctx, const char* fn,
                          const char* counter) {
    char body[1024];
    snprintf(body, sizeof(body),
             "%s = If[IntegerQ[%s], %s + 1, 1];\n"
             "BeginPackage[\"%s\"];\n"
             "%s::usage = \"%s[x] squares x.\";\n"
             "Begin[\"`Private`\"];\n"
             "%s[x_] := x^2;\n"
             "End[];\n"
             "EndPackage[];\n",
             counter, counter, counter, ctx, fn, fn, fn);
    write_file(rel, body);
}

static void eval_discard(const char* input) {
    Expr* parsed = parse_expression(input);
    assert(parsed != NULL);
    Expr* r = evaluate(parsed);
    expr_free(parsed);
    expr_free(r);
}

static void setup_scratch(void) {
    snprintf(g_dir, sizeof(g_dir), "/tmp/Mathilda_needs_%d", (int)getpid());
    mkdir(g_dir, 0755);

    write_package("NdA.m", "NdA`", "ndaF", "ndCountA");
    write_package("NdDir/init.m", "NdDir`", "nddirF", "ndCountDir");
    write_package("NdKer/Kernel/init.m", "NdKer`", "ndkerF", "ndCountKer");
    write_package("NdWl.wl", "NdWl`", "ndwlF", "ndCountWl");
    write_package("NdOuter/NdInner.m", "NdOuter`NdInner`", "ndinF", "ndCountIn");
    write_package("NdAl.m", "NdAl`", "ndalF", "ndCountAl");
    write_package("NdNone.m", "NdNone`", "ndnoneF", "ndCountNone");
    write_package("custom_file.m", "NdCustom`", "ndcustF", "ndCountCust");
    write_package("NdTwice.m", "NdTwice`", "ndtwF", "ndCountTw");
    write_file("NdNoCtx.m", "ndNoCtxRead = True; 17\n");
    write_file("plain_file.m", "ndPlain = 42\n");

    char set_path[512];
    snprintf(set_path, sizeof(set_path), "$Path = {\"%s\"}", g_dir);
    eval_discard(set_path);
}

/* ---------------------------------------------------------------- state */

static void test_initial_state(void) {
    /* Run before setup_scratch, so $Path still has its default. */
    assert_eval_eq("$Path", "{\".\"}", 0);
    assert_eval_eq("$ContextAliases", "<||>", 0);
    assert_eval_eq("Head[$ContextAliases]", "Association", 0);
    /* The bootstrap's FullSimplify package is a kernel package: loaded (so
     * Needs never re-reads it) but not listed in $Packages. */
    assert_eval_eq("SubsetQ[$Packages, {\"System`\", \"Global`\"}]", "True", 0);
    assert_eval_eq("MemberQ[$Packages, \"FullSimplify`\"]", "False", 0);
    ASSERT(context_packages_contains("System`"));
    ASSERT(context_packages_contains("FullSimplify`"));
    ASSERT(!context_packages_contains("NdA`"));
}

static void test_attributes_and_docs(void) {
    assert_eval_eq("Attributes[Needs]", "{Protected}", 0);
    const char* doc = symtab_get_docstring("Needs");
    ASSERT(doc != NULL && strstr(doc, "$ContextAliases") != NULL);
    ASSERT(symtab_get_docstring("$Packages") != NULL);
    ASSERT(symtab_get_docstring("$ContextAliases") != NULL);
    ASSERT(symtab_get_docstring("$Path") != NULL);
    ASSERT(strstr(symtab_get_docstring("Get"), "$Path") != NULL);
    assert_eval_eq("Context[Needs]", "\"System`\"", 0);
    /* Protected: no user definitions. */
    assert_eval_eq("Check[Needs[\"x`\"] := 1, \"msg\"]", "\"msg\"", 0);
}

static void test_begin_package_records_packages(void) {
    assert_eval_eq("BeginPackage[\"NdBp`\"]; EndPackage[]; First[$Packages]",
                   "\"NdBp`\"", 0);
    ASSERT(context_packages_contains("NdBp`"));
    /* Re-entering moves it to the front without duplicating it. */
    assert_eval_eq("BeginPackage[\"NdBp2`\"]; EndPackage[]; "
                   "BeginPackage[\"NdBp`\"]; EndPackage[]; "
                   "{Take[$Packages, 2], Count[$Packages, \"NdBp`\"]}",
                   "{{\"NdBp`\", \"NdBp2`\"}, 1}", 0);
}

/* ------------------------------------------------------------------ Get */

static void test_get_context_file_forms(void) {
    /* Pkg.m; Get returns the file's last value (EndPackage[] -> Null). */
    assert_eval_eq("Get[\"NdTwice`\"]", "Null", 0);
    assert_eval_eq("{ndCountTw, ndtwF[3], First[$ContextPath]}", "{1, 9, \"NdTwice`\"}", 0);
    /* Get always re-reads (only Needs is load-once). */
    assert_eval_eq("Get[\"NdTwice`\"]; ndCountTw", "2", 0);
    /* Pkg/init.m */
    eval_discard("Get[\"NdDir`\"]");
    assert_eval_eq("{ndCountDir, nddirF[4]}", "{1, 16}", 0);
    /* Pkg/Kernel/init.m */
    eval_discard("Get[\"NdKer`\"]");
    assert_eval_eq("{ndCountKer, ndkerF[5]}", "{1, 25}", 0);
    /* Pkg.wl */
    eval_discard("Get[\"NdWl`\"]");
    assert_eval_eq("{ndCountWl, ndwlF[6]}", "{1, 36}", 0);
}

static void test_get_resolve_api(void) {
    char* p = readwrite_resolve("NdOuter`NdInner`");
    ASSERT(p != NULL && strstr(p, "NdOuter/NdInner.m") != NULL);
    free(p);
    p = readwrite_resolve("NdKer`");
    ASSERT(p != NULL && strstr(p, "NdKer/Kernel/init.m") != NULL);
    free(p);
    ASSERT(readwrite_resolve("NoSuchNd`") == NULL);
    ASSERT(readwrite_resolve("`NdA`") == NULL);      /* relative: not loadable */
    ASSERT(readwrite_resolve("") == NULL);
    ASSERT(readwrite_resolve("/no/such/abs/file.m") == NULL);
}

static void test_get_relative_file_along_path(void) {
    assert_eval_eq("Get[\"plain_file.m\"]", "42", 0);
    assert_eval_eq("ndPlain", "42", 0);
    /* Searched in $Path order: a missing first directory is skipped. */
    assert_eval_eq("Block[{$Path = Prepend[$Path, \"/no/such/dir\"]}, Get[\"plain_file.m\"]]",
                   "42", 0);
    /* Non-string and empty entries of $Path are ignored. */
    assert_eval_eq("Block[{$Path = Join[{3, \"\"}, $Path]}, Get[\"plain_file.m\"]]", "42", 0);
}

static void test_get_failures(void) {
    assert_eval_eq("Quiet[Get[\"NoSuchNd`\"]]", "$Failed", 0);
    assert_eval_eq("Check[Get[\"NoSuchNd`\"], \"msg\"]", "\"msg\"", 0);
    assert_eval_eq("Check[Get[\"no_such_file.m\"], \"msg\"]", "\"msg\"", 0);
    /* Not found when the directory holding it leaves $Path. */
    assert_eval_eq("Block[{$Path = {}}, Quiet[Get[\"plain_file.m\"]]]", "$Failed", 0);
    /* Non-string argument stays unevaluated. */
    assert_eval_eq("Get[3]", "Get[3]", 0);
}

/* ---------------------------------------------------------------- Needs */

static void test_needs_basic_loads_once(void) {
    assert_eval_eq("Needs[\"NdA`\"]", "Null", 0);
    ASSERT(context_packages_contains("NdA`"));
    assert_eval_eq("{ndCountA, ndaF[7], First[$ContextPath]}", "{1, 49, \"NdA`\"}", 0);
    /* Second Needs does not re-read the file. */
    assert_eval_eq("Needs[\"NdA`\"]; ndCountA", "1", 0);
    assert_eval_eq("Count[$ContextPath, \"NdA`\"]", "1", 0);
    assert_eval_eq("Count[$Packages, \"NdA`\"]", "1", 0);
    /* Off $ContextPath: Needs puts it back without reloading. */
    assert_eval_eq("$ContextPath = DeleteCases[$ContextPath, \"NdA`\"]; "
                   "Needs[\"NdA`\"]; {ndCountA, First[$ContextPath]}",
                   "{1, \"NdA`\"}", 0);
}

static void test_needs_already_in_packages(void) {
    /* System` is in $Packages: nothing is read, nothing fails. */
    assert_eval_eq("Check[Needs[\"System`\"], \"msg\"]", "Null", 0);
    /* A context created by BeginPackage in the session counts as loaded even
     * with no file anywhere. */
    assert_eval_eq("BeginPackage[\"NdInline`\"]; EndPackage[]; "
                   "Check[Needs[\"NdInline`\"], \"msg\"]", "Null", 0);
}

static void test_needs_nested_context(void) {
    assert_eval_eq("Needs[\"NdOuter`NdInner`\"]", "Null", 0);
    assert_eval_eq("{ndCountIn, ndinF[2]}", "{1, 4}", 0);
    assert_eval_eq("MemberQ[$ContextPath, \"NdOuter`NdInner`\"]", "True", 0);
}

static void test_needs_two_arg_file(void) {
    assert_eval_eq("Needs[\"NdCustom`\", \"custom_file.m\"]", "Null", 0);
    assert_eval_eq("{ndCountCust, ndcustF[9], MemberQ[$ContextPath, \"NdCustom`\"]}",
                   "{1, 81, True}", 0);
    /* Already loaded: the file is not read again. */
    assert_eval_eq("Needs[\"NdCustom`\", \"custom_file.m\"]; ndCountCust", "1", 0);
}

static void test_needs_nocont(void) {
    /* The file is read but never creates the context: Needs::nocont, and the
     * value of the Get comes back. */
    assert_eval_eq("Check[Needs[\"NdNever`\", \"NdNoCtx.m\"], \"msg\"]", "\"msg\"", 0);
    assert_eval_eq("Quiet[Needs[\"NdNever`\", \"NdNoCtx.m\"]]", "17", 0);
    assert_eval_eq("ndNoCtxRead", "True", 0);
    ASSERT(!context_packages_contains("NdNever`"));
    assert_eval_eq("MemberQ[$ContextPath, \"NdNever`\"]", "False", 0);
    /* No file at all: Get::noopen, Needs::nocont, $Failed. */
    assert_eval_eq("Quiet[Needs[\"NdMissing`\"]]", "$Failed", 0);
    assert_eval_eq("Check[Needs[\"NdMissing`\"], \"msg\"]", "\"msg\"", 0);
    assert_eval_eq("Quiet[Needs[\"NdMissing`\", \"no_such_file.m\"]]", "$Failed", 0);
    /* A failed load leaves $ContextPath alone. */
    assert_eval_eq("p0 = $ContextPath; Quiet[Needs[\"NdMissing`\"]]; $ContextPath === p0",
                   "True", 0);
}

static void test_needs_alias(void) {
    assert_eval_eq("p0 = $ContextPath; Needs[\"NdAl`\" -> \"NdA1`\"]", "Null", 0);
    ASSERT(context_packages_contains("NdAl`"));
    /* Not added to $ContextPath, although the package's EndPackage did. */
    assert_eval_eq("{$ContextPath === p0, MemberQ[$ContextPath, \"NdAl`\"]}",
                   "{True, False}", 0);
    assert_eval_eq("$ContextAliases[\"NdA1`\"]", "\"NdAl`\"", 0);
    assert_eval_eq("ndCountAl", "1", 0);
    /* The alias resolves at read time: NdA1`ndalF is NdAl`ndalF. */
    assert_eval_eq("NdA1`ndalF[5]", "25", 0);
    assert_eval_eq("Context[NdA1`ndalF]", "\"NdAl`\"", 0);
    assert_eval_eq("NdAl`ndalF === NdA1`ndalF", "True", 0);
    /* Not on the path: the short name is a fresh Global` symbol. */
    assert_eval_eq("Context[ndalF]", "\"Global`\"", 0);
    /* Private symbols are reachable through the alias too. */
    assert_eval_eq("Context[NdA1`Private`x]", "\"NdAl`Private`\"", 0);
    /* A second alias for the same, already-loaded context: no reload. */
    eval_discard("Needs[\"NdAl`\" -> \"NdA2`\"]");
    assert_eval_eq("{ndCountAl, NdA2`ndalF[3]}", "{1, 9}", 0);
    assert_eval_eq("Length[$ContextAliases]", "2", 0);
    /* Re-aliasing a key replaces its target. */
    eval_discard("Needs[\"NdA`\" -> \"NdA2`\"]");
    assert_eval_eq("{$ContextAliases[\"NdA2`\"], NdA2`ndaF[2]}", "{\"NdA`\", 4}", 0);
    assert_eval_eq("Length[$ContextAliases]", "2", 0);
}

static void test_needs_alias_assignable(void) {
    /* $ContextAliases may be set directly, as an Association or List of rules. */
    eval_discard("$ContextAliases = <|\"NdQ`\" -> \"NdA`\"|>");
    assert_eval_eq("NdQ`ndaF[3]", "9", 0);
    eval_discard("$ContextAliases = {\"NdR`\" -> \"NdA`\"}");
    assert_eval_eq("NdR`ndaF[4]", "16", 0);
    /* Needs on top of a List of rules rebuilds an Association keeping them. */
    assert_eval_eq("Needs[\"NdA`\" -> \"NdS`\"]; {Head[$ContextAliases], Keys[$ContextAliases]}",
                   "{Association, {\"NdR`\", \"NdS`\"}}", 0);
    /* Garbage is read as no aliases. */
    eval_discard("$ContextAliases = 7");
    assert_eval_eq("Context[NdS`zz]", "\"NdS`\"", 0);
    eval_discard("$ContextAliases = <||>");
    assert_eval_eq("Context[NdQ`zz]", "\"NdQ`\"", 0);
}

static void test_needs_none(void) {
    assert_eval_eq("p0 = $ContextPath; a0 = $ContextAliases; Needs[\"NdNone`\" -> None]",
                   "Null", 0);
    assert_eval_eq("{$ContextPath === p0, $ContextAliases === a0, ndCountNone}",
                   "{True, True, 1}", 0);
    ASSERT(context_packages_contains("NdNone`"));
    assert_eval_eq("NdNone`ndnoneF[6]", "36", 0);
    /* Already loaded via -> None: plain Needs now adds it to the path. */
    assert_eval_eq("Needs[\"NdNone`\"]; {ndCountNone, First[$ContextPath]}",
                   "{1, \"NdNone`\"}", 0);
}

static void test_needs_bad_arguments(void) {
    assert_eval_eq("Quiet[Needs[1]]", "Needs[1]", 0);
    assert_eval_eq("Check[Needs[1], \"msg\"]", "\"msg\"", 0);
    assert_eval_eq("Quiet[Needs[\"NoTick\"]]", "Needs[\"NoTick\"]", 0);
    assert_eval_eq("Check[Needs[\"NoTick\"], \"msg\"]", "\"msg\"", 0);
    assert_eval_eq("Check[Needs[\"`Rel`\"], \"msg\"]", "\"msg\"", 0);
    assert_eval_eq("Check[Needs[\"A``B`\"], \"msg\"]", "\"msg\"", 0);
    assert_eval_eq("Check[Needs[\"1A`\"], \"msg\"]", "\"msg\"", 0);
    assert_eval_eq("Check[Needs[\"Ok`\" -> \"bad\"], \"msg\"]", "\"msg\"", 0);
    assert_eval_eq("Check[Needs[\"Ok`\" -> 3], \"msg\"]", "\"msg\"", 0);
    assert_eval_eq("Check[Needs[x -> \"Al`\"], \"msg\"]", "\"msg\"", 0);
    assert_eval_eq("Check[Needs[\"NdA`\", 3], \"msg\"]", "\"msg\"", 0);
    assert_eval_eq("Quiet[Needs[\"NdA`\", 3]]", "Needs[\"NdA`\", 3]", 0);
    assert_eval_eq("Check[Needs[], \"msg\"]", "\"msg\"", 0);
    assert_eval_eq("Quiet[Needs[\"a`\", \"b\", \"c\"]]", "Needs[\"a`\", \"b\", \"c\"]", 0);
    /* Bad arguments never touch the state. */
    assert_eval_eq("p0 = {$ContextPath, $Packages, $ContextAliases}; "
                   "Quiet[Needs[\"NoTick\"]; Needs[1]; Needs[\"Ok`\" -> 3]]; "
                   "p0 === {$ContextPath, $Packages, $ContextAliases}", "True", 0);
}

static void test_needs_inside_package(void) {
    /* Needs inside BeginPackage extends the package's restricted path; the
     * caller's $ContextPath comes back from EndPackage. */
    assert_eval_eq("BeginPackage[\"NdHost`\"]; Needs[\"NdA`\"]; "
                   "inPath = $ContextPath; EndPackage[]; "
                   "{MemberQ[inPath, \"NdA`\"], First[$ContextPath]}",
                   "{True, \"NdHost`\"}", 0);
}

static void test_alias_resolver_direct(void) {
    eval_discard("$ContextAliases = <|\"Zq`\" -> \"NdA`\", \"Zq`Sub`\" -> \"NdDir`\", "
                 "\"Zs`\" -> \"System`\"|>");
    char* r = context_resolve_name("Zq`f");
    ASSERT_STR_EQ(r, "NdA`f");
    free(r);
    /* Longest alias prefix wins. */
    r = context_resolve_name("Zq`Sub`g");
    ASSERT_STR_EQ(r, "NdDir`g");
    free(r);
    r = context_resolve_name("Zq`Other`g");
    ASSERT_STR_EQ(r, "NdA`Other`g");
    free(r);
    /* An alias for System` still canonicalises to the bare name. */
    r = context_resolve_name("Zs`Plus");
    ASSERT_STR_EQ(r, "Plus");
    free(r);
    /* Only whole leading contexts match. */
    r = context_resolve_name("ZqZ`f");
    ASSERT_STR_EQ(r, "ZqZ`f");
    free(r);
    r = context_resolve_name("Abc`Zq`f");
    ASSERT_STR_EQ(r, "Abc`Zq`f");
    free(r);
    /* Bare names are unaffected. */
    r = context_resolve_name("Plus");
    ASSERT_STR_EQ(r, "Plus");
    free(r);
    eval_discard("$ContextAliases = <||>");
}

static void cleanup_scratch(void) {
    char cmd[512];
    snprintf(cmd, sizeof(cmd), "rm -rf '%s'", g_dir);
    if (system(cmd) != 0) fprintf(stderr, "warning: could not remove %s\n", g_dir);
}

int main(void) {
    symtab_init();
    core_init();

    TEST(test_initial_state);
    setup_scratch();
    TEST(test_attributes_and_docs);
    TEST(test_begin_package_records_packages);
    TEST(test_get_context_file_forms);
    TEST(test_get_resolve_api);
    TEST(test_get_relative_file_along_path);
    TEST(test_get_failures);
    TEST(test_needs_basic_loads_once);
    TEST(test_needs_already_in_packages);
    TEST(test_needs_nested_context);
    TEST(test_needs_two_arg_file);
    TEST(test_needs_nocont);
    TEST(test_needs_alias);
    TEST(test_needs_alias_assignable);
    TEST(test_needs_none);
    TEST(test_needs_bad_arguments);
    TEST(test_needs_inside_package);
    TEST(test_alias_resolver_direct);
    cleanup_scratch();

    printf("All Needs tests passed!\n");
    return 0;
}
