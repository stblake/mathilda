/*
 * test_message.c -- Message formatting and Off / On.
 *
 * Message[sym::tag, args] prints "sym::tag: text" through the mth_message
 * funnel, with the `n` / `` slots of the text filled from args (OutputForm).
 * The text is sym::tag's, else General::tag's (a built-in table of the
 * standard texts), else "-- Message text not found --".  Off[sym::tag]
 * removes the message entirely -- not printed, not seen by Check.
 *
 * Printed text is checked by capturing stderr (fd 2) into a scratch file.
 */
#include "test_utils.h"
#include "symtab.h"
#include "core.h"

#include <fcntl.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

/* Evaluate `src` with stderr captured; returns what was written (caller frees). */
static char* stderr_of(const char* src) {
    char path[128];
    snprintf(path, sizeof path, "/tmp/Mathilda_message_%d", (int)getpid());
    fflush(stderr);
    int saved = dup(2);
    int fd = open(path, O_WRONLY | O_CREAT | O_TRUNC, 0600);
    ASSERT(saved >= 0 && fd >= 0);
    dup2(fd, 2);
    close(fd);
    Expr* e = parse_expression(src);
    ASSERT(e != NULL);
    Expr* v = evaluate(e);
    expr_free(e);
    expr_free(v);
    fflush(stderr);
    dup2(saved, 2);
    close(saved);
    FILE* fp = fopen(path, "r");
    ASSERT(fp != NULL);
    char* buf = calloc(1, 8192);
    size_t n = fread(buf, 1, 8191, fp);
    buf[n] = '\0';
    fclose(fp);
    unlink(path);
    return buf;
}

static void assert_stderr(const char* src, const char* expected) {
    char* got = stderr_of(src);
    if (strcmp(got, expected) != 0) {
        printf("FAIL: %s\n  Expected: [%s]\n  Actual:   [%s]\n", src, expected, got);
        ASSERT(0);
    }
    free(got);
}

static void test_message_formatting(void) {
    assert_eval_eq("tm::m = \"bad `1` and `2`\"", "\"bad `1` and `2`\"", 0);
    assert_stderr("Message[tm::m, 7, x^2]", "tm::m: bad 7 and x^2\n");
    assert_stderr("Message[tm::m, \"str\", {1, 2}]", "tm::m: bad str and {1, 2}\n");
    /* A missing argument: StringForm::sfr, the slot left as written. */
    assert_stderr("Message[tm::m, 7]",
                  "StringForm::sfr: Item 2 requested in \"bad `1` and `2`\" out of range; 1 items available.\n"
                  "tm::m: bad 7 and `2`\n");
    /* Sequential `` slots; extra arguments are ignored. */
    assert_eval_eq("tm::t = \"`` then ``\"", "\"`` then ``\"", 0);
    assert_stderr("Message[tm::t, a, b]", "tm::t: a then b\n");
    assert_eval_eq("tm::n = \"no slots\"", "\"no slots\"", 0);
    assert_stderr("Message[tm::n, 3]", "tm::n: no slots\n");
    /* No text of its own: General::tag, else not found. */
    assert_stderr("Message[tm::argx, tm, 3]", "tm::argx: tm called with 3 arguments; 1 argument is expected.\n");
    assert_stderr("Message[tm::undefinedtag]", "tm::undefinedtag: -- Message text not found --\n");
    assert_stderr("Message[tm::undefinedtag, 1, 2]",
                  "tm::undefinedtag: -- Message text not found -- (1) (2)\n");
    /* A General::tag definition takes precedence over the built-in text. */
    assert_eval_eq("General::tmown = \"own `1`\"", "\"own `1`\"", 0);
    assert_stderr("Message[tm::tmown, 5]", "tm::tmown: own 5\n");
    assert_eval_eq("Message[tm::n]", "Null", 0);
    /* Quiet silences the print; Check still sees it. */
    assert_stderr("Quiet[Message[tm::n]]", "");
    assert_eval_eq("Check[Message[tm::n]; 1, \"caught\"]", "\"caught\"", 0);
    assert_eval_eq("Quiet[Check[Message[tm::n]; 1, \"caught\"]]", "\"caught\"", 0);
    assert_eval_eq("{Attributes[Message], Attributes[Off], Attributes[On]}",
                   "{{HoldFirst, Protected}, {HoldAll, Protected}, {HoldAll, Protected}}", 0);
}

static void test_off_on(void) {
    /* Off: neither printed nor seen by Check. */
    assert_stderr("Off[tm::n]; Message[tm::n]", "");
    assert_eval_eq("Check[Message[tm::n]; 1, \"caught\"]", "1", 0);
    assert_stderr("On[tm::n]; Message[tm::n]", "tm::n: no slots\n");
    /* Several at once, and a list. */
    assert_stderr("Off[tm::n, tm::m]; Message[tm::n]; Message[tm::m, 1, 2]", "");
    assert_stderr("On[{tm::n, tm::m}]; Message[tm::m, 1, 2]", "tm::m: bad 1 and 2\n");
    /* General::tag switches tag off for every symbol -- built-in messages too. */
    assert_stderr("Off[General::argx]; Message[tm::argx, tm, 3]; Message[g::argx, g, 3]; Sin[1, 2]", "");
    assert_stderr("On[General::argx]; Message[tm::argx, tm, 3]",
                  "tm::argx: tm called with 3 arguments; 1 argument is expected.\n");
    /* A built-in diagnostic through the funnel. */
    assert_stderr("Off[Environment::string]; Environment[3]", "");
    assert_eval_eq("Check[Environment[3], \"caught\"]", "Environment[3]", 0);
    assert_stderr("On[Environment::string]; Environment[3]",
                  "Environment::string: String expected at position 1 in Environment[3].\n");
    /* Off[sym]: every message of sym. */
    assert_stderr("Off[tm]; Message[tm::n]; Message[tm::m, 1, 2]", "");
    assert_stderr("On[tm]; Message[tm::n]", "tm::n: no slots\n");
    /* A message group is not supported: Off::mgre. */
    assert_eval_eq("Check[Off[\"Global`*\"], \"caught\"]", "\"caught\"", 0);
    assert_eval_eq("Off[tm::n]", "Null", 0);
    assert_eval_eq("On[tm::n]", "Null", 0);
}

static void test_message_list(void) {
    /* Each PRINTED message, as HoldForm[sym::tag]; Quiet / Off ones are not listed. */
    assert_eval_eq("FullForm[Block[{$MessageList = {}}, Quiet[Message[tm::n]]; Off[tm::n]; Message[tm::n]; "
                   "On[tm::n]; Quiet[1/0]; $MessageList]]", "List[]", 0);
    assert_stderr("ml = Block[{$MessageList = {}}, Message[tm::n]; 1/0; Environment[3]; $MessageList]",
                  "tm::n: no slots\nPower::infy: Infinite expression 1/0 encountered.\n"
                  "Environment::string: String expected at position 1 in Environment[3].\n");
    assert_eval_eq("FullForm[ml]", "List[HoldForm[MessageName[tm, \"n\"]], HoldForm[MessageName[Power, \"infy\"]], "
                   "HoldForm[MessageName[Environment, \"string\"]]]", 0);
}

/* The arithmetic, NDArray and linear-algebra warning helpers emit through the
 * funnel: Off, Quiet and Check all apply to them. */
static void test_helper_warnings_routed(void) {
    assert_stderr("Off[Power::infy]; r = 1/0", "");
    assert_eval_eq("Check[1/0, \"caught\"]", "ComplexInfinity", 0);
    assert_eval_eq("On[Power::infy]; Quiet[Check[1/0, \"caught\"]]", "\"caught\"", 0);
    assert_eval_eq("Quiet[Check[NDArray[{1, 2}, DataType -> \"bogus\"], \"caught\"]]", "\"caught\"", 0);
    assert_stderr("Quiet[Clip[{x}, {0, 1}]]", "");
}

/* A site guarded by arith_warnings_muted() returned before the funnel under
 * Quiet, so an enclosing Check never saw it (Check took the wrong branch).
 * These sites now gate on the raw internal mute; tools/check_message_routing.py
 * rejects the old guard.  Check also takes several message names. */
static void test_check_sees_quiet_guarded_sites(void) {
    assert_eval_eq("Quiet[Check[Product[2, {k, 1, Infinity}], \"caught\"]]", "\"caught\"", 0);
    assert_eval_eq("Quiet[Check[Sum[1/k, {k, 1, Infinity}], \"caught\"]]", "\"caught\"", 0);
    assert_eval_eq("Quiet[Check[Prime[0], \"caught\"]]", "\"caught\"", 0);
    assert_eval_eq("Quiet[Check[ExtendedGCD[3/2, 2], \"caught\"]]", "\"caught\"", 0);
    assert_eval_eq("Quiet[Check[MoebiusMu[1, 2], \"caught\"]]", "\"caught\"", 0);
    assert_eval_eq("Quiet[Check[1/0, \"caught\", Sum::div, Power::infy]]", "\"caught\"", 0);
    assert_eval_eq("Quiet[Check[1/0, \"caught\", Sum::div, Product::div]]", "ComplexInfinity", 0);
    assert_eval_eq("Quiet[Check[1/0, \"caught\", {Sum::div}, {Power::infy}]]", "\"caught\"", 0);
}

static void test_printtemporary(void) {
    assert_eval_eq("Quiet[PrintTemporary[]]; Attributes[PrintTemporary]", "{Protected}", 0);
}

/* Style is presentation only: OutputForm (Print, ToString) shows the content. */
static void test_style_outputform(void) {
    assert_stderr("Print[1]", "");      /* Print writes stdout; nothing on stderr */
    assert_eval_eq("{ToString[Style[\"s\", Bold]], ToString[Style[x^2, Red, 12]], ToString[{Style[a, Italic], b}]}",
                   "{\"s\", \"x^2\", \"{a, b}\"}", 0);
    assert_eval_eq("FullForm[Style[\"s\", Bold]]", "Style[\"s\", Bold]", 0);
    assert_eval_eq("{Context[DateString], Context[FileDate]}", "{\"System`\", \"System`\"}", 0);
}

/* Quiet / Check nesting (src/message.c frame model).  Every row was checked
 * against the reference: a message is silenced by the innermost Quiet whose
 * spec names it and is seen only by the Checks inside that Quiet. */
static void test_check_quiet_nesting(void) {
    assert_eval_eq("Check[Quiet[1/0]; 5, 7]", "5", 0);
    assert_eval_eq("Quiet[Check[1/0; 5, 7]]", "7", 0);
    assert_eval_eq("Check[Quiet[Check[1/0; 5, 7]], 9]", "7", 0);
    assert_eval_eq("Quiet[Check[Check[1/0, a], b]]", "b", 0);
    assert_eval_eq("Quiet[Check[Quiet[1/0]; 5, 7]]", "5", 0);
    assert_eval_eq("Check[Quiet[Message[cqn::x]]; 1, 2]", "1", 0);
    assert_eval_eq("Check[Off[Power::infy]; 1/0; On[Power::infy]; 5, 7]", "5", 0);
    /* Quiet with a message spec silences only what it names. */
    assert_eval_eq("Check[Quiet[1/0, Power::infy]; 5, 7]", "5", 0);
    assert_eval_eq("Quiet[Check[Quiet[1/0, General::foo]; 5, 7]]", "7", 0);
    assert_eval_eq("Check[Quiet[1/0, {General::foo, Power::infy}]; 5, 7]", "5", 0);
    assert_eval_eq("Check[Quiet[1/0, All]; 5, 7]", "5", 0);
    assert_eval_eq("Quiet[Check[Quiet[1/0, None]; 5, 7]]", "7", 0);
    assert_eval_eq("Check[Quiet[Message[cqn::argx, cqn, 2], General::argx]; 5, 7]", "5", 0);
    /* The spec is held: a symbol whose value is a list names nothing. */
    assert_eval_eq("cqnms = {Power::infy}; Quiet[Check[Quiet[1/0, cqnms]; 5, 7]]", "7", 0);
    assert_stderr("Quiet[1/0, General::foo]", "Power::infy: Infinite expression 1/0 encountered.\n");
    /* Check with a spec counts only what it names. */
    assert_eval_eq("Quiet[Check[1/0; 5, 7, Power::infy]]", "7", 0);
    assert_eval_eq("Quiet[Check[1/0; 5, 7, General::foo]]", "5", 0);
    assert_eval_eq("Quiet[Check[1/0; 5, 7, {Power::infy}]]", "7", 0);
    /* $MessageList holds only the messages that were not silenced. */
    assert_eval_eq("Block[{$MessageList = {}}, Quiet[1/0]; $MessageList]", "{}", 0);
    assert_eval_eq("Block[{$MessageList = {}}, Quiet[Check[1/0, 0]]; $MessageList]", "{}", 0);
    /* A TimeConstrained timeout inside Check / Quiet leaves no stale frame. */
    assert_eval_eq("TimeConstrained[Check[Quiet[While[True, 1]], 0], 0.2, t]; Check[1/0; 5, 7]", "7", 1);
}

int main(void) {
    symtab_init();
    core_init();
    TEST(test_message_formatting);
    TEST(test_off_on);
    TEST(test_message_list);
    TEST(test_helper_warnings_routed);
    TEST(test_check_sees_quiet_guarded_sites);
    TEST(test_printtemporary);
    TEST(test_style_outputform);
    TEST(test_check_quiet_nesting);
    printf("All message tests passed.\n");
    return 0;
}
