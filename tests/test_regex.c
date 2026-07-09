/*
 * test_regex.c - unit tests for the PCRE2 wrapper (regex_engine.h), the
 * $n template expander (regex_common.h), and the inert RegularExpression[]
 * head.  When Mathilda is built without PCRE2 (USE_REGEX undefined) the
 * engine-level checks are skipped; the head still stays inert.
 */

#include "test_utils.h"
#include "expr.h"
#include "eval.h"
#include "core.h"
#include "symtab.h"
#include "regex_engine.h"
#include "regex_common.h"

#include <stdlib.h>
#include <string.h>

/* ============================ engine wrapper ============================ */

#ifdef USE_REGEX

static void test_engine_available(void) {
    ASSERT(regex_available() == 1);
}

static void test_engine_compile_match(void) {
    char err[256];
    RegexProgram* p = regex_compile("a(\\d+)b", err, sizeof err);
    ASSERT(p != NULL);
    ASSERT(regex_group_count(p) == 1);

    const char* s = "xxa123byy";
    size_t ov[4];
    int r = regex_match(p, s, strlen(s), 0, ov, 2);
    ASSERT(r == 1);
    ASSERT(ov[0] == 2 && ov[1] == 7);   /* whole match "a123b" */
    ASSERT(ov[2] == 3 && ov[3] == 6);   /* group 1 "123"       */

    /* No further match past the first one. */
    r = regex_match(p, s, strlen(s), 7, ov, 2);
    ASSERT(r == 0);

    regex_free(p);
}

static void test_engine_unset_group(void) {
    char err[256];
    RegexProgram* p = regex_compile("(a)|(b)", err, sizeof err);
    ASSERT(p != NULL);
    ASSERT(regex_group_count(p) == 2);

    const char* s = "b";
    size_t ov[6];
    int r = regex_match(p, s, strlen(s), 0, ov, 3);
    ASSERT(r == 1);
    ASSERT(ov[2] == REGEX_UNSET);       /* group 1 (a) did not participate */
    ASSERT(ov[4] == 0 && ov[5] == 1);   /* group 2 (b) matched            */
    regex_free(p);
}

static void test_engine_bad_pattern(void) {
    char err[256];
    err[0] = '\0';
    RegexProgram* p = regex_compile("a(", err, sizeof err);
    ASSERT(p == NULL);
    ASSERT(err[0] != '\0');             /* a diagnostic was written */
}

/* ============================ $n expansion ============================= */

static void test_template_expand(void) {
    const char* s = "xxa123byy";
    size_t ov[4] = {2, 7, 3, 6};        /* whole "a123b", group1 "123" */
    char* out = regex_expand_template("[$1]<$0>", s, ov, 2);
    ASSERT_STR_EQ(out, "[123]<a123b>");
    free(out);

    out = regex_expand_template("a$$b", s, ov, 2);   /* $$ -> literal $ */
    ASSERT_STR_EQ(out, "a$b");
    free(out);

    out = regex_expand_template("$9", s, ov, 2);      /* absent group -> "" */
    ASSERT_STR_EQ(out, "");
    free(out);
}

#endif /* USE_REGEX */

/* ===================== RegularExpression head (inert) ================== */

static void test_regularexpression_inert(void) {
    /* RegularExpression["re"] evaluates to itself (an inert data head). */
    assert_eval_eq("RegularExpression[\"a+\"]", "RegularExpression[\"a+\"]", 0);
    assert_eval_eq("Head[RegularExpression[\"a+\"]]", "RegularExpression", 0);
}

int main(void) {
    symtab_init();
    core_init();

#ifdef USE_REGEX
    TEST(test_engine_available);
    TEST(test_engine_compile_match);
    TEST(test_engine_unset_group);
    TEST(test_engine_bad_pattern);
    TEST(test_template_expand);
#else
    printf("USE_REGEX not defined; skipping engine-level regex tests\n");
#endif

    TEST(test_regularexpression_inert);

    printf("All regex tests passed!\n");
    return 0;
}
