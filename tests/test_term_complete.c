/*
 * test_term_complete.c -- Tab completion candidates (src/term/term_complete.c,
 * names_complete in src/names.c): symbols, the options of the enclosing call
 * head, strings and comments, context-qualified names, Module temporaries,
 * and the guarantee that completing never grows the symbol table.
 */
#include "test_utils.h"
#include "symtab.h"
#include "core.h"
#include "term/term_complete.h"

#include <regex.h>
#include <stdio.h>
#include <string.h>

static void run(const char* src) {
    Expr* parsed = parse_expression(src);
    ASSERT(parsed != NULL);
    Expr* result = evaluate(parsed);
    expr_free(parsed);
    expr_free(result);
}

static TermCompletions comp(const char* buf) {
    TermCompletions c;
    term_complete(buf, strlen(buf), &c);
    return c;
}

static bool has(const TermCompletions* c, const char* s) {
    for (size_t i = 0; i < c->count; i++) if (strcmp(c->items[i], s) == 0) return true;
    return false;
}

static void test_symbols(void) {
    /* 1. A symbol prefix. */
    TermCompletions c = comp("Integr");
    ASSERT(has(&c, "Integrate"));
    for (size_t i = 0; i < c.count; i++)
        ASSERT_MSG(strncmp(c.items[i], "Integr", 6) == 0, "%s", c.items[i]);
    /* Sorted, no duplicates. */
    for (size_t i = 1; i < c.count; i++) ASSERT(strcmp(c.items[i - 1], c.items[i]) < 0);
    ASSERT(c.word_start == 0);
    term_completions_free(&c);

    /* The word after an operator. */
    c = comp("1 + Integr");
    ASSERT(c.word_start == 4 && has(&c, "Integrate"));
    term_completions_free(&c);

    /* 10. A defined user symbol is offered, a bare interned one is not. */
    run("myvar = 3");
    run("qqqUnbound");
    c = comp("myv");
    ASSERT(has(&c, "myvar"));
    term_completions_free(&c);
    c = comp("qqqUnb");
    ASSERT(c.count == 0);
    term_completions_free(&c);
}

static void test_options(void) {
    /* 2. A builtin head's registered options come first. */
    TermCompletions c = comp("ReplaceRepeated[x, {}, MaxI");
    ASSERT(c.count >= 1);
    ASSERT_STR_EQ(c.items[0], "MaxIterations -> ");
    ASSERT(!has(&c, "MaxIterations"));      /* the option covers the symbol */
    term_completions_free(&c);

    /* 3. User options, symbol and string named. */
    run("Options[myf] = {Alpha -> 1, \"Beta\" -> 2}");
    c = comp("myf[Al");
    ASSERT(c.count >= 1);
    ASSERT_STR_EQ(c.items[0], "Alpha -> ");
    term_completions_free(&c);
    c = comp("myf[\"Be");
    ASSERT(c.count == 1);
    ASSERT_STR_EQ(c.items[0], "Beta\" -> ");
    ASSERT(c.word_start == 5);
    term_completions_free(&c);
    c = comp("myf[1, Al");
    ASSERT(c.count >= 1);
    ASSERT_STR_EQ(c.items[0], "Alpha -> ");
    term_completions_free(&c);
    c = comp("myf[x -> Al");               /* not an argument position */
    ASSERT(!has(&c, "Alpha -> "));
    term_completions_free(&c);
    /* An empty word in argument position: options only. */
    c = comp("myf[1, ");
    ASSERT(c.count == 1);
    ASSERT_STR_EQ(c.items[0], "Alpha -> ");
    term_completions_free(&c);

    /* 6. A Part is not a call. */
    c = comp("a[[Ma");
    for (size_t i = 0; i < c.count; i++) ASSERT(strstr(c.items[i], " -> ") == NULL);
    term_completions_free(&c);
}

static void test_temporaries(void) {
    /* 4. Module temporaries are hidden. */
    run("Module[{zz}, zz]");
    run("Module[{zz}, zz = 1; zz]");
    TermCompletions c = comp("zz");
    regex_t rx;
    ASSERT(regcomp(&rx, "\\$[0-9]+$", REG_EXTENDED | REG_NOSUB) == 0);
    for (size_t i = 0; i < c.count; i++)
        ASSERT_MSG(regexec(&rx, c.items[i], 0, NULL, 0) != 0, "temporary %s offered", c.items[i]);
    regfree(&rx);
    term_completions_free(&c);
}

static void test_strings_and_comments(void) {
    /* 5. Nothing inside a comment or a plain string. */
    TermCompletions c = comp("(* Integr");
    ASSERT(c.count == 0);
    term_completions_free(&c);
    c = comp("\"Integr");
    ASSERT(c.count == 0);
    term_completions_free(&c);
    /* A closed comment does not block what follows. */
    c = comp("(* x *) Integr");
    ASSERT(has(&c, "Integrate"));
    term_completions_free(&c);
}

static void test_context(void) {
    /* 7. A backticked prefix matches fully qualified names. */
    TermCompletions c = comp("System`Integr");
    ASSERT(has(&c, "System`Integrate"));
    ASSERT(c.word_start == 0);
    term_completions_free(&c);
}

static void test_empty(void) {
    /* 8. An empty word outside argument position completes nothing. */
    TermCompletions c = comp("");
    ASSERT(c.count == 0);
    term_completions_free(&c);
    c = comp("1 + ");
    ASSERT(c.count == 0);
    term_completions_free(&c);
    /* A number is not a symbol prefix. */
    c = comp("12");
    ASSERT(c.count == 0);
    term_completions_free(&c);
    term_complete(NULL, 0, &c);
    ASSERT(c.count == 0);
    term_completions_free(&c);
}

static size_t names_length(void) {
    Expr* p = parse_expression("Length[Names[]]");
    Expr* r = evaluate(p);
    expr_free(p);
    ASSERT(r && r->type == EXPR_INTEGER);
    size_t n = (size_t)r->data.integer;
    expr_free(r);
    return n;
}

static void test_no_growth(void) {
    /* 9. Completing never interns a symbol. */
    size_t before = names_length();
    for (int i = 0; i < 100; i++) {
        char buf[64];
        snprintf(buf, sizeof(buf), "neverDefined%dHead[Opt%d", i, i);
        TermCompletions c = comp(buf);
        term_completions_free(&c);
        c = comp("Integr");
        term_completions_free(&c);
    }
    ASSERT(names_length() == before);
}

int main(void) {
    symtab_init();
    core_init();
    TEST(test_symbols);
    TEST(test_options);
    TEST(test_temporaries);
    TEST(test_strings_and_comments);
    TEST(test_context);
    TEST(test_empty);
    TEST(test_no_growth);
    printf("All term_complete tests passed.\n");
    return 0;
}
