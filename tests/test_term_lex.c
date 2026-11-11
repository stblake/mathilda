/*
 * test_term_lex.c -- the interactive terminal's input lexer (src/term/term_lex.c):
 * the completeness predicate behind the smart Return key, and the cursor
 * context used by Tab completion.
 */
#include "test_utils.h"
#include "term/term_lex.h"

#include <stdio.h>
#include <string.h>

static void test_completeness(void) {
    struct { const char* s; bool complete; } rows[] = {
        { "f[",            false },
        { "\"a",           false },
        { "(* x",          false },
        { "f[1]]",         true  },   /* over-closed submits */
        { "",              true  },
        { "\"(*\"",        true  },   /* a comment opener inside a string */
        { "(* \" *)",      true  },   /* a quote inside a comment */
        { "f[\n1]",        true  },
        { "{1, (2",        false },
        { "(* (* *) *)",   true  },   /* nested comments */
        { "(* (* *)",      false },
        { "\"a\\\"b\"",    true  },   /* escaped quote */
        { "\"a\\\"b",      false },
    };
    for (size_t i = 0; i < sizeof(rows) / sizeof(rows[0]); i++) {
        bool got = term_lex_input_complete(rows[i].s);
        ASSERT_MSG(got == rows[i].complete, "complete(%s) = %d, want %d",
                   rows[i].s, got, rows[i].complete);
    }
}

static TermLexContext ctx(const char* s) {
    TermLexContext c;
    term_lex_context(s, strlen(s), &c);
    return c;
}

static void assert_head(const char* s, const TermLexContext* c, const char* head) {
    ASSERT_MSG(c->has_call_head, "%s: expected call head %s", s, head);
    ASSERT_MSG(c->head_len == strlen(head) && strncmp(s + c->head_start, head, c->head_len) == 0,
               "%s: head is %.*s, want %s", s, (int)c->head_len, s + c->head_start, head);
}

static void assert_word(const char* s, const TermLexContext* c, const char* word) {
    ASSERT_MSG(strcmp(s + c->word_start, word) == 0, "%s: word is %s, want %s",
               s, s + c->word_start, word);
}

static void test_context(void) {
    const char* s;
    TermLexContext c;

    s = "Plot[x, Pl"; c = ctx(s);
    assert_head(s, &c, "Plot"); ASSERT(c.arg_position); assert_word(s, &c, "Pl");
    ASSERT(!c.in_string && !c.in_comment);

    s = "a[[Pl"; c = ctx(s);
    ASSERT(!c.has_call_head); assert_word(s, &c, "Pl");

    s = "f[g[1], "; c = ctx(s);
    assert_head(s, &c, "f"); ASSERT(c.arg_position); assert_word(s, &c, "");

    s = "f[\"Ve"; c = ctx(s);
    ASSERT(c.in_string); ASSERT(c.string_start == 2); assert_word(s, &c, "Ve");
    assert_head(s, &c, "f"); ASSERT(c.arg_position);

    s = "(* f[ *) Int"; c = ctx(s);
    ASSERT(!c.has_call_head); ASSERT(!c.in_comment); assert_word(s, &c, "Int");

    s = "f[x -> Aut"; c = ctx(s);
    assert_head(s, &c, "f"); ASSERT(!c.arg_position); assert_word(s, &c, "Aut");

    /* Inside a comment. */
    s = "(* Integr"; c = ctx(s);
    ASSERT(c.in_comment); ASSERT(!c.in_string);

    /* Whitespace between the head and '['. */
    s = "Plot  [x, Pl"; c = ctx(s);
    assert_head(s, &c, "Plot"); ASSERT(c.arg_position);

    /* A Part closes; the enclosing call is the context again. */
    s = "f[a[[1]], Ma"; c = ctx(s);
    assert_head(s, &c, "f"); ASSERT(c.arg_position);

    /* The innermost call wins; a closed call does not. */
    s = "f[g[x], h[y, Op"; c = ctx(s);
    assert_head(s, &c, "h"); ASSERT(c.arg_position);
    s = "f[g[x"; c = ctx(s);
    assert_head(s, &c, "g");

    /* Inside a list inside a call: no call head at that depth. */
    s = "f[{Ma"; c = ctx(s);
    ASSERT(!c.has_call_head);

    /* A backticked word. */
    s = "System`Integr"; c = ctx(s);
    assert_word(s, &c, "System`Integr"); ASSERT(!c.has_call_head);

    /* A number before '[' is not a call head. */
    s = "2[Ma"; c = ctx(s);
    ASSERT(!c.has_call_head);

    /* The cursor in the middle of the buffer. */
    s = "Plot[x, Pl] + 1";
    term_lex_context(s, 10, &c);
    assert_head(s, &c, "Plot"); ASSERT(c.word_start == 8);

    /* Deeper than the frame stack: no head, no crash. */
    char deep[200];
    size_t k = 0;
    for (int i = 0; i < 70; i++) { deep[k++] = 'f'; deep[k++] = '['; }
    deep[k] = '\0';
    c = ctx(deep);
    ASSERT(!c.has_call_head);

    /* NULL and an empty buffer. */
    term_lex_context(NULL, 0, &c);
    ASSERT(!c.has_call_head && c.word_start == 0);
    c = ctx("");
    ASSERT(!c.has_call_head && c.word_start == 0);
}

int main(void) {
    TEST(test_completeness);
    TEST(test_context);
    printf("All term_lex tests passed.\n");
    return 0;
}
