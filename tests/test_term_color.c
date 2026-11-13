/*
 * test_term_color.c -- the REPL colour rules (src/term/term_color.c): the
 * decision table, prompt wrapping for readline, the elapsed-time format, and
 * the message head style hook in src/message.c.
 */
#ifndef _POSIX_C_SOURCE
#define _POSIX_C_SOURCE 200809L
#endif

#include "test_utils.h"
#include "symtab.h"
#include "core.h"
#include "message.h"
#include "term/term_color.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

static TermColorDecision dec(const char* no_color, const char* force, const char* term,
                             bool interactive, bool out_tty, bool err_tty) {
    TermColorEnv e = { no_color, force, term, interactive, out_tty, err_tty };
    return term_color_decide(&e);
}

#define EXPECT(d, o, e) ASSERT_MSG((d).out == (o) && (d).err == (e), \
    "got out=%d err=%d, want out=%d err=%d", (d).out, (d).err, (o), (e))

static void test_decision(void) {
    TermColorDecision d;
    /* 1. Rule rows, in order. */
    d = dec(NULL, NULL, "xterm", false, true, true);     EXPECT(d, false, false);
    d = dec(NULL, "1", "xterm", false, true, true);      EXPECT(d, false, false); /* mode beats force */
    d = dec("1", NULL, "xterm", true, true, true);       EXPECT(d, false, false);
    d = dec("1", "1", "xterm", true, true, true);        EXPECT(d, false, false); /* NO_COLOR beats force */
    d = dec("", NULL, "xterm", true, true, true);        EXPECT(d, true, true);   /* empty NO_COLOR */
    d = dec(NULL, "1", NULL, true, false, false);        EXPECT(d, true, true);   /* force */
    d = dec(NULL, "1", "dumb", true, false, false);      EXPECT(d, true, true);
    d = dec(NULL, "0", "xterm", true, false, false);     EXPECT(d, false, false); /* "0" does not force */
    d = dec(NULL, "", "xterm", true, false, true);       EXPECT(d, false, true);
    d = dec(NULL, NULL, NULL, true, true, true);         EXPECT(d, false, false); /* TERM unset */
    d = dec(NULL, NULL, "", true, true, true);           EXPECT(d, false, false);
    d = dec(NULL, NULL, "dumb", true, true, true);       EXPECT(d, false, false);
    /* Per stream. */
    d = dec(NULL, NULL, "xterm-256color", true, true, false); EXPECT(d, true, false);
    d = dec(NULL, NULL, "xterm-256color", true, false, true); EXPECT(d, false, true);
    d = dec(NULL, NULL, "xterm-256color", true, true, true);  EXPECT(d, true, true);
    ASSERT(term_color_decide(NULL).out == false);
    /* Nothing is on before term_color_init. */
    ASSERT(!term_color_stdout() && !term_color_stderr() && !term_tty_stdout());
    term_color_init(false);
    ASSERT(!term_color_stdout() && !term_color_stderr() && !term_tty_stdout());
}

static void test_wrap(void) {
    /* 2. Readline wrapping. */
    char buf[128];
    term_style_wrap(buf, sizeof(buf), TERM_STYLE_IN, "In[1]:= ", true, true);
    ASSERT_STR_EQ(buf, "\001\033[1;34m\002In[1]:= \001\033[0m\002");
    term_style_wrap(buf, sizeof(buf), TERM_STYLE_IN, "  In[3]:= ", true, true);
    ASSERT_STR_EQ(buf, "  \001\033[1;34m\002In[3]:= \001\033[0m\002");
    term_style_wrap(buf, sizeof(buf), TERM_STYLE_OUT, "Out", true, false);
    ASSERT_STR_EQ(buf, "\033[1;35mOut\033[0m");
    term_style_wrap(buf, sizeof(buf), TERM_STYLE_IN, "  In[3]:= ", false, true);
    ASSERT_STR_EQ(buf, "  In[3]:= ");
    /* snprintf semantics. */
    int need = term_style_wrap(buf, 4, TERM_STYLE_IN, "In[1]:= ", true, true);
    ASSERT(need > 4 && strlen(buf) == 3);
    ASSERT_STR_EQ(term_style_on(TERM_STYLE_DIM, true), "\033[2m");
    ASSERT_STR_EQ(term_style_on(TERM_STYLE_MSG_HEAD, true), "\033[1;31m");
    ASSERT_STR_EQ(term_style_on(TERM_STYLE_DIM, false), "");
    ASSERT_STR_EQ(term_style_off(false), "");
    ASSERT_STR_EQ(term_style_off(true), "\033[0m");
}

static void test_elapsed(void) {
    /* 3. The "took" format. */
    char buf[64];
    term_format_elapsed(buf, sizeof(buf), 1.234);  ASSERT_STR_EQ(buf, "1.23 s");
    term_format_elapsed(buf, sizeof(buf), 75.2);   ASSERT_STR_EQ(buf, "1 min 15 s");
    term_format_elapsed(buf, sizeof(buf), 3725);   ASSERT_STR_EQ(buf, "1 h 2 min");
    term_format_elapsed(buf, sizeof(buf), 59.994); ASSERT_STR_EQ(buf, "59.99 s");
}

static char* read_file(const char* path) {
    FILE* f = fopen(path, "rb");
    if (!f) return NULL;
    char* b = calloc(1, 4096);
    size_t n = fread(b, 1, 4095, f);
    b[n] = '\0';
    fclose(f);
    return b;
}

static void eval_discard(const char* src) {
    Expr* p = parse_expression(src);
    Expr* r = evaluate(p);
    expr_free(p);
    expr_free(r);
}

static void test_message_style(void) {
    /* 4. The message head style hook. */
    char path[64];
    snprintf(path, sizeof(path), "/tmp/mth_color_msg_%ld.txt", (long)getpid());
    mth_message_set_style("<", ">");
    int saved = test_stderr_redirect(path);
    ASSERT(saved >= 0);
    eval_discard("f::x = \"boom\"; Message[f::x]");
    test_stderr_restore(saved);
    char* text = read_file(path);
    ASSERT_MSG(text && strstr(text, "<f::x>: boom"), "stderr was %s", text ? text : "(null)");
    free(text);

    mth_message_set_style(NULL, NULL);
    saved = test_stderr_redirect(path);
    eval_discard("Message[f::x]");
    test_stderr_restore(saved);
    text = read_file(path);
    ASSERT_MSG(text && strstr(text, "f::x: boom") && !strchr(text, '<'), "stderr was %s", text);
    free(text);
    remove(path);
}

int main(void) {
    symtab_init();
    core_init();
    TEST(test_decision);
    TEST(test_wrap);
    TEST(test_elapsed);
    TEST(test_message_style);
    printf("All term_color tests passed.\n");
    return 0;
}
