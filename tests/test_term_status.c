/*
 * test_term_status.c -- the REPL's status thread (src/term/term_status.c):
 * the enable rule, the locale and label helpers, and the forwarder itself,
 * driven through its test seam with a pipe standing in for the terminal:
 * pass-through bytes and their order, no output at all for a quick
 * evaluation, the delayed spinner, a partial line never overwritten,
 * suspend/resume, and descriptors restored afterwards.
 */
#ifndef _POSIX_C_SOURCE
#define _POSIX_C_SOURCE 200809L
#endif

#include "test_utils.h"
#include "symtab.h"
#include "core.h"
#include "parse.h"
#include "expr.h"
#include "term/term_status.h"
#include "term/term_progress.h"

#include <fcntl.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>

static void sleep_ms(long ms) {
    struct timespec ts = { ms / 1000, (ms % 1000) * 1000000L };
    nanosleep(&ts, NULL);
}

static void test_decide(void) {
    TermStatusEnv e = { true, true, true, "xterm-256color", NULL };
    ASSERT(term_status_decide(&e));
    e.spinner = "1";           ASSERT(term_status_decide(&e));
    e.spinner = "0";           ASSERT(!term_status_decide(&e));
    e.spinner = NULL;
    e.term = "dumb";           ASSERT(!term_status_decide(&e));
    e.term = "";               ASSERT(!term_status_decide(&e));
    e.term = NULL;             ASSERT(!term_status_decide(&e));
    e.term = "xterm";
    e.ttys = false;            ASSERT(!term_status_decide(&e));
    e.ttys = true; e.threads = false; ASSERT(!term_status_decide(&e));
    e.threads = true; e.interactive = false; ASSERT(!term_status_decide(&e));
    ASSERT(!term_status_decide(NULL));

    ASSERT(term_status_locale_utf8(NULL, NULL, "en_US.UTF-8"));
    ASSERT(term_status_locale_utf8("", "C.utf8", "C"));
    ASSERT(!term_status_locale_utf8("C", "en_US.UTF-8", "en_US.UTF-8"));   /* LC_ALL wins */
    ASSERT(!term_status_locale_utf8(NULL, NULL, NULL));
    ASSERT(!term_status_locale_utf8(NULL, NULL, "POSIX"));
}

static void test_label(void) {
    struct { const char* src; const char* want; } c[] = {
        { "Integrate[x, x]", "Integrate" },
        { "a = Table[i, {i, 3}]", "Table" },
        { "f[x_] := Expand[x]", "Expand" },
        { "a = 1; Do[i, {i, 3}]; b", "Do" },
        { "x", NULL },
        { "3", NULL },
    };
    for (size_t i = 0; i < sizeof c / sizeof c[0]; i++) {
        Expr* p = parse_expression(c[i].src);
        const char* got = term_status_label_of(p);
        ASSERT_MSG((got == NULL) == (c[i].want == NULL) && (!got || strcmp(got, c[i].want) == 0),
                   "%s: got %s", c[i].src, got ? got : "(null)");
        expr_free(p);
    }
    ASSERT(term_status_label_of(NULL) == NULL);
}

/* Everything readable from fd within ms milliseconds. */
static char* slurp(int fd, size_t* len) {
    size_t cap = 1 << 16, n = 0;
    char* b = malloc(cap + 1);
    for (;;) {
        long k = (long)read(fd, b + n, cap - n);
        if (k <= 0) break;
        n += (size_t)k;
        if (n == cap) { cap *= 2; b = realloc(b, cap + 1); }
    }
    b[n] = '\0';
    *len = n;
    return b;
}

static int g_tty[2];       /* the fake terminal: [0] the test reads, [1] the thread writes */
static int g_out, g_err;   /* the routed descriptors (files) */
static char g_out_path[64], g_err_path[64];

static void open_seam(void) {
    ASSERT(pipe(g_tty) == 0);
    int fl = fcntl(g_tty[0], F_GETFL);
    fcntl(g_tty[0], F_SETFL, fl | O_NONBLOCK);
    snprintf(g_out_path, sizeof g_out_path, "/tmp/mth_status_out_%ld", (long)getpid());
    snprintf(g_err_path, sizeof g_err_path, "/tmp/mth_status_err_%ld", (long)getpid());
    g_out = open(g_out_path, O_RDWR | O_CREAT | O_TRUNC, 0600);
    g_err = open(g_err_path, O_RDWR | O_CREAT | O_TRUNC, 0600);
    ASSERT(g_out >= 0 && g_err >= 0);
    TermStatusConfig c;
    memset(&c, 0, sizeof c);
    c.tty_fd = g_tty[1];
    c.route[0] = g_out; c.route[1] = g_err; c.n_route = 2;
    c.echo_fd = -1;
    c.color = false; c.unicode = true;
    c.delay_ms = 150; c.period_ms = 30; c.cols = 80;
    ASSERT(term_status_open(&c));
    ASSERT(term_status_enabled());
}

static void wr(int fd, const char* s) { ASSERT(write(fd, s, strlen(s)) == (long)strlen(s)); }

static void test_quick_is_silent_passthrough(void) {
    term_status_begin("Plus");
    wr(g_out, "one\n");
    wr(g_err, "two\n");
    wr(g_out, "three");
    TermStatusStats st;
    term_status_end(&st);
    size_t n;
    char* got = slurp(g_tty[0], &n);
    /* byte for byte, in order, and no spinner for a quick evaluation */
    ASSERT_MSG(strcmp(got, "one\ntwo\nthree") == 0, "got %s", got);
    free(got);
    ASSERT(st.elapsed >= 0 && st.elapsed < 0.15);
    /* the descriptors point at their files again */
    wr(g_out, "after");
    ASSERT(lseek(g_out, 0, SEEK_CUR) == 5);
    got = slurp(g_tty[0], &n);
    ASSERT(n == 0);
    free(got);
}

static void test_spinner_and_output(void) {
    term_status_begin("Table");
    int tok = mth_progress_push("Table", 10);
    mth_progress_step(tok, 3);
    sleep_ms(300);                       /* the spinner appears */
    wr(g_out, "printed\n");              /* erased, printed, redrawn */
    sleep_ms(120);
    wr(g_out, "partial");                /* never overwritten */
    sleep_ms(150);
    mth_progress_pop(tok);
    TermStatusStats st;
    term_status_end(&st);
    size_t n;
    char* got = slurp(g_tty[0], &n);
    ASSERT_MSG(strstr(got, "Table 3/10"), "no progress in %s", got);
    ASSERT(!strchr(got, '\033'));        /* colour off: no escapes at all */
    char* p = strstr(got, "printed\n");
    ASSERT(p);
    /* before "printed": a frame, then the erase (\r spaces \r) */
    ASSERT(p > got && p[-1] == '\r');
    /* after "partial": nothing drawn until the end */
    char* q = strstr(got, "partial");
    ASSERT(q && strcmp(q, "partial") == 0);
    /* something was drawn between the two outputs */
    ASSERT(strstr(p + 8, "Table 3/10") && strstr(p + 8, "Table 3/10") < q);
    free(got);
    ASSERT(st.elapsed >= 0.5);
    ASSERT(st.peak > 0);
}

static void test_suspend_resume(void) {
    term_status_begin("Dialog");
    wr(g_out, "a\n");
    term_status_suspend();
    wr(g_out, "direct");                 /* goes to the file, not the tty */
    term_status_resume();
    wr(g_out, "b\n");
    /* nested begin/end only count */
    term_status_begin("Inner");
    TermStatusStats in;
    term_status_end(&in);
    wr(g_out, "c\n");
    TermStatusStats st;
    term_status_end(&st);
    size_t n;
    char* got = slurp(g_tty[0], &n);
    ASSERT_MSG(strcmp(got, "a\nb\nc\n") == 0, "got %s", got);
    free(got);
    char buf[64] = { 0 };
    ASSERT(pread(g_out, buf, sizeof buf - 1, 5) == 6);
    ASSERT_STR_EQ(buf, "direct");
}

int main(void) {
    symtab_init();
    core_init();
    TEST(test_decide);
    TEST(test_label);
    open_seam();
    TEST(test_quick_is_silent_passthrough);
    TEST(test_spinner_and_output);
    TEST(test_suspend_resume);
    term_status_shutdown();
    ASSERT(!term_status_enabled());
    close(g_out); close(g_err); close(g_tty[0]);
    remove(g_out_path); remove(g_err_path);
    printf("All term_status tests passed.\n");
    return 0;
}
