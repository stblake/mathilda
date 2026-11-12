/*
 * test_term_spinner.c -- layout of the REPL's live status line, the took note
 * and the status line (src/term/term_spinner.c): formatters, width fitting
 * and segment dropping, colour and ASCII fallbacks, the fraction and ETA.
 */
#include "test_utils.h"
#include "term/term_spinner.h"

#include <stdio.h>
#include <string.h>

static TermSpinnerSnap snap(void) {
    TermSpinnerSnap s;
    memset(&s, 0, sizeof s);
    s.label = "Integrate";
    s.live = "Together";
    s.elapsed = 3.21;
    s.cpu = 3.18;
    s.frac = -1;
    s.eta = -1;
    s.mem = 412ull << 20;
    s.peak = 1200ull << 20;
    s.budget = 16ull << 30;
    s.workers = 4;
    return s;
}

static bool has_esc(const char* s) { return strchr(s, '\033') != NULL; }
static bool all_ascii(const char* s) {
    for (; *s; s++) if ((unsigned char)*s >= 0x80) return false;
    return true;
}

static void test_formatters(void) {
    char b[64];
    term_fmt_bytes(b, sizeof b, 512ull << 10);   ASSERT_STR_EQ(b, "512 KB");
    term_fmt_bytes(b, sizeof b, 412ull << 20);   ASSERT_STR_EQ(b, "412 MB");
    term_fmt_bytes(b, sizeof b, 1536ull << 20);  ASSERT_STR_EQ(b, "1.50 GB");
    term_fmt_bytes(b, sizeof b, 24ull << 30);    ASSERT_STR_EQ(b, "24.0 GB");
    term_fmt_secs(b, sizeof b, 3.21);            ASSERT_STR_EQ(b, "3.2s");
    term_fmt_secs(b, sizeof b, 65);              ASSERT_STR_EQ(b, "1m05s");
    term_fmt_secs(b, sizeof b, 3725);            ASSERT_STR_EQ(b, "1h02m");
    term_fmt_secs(b, sizeof b, -1);              ASSERT_STR_EQ(b, "0.0s");
    term_fmt_count(b, sizeof b, 3400);           ASSERT_STR_EQ(b, "3400");
    term_fmt_count(b, sizeof b, 1250000);        ASSERT_STR_EQ(b, "1.25M");
    term_fmt_count(b, sizeof b, 22400000);       ASSERT_STR_EQ(b, "22.4M");
    term_fmt_count(b, sizeof b, 1000000000L);    ASSERT_STR_EQ(b, "1G");
    term_fmt_bar(b, sizeof b, 0.5, 10, false);   ASSERT_STR_EQ(b, "[#####.....]");
    term_fmt_bar(b, sizeof b, 1.5, 4, false);    ASSERT_STR_EQ(b, "[####]");
    term_fmt_bar(b, sizeof b, -1, 4, false);     ASSERT_STR_EQ(b, "[....]");
    term_fmt_bar(b, sizeof b, 0.34, 10, true);
    ASSERT(term_display_width(b) == 12);         /* edges + 10 cells */
    ASSERT_STR_EQ(b, "\xe2\x96\x95\xe2\x96\x88\xe2\x96\x88\xe2\x96\x88\xe2\x96\x8d      \xe2\x96\x8f");
    ASSERT(term_display_width("\033[1;31mab\xc2\xb7\033[0m") == 3);
}

static void test_render_full(void) {
    TermSpinnerSnap s = snap();
    char b[1024];
    TermSpinnerStyle plain = { false, true };
    int w = term_spinner_render(&s, 200, plain, b, sizeof b);
    ASSERT_STR_EQ(b, "\xe2\xa0\x8b Integrate \xc2\xb7 Together \xc2\xb7 3.2s \xc2\xb7 cpu 3.2s"
                     " \xc2\xb7 mem 412 MB (2.5%) \xc2\xb7 peak 1.17 GB \xc2\xb7 4 threads");
    ASSERT(w == term_display_width(b));
    ASSERT(!has_esc(b));
    /* ASCII */
    TermSpinnerStyle ascii = { false, false };
    term_spinner_render(&s, 200, ascii, b, sizeof b);
    ASSERT(all_ascii(b));
    ASSERT(strncmp(b, "| Integrate | Together | 3.2s", 29) == 0);
    s.frame = 1;
    term_spinner_render(&s, 200, ascii, b, sizeof b);
    ASSERT(b[0] == '/');
    /* colour */
    TermSpinnerStyle col = { true, true };
    w = term_spinner_render(&s, 200, col, b, sizeof b);
    ASSERT(has_esc(b));
    ASSERT(w == term_display_width(b));
    ASSERT(strstr(b, "\033[36m") && strstr(b, "\033[1mIntegrate"));
}

static void test_render_fit(void) {
    TermSpinnerSnap s = snap();
    char b[1024];
    TermSpinnerStyle st = { true, true };
    int full = term_spinner_render(&s, 400, st, b, sizeof b);
    int prev = full;
    for (int cols = 120; cols >= 1; cols--) {
        int w = term_spinner_render(&s, cols, st, b, sizeof b);
        ASSERT_MSG(w <= (cols > 1 ? cols - 1 : 0), "cols %d gave width %d", cols, w);
        ASSERT(w == term_display_width(b));
        ASSERT(w <= prev);
        prev = w;
    }
    /* drop order: the live head goes first, elapsed last */
    term_spinner_render(&s, full, st, b, sizeof b);   /* one short */
    ASSERT(!strstr(b, "Together") && strstr(b, "4 threads"));
    term_spinner_render(&s, 22, st, b, sizeof b);
    ASSERT(strstr(b, "Integrate") && strstr(b, "3.2s") && !strstr(b, "mem"));
    /* tiny: the label is cut, the glyph stays */
    TermSpinnerStyle pl = { false, true };
    term_spinner_render(&s, 6, pl, b, sizeof b);
    ASSERT_STR_EQ(b, "\xe2\xa0\x8b Int");      /* 5 cells */
    ASSERT(term_spinner_render(&s, 1, pl, b, sizeof b) == 0 && b[0] == '\0');
}

static void test_render_progress(void) {
    TermSpinnerSnap s = snap();
    s.label = "Table"; s.live = NULL; s.workers = 0; s.peak = 0;
    s.outer.head = "Table"; s.outer.done = 3400; s.outer.total = 10000;
    s.frac = term_spinner_fraction(&s.outer, &s.inner, 1);
    s.eta = 6.4;
    char b[1024];
    TermSpinnerStyle st = { false, true };
    term_spinner_render(&s, 200, st, b, sizeof b);
    /* the loop is the input: it takes the label's place */
    ASSERT_STR_EQ(b, "\xe2\xa0\x8b Table 3400/10000 \xe2\x96\x95\xe2\x96\x88\xe2\x96\x88\xe2\x96\x88"
                     "\xe2\x96\x8d      \xe2\x96\x8f 34% \xc2\xb7 eta 6.4s \xc2\xb7 3.2s \xc2\xb7 cpu 3.2s"
                     " \xc2\xb7 mem 412 MB (2.5%)");
    /* a loop inside another input, nested */
    s.label = "Integrate";
    s.inner.head = "Do"; s.inner.done = 5; s.inner.total = 10;
    s.frac = term_spinner_fraction(&s.outer, &s.inner, 2);
    ASSERT(s.frac > 0.34004 && s.frac < 0.34006);
    TermSpinnerStyle as = { false, false };
    term_spinner_render(&s, 200, as, b, sizeof b);
    ASSERT_MSG(strstr(b, "| Integrate | Table 3400/10000 > Do 5/10 [###.......] 34% | eta 6.4s"), "%s", b);
    /* unknown total: a count and no bar */
    s.outer.total = 0; s.inner.head = NULL;
    s.frac = term_spinner_fraction(&s.outer, &s.inner, 1);
    ASSERT(s.frac < 0);
    term_spinner_render(&s, 200, as, b, sizeof b);
    ASSERT(strstr(b, "Table 3400 |") && !strchr(b, '['));
}

static void test_memory_tones(void) {
    TermSpinnerSnap s = snap();
    s.peak = 0;
    char b[1024];
    TermSpinnerStyle st = { true, true };
    s.mem = 5ull << 30;  s.budget = 16ull << 30;          /* 31 % */
    term_spinner_render(&s, 200, st, b, sizeof b);
    ASSERT(!strstr(b, "\033[33m") && !strstr(b, "\033[31m"));
    s.mem = 10ull << 30;                                   /* 62 % */
    term_spinner_render(&s, 200, st, b, sizeof b);
    ASSERT(strstr(b, "\033[33mmem 10.0 GB (62%)"));
    s.mem = 14ull << 30;                                   /* 88 % */
    term_spinner_render(&s, 200, st, b, sizeof b);
    ASSERT(strstr(b, "\033[31mmem 14.0 GB (88%)"));
    s.budget = 0;                                          /* unknown budget */
    term_spinner_render(&s, 200, st, b, sizeof b);
    ASSERT(strstr(b, "mem 14.0 GB\033[0m"));
    s.mem = 1ull << 20; s.budget = 64ull << 30;
    TermSpinnerStyle pl = { false, true };
    term_spinner_render(&s, 200, pl, b, sizeof b);
    ASSERT(strstr(b, "mem 1 MB (<0.1%)"));
}

static void test_eta(void) {
    TermEta e = { 0, 0, 0, false };
    ASSERT(term_eta_update(&e, 7, 0.10, 100.0) < 0);   /* first sight */
    ASSERT(term_eta_update(&e, 7, 0.15, 100.5) < 0);   /* < 1 s observed */
    double eta = term_eta_update(&e, 7, 0.30, 102.0);   /* 0.2 in 2 s: 7 s for 0.7 */
    ASSERT_MSG(eta > 6.99 && eta < 7.01, "eta %g", eta);
    ASSERT(term_eta_update(&e, 8, 0.50, 103.0) < 0);   /* new loop: reset */
    ASSERT(term_eta_update(&e, 8, -1, 104.0) < 0);
    TermSpinnerLoop o = { "Do", 1, 4 }, in = { "Do", 3, 0 };
    ASSERT(term_spinner_fraction(&o, &in, 2) == 0.25);  /* inner total unknown */
    TermSpinnerLoop none = { NULL, 0, 0 };
    ASSERT(term_spinner_fraction(&none, NULL, 0) < 0);
}

static void test_took_and_status(void) {
    char b[256];
    TermSpinnerStyle pl = { false, true }, col = { true, false };
    term_took_line(b, sizeof b, 3.214, 3.18, 412ull << 20, pl);
    ASSERT_STR_EQ(b, "took 3.21 s \xc2\xb7 cpu 3.18 s \xc2\xb7 peak 412 MB\n");
    term_took_line(b, sizeof b, 75, -1, 0, pl);
    ASSERT_STR_EQ(b, "took 1 min 15 s\n");
    term_took_line(b, sizeof b, 1.2, 0.5, 0, col);
    ASSERT_STR_EQ(b, "\033[2mtook 1.20 s | cpu 0.50 s\033[0m\n");
    term_status_line(b, sizeof b, 3.21, 1536ull << 20, "Global`", pl);
    ASSERT_STR_EQ(b, "last 3.21 s \xc2\xb7 mem 1.50 GB \xc2\xb7 Global`\n");
    term_status_line(b, sizeof b, 0.004, 0, NULL, col);
    ASSERT_STR_EQ(b, "\033[2mlast 0.00 s\033[0m\n");
}

int main(void) {
    TEST(test_formatters);
    TEST(test_render_full);
    TEST(test_render_fit);
    TEST(test_render_progress);
    TEST(test_memory_tones);
    TEST(test_eta);
    TEST(test_took_and_status);
    printf("All term_spinner tests passed.\n");
    return 0;
}
