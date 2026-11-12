/*
 * test_term_editor.c -- the in-tree line editor (src/term/term_editor.c and
 * term_render.c), driven through an in-memory VT screen (tests/term_vt.h):
 * what is drawn and where the cursor is, multi-line editing, history at the
 * buffer edges, smart and forced Return, emacs keys, kill/yank, undo,
 * bracketed paste, typeahead, wide characters, wrapping, bracket highlight,
 * the selection rules of term_editor_select, and one write per fed batch.
 */
#include "test_utils.h"
#include "symtab.h"
#include "core.h"
#include "term/term_editor.h"
#include "term/term_editor_tty.h"
#include "term/term_history.h"
#include "term_vt.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static Vt          g_vt;
static TermEditor* g_ed;

static TermEditor* open_ed(int cols, int rows, bool color, TermHistory* h) {
    vt_init(&g_vt, cols, rows);
    TermEditorConfig cfg;
    memset(&cfg, 0, sizeof(cfg));
    cfg.cols = cols;
    cfg.rows = rows;
    cfg.color = color;
    cfg.history = h;
    cfg.classify = term_hl_classify_symtab;
    cfg.write = vt_write_cb;
    cfg.write_ctx = &g_vt;
    TermEditor* ed = term_editor_new(&cfg);
    ASSERT(ed != NULL);
    return ed;
}

static TermEdStatus feed(const char* s) { return term_editor_feed(g_ed, s, strlen(s)); }

static const char* text(void) { return term_editor_text(g_ed, NULL); }

#define ROW(r, want) ASSERT_MSG(strcmp(vt_row(&g_vt, (r)), (want)) == 0, \
    "row %d is \"%s\", want \"%s\"", (r), vt_row(&g_vt, (r)), (want))

static void test_typing_and_multiline(void) {
    g_ed = open_ed(40, 10, false, NULL);
    term_editor_begin(g_ed, "In[1]:= ", 0, 0);
    ROW(0, "In[1]:=");
    ASSERT(g_vt.col == 8);
    ASSERT(feed("f[") == TERM_ED_CONTINUE);
    ROW(0, "In[1]:= f[");
    ASSERT(feed("\r") == TERM_ED_CONTINUE);          /* incomplete: a new line */
    ASSERT(strcmp(text(), "f[\n") == 0);
    ASSERT(g_vt.row == 1 && g_vt.col == 8);           /* indented under the prompt */
    ASSERT(feed("1]") == TERM_ED_CONTINUE);
    ROW(1, "        1]");
    ASSERT(feed("\r") == TERM_ED_ACCEPT);
    ASSERT(strcmp(text(), "f[\n1]") == 0);
    ASSERT(g_vt.row == 2 && g_vt.col == 0);           /* below the buffer */
    term_editor_free(g_ed);
}

static void test_force_return_and_typeahead(void) {
    g_ed = open_ed(40, 10, false, NULL);
    term_editor_begin(g_ed, "> ", 0, 0);
    ASSERT(feed("f[") == TERM_ED_CONTINUE);
    ASSERT(feed("\x1b\r") == TERM_ED_ACCEPT);         /* Esc Return */
    ASSERT(strcmp(text(), "f[") == 0);

    /* Typeahead: the second line waits for the next read. */
    term_editor_begin(g_ed, "> ", -1, 0);
    ASSERT(feed("1+1\r2+2\r") == TERM_ED_ACCEPT);
    ASSERT(strcmp(text(), "1+1") == 0);
    term_editor_begin(g_ed, "> ", -1, 0);
    ASSERT(term_editor_feed(g_ed, NULL, 0) == TERM_ED_ACCEPT);
    ASSERT(strcmp(text(), "2+2") == 0);

    /* A lone ESC waits for its sequence, then times out harmlessly. */
    term_editor_begin(g_ed, "> ", -1, 0);
    ASSERT(feed("ab\x1b") == TERM_ED_CONTINUE);
    ASSERT(term_editor_pending(g_ed));
    ASSERT(term_editor_timeout(g_ed) == TERM_ED_CONTINUE);
    ASSERT(!term_editor_pending(g_ed));
    ASSERT(strcmp(text(), "ab") == 0);
    /* A sequence split across feeds. */
    ASSERT(feed("\x1b[") == TERM_ED_CONTINUE);
    ASSERT(feed("D") == TERM_ED_CONTINUE);            /* Left */
    ASSERT(term_editor_cursor(g_ed) == 1);
    term_editor_free(g_ed);
}

static void test_vertical_and_history(void) {
    TermHistory* h = term_history_open(NULL, 100);
    ASSERT(h != NULL);
    ASSERT(term_history_add(h, "first"));
    ASSERT(term_history_add(h, "multi[\nline]"));
    g_ed = open_ed(40, 10, false, h);
    term_editor_begin(g_ed, "In[1]:= ", 0, 0);
    feed("draft");
    feed("\x1b[A");                                   /* Up at the top: history */
    ASSERT(strcmp(text(), "multi[\nline]") == 0);
    ASSERT(term_editor_cursor(g_ed) == 6);            /* end of its first line */
    ROW(0, "In[1]:= multi[");
    ROW(1, "        line]");
    feed("\x1b[B");                                   /* Down inside the entry */
    ASSERT(term_editor_cursor(g_ed) == 7 + 5);        /* goal column 6, line len 5 */
    feed("\x1b[A");                                   /* back up a line */
    ASSERT(term_editor_cursor(g_ed) == 6);
    feed("\x1b[A");                                   /* top edge: older entry */
    ASSERT(strcmp(text(), "first") == 0);
    ROW(0, "In[1]:= first");
    ROW(1, "");                                       /* the second row is cleared */
    feed("\x1b[A");                                   /* oldest: stays */
    ASSERT(strcmp(text(), "first") == 0);
    feed("\x1b[B");
    ASSERT(strcmp(text(), "multi[\nline]") == 0);
    ASSERT(term_editor_cursor(g_ed) == strlen("multi[\nline]"));   /* end of last line */
    feed("\x1b[B");                                   /* bottom edge: the draft */
    ASSERT(strcmp(text(), "draft") == 0);
    feed("\x1b<");                                    /* Alt-<: the oldest */
    ASSERT(strcmp(text(), "first") == 0);
    feed("\x1b>");                                    /* Alt->: back to the draft */
    ASSERT(strcmp(text(), "draft") == 0);
    term_editor_free(g_ed);
    term_history_close(h);
}

static void test_emacs_keys(void) {
    g_ed = open_ed(60, 10, false, NULL);
    term_editor_begin(g_ed, "> ", 0, 0);
    feed("hello big world");
    feed("\x17");                                     /* Ctrl-W */
    ASSERT(strcmp(text(), "hello big ") == 0);
    feed("\x19");                                     /* Ctrl-Y */
    ASSERT(strcmp(text(), "hello big world") == 0);
    feed("\x01");                                     /* Ctrl-A */
    ASSERT(term_editor_cursor(g_ed) == 0);
    feed("\x1b" "f");                                 /* Alt-F */
    ASSERT(term_editor_cursor(g_ed) == 5);
    feed("\x1b[1;5C");                                /* Ctrl-Right */
    ASSERT(term_editor_cursor(g_ed) == 9);
    feed("\x1b" "b");                                 /* Alt-B */
    ASSERT(term_editor_cursor(g_ed) == 6);
    feed("\x0b");                                     /* Ctrl-K */
    ASSERT(strcmp(text(), "hello ") == 0);
    feed("\x05");                                     /* Ctrl-E */
    feed("\x19");
    ASSERT(strcmp(text(), "hello big world") == 0);
    feed("\x15");                                     /* Ctrl-U: to line start */
    ASSERT(strcmp(text(), "") == 0);
    feed("ab");
    feed("\x14");                                     /* Ctrl-T at the end: swap */
    ASSERT(strcmp(text(), "ba") == 0);
    feed("\x02\x02");                                 /* Ctrl-B twice */
    feed("\x04");                                     /* Ctrl-D: delete under cursor */
    ASSERT(strcmp(text(), "a") == 0);
    feed("\x1b[3~");                                  /* Delete */
    ASSERT(strcmp(text(), "") == 0);
    ASSERT(feed("\x04") == TERM_ED_EOF);              /* Ctrl-D on empty: EOF */
    term_editor_free(g_ed);

    /* Ctrl-C clears a multi-line buffer and shows ^C. */
    g_ed = open_ed(60, 10, false, NULL);
    term_editor_begin(g_ed, "> ", 0, 0);
    feed("f[\rx");
    ASSERT(feed("\x03") == TERM_ED_INTERRUPT);
    ROW(1, "  x^C");
    ASSERT(g_vt.row == 2 && g_vt.col == 0);
    term_editor_free(g_ed);
}

static void test_undo(void) {
    g_ed = open_ed(60, 10, false, NULL);
    term_editor_begin(g_ed, "> ", 0, 0);
    feed("abc");
    feed(" def");
    feed("\x1f");                                     /* Ctrl-_: one typing run */
    ASSERT(strcmp(text(), "") == 0);
    feed("abc");
    feed("\x7f\x7f");                                 /* two Backspaces */
    ASSERT(strcmp(text(), "a") == 0);
    feed("\x1f");                                     /* undo both deletions */
    ASSERT(strcmp(text(), "abc") == 0);
    ASSERT(term_editor_cursor(g_ed) == 3);
    feed("\x01\x0b");                                 /* Ctrl-A Ctrl-K */
    ASSERT(strcmp(text(), "") == 0);
    feed("\x1f");
    ASSERT(strcmp(text(), "abc") == 0);
    feed("\x1f\x1f\x1f");                             /* past the start: harmless */
    ASSERT(strcmp(text(), "") == 0);
    term_editor_free(g_ed);
}

static void test_paste(void) {
    g_ed = open_ed(60, 10, false, NULL);
    term_editor_begin(g_ed, "> ", 0, 0);
    /* A pasted Return never submits; CR LF becomes LF; tabs become spaces. */
    ASSERT(feed("\x1b[200~a = 1;\r\nb = 2\r\tc\x1b[201~") == TERM_ED_CONTINUE);
    ASSERT_MSG(strcmp(text(), "a = 1;\nb = 2\n    c") == 0, "%s", text());
    ROW(2, "      c");
    /* The end marker split across feeds. */
    feed("\x1b[200~xy\x1b[20");
    feed("1~z");
    ASSERT_MSG(strcmp(text(), "a = 1;\nb = 2\n    cxyz") == 0, "%s", text());
    /* Up moves inside the pasted lines, keeping the column. */
    feed("\x1b[A");
    ASSERT(term_editor_cursor(g_ed) == strlen("a = 1;\nb = 2"));
    ASSERT(feed("\x1f") == TERM_ED_CONTINUE);         /* undo the "z" */
    term_editor_free(g_ed);
}

static void test_wide_and_wrap(void) {
    g_ed = open_ed(20, 10, false, NULL);
    term_editor_begin(g_ed, "> ", 0, 0);
    feed("\xE4\xB8\xAD" "a");                         /* 中a */
    ASSERT(g_vt.row == 0 && g_vt.col == 2 + 3);
    feed("\x1b[D");                                   /* Left over 'a' */
    ASSERT(g_vt.col == 4);
    feed("\x1b[D");                                   /* Left over the wide char */
    ASSERT(g_vt.col == 2);
    feed("\x05");
    /* 28 more columns: 15 fit after the prompt and "中a"; the rest wrap to a
     * row without a gutter (20 - 2 < 20 drops the continuation indent). */
    feed("0123456789012345678901234567");
    ROW(0, "> \xE4\xB8\xAD" "a012345678901234");
    ROW(1, "5678901234567");
    ASSERT(g_vt.row == 1 && g_vt.col == 13);
    feed("\x01");                                     /* Ctrl-A: back to row 0 */
    ASSERT(g_vt.row == 0 && g_vt.col == 2);
    term_editor_free(g_ed);
}

static void test_colour_and_brackets(void) {
    g_ed = open_ed(60, 10, true, NULL);
    term_editor_begin(g_ed, "In[1]:= ", 0, 0);
    ASSERT(strstr(g_vt.sgr_log, "[1;34m") != NULL);   /* the In style */
    g_vt.sgr_len = 0; g_vt.sgr_log[0] = '\0';
    feed("Sin[x)");
    ASSERT(strstr(g_vt.sgr_log, "[36m") != NULL);     /* Sin: builtin */
    ASSERT(strstr(g_vt.sgr_log, "[1;31m") != NULL);   /* the stray ')' */
    feed("\x7f]");
    g_vt.sgr_len = 0; g_vt.sgr_log[0] = '\0';
    feed("\x1b[D\x1b[C");                             /* redraw with the cursor after ']' */
    term_editor_free(g_ed);

    /* The match highlight: cursor after a ']' whose '[' is on the line. */
    g_ed = open_ed(60, 10, true, NULL);
    term_editor_begin(g_ed, "> ", 0, 0);
    g_vt.sgr_len = 0; g_vt.sgr_log[0] = '\0';
    feed("f[x]");
    ASSERT(strstr(g_vt.sgr_log, "[1;4m") != NULL);
    /* The accepted line is redrawn without the highlight. */
    g_vt.sgr_len = 0; g_vt.sgr_log[0] = '\0';
    ASSERT(feed("\r") == TERM_ED_ACCEPT);
    ASSERT(strstr(g_vt.sgr_log, "[1;4m") == NULL);
    term_editor_free(g_ed);

    /* No colour: no SGR at all. */
    g_ed = open_ed(60, 10, false, NULL);
    term_editor_begin(g_ed, "In[1]:= ", 0, 0);
    feed("Sin[x)] \"s\" (* c *) x_");
    ASSERT(g_vt.sgr_len == 0);
    term_editor_free(g_ed);
}

static void test_one_write_per_batch(void) {
    g_ed = open_ed(60, 10, false, NULL);
    term_editor_begin(g_ed, "> ", 0, 0);
    size_t w0 = term_editor_write_count(g_ed);
    feed("abc");
    ASSERT(term_editor_write_count(g_ed) == w0 + 1);
    feed("x");
    ASSERT(term_editor_write_count(g_ed) == w0 + 2);
    feed("\x1b[D");
    ASSERT(term_editor_write_count(g_ed) == w0 + 3);
    feed("\x0c");                                     /* Ctrl-L is one write too */
    ASSERT(term_editor_write_count(g_ed) == w0 + 4);
    ROW(0, "> abcx");
    ASSERT(feed("\r") == TERM_ED_ACCEPT);
    ASSERT(term_editor_write_count(g_ed) == w0 + 5);
    term_editor_free(g_ed);
}

static void test_resize_and_tall(void) {
    g_ed = open_ed(40, 6, false, NULL);
    term_editor_begin(g_ed, "> ", 5, 0);              /* on the last screen row */
    feed("\x1b[200~a1\na2\na3\na4\na5\na6\na7\na8\x1b[201~");
    /* Eight lines on a six-row screen: the viewport follows the cursor. */
    ROW(5, "  a8");
    ROW(0, "  a3");
    feed("\x1b[5~");                                  /* PageUp: buffer start */
    ROW(0, "> a1");
    term_editor_resize(g_ed, 30, 6);
    ROW(0, "> a1");
    ASSERT(feed("\x1b\r") == TERM_ED_ACCEPT);         /* the whole buffer is written out */
    ASSERT(g_vt.scrolled > 0);
    term_editor_free(g_ed);
}

static TermEditor* open_pop(int cols, int rows) {
    vt_init(&g_vt, cols, rows);
    TermEditorConfig cfg;
    memset(&cfg, 0, sizeof(cfg));
    cfg.cols = cols;
    cfg.rows = rows;
    cfg.classify = term_hl_classify_symtab;
    cfg.doc = term_hl_doc_symtab;
    cfg.write = vt_write_cb;
    cfg.write_ctx = &g_vt;
    TermEditor* ed = term_editor_new(&cfg);
    ASSERT(ed != NULL);
    return ed;
}

/* The screen row holding `needle`, or -1. */
static int find_row(const char* needle) {
    for (int r = 0; r < g_vt.rows; r++) if (strstr(vt_row(&g_vt, r), needle)) return r;
    return -1;
}

/* The row whose text, leading spaces dropped, is exactly `item`, or -1. */
static int find_item(const char* item) {
    for (int r = 0; r < g_vt.rows; r++) {
        const char* t = vt_row(&g_vt, r);
        while (*t == ' ') t++;
        if (strcmp(t, item) == 0) return r;
    }
    return -1;
}

static void test_popup(void) {
    g_ed = open_pop(80, 20);
    term_editor_begin(g_ed, "In[1]:= ", 0, 0);
    feed("Si");
    feed("n\t");                                     /* Sin, SinIntegral, Sinc, ... */
    ASSERT_MSG(strcmp(text(), "Sin") == 0, "%s", text());
    int r = find_item(">Sin");                        /* selected, marked without colour */
    ASSERT_MSG(r == 1, "popup row %d", r);
    ASSERT(find_row(" SinIntegral") == 2);
    /* The description row: Sin's docstring, first line. */
    const char* doc = term_hl_doc_symtab("Sin", 3, NULL);
    ASSERT(doc != NULL);
    char first[24];
    size_t fl = strcspn(doc, "\n");
    if (fl > 20) fl = 20;
    memcpy(first, doc, fl); first[fl] = '\0';
    ASSERT_MSG(find_row(first) > r, "no description row for \"%s\"", first);
    feed("\t");                                      /* Tab moves the selection */
    ASSERT(find_row(">SinIntegral") == 2);
    feed("\x1b[A");                                  /* Up moves it back */
    ASSERT(find_item(">Sin") == 1);
    feed("\x1b[B\r");                               /* Down, Enter accepts */
    ASSERT_MSG(strcmp(text(), "SinIntegral") == 0, "%s", text());
    ASSERT(find_row(">") == -1);                      /* popup gone */
    ASSERT(find_row("Sinc") == -1);
    feed("[");
    ASSERT(strcmp(text(), "SinIntegral[") == 0);      /* Enter did not submit */

    /* Typing filters; a non-word key closes without accepting. */
    feed("Sin\t");
    ASSERT(find_row(" Sinc") >= 1);
    feed("h");
    ASSERT(find_row(">Sinh") >= 1 && find_row(" SinhIntegral") >= 1 && find_row(" Sinc") == -1);
    feed(" ");
    ASSERT_MSG(strcmp(text(), "SinIntegral[Sinh ") == 0, "%s", text());
    ASSERT(find_row(">") == -1);
    /* Backspace past the queried word re-queries. */
    feed("Sinh\t");
    feed("\x7f");
    ASSERT(find_row(" Sinc") >= 1);
    feed("\x1b");
    term_editor_timeout(g_ed);
    term_editor_free(g_ed);

    /* Option names first inside a call's argument position. */
    g_ed = open_pop(80, 20);
    term_editor_begin(g_ed, "In[1]:= ", 0, 0);
    feed("NIntegrate[f, {x, 0, 1}, Work\t");
    ASSERT_MSG(strcmp(text(), "NIntegrate[f, {x, 0, 1}, WorkingPrecision -> ") == 0, "%s", text());
    feed("\x15");
    feed("NIntegrate[f, {x, 0, 1}, \t");             /* empty word: the options */
    ASSERT(find_row("AccuracyGoal ->") >= 1);
    feed("\x1b");
    ASSERT(term_editor_timeout(g_ed) == TERM_ED_CONTINUE);   /* Esc closes */
    ASSERT(find_row("AccuracyGoal") == -1);
    term_editor_free(g_ed);

    /* At the bottom of the screen with room above: the popup goes above. */
    g_ed = open_pop(80, 14);
    g_vt.row = 3;                                     /* the prompt is on row 3 */
    term_editor_begin(g_ed, "> ", 3, 0);
    feed("\x1b[200~a\nb\nc\nd\ne\nf\ng\nh\ni\n\x1b[201~");
    feed("Sin\t");
    int cur = g_vt.row;
    ASSERT(cur == 12);                                /* one row left below */
    r = find_item(">Sin");
    ASSERT_MSG(r >= 0 && r < cur, "popup row %d, cursor row %d", r, cur);
    term_editor_free(g_ed);

    /* Too small for a popup: Tab cycles inline. */
    g_ed = open_pop(14, 4);
    term_editor_begin(g_ed, "> ", 0, 0);
    feed("Sin\t");
    ASSERT_MSG(strcmp(text(), "Sin") == 0, "%s", text());
    feed("\t");
    ASSERT_MSG(strcmp(text(), "SinIntegral") == 0, "%s", text());
    feed("x");
    ASSERT(strcmp(text(), "SinIntegralx") == 0);
    term_editor_free(g_ed);
}

static void test_search(void) {
    TermHistory* h = term_history_open(NULL, 100);
    ASSERT(term_history_add(h, "a = 1"));
    ASSERT(term_history_add(h, "Plot[Sin[x],\n {x, 0, 1}]"));
    ASSERT(term_history_add(h, "b = 2"));
    ASSERT(term_history_add(h, "a = 1"));            /* not consecutive: kept */
    g_ed = open_ed(60, 20, false, h);
    term_editor_begin(g_ed, "In[5]:= ", 0, 0);
    feed("draft");
    feed("\x12");                                    /* Ctrl-R */
    ROW(0, "search:");
    ROW(1, "> a = 1");                                /* newest first, once */
    ROW(2, "  b = 2");
    ROW(3, "  Plot[Sin[x],");
    ROW(4, "     {x, 0, 1}]");                        /* multi-line: whole */
    ROW(5, "");
    feed("Sin");                                      /* filter */
    ROW(0, "search: Sin");
    ROW(1, "> Plot[Sin[x],");
    ROW(2, "     {x, 0, 1}]");
    ROW(3, "");
    ASSERT(g_vt.row == 0 && g_vt.col == 11);
    ASSERT(feed("\r") == TERM_ED_CONTINUE);          /* into the buffer, not run */
    ASSERT(strcmp(text(), "Plot[Sin[x],\n {x, 0, 1}]") == 0);
    ROW(0, "In[5]:= Plot[Sin[x],");
    ROW(1, "         {x, 0, 1}]");
    ROW(2, "");

    feed("\x12zzz");                                 /* no match, Esc restores */
    ROW(1, "  (no match)");
    feed("\x1b");
    term_editor_timeout(g_ed);
    ASSERT(strcmp(text(), "Plot[Sin[x],\n {x, 0, 1}]") == 0);

    feed("\x15\x15\x08");                           /* clear (two lines) */
    feed("\x01\x0b\x0b\x0b");
    ASSERT_MSG(strcmp(text(), "") == 0, "%s", text());
    feed("\x12\x12");                               /* Ctrl-R again: older */
    ROW(2, "> b = 2");
    feed("\x1b[D");                                  /* another key: take it, then act */
    ASSERT(strcmp(text(), "b = 2") == 0);
    ASSERT(term_editor_cursor(g_ed) == 4);
    ASSERT(feed("\r") == TERM_ED_ACCEPT);
    term_editor_free(g_ed);
    term_history_close(h);
}

static void test_select(void) {
    ASSERT(term_editor_select(NULL, "xterm-256color", true, true));
    ASSERT(term_editor_select("", "xterm", true, true));
    ASSERT(term_editor_select("builtin", "xterm", true, true));
    ASSERT(!term_editor_select("readline", "xterm", true, true));
    ASSERT(!term_editor_select("libedit", "xterm", true, true));
    ASSERT(!term_editor_select(NULL, "dumb", true, true));
    ASSERT(!term_editor_select(NULL, NULL, true, true));
    ASSERT(!term_editor_select(NULL, "", true, true));
    ASSERT(!term_editor_select(NULL, "xterm", false, true));
    ASSERT(!term_editor_select(NULL, "xterm", true, false));
}

int main(void) {
    symtab_init();
    core_init();
    TEST(test_typing_and_multiline);
    TEST(test_force_return_and_typeahead);
    TEST(test_vertical_and_history);
    TEST(test_emacs_keys);
    TEST(test_undo);
    TEST(test_paste);
    TEST(test_wide_and_wrap);
    TEST(test_colour_and_brackets);
    TEST(test_one_write_per_batch);
    TEST(test_resize_and_tall);
    TEST(test_popup);
    TEST(test_search);
    TEST(test_select);
    printf("All term_editor tests passed.\n");
    return 0;
}
