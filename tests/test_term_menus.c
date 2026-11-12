/*
 * test_term_menus.c -- REPL PR 4: themes (src/term/term_theme.c), the config
 * file (term_config.c), the slash-command registry (term_slash.c), the
 * editor's slash menu, picker, redo, wrapped-row Up/Down and cursor-report
 * handling (term_editor.c, through the VT emulator), the status line's
 * Monitor segment (term_spinner.c), and the SlashCommand / Monitor builtins
 * (frontend.c).
 */
#include "test_utils.h"
#include "symtab.h"
#include "core.h"
#include "parse.h"
#include "eval.h"
#include "expr.h"
#include "print.h"
#include "sym_names.h"
#include "term/term_theme.h"
#include "term/term_config.h"
#include "term/term_slash.h"
#include "term/term_color.h"
#include "term/term_editor.h"
#include "term/term_highlight.h"
#include "term/term_spinner.h"
#include "term/term_progress.h"
#include "term_vt.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* ---- themes ------------------------------------------------------------------ */

static void test_theme_depth(void) {
    ASSERT(term_theme_depth_decide(NULL, NULL, "xterm") == TERM_DEPTH_16);
    ASSERT(term_theme_depth_decide(NULL, NULL, "xterm-256color") == TERM_DEPTH_256);
    ASSERT(term_theme_depth_decide(NULL, "truecolor", "xterm-256color") == TERM_DEPTH_24BIT);
    ASSERT(term_theme_depth_decide(NULL, "24bit", "screen") == TERM_DEPTH_24BIT);
    ASSERT(term_theme_depth_decide(NULL, NULL, "xterm-direct") == TERM_DEPTH_24BIT);
    ASSERT(term_theme_depth_decide(NULL, "truecolor", "dumb") == TERM_DEPTH_NONE);
    ASSERT(term_theme_depth_decide(NULL, NULL, NULL) == TERM_DEPTH_NONE);
    ASSERT(term_theme_depth_decide("16", "truecolor", "xterm") == TERM_DEPTH_16);
    ASSERT(term_theme_depth_decide("24", NULL, "xterm") == TERM_DEPTH_24BIT);
    ASSERT(term_theme_depth_decide("256", NULL, "xterm") == TERM_DEPTH_256);
    ASSERT(term_theme_depth_decide("0", NULL, "xterm") == TERM_DEPTH_NONE);
    ASSERT(term_theme_depth_decide("bogus", NULL, "xterm") == TERM_DEPTH_16);
    ASSERT(term_theme_light_background("0;15"));
    ASSERT(term_theme_light_background("0;default;7"));
    ASSERT(!term_theme_light_background("15;0"));
    ASSERT(!term_theme_light_background("default;default"));
    ASSERT(!term_theme_light_background(NULL));
}

/* At 16 colours the default theme is the REPL's original palette. */
static void test_theme_default_is_legacy(void) {
    static const char* const WANT[TERM_ROLE_COUNT] = {
        "\033[1;34m", "\033[1;35m", "\033[1;31m", "\033[2m",
        "\033[36m", "\033[1m", "\033[2m", "\033[36m", "\033[33m", "\033[31m",
        "\033[36m", "\033[34m", "\033[35m", "\033[32m", "\033[90m", "\033[33m", "",
        "\033[1;4m", "\033[1;31m", "\033[97;100m", "\033[1;97;44m", "\033[37;100m", "\033[2m",
    };
    size_t dark = (size_t)term_theme_find("dark"), ansi = (size_t)term_theme_find("ansi");
    ASSERT(dark == 0);
    char buf[64];
    for (int r = 0; r < TERM_ROLE_COUNT; r++) {
        term_theme_sgr_for(dark, (TermRole)r, TERM_DEPTH_16, buf, sizeof buf);
        ASSERT_MSG(strcmp(buf, WANT[r]) == 0, "dark role %d: %s", r, buf + 1);
        /* ansi is the same at every depth */
        term_theme_sgr_for(ansi, (TermRole)r, TERM_DEPTH_24BIT, buf, sizeof buf);
        ASSERT_MSG(strcmp(buf, WANT[r]) == 0, "ansi role %d: %s", r, buf + 1);
    }
    /* before any apply, the live table is the same */
    ASSERT(strcmp(term_style_on(TERM_STYLE_OUT, true), "\033[1;35m") == 0);
    ASSERT(strcmp(term_tone_on(TERM_TONE_WARN, true), "\033[33m") == 0);
    ASSERT(strcmp(term_hl_sgr(TERM_HL_POPUP_SEL), "\033[1;97;44m") == 0);
    ASSERT(strcmp(term_style_on(TERM_STYLE_OUT, false), "") == 0);
}

static void test_theme_depths_and_apply(void) {
    char buf[64];
    size_t dark = 0;
    term_theme_sgr_for(dark, TERM_ROLE_OUT, TERM_DEPTH_24BIT, buf, sizeof buf);
    ASSERT_MSG(strcmp(buf, "\033[1;38;2;198;120;221m") == 0, "%s", buf + 1);
    term_theme_sgr_for(dark, TERM_ROLE_POPUP_SEL, TERM_DEPTH_24BIT, buf, sizeof buf);
    ASSERT_MSG(strcmp(buf, "\033[1;38;2;255;255;255;48;2;61;111;180m") == 0, "%s", buf + 1);
    term_theme_sgr_for(dark, TERM_ROLE_OUT, TERM_DEPTH_256, buf, sizeof buf);
    ASSERT_MSG(strncmp(buf, "\033[1;38;5;", 9) == 0, "%s", buf + 1);
    term_theme_sgr_for(dark, TERM_ROLE_OUT, TERM_DEPTH_NONE, buf, sizeof buf);
    ASSERT(buf[0] == '\0');
    /* 256-colour mapping: cube corners and the grey ramp */
    ASSERT(term_theme_rgb_to_256(0x000000) == 16);
    ASSERT(term_theme_rgb_to_256(0xFFFFFF) == 231);
    ASSERT(term_theme_rgb_to_256(0xFF0000) == 196);
    ASSERT(term_theme_rgb_to_256(0x808080) == 244);
    /* mono: attributes only, at every depth */
    size_t mono = (size_t)term_theme_find("MONO");
    ASSERT(mono == 3);
    for (int r = 0; r < TERM_ROLE_COUNT; r++) {
        term_theme_sgr_for(mono, (TermRole)r, TERM_DEPTH_24BIT, buf, sizeof buf);
        ASSERT_MSG(strstr(buf, "38;") == NULL && strstr(buf, "48;") == NULL, "mono role %d", r);
    }
    ASSERT(term_theme_find("nope") == -1);
    ASSERT(term_theme_find("High-Contrast") == 2);
    /* apply changes the live strings in place (pointers stay valid) */
    const char* out = term_style_on(TERM_STYLE_OUT, true);
    term_theme_apply((size_t)term_theme_find("light"), TERM_DEPTH_24BIT);
    ASSERT_MSG(strcmp(out, "\033[1;38;2;130;80;223m") == 0, "%s", out + 1);
    ASSERT(term_theme_current() == 1 && term_theme_depth() == TERM_DEPTH_24BIT);
    term_theme_apply(0, TERM_DEPTH_16);
    ASSERT(strcmp(out, "\033[1;35m") == 0);
    /* every theme has a name and a summary */
    for (size_t i = 0; i < term_theme_count(); i++)
        ASSERT(term_theme_name(i) && term_theme_summary(i) && *term_theme_summary(i));
}

/* ---- config ------------------------------------------------------------------- */

static void test_config(void) {
    char* d = term_config_dir(NULL, "/home/u");
    ASSERT_STR_EQ(d, "/home/u/.mathilda");
    free(d);
    d = term_config_dir("/x/y/", "/home/u");
    ASSERT_STR_EQ(d, "/x/y/");
    char* j = term_config_join(d, "config");
    ASSERT_STR_EQ(j, "/x/y/config");
    free(j);
    free(d);
    ASSERT(term_config_dir("", "/home/u") == NULL);
    ASSERT(term_config_dir(NULL, NULL) == NULL);
    ASSERT(term_config_join(NULL, "config") == NULL);

    char* t = term_config_edit("# mine\ntheme = dark\nother=1\ntheme=light\n", "theme", "mono");
    ASSERT_STR_EQ(t, "# mine\ntheme = mono\nother=1\n");
    free(t);
    t = term_config_edit("other = 2", "theme", "ansi");
    ASSERT_STR_EQ(t, "other = 2\ntheme = ansi\n");
    free(t);
    t = term_config_edit(NULL, "theme", "dark");
    ASSERT_STR_EQ(t, "theme = dark\n");
    free(t);

    char dir[256], path[300], buf[64];
    snprintf(dir, sizeof dir, "/tmp/mth_cfg_test_%ld", (long)(size_t)&dir % 100000);
    snprintf(path, sizeof path, "%s/config", dir);
    remove(path);
    remove(dir);
    ASSERT(!term_config_get(path, "theme", buf, sizeof buf));
    ASSERT(term_config_set(path, "theme", "light"));      /* creates the directory */
    ASSERT(term_config_get(path, "theme", buf, sizeof buf));
    ASSERT_STR_EQ(buf, "light");
    ASSERT(term_config_set(path, "theme", "mono"));
    ASSERT(term_config_get(path, "theme", buf, sizeof buf));
    ASSERT_STR_EQ(buf, "mono");
    ASSERT(!term_config_get(path, "missing", buf, sizeof buf));
    FILE* f = fopen(path, "a");
    ASSERT(f != NULL);
    fputs("  pad   =  a b   # comment\n", f);
    fclose(f);
    ASSERT(term_config_get(path, "pad", buf, sizeof buf));
    ASSERT_STR_EQ(buf, "a b");
    remove(path);
    remove(dir);
}

/* ---- the registry ----------------------------------------------------------------- */

static int g_freed = 0;
static void count_free(void* p) { (void)p; g_freed++; }

static void register_builtins(void) {
    term_slash_clear();
    ASSERT(term_slash_add_builtin("help", "the help", 1));
    ASSERT(term_slash_add_builtin("theme", "themes", 2));
    ASSERT(term_slash_add_builtin("history", "past inputs", 3));
}

static void test_slash_registry(void) {
    ASSERT(term_slash_valid_name("a1", 2));
    ASSERT(!term_slash_valid_name("1a", 2));
    ASSERT(!term_slash_valid_name("a-b", 3));
    ASSERT(!term_slash_valid_name("", 0));
    register_builtins();
    ASSERT(!term_slash_add_builtin("help", "dup", 9));
    static int payload1, payload2;
    ASSERT(!term_slash_add_mod("help", "no", NULL, &payload1, count_free));   /* built-in */
    ASSERT(g_freed == 0);
    ASSERT(term_slash_add_mod("zeta", "z", "init.m", &payload1, count_free));
    ASSERT(term_slash_add_mod("alpha", "a", NULL, &payload2, count_free));
    ASSERT(term_slash_count() == 5);
    /* built-ins first (registration order), then mods by name */
    ASSERT_STR_EQ(term_slash_at(0)->name, "help");
    ASSERT_STR_EQ(term_slash_at(2)->name, "history");
    ASSERT_STR_EQ(term_slash_at(3)->name, "alpha");
    ASSERT_STR_EQ(term_slash_at(4)->name, "zeta");
    ASSERT_STR_EQ(term_slash_at(4)->owner, "init.m");
    /* replacing a mod command frees the old payload */
    ASSERT(term_slash_add_mod("zeta", "z2", NULL, &payload2, count_free));
    ASSERT(g_freed == 1 && term_slash_count() == 5);
    ASSERT_STR_EQ(term_slash_find("/zeta", 5)->summary, "z2");
    ASSERT_STR_EQ(term_slash_summary("/theme", 6), "themes");
    ASSERT(term_slash_summary("/nope", 5) == NULL);

    TermCompletions c;
    term_slash_complete("/", 1, &c);
    ASSERT(c.count == 5 && c.word_start == 0);
    ASSERT_STR_EQ(c.items[0], "/help");
    term_completions_free(&c);
    term_slash_complete("/H", 2, &c);
    ASSERT(c.count == 2);
    ASSERT_STR_EQ(c.items[1], "/history");
    term_completions_free(&c);
    term_slash_complete("/x", 2, &c);
    ASSERT(c.count == 0 && c.items == NULL);
    term_slash_complete("/he lp", 6, &c);
    ASSERT(c.count == 0);
    term_slash_complete("x/", 2, &c);
    ASSERT(c.count == 0);

    const char *nm, *ar;
    size_t nl, al;
    ASSERT(term_slash_parse("/theme  light \n", &nm, &nl, &ar, &al));
    ASSERT(nl == 5 && strncmp(nm, "theme", 5) == 0 && al == 5 && strncmp(ar, "light", 5) == 0);
    ASSERT(term_slash_parse("/", &nm, &nl, &ar, &al) && nl == 0 && al == 0);
    ASSERT(!term_slash_parse("x/y", &nm, &nl, &ar, &al));

    ASSERT(term_slash_remove("alpha"));
    ASSERT(!term_slash_remove("alpha"));
    ASSERT(!term_slash_remove("help"));                /* built-ins stay */
    ASSERT(g_freed == 2);
    term_slash_clear();
    ASSERT(g_freed == 3 && term_slash_count() == 0);
}

/* ---- the editor ----------------------------------------------------------------- */

static Vt          g_vt;
static TermEditor* g_ed;

static const char* slash_doc(const char* n, size_t len, void* ctx) {
    (void)ctx;
    return term_slash_summary(n, len);
}

static TermEditor* open_ed(int cols, int rows, bool slash) {
    vt_init(&g_vt, cols, rows);
    TermEditorConfig cfg;
    memset(&cfg, 0, sizeof(cfg));
    cfg.cols = cols;
    cfg.rows = rows;
    cfg.color = false;
    cfg.write = vt_write_cb;
    cfg.write_ctx = &g_vt;
    if (slash) {
        cfg.slash_complete = term_slash_complete;
        cfg.slash_doc = slash_doc;
    }
    TermEditor* ed = term_editor_new(&cfg);
    ASSERT(ed != NULL);
    return ed;
}

static TermEdStatus feed(const char* s) { return term_editor_feed(g_ed, s, strlen(s)); }
static const char* text(void) { return term_editor_text(g_ed, NULL); }

static int find_row(const char* needle) {
    for (int r = 0; r < g_vt.rows; r++)
        if (strstr(vt_row(&g_vt, r), needle)) return r;
    return -1;
}

static void test_slash_menu(void) {
    register_builtins();
    g_ed = open_ed(80, 20, true);
    term_editor_begin(g_ed, "In[1]:= ", 0, 0);
    ASSERT(feed("/") == TERM_ED_CONTINUE);
    ASSERT(find_row(">/help") == 1);                  /* menu under the prompt */
    ASSERT(find_row("the help") == 1);                /* summary beside the name */
    ASSERT(find_row("/history") == 3);
    feed("h");                                        /* filter: help, history */
    ASSERT(find_row("/theme") < 0);
    feed("\x1b[B");                                   /* Down: history */
    ASSERT(feed("\r") == TERM_ED_ACCEPT);             /* Return runs the selection */
    ASSERT_STR_EQ(text(), "/history");
    ASSERT(find_row(">/") < 0);                       /* the menu is gone */

    /* Tab completes the name and leaves room for an argument */
    term_editor_begin(g_ed, "In[1]:= ", -1, 0);
    feed("/th\t");
    ASSERT_STR_EQ(text(), "/theme ");
    feed("light");
    ASSERT(feed("\r") == TERM_ED_ACCEPT);             /* one line: always accepted */
    ASSERT_STR_EQ(text(), "/theme light");

    /* no match closes the menu; Backspace reopens it */
    term_editor_begin(g_ed, "In[1]:= ", -1, 0);
    feed("/zz");
    ASSERT(find_row(">/") < 0);
    feed("\x7f\x7f");
    ASSERT(find_row(">/help") >= 0);
    feed("\x1b");                                      /* Esc closes it */
    term_editor_timeout(g_ed);
    ASSERT(find_row(">/help") < 0);
    ASSERT_STR_EQ(text(), "/");
    ASSERT(feed("\r") == TERM_ED_ACCEPT);

    /* '/' anywhere but the start of an empty buffer is ordinary input */
    term_editor_begin(g_ed, "In[1]:= ", -1, 0);
    feed("x /. x -> 2");
    ASSERT(find_row(">/") < 0);
    feed("\x01");                                      /* Ctrl-A */
    feed("/");                                         /* "/x /. ..." not at the end */
    ASSERT(find_row(">/") < 0);
    feed("\x05\x15\x15");                              /* clear */
    feed("a //. b");
    ASSERT(find_row(">/") < 0);
    ASSERT(feed("\r") == TERM_ED_ACCEPT);
    ASSERT_STR_EQ(text(), "a //. b");
    /* an incomplete expression still continues on Return */
    term_editor_begin(g_ed, "In[1]:= ", -1, 0);
    feed("f[1/");
    ASSERT(feed("\r") == TERM_ED_CONTINUE);
    term_editor_free(g_ed);

    /* without a command source, '/' is just a character */
    g_ed = open_ed(80, 20, false);
    term_editor_begin(g_ed, "In[1]:= ", 0, 0);
    feed("/");
    ASSERT(find_row(">/") < 0);
    term_editor_free(g_ed);
    term_slash_clear();
}

static size_t g_moves[16];
static int    g_nmoves = 0;
static void on_move(size_t i, void* ctx) { (void)ctx; if (g_nmoves < 16) g_moves[g_nmoves++] = i; }

static void test_picker(void) {
    g_ed = open_ed(60, 20, false);
    const char* items[3] = { "dark", "light", "mono" };
    const char* docs[3] = { "for dark", "for light", "plain" };
    TermPreviewLine pv[2] = {
        { TERM_PV_INPUT, "In[1]:= ", "Sin[x] + 1" },
        { TERM_PV_MESSAGE, NULL, "Power::infy: Infinite expression." },
    };
    TermPick p;
    memset(&p, 0, sizeof p);
    p.title = "Theme";
    p.items = items;
    p.docs = docs;
    p.count = 3;
    p.selected = 1;
    p.preview = pv;
    p.npreview = 2;
    p.on_move = on_move;
    g_nmoves = 0;
    term_editor_pick_begin(g_ed, &p, 0);
    ASSERT(strncmp(vt_row(&g_vt, 0), "Theme", 5) == 0);
    ASSERT(find_row(">") == 2);                        /* light selected, plain marker */
    ASSERT(find_row("for light") == 2);
    ASSERT(find_row("In[1]:= Sin[x] + 1") == 5);
    ASSERT(find_row("Power::infy") == 6);
    feed("\x1b[B");
    ASSERT(g_nmoves == 1 && g_moves[0] == 2 && find_row(">") == 3);
    feed("\x1b[B");                                    /* wraps to the top */
    ASSERT(g_nmoves == 2 && g_moves[1] == 0);
    feed("3");                                         /* digits jump */
    ASSERT(g_moves[2] == 2);
    ASSERT(feed("\r") == TERM_ED_ACCEPT);
    ASSERT(term_editor_pick_index(g_ed) == 2);
    for (int r = 0; r < 8; r++) ASSERT_MSG(vt_row(&g_vt, r)[0] == '\0', "row %d left: %s", r, vt_row(&g_vt, r));
    ASSERT(g_vt.row == 0 && g_vt.col == 0);

    /* Esc cancels */
    term_editor_pick_begin(g_ed, &p, 0);
    feed("\x1b[A");
    feed("\x1b");
    ASSERT(term_editor_timeout(g_ed) == TERM_ED_INTERRUPT);
    /* a normal read afterwards works */
    term_editor_begin(g_ed, "> ", 0, 0);
    ASSERT(feed("1+1\r") == TERM_ED_ACCEPT);
    ASSERT_STR_EQ(text(), "1+1");
    term_editor_free(g_ed);
}

static void test_redo(void) {
    g_ed = open_ed(40, 10, false);
    term_editor_begin(g_ed, "> ", 0, 0);
    feed("abc");
    feed("\x1b[D\x1b[C");                              /* a move ends the undo group */
    feed(" def");
    feed("\x1f");                                      /* undo " def" */
    ASSERT_STR_EQ(text(), "abc");
    feed("\x1f");
    ASSERT_STR_EQ(text(), "");
    feed("\x1b_");                                     /* Alt-_: redo */
    ASSERT_STR_EQ(text(), "abc");
    feed("\x1b_");
    ASSERT_STR_EQ(text(), "abc def");
    ASSERT(term_editor_cursor(g_ed) == 7);
    feed("\x1b_");                                     /* nothing left */
    ASSERT_STR_EQ(text(), "abc def");
    feed("\x1f");
    ASSERT_STR_EQ(text(), "abc");
    feed("X");                                         /* a new edit ends the chain */
    feed("\x1b_");
    ASSERT_STR_EQ(text(), "abcX");
    term_editor_free(g_ed);
}

static void test_wrapped_updown(void) {
    g_ed = open_ed(40, 10, false);
    term_editor_begin(g_ed, "> ", 0, 0);
    /* 80 characters: rows of 38 after the 2-column gutter */
    char line[128];
    memset(line, 'a', 80);
    line[80] = '\0';
    feed(line);
    ASSERT(term_editor_cursor(g_ed) == 80);            /* row 2, screen column 6 */
    feed("\x1b[A");                                    /* up one screen row */
    ASSERT_MSG(term_editor_cursor(g_ed) == 42, "%zu", term_editor_cursor(g_ed));
    feed("\x1b[A");
    ASSERT(term_editor_cursor(g_ed) == 4);
    feed("\x1b[A");                                    /* the top row: history (none) */
    ASSERT(term_editor_cursor(g_ed) == 4);
    feed("\x1b[B");
    ASSERT(term_editor_cursor(g_ed) == 42);
    feed("\x1b[B");
    ASSERT(term_editor_cursor(g_ed) == 80);
    /* across lines: from the start of line 2 up to the last row of line 1 */
    term_editor_begin(g_ed, "> ", 0, 0);
    memset(line, 'b', 52);
    memcpy(line, "f[", 2);
    line[52] = ',';
    line[53] = '\0';
    feed(line);
    feed("\r1]");
    ASSERT(strchr(text(), '\n') != NULL);
    feed("\x01");                                      /* Ctrl-A: line 2, screen col 2 */
    feed("\x1b[A");
    ASSERT_MSG(term_editor_cursor(g_ed) == 38, "%zu", term_editor_cursor(g_ed));
    feed("\x1b[A");                                    /* the first row of line 1 */
    ASSERT_MSG(term_editor_cursor(g_ed) == 0, "%zu", term_editor_cursor(g_ed));
    term_editor_free(g_ed);
}

static void test_cpr(void) {
    g_ed = open_ed(40, 10, false);
    term_editor_begin(g_ed, "> ", -1, 0);
    term_editor_expect_cpr(g_ed);
    ASSERT(term_editor_cpr_pending(g_ed) && !term_editor_cpr_seen(g_ed));
    feed("ab\x1b[5;1Rc");                              /* the reply among the keys */
    ASSERT_STR_EQ(text(), "abc");
    ASSERT(term_editor_cpr_seen(g_ed) && !term_editor_cpr_pending(g_ed));
    feed("\x1b[1;2R");                                 /* unexpected: ignored (Shift-F3) */
    ASSERT_STR_EQ(text(), "abc");
    term_editor_free(g_ed);
}

/* ---- the Monitor segment ------------------------------------------------------- */

static void test_spinner_monitor(void) {
    TermSpinnerSnap s;
    memset(&s, 0, sizeof s);
    s.label = "Do";
    s.elapsed = 1.5;
    s.cpu = -1;
    s.frac = -1;
    s.eta = -1;
    s.monitor = "i = 42";
    TermSpinnerStyle st = { false, false };
    char buf[512];
    term_spinner_render(&s, 80, st, buf, sizeof buf);
    ASSERT_MSG(strstr(buf, "Do | i = 42 | 1.5s") != NULL, "%s", buf);
    /* a narrow line keeps the elapsed time longest */
    int w = term_spinner_render(&s, 16, st, buf, sizeof buf);
    ASSERT_MSG(strstr(buf, "1.5s") != NULL && strstr(buf, "i = 42") == NULL && w <= 15, "%s", buf);
    w = term_spinner_render(&s, 24, st, buf, sizeof buf);
    ASSERT_MSG(strstr(buf, "i = 42") != NULL && w <= 23, "%s", buf);
    /* a long monitor text is cut to half the line */
    char longt[200];
    memset(longt, 'x', sizeof longt - 1);
    longt[sizeof longt - 1] = '\0';
    s.monitor = longt;
    w = term_spinner_render(&s, 80, st, buf, sizeof buf);
    ASSERT_MSG(w <= 79 && strstr(buf, "1.5s") != NULL, "%s", buf);
}

/* ---- the builtins -------------------------------------------------------------- */

static char g_grabbed[MTH_MONITOR_TEXT];
static Expr* builtin_set_due(Expr* res) { (void)res; g_mth_monitor_due = 1; return expr_new_symbol(SYM_Null); }
static Expr* builtin_grab(Expr* res) {
    (void)res;
    mth_monitor_text(g_grabbed, sizeof g_grabbed);
    return expr_new_symbol(SYM_Null);
}

static char* eval_str(const char* src) {
    Expr* p = parse_expression(src);
    ASSERT(p != NULL);
    Expr* r = evaluate(p);
    char* s = expr_to_string_display(r);
    expr_free(r);
    expr_free(p);
    return s;
}

#define EVAL_IS(src, want) do { char* s_ = eval_str(src); \
    ASSERT_MSG(strcmp(s_, (want)) == 0, "%s gave %s, want %s", (src), s_, (want)); free(s_); } while (0)

static void test_builtins(void) {
    register_builtins();
    EVAL_IS("SlashCommand[\"greet\", (\"hi \" <> #)&, \"say hi\"]", "\"/greet\"");
    EVAL_IS("SlashCommand[\"/greet\"]", "\"say hi\"");
    EVAL_IS("SlashCommand[\"nope\"]", "None");
    EVAL_IS("SlashCommand[\"help\", f]", "$Failed");           /* built-in */
    EVAL_IS("SlashCommand[\"bad name\", f]", "$Failed");
    EVAL_IS("SlashCommand[1, f]", "SlashCommand[1, f]");
    EVAL_IS("Length[SlashCommand[]]", "4");
    EVAL_IS("SlashCommand[][[4]]", "\"/greet\" -> \"say hi\"");
    EVAL_IS("SlashCommand[\"greet\", None]", "Null");
    EVAL_IS("Length[SlashCommand[]]", "3");
    EVAL_IS("MemberQ[Attributes[SlashCommand], Protected]", "True");
    /* Monitor: HoldAll, evaluates expr, returns its value */
    EVAL_IS("Monitor[Table[k^2, {k, 4}], k]", "{1, 4, 9, 16}");
    EVAL_IS("Attributes[Monitor]", "{HoldAll, Protected}");
    EVAL_IS("Monitor[1 + 1, Print[\"never\"]]", "2");
    /* a due refresh (the status thread's request, here a probe) evaluates
     * mon on this thread at the next loop step and publishes its text */
    symtab_add_builtin("mthSetDue", builtin_set_due);
    symtab_add_builtin("mthGrab", builtin_grab);
    EVAL_IS("Monitor[(mx = 7; mthSetDue[]; Do[Null, {2}]; mthGrab[]), {mx, \"of\", 9}]", "Null");
    ASSERT_MSG(strcmp(g_grabbed, "{7, \"of\", 9}") == 0, "%s", g_grabbed);
    ASSERT(!mth_monitor_active());
    term_slash_clear();
}

int main(void) {
    symtab_init();
    core_init();
    TEST(test_theme_depth);
    TEST(test_theme_default_is_legacy);
    TEST(test_theme_depths_and_apply);
    TEST(test_config);
    TEST(test_slash_registry);
    TEST(test_slash_menu);
    TEST(test_picker);
    TEST(test_redo);
    TEST(test_wrapped_updown);
    TEST(test_cpr);
    TEST(test_spinner_monitor);
    TEST(test_builtins);
    printf("All term_menus tests passed.\n");
    return 0;
}
