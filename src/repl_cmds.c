/*
 * repl_cmds.c -- see repl_cmds.h.
 */
/* opendir/readdir are POSIX, hidden by glibc under -std=c99. */
#ifndef _POSIX_C_SOURCE
#define _POSIX_C_SOURCE 200809L
#endif

#include "repl_cmds.h"

#include "eval.h"
#include "expr.h"
#include "frontend.h"
#include "loadmodule.h"
#include "message.h"
#include "print.h"
#include "quit.h"
#include "sym_names.h"
#include "symtab.h"
#include "term/term_color.h"
#include "term/term_config.h"
#include "term/term_editor.h"
#include "term/term_editor_tty.h"
#include "term/term_history.h"
#include "term/term_interrupt.h"
#include "term/term_slash.h"
#include "term/term_status.h"
#include "term/term_theme.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifndef _WIN32
#include <dirent.h>
#endif

enum { CMD_HELP = 1, CMD_THEME, CMD_HISTORY, CMD_TIMING, CMD_STATUSLINE, CMD_CLEAR, CMD_MODS,
       CMD_QUIT };

static bool         g_on = false;
static TermReader*  g_rd = NULL;
static TermHistory* g_hist = NULL;
static char*        g_dir = NULL;       /* the per-user directory, or NULL */
static int          g_timing = 1;       /* 0 off, 1 slow results, 2 always */

/* ---- output helpers --------------------------------------------------------- */

static bool colour(void) { return term_color_stdout(); }
static const char* dim(void)  { return term_style_on(TERM_STYLE_DIM, colour()); }
static const char* bold(void) { return colour() ? term_theme_sgr(TERM_ROLE_LABEL) : ""; }
static const char* off(void)  { return term_style_off(colour()); }

/* `path` with $HOME shown as "~". */
static void show_path(char* buf, size_t cap, const char* path) {
    const char* home = getenv("HOME");
    size_t hl = home ? strlen(home) : 0;
    if (hl > 1 && path && strncmp(path, home, hl) == 0 && (path[hl] == '/' || !path[hl]))
        snprintf(buf, cap, "~%s", path + hl);
    else
        snprintf(buf, cap, "%s", path ? path : "");
}

/* ---- mods --------------------------------------------------------------------- */

#define MAX_MODS 64
enum { MOD_OK, MOD_MESSAGES, MOD_ABORTED, MOD_QUIT };

typedef struct {
    char*         name;       /* "init.m", "mods/units.m" */
    char*         path;
    int           status;
    unsigned long messages;
    char          hooks[64];  /* "$Pre $Post" */
} ModInfo;

static ModInfo g_mods[MAX_MODS];
static int     g_nmods = 0;
static bool    g_mods_skipped = false;

static const char* const HOOKS[] = { "$PreRead", "$Pre", "$Post", "$PrePrint" };
#define N_HOOKS 4

static const void* hook_value(const char* name) {
    SymbolDef* d = symtab_lookup(name);
    return d && d->own_values ? (const void*)d->own_values->replacement : NULL;
}

static void load_one(const char* name, const char* path) {
    FILE* f = fopen(path, "r");
    if (!f) return;
    fclose(f);
    if (g_nmods >= MAX_MODS) return;
    ModInfo* m = &g_mods[g_nmods];
    memset(m, 0, sizeof *m);
    size_t nl = strlen(name), pl = strlen(path);
    m->name = malloc(nl + 1);
    m->path = malloc(pl + 1);
    if (!m->name || !m->path) { free(m->name); free(m->path); return; }
    memcpy(m->name, name, nl + 1);
    memcpy(m->path, path, pl + 1);
    g_nmods++;

    const void* before[N_HOOKS];
    for (int i = 0; i < N_HOOKS; i++) before[i] = hook_value(HOOKS[i]);
    unsigned long msgs0 = mth_msg_fired_count();
    frontend_set_mod_owner(name);
    TermState prev = term_interrupt_set_state(TERM_STATE_EVALUATING);
    mth_interrupt_clear();
    (void)eval_take_top_abort();
    int opened = 0;
    Expr* r = mathilda_run_file(path, &opened);
    bool aborted = eval_take_top_abort()
        || (r && r->type == EXPR_SYMBOL && r->data.symbol.name == SYM_DollarAborted);
    expr_free(r);
    mth_interrupt_clear();
    term_interrupt_set_state(prev);
    frontend_set_mod_owner(NULL);

    m->messages = mth_msg_fired_count() - msgs0;
    m->status = m->messages ? MOD_MESSAGES : MOD_OK;
    if (aborted) m->status = MOD_ABORTED;
    if (mth_quit_requested()) { mth_quit_cancel(); m->status = MOD_QUIT; }
    size_t o = 0;
    for (int i = 0; i < N_HOOKS; i++) {
        if (hook_value(HOOKS[i]) == before[i]) continue;
        int k = snprintf(m->hooks + o, sizeof m->hooks - o, "%s%s", o ? " " : "", HOOKS[i]);
        if (k > 0 && (size_t)k < sizeof m->hooks - o) o += (size_t)k;
    }
    if (m->status != MOD_OK) {
        static const char* const WHY[] = { "", "reported messages", "was aborted", "called Quit[], ignored" };
        fprintf(stderr, "%s%s %s; the session continues (/mods).%s\n", dim(), name, WHY[m->status], off());
    }
}

#ifndef _WIN32
/* qsort comparator for the mods/ directory listing (POSIX branch only). */
static int cmp_str(const void* a, const void* b) {
    return strcmp(*(const char* const*)a, *(const char* const*)b);
}
#endif

void repl_cmds_load_mods(void) {
    if (!g_on || !g_dir) return;
    const char* e = getenv("MATHILDA_MODS");
    if (e && strcmp(e, "0") == 0) { g_mods_skipped = true; return; }
    char* p = term_config_join(g_dir, "init.m");
    if (p) load_one("init.m", p);
    free(p);
#ifndef _WIN32
    char* md = term_config_join(g_dir, "mods");
    DIR* d = md ? opendir(md) : NULL;
    if (d) {
        char* names[MAX_MODS];
        int n = 0;
        struct dirent* de;
        while ((de = readdir(d)) != NULL && n < MAX_MODS) {
            size_t l = strlen(de->d_name);
            if (l < 3 || de->d_name[0] == '.' || strcmp(de->d_name + l - 2, ".m") != 0) continue;
            names[n] = malloc(l + 1);
            if (names[n]) { memcpy(names[n], de->d_name, l + 1); n++; }
        }
        closedir(d);
        qsort(names, (size_t)n, sizeof names[0], cmp_str);
        for (int i = 0; i < n; i++) {
            char label[300];
            snprintf(label, sizeof label, "mods/%s", names[i]);
            char* path = term_config_join(md, names[i]);
            if (path) load_one(label, path);
            free(path);
            free(names[i]);
        }
    }
    free(md);
#endif
}

/* ---- the commands --------------------------------------------------------------- */

static void cmd_help(void) {
    printf("%sCommands%s\n", bold(), off());
    size_t n = term_slash_count();
    for (size_t i = 0; i < n; i++) {
        const TermSlashCmd* c = term_slash_at(i);
        printf("  /%-12s %s%s", c->name, dim(), c->summary);
        if (c->owner) printf("  (%s)", c->owner);
        printf("%s\n", off());
    }
    printf("%sKeys%s\n", bold(), off());
    if (term_reader_kind(g_rd) == TERM_BACKEND_EDITOR) {
        printf("  %sReturn evaluates a complete input, Esc Return forces it; Up/Down move, then recall.%s\n", dim(), off());
        printf("  %sTab completes, Ctrl-R searches the history, Ctrl-_ undoes, Alt-_ redoes.%s\n", dim(), off());
        printf("  %sCtrl-C clears the input or interrupts (twice quits), Ctrl-Z suspends, Ctrl-D exits.%s\n", dim(), off());
        printf("  %sType / at an empty prompt for this menu.%s\n", dim(), off());
    } else {
        printf("  %sReturn evaluates a complete input; Tab completes; Ctrl-C interrupts (twice quits); Ctrl-D exits.%s\n",
               dim(), off());
    }
    if (g_dir) {
        char sp[512];
        show_path(sp, sizeof sp, g_dir);
        printf("%sFiles%s\n  %s%s/config (theme), %s/init.m and %s/mods/*.m (mods)%s\n",
               bold(), off(), dim(), sp, sp, sp, off());
    }
}

static void cmd_mods(void) {
    char sp[512];
    if (!g_dir) {
        printf("Mods are off: MATHILDA_CONFIG_DIR is empty (or HOME is unset).\n");
        return;
    }
    show_path(sp, sizeof sp, g_dir);
    if (g_mods_skipped) printf("%sMATHILDA_MODS=0: no mods were loaded.%s\n", dim(), off());
    if (g_nmods == 0) {
        if (!g_mods_skipped)
            printf("No mods loaded. Mathilda code in %s/init.m or %s/mods/NAME.m runs at start-up;\n"
                   "SlashCommand[\"name\", f, \"summary\"] there adds a /name command.\n", sp, sp);
    } else {
        printf("%sMods%s %s(%s)%s\n", bold(), off(), dim(), sp, off());
        for (int i = 0; i < g_nmods; i++) {
            const ModInfo* m = &g_mods[i];
            static const char* const ST[] = { "loaded", "messages", "aborted", "Quit[] ignored" };
            printf("  %-18s %s", m->name, m->status == MOD_OK ? "" : term_style_on(TERM_STYLE_MSG_HEAD, colour()));
            if (m->status == MOD_MESSAGES) printf("%lu message%s", m->messages, m->messages == 1 ? "" : "s");
            else printf("%s", ST[m->status]);
            printf("%s", m->status == MOD_OK ? "" : off());
            /* the commands it defined */
            size_t n = term_slash_count();
            for (size_t k = 0; k < n; k++) {
                const TermSlashCmd* c = term_slash_at(k);
                if (c->owner && strcmp(c->owner, m->name) == 0) printf("%s \xC2\xB7 /%s%s", dim(), c->name, off());
            }
            if (m->hooks[0]) printf("%s \xC2\xB7 sets %s%s", dim(), m->hooks, off());
            printf("\n");
        }
    }
    /* commands defined at the prompt (no owner) */
    size_t n = term_slash_count();
    bool head = false;
    for (size_t k = 0; k < n; k++) {
        const TermSlashCmd* c = term_slash_at(k);
        if (c->kind != TERM_SLASH_MOD || c->owner) continue;
        if (!head) { printf("%sDefined this session%s\n", bold(), off()); head = true; }
        printf("  /%-12s %s%s%s\n", c->name, dim(), c->summary, off());
    }
}

static void cmd_history(const char* args) {
    size_t want = 20;
    if (*args) {
        long v = strtol(args, NULL, 10);
        if (v > 0) want = (size_t)v;
    }
    size_t n = g_hist ? term_history_count(g_hist) : 0;
    if (n == 0) { printf("No history yet.\n"); return; }
    size_t from = n > want ? n - want : 0;
    for (size_t i = from; i < n; i++) {
        const char* e = term_history_get(g_hist, i);
        printf("%s%5zu%s  ", dim(), i + 1, off());
        for (const char* p = e; *p; p++) {
            putchar(*p);
            if (*p == '\n') fputs("       ", stdout);
        }
        putchar('\n');
    }
    if (term_reader_kind(g_rd) == TERM_BACKEND_EDITOR)
        printf("%sCtrl-R searches it; Up recalls.%s\n", dim(), off());
}

static bool on_off(const char* a, bool cur, bool* ok) {
    *ok = true;
    if (!*a) return !cur;
    if (strcmp(a, "on") == 0 || strcmp(a, "1") == 0) return true;
    if (strcmp(a, "off") == 0 || strcmp(a, "0") == 0) return false;
    *ok = false;
    return cur;
}

static void cmd_timing(const char* a) {
    if (strcmp(a, "always") == 0) g_timing = 2;
    else {
        bool ok;
        bool v = on_off(a, g_timing != 0, &ok);
        if (!ok) { printf("Usage: /timing [on|off|always]\n"); return; }
        g_timing = v ? 1 : 0;
    }
    static const char* const MSG[] = {
        "Timing notes off.",
        "Timing notes on: results that take 1 s or more show the time taken.",
        "Timing notes on for every result."
    };
    printf("%s\n", MSG[g_timing]);
}

static void cmd_statusline(const char* a) {
    bool ok;
    bool v = on_off(a, term_status_statusline(), &ok);
    if (!ok) { printf("Usage: /statusline [on|off]\n"); return; }
    term_status_set_statusline(v);
    fputs(v ? "Status line on: a summary line follows each result.\n" : "Status line off.\n", stdout);
}

/* ---- /theme ---------------------------------------------------------------------- */

static void theme_apply(size_t i) { term_theme_apply(i, term_theme_depth()); }

static void theme_move(size_t i, void* ctx) { (void)ctx; theme_apply(i); }

static void theme_saved(size_t i) {
    char* cfg = term_config_join(g_dir, "config");
    char sp[512];
    if (cfg && term_config_set(cfg, "theme", term_theme_name(i))) {
        show_path(sp, sizeof sp, cfg);
        printf("Theme %s%s%s %s(saved in %s)%s\n", bold(), term_theme_name(i), off(), dim(), sp, off());
    } else {
        printf("Theme %s%s%s %s(for this session: %s)%s\n", bold(), term_theme_name(i), off(), dim(),
               cfg ? "the config file could not be written" : "no config file", off());
    }
    free(cfg);
    if (getenv("MATHILDA_THEME") && *getenv("MATHILDA_THEME"))
        printf("%sMATHILDA_THEME is set and chooses the theme at start-up.%s\n", dim(), off());
    if (!colour()) printf("%sColour is off in this session (NO_COLOR, TERM or not a terminal).%s\n", dim(), off());
}

static void theme_list(void) {
    printf("%sThemes%s\n", bold(), off());
    for (size_t i = 0; i < term_theme_count(); i++)
        printf("  %s %-14s %s%s%s\n", i == term_theme_current() ? "*" : " ", term_theme_name(i), dim(),
               term_theme_summary(i), off());
    printf("%s/theme NAME chooses one.%s\n", dim(), off());
}

static void cmd_theme(const char* a) {
    if (*a) {
        int t = term_theme_find(a);
        if (t < 0) {
            printf("Unknown theme \"%s\".\n", a);
            theme_list();
            return;
        }
        theme_apply((size_t)t);
        theme_saved((size_t)t);
        return;
    }
    size_t n = term_theme_count();
    const char* items[16];
    char docs[16][96];
    const char* dp[16];
    size_t cur = term_theme_current();
    if (n > 16) n = 16;
    for (size_t i = 0; i < n; i++) {
        items[i] = term_theme_name(i);
        snprintf(docs[i], sizeof docs[i], "%s%s", term_theme_summary(i), i == cur ? " (current)" : "");
        dp[i] = docs[i];
    }
    bool uni = term_status_unicode();
    TermPreviewLine pv[5] = {
        { TERM_PV_INPUT,   "In[1]:= ", "f[x_] := Expand[(x + 1)^3] (* a rule *)" },
        { TERM_PV_OUTPUT,  "Out[1]= ", "{f[2], \"text\", 1.5 y^2}" },
        { TERM_PV_MESSAGE, NULL,       "Power::infy: Infinite expression 1/0 encountered." },
        { TERM_PV_STATUS,  uni ? "\xE2\xA0\xB9 Table 340/1000 \xE2\x96\x95\xE2\x96\x88\xE2\x96\x88\xE2\x96\x88\xE2\x96\x8D      \xE2\x96\x8F 34%"
                               : "| Table 340/1000 [###.......] 34%",
                           uni ? " \xC2\xB7 1.2s \xC2\xB7 cpu 1.2s \xC2\xB7 mem 96 MB" : " | 1.2s | cpu 1.2s | mem 96 MB" },
        { TERM_PV_NOTE,    NULL,       uni ? "took 3.21 s \xC2\xB7 cpu 3.18 s \xC2\xB7 peak 412 MB"
                                           : "took 3.21 s | cpu 3.18 s | peak 412 MB" },
    };
    TermPick p;
    memset(&p, 0, sizeof p);
    p.title = "Theme";
    p.items = items;
    p.docs = dp;
    p.count = n;
    p.selected = cur;
    p.preview = colour() ? pv : NULL;
    p.npreview = colour() ? 5 : 0;
    p.on_move = theme_move;
    size_t pick = cur;
    TermReadStatus st = term_reader_pick(g_rd, &p, &pick);
    if (st == TERM_READ_OK) {
        theme_apply(pick);
        theme_saved(pick);
    } else {
        theme_apply(cur);                      /* Esc: back to what it was */
        if (st == TERM_READ_EOF) theme_list(); /* no picker on this reader */
    }
}

/* ---- running a line -------------------------------------------------------------- */

static void run_mod(const TermSlashCmd* c, const char* args, size_t alen) {
    static char label[40];
    snprintf(label, sizeof label, "/%s", c->name);
    Expr* arg = expr_new_string_n(args, alen);
    Expr* call = expr_new_function(expr_copy((Expr*)c->payload), &arg, 1);
    TermState prev = term_interrupt_set_state(TERM_STATE_EVALUATING);
    mth_interrupt_clear();
    (void)eval_take_top_abort();
    term_status_begin(label);
    Expr* r = evaluate(call);
    TermStatusStats st;
    term_status_end(&st);
    bool aborted = eval_take_top_abort()
        || (r && r->type == EXPR_SYMBOL && r->data.symbol.name == SYM_DollarAborted);
    mth_interrupt_clear();
    term_interrupt_set_state(prev);
    if (aborted) {
        printf("%s%s: aborted%s\n", dim(), label, off());
    } else if (r && r->type == EXPR_STRING) {
        printf("%s\n", r->data.string);             /* text, as Print shows it */
    } else if (r && !(r->type == EXPR_SYMBOL && r->data.symbol.name == SYM_Null)) {
        expr_print(r);
    }
    expr_free(r);
    expr_free(call);
    fflush(stdout);
}

void repl_cmds_run(const char* line) {
    const char* name;
    const char* args;
    size_t nl, al;
    if (!term_slash_parse(line, &name, &nl, &args, &al)) return;
    char a[256];
    size_t cl = al < sizeof a - 1 ? al : sizeof a - 1;
    memcpy(a, args, cl);
    a[cl] = '\0';
    const TermSlashCmd* c = nl ? term_slash_find(name, nl) : NULL;
    if (!c) {
        if (nl == 0) cmd_help();
        else printf("Unknown command /%.*s %s(/help lists them)%s\n", (int)nl, name, dim(), off());
        fflush(stdout);
        return;
    }
    if (c->kind == TERM_SLASH_MOD) { run_mod(c, args, al); return; }
    switch (c->id) {
    case CMD_HELP:       cmd_help(); break;
    case CMD_THEME:      cmd_theme(a); break;
    case CMD_HISTORY:    cmd_history(a); break;
    case CMD_TIMING:     cmd_timing(a); break;
    case CMD_STATUSLINE: cmd_statusline(a); break;
    case CMD_CLEAR:
        if (term_tty_stdout()) fputs("\033[H\033[2J\033[3J", stdout);
        break;
    case CMD_MODS:       cmd_mods(); break;
    case CMD_QUIT: {
        Expr* q = expr_new_function(expr_new_symbol(SYM_Quit), NULL, 0);
        Expr* r = evaluate(q);
        expr_free(r);
        expr_free(q);
        break;
    }
    default: break;
    }
    fflush(stdout);
}

bool repl_cmds_show_took(double seconds) {
    return g_timing == 2 || (g_timing == 1 && seconds >= 1.0);
}

/* ---- lifecycle -------------------------------------------------------------------- */

void repl_cmds_init(TermReader* rd, TermHistory* hist, const char* config_dir) {
    g_rd = rd;
    g_hist = hist;
    free(g_dir);
    g_dir = NULL;
    if (config_dir) {
        size_t n = strlen(config_dir);
        g_dir = malloc(n + 1);
        if (g_dir) memcpy(g_dir, config_dir, n + 1);
    }
    term_slash_clear();
    term_slash_add_builtin("help",       "list the commands and keys", CMD_HELP);
    term_slash_add_builtin("theme",      "choose a colour theme, with a live preview", CMD_THEME);
    term_slash_add_builtin("history",    "show recent inputs: /history [n]", CMD_HISTORY);
    term_slash_add_builtin("timing",     "the time-taken note: /timing [on|off|always]", CMD_TIMING);
    term_slash_add_builtin("statusline", "a status line after each result: /statusline [on|off]", CMD_STATUSLINE);
    term_slash_add_builtin("clear",      "clear the screen", CMD_CLEAR);
    term_slash_add_builtin("mods",       "the loaded mods, their commands and hooks", CMD_MODS);
    term_slash_add_builtin("quit",       "end the session (as Quit[])", CMD_QUIT);
    g_on = true;
}

bool repl_cmds_enabled(void) { return g_on; }

void repl_cmds_shutdown(void) {
    term_slash_clear();
    for (int i = 0; i < g_nmods; i++) { free(g_mods[i].name); free(g_mods[i].path); }
    g_nmods = 0;
    free(g_dir);
    g_dir = NULL;
    frontend_set_mod_owner(NULL);
    g_on = false;
}
