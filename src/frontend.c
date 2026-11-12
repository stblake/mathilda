/*
 * frontend.c -- SlashCommand and Monitor; see frontend.h.
 */
#include "frontend.h"

#include "attr.h"
#include "core.h"
#include "eval.h"
#include "message.h"
#include "print.h"
#include "sym_names.h"
#include "symtab.h"
#include "term/term_progress.h"
#include "term/term_slash.h"

#include <stdlib.h>
#include <string.h>

/* ---- SlashCommand ---------------------------------------------------------- */

static char* g_owner = NULL;

void frontend_set_mod_owner(const char* name) {
    free(g_owner);
    g_owner = NULL;
    if (name) {
        size_t n = strlen(name);
        g_owner = malloc(n + 1);
        if (g_owner) memcpy(g_owner, name, n + 1);
    }
}

static void free_expr_payload(void* p) { expr_free((Expr*)p); }

static Expr* failed(void) { return expr_new_symbol(SYM_DollarFailed); }

/* SlashCommand[]: {"/help" -> "summary", ...} in menu order. */
static Expr* slash_list(void) {
    size_t n = term_slash_count();
    Expr** items = n ? malloc(n * sizeof(Expr*)) : NULL;
    if (n && !items) return NULL;
    for (size_t i = 0; i < n; i++) {
        const TermSlashCmd* c = term_slash_at(i);
        size_t ln = strlen(c->name);
        char* nm = malloc(ln + 2);
        if (!nm) { for (size_t k = 0; k < i; k++) expr_free(items[k]); free(items); return NULL; }
        nm[0] = '/';
        memcpy(nm + 1, c->name, ln + 1);
        Expr* kv[2] = { expr_new_string_adopt(nm), expr_new_string(c->summary) };
        items[i] = expr_new_function(expr_new_symbol(SYM_Rule), kv, 2);
    }
    Expr* r = expr_new_function(expr_new_symbol(SYM_List), items, n);
    free(items);
    return r;
}

Expr* builtin_slashcommand(Expr* res) {
    size_t argc = res->data.function.arg_count;
    Expr** a = res->data.function.args;
    if (argc == 0) return slash_list();
    if (argc > 3) {
        mth_message("SlashCommand", "argb",
                    "SlashCommand called with %zu arguments; between 0 and 3 arguments are expected.", argc);
        return NULL;
    }
    if (a[0]->type != EXPR_STRING) {
        mth_message("SlashCommand", "string", "String expected at position 1 in SlashCommand.");
        return NULL;
    }
    const char* name = a[0]->data.string;
    if (name[0] == '/') name++;
    size_t nl = strlen(name);
    if (!term_slash_valid_name(name, nl)) {
        mth_message("SlashCommand", "name",
                    "\"%s\" is not a command name: a letter followed by letters and digits.", name);
        return failed();
    }
    const TermSlashCmd* old = term_slash_find(name, nl);
    if (old && old->kind == TERM_SLASH_BUILTIN) {
        mth_message("SlashCommand", "builtin", "/%s is a built-in command and cannot be redefined.", name);
        return failed();
    }
    if (argc == 1) {   /* SlashCommand["name"]: its summary, or None */
        if (!old) return expr_new_symbol(SYM_None);
        return expr_new_string(old->summary);
    }
    if (a[1]->type == EXPR_SYMBOL && a[1]->data.symbol.name == SYM_None) {
        term_slash_remove(name);
        return expr_new_symbol(SYM_Null);
    }
    const char* summary = "";
    if (argc == 3) {
        if (a[2]->type != EXPR_STRING) {
            mth_message("SlashCommand", "string", "String expected at position 3 in SlashCommand.");
            return NULL;
        }
        summary = a[2]->data.string;
    }
    Expr* f = expr_copy(a[1]);
    if (!term_slash_add_mod(name, summary, g_owner, f, free_expr_payload)) {
        expr_free(f);
        return failed();
    }
    char* out = malloc(nl + 2);
    if (!out) return NULL;
    out[0] = '/';
    memcpy(out + 1, name, nl + 1);
    return expr_new_string_adopt(out);
}

/* ---- Monitor ------------------------------------------------------------------ */

#define MON_MAX 16
static Expr* g_mon[MON_MAX];
static int   g_mon_depth = 0;     /* registered with the TimeConstrained unwind */
static int   g_mon_ticking = 0;   /* likewise: a timeout inside a tick */

Expr* builtin_monitor(Expr* res) {
    if (res->data.function.arg_count != 2) {
        mth_message("Monitor", "argrx", "Monitor called with %zu arguments; 2 arguments are expected.",
                    res->data.function.arg_count);
        return NULL;
    }
    Expr** a = res->data.function.args;
    int d = g_mon_depth;
    if (d >= 0 && d < MON_MAX) g_mon[d] = expr_copy(a[1]);
    g_mon_depth = d + 1;
    mth_monitor_set_active(g_mon_depth);
    Expr* r = evaluate(a[0]);
    /* pop to the depth this call saw (a missed inner pop is repaired too) */
    for (int i = g_mon_depth; i > d; i--)
        if (i - 1 < MON_MAX && g_mon[i - 1]) { expr_free(g_mon[i - 1]); g_mon[i - 1] = NULL; }
    g_mon_depth = d;
    mth_monitor_set_active(g_mon_depth);
    return r;
}

void mth_monitor_tick(void) {
    g_mth_monitor_due = 0;
    if (g_mon_ticking) return;
    int d = g_mon_depth;
    if (d <= 0) { mth_monitor_set_active(0); return; }
    if (d > MON_MAX) d = MON_MAX;
    Expr* mon = g_mon[d - 1];
    if (!mon) return;
    g_mon_ticking = 1;
    mth_msg_suppress_push();               /* the line shows a value, not messages */
    Expr* v = evaluate(mon);
    mth_msg_suppress_pop();
    char* s = v ? expr_to_string_display(v) : NULL;
    if (s) {
        for (char* p = s; *p; p++) if (*p == '\n' || *p == '\t') *p = ' ';
        mth_monitor_publish(s);
        free(s);
    }
    expr_free(v);
    g_mon_ticking = 0;
}

/* ---- registration ----------------------------------------------------------- */

void frontend_init(void) {
    symtab_add_builtin("SlashCommand", builtin_slashcommand);
    symtab_get_def("SlashCommand")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("SlashCommand",
        "SlashCommand[\"name\", f, \"summary\"]\n"
        "\tadds the interactive REPL command /name: typing \"/name text\" evaluates\n"
        "\tf[\"text\"] (the trimmed rest of the line) and prints a non-Null result.\n"
        "SlashCommand[\"name\", f]\n"
        "\tadds it without a summary.\n"
        "SlashCommand[\"name\", None]\n"
        "\tremoves it.\n"
        "SlashCommand[\"name\"]\n"
        "\tgives its summary, or None.\n"
        "SlashCommand[]\n"
        "\tlists every command as \"/name\" -> \"summary\".\n"
        "Names are a letter followed by letters and digits; the REPL's own commands\n"
        "(/help, /theme, ...) cannot be redefined. Usually called from\n"
        "~/.mathilda/init.m or a file in ~/.mathilda/mods/. Attributes: Protected.");

    symtab_add_builtin("Monitor", builtin_monitor);
    symtab_get_def("Monitor")->attributes |= ATTR_HOLDALL | ATTR_PROTECTED;
    symtab_set_docstring("Monitor",
        "Monitor[expr, mon]\n"
        "\tevaluates expr and returns its value, showing the current value of mon\n"
        "\tin the interactive REPL's status line while it runs (refreshed about ten\n"
        "\ttimes a second; messages from mon are suppressed). Elsewhere it simply\n"
        "\tevaluates expr. Attributes: HoldAll, Protected.");

    tc_register_unwind_counter(&g_mon_depth);
    tc_register_unwind_counter(&g_mon_ticking);
}
