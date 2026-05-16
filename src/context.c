#include "context.h"
#include "symtab.h"
#include "attr.h"
#include "expr.h"
#include "sym_names.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* -------------------- State -------------------- */

static char*  g_current     = NULL;   /* always ends with a backtick */
static char** g_path        = NULL;   /* ordered; first entry is most-specific */
static size_t g_path_count  = 0;
static size_t g_path_cap    = 0;

typedef enum {
    FRAME_BEGIN,
    FRAME_PACKAGE
} FrameKind;

typedef struct CtxFrame {
    FrameKind kind;
    char*     saved_context;
    char**    saved_path;
    size_t    saved_path_count;
    struct CtxFrame* next;
} CtxFrame;

static CtxFrame* g_stack = NULL;

/* -------------------- Helpers -------------------- */

static char* ctx_strdup(const char* s) {
    if (!s) return NULL;
    size_t n = strlen(s);
    char* r = (char*)malloc(n + 1);
    memcpy(r, s, n + 1);
    return r;
}

static char* ctx_concat(const char* a, const char* b) {
    size_t na = strlen(a);
    size_t nb = strlen(b);
    char* r = (char*)malloc(na + nb + 1);
    memcpy(r, a, na);
    memcpy(r + na, b, nb + 1);
    return r;
}

static bool ends_with_tick(const char* s) {
    size_t n = s ? strlen(s) : 0;
    return n > 0 && s[n - 1] == '`';
}

static bool has_backtick(const char* s) {
    for (; *s; s++) if (*s == '`') return true;
    return false;
}

/* ----- path manipulation ----- */

static void path_ensure(size_t need) {
    if (need <= g_path_cap) return;
    size_t cap = g_path_cap ? g_path_cap : 4;
    while (cap < need) cap *= 2;
    g_path = (char**)realloc(g_path, cap * sizeof(char*));
    g_path_cap = cap;
}

static void path_append(const char* s) {
    path_ensure(g_path_count + 1);
    g_path[g_path_count++] = ctx_strdup(s);
}

static void path_prepend(const char* s) {
    path_ensure(g_path_count + 1);
    for (size_t i = g_path_count; i > 0; i--) g_path[i] = g_path[i - 1];
    g_path[0] = ctx_strdup(s);
    g_path_count++;
}

static bool path_contains(const char* s) {
    for (size_t i = 0; i < g_path_count; i++) {
        if (strcmp(g_path[i], s) == 0) return true;
    }
    return false;
}

static void path_clear(void) {
    for (size_t i = 0; i < g_path_count; i++) free(g_path[i]);
    g_path_count = 0;
}

static char** path_snapshot(size_t* out_count) {
    *out_count = g_path_count;
    if (g_path_count == 0) return NULL;
    char** copy = (char**)malloc(g_path_count * sizeof(char*));
    for (size_t i = 0; i < g_path_count; i++) copy[i] = ctx_strdup(g_path[i]);
    return copy;
}

static void path_restore(char** snap, size_t count) {
    path_clear();
    path_ensure(count);
    for (size_t i = 0; i < count; i++) g_path[i] = snap[i]; /* ownership transferred */
    g_path_count = count;
    free(snap);
}

/* -------------------- $Context / $ContextPath publishing -------------------- */

/* Update the OwnValue of the named symbol to `replacement`. Replaces any
 * existing own-value by clearing them first. Used to republish $Context /
 * $ContextPath whenever state changes. */
static void publish_own(const char* name, Expr* replacement) {
    symtab_clear_symbol(name);
    Expr* pat = expr_new_symbol(name);
    symtab_add_own_value(name, pat, replacement);
    expr_free(pat);
    expr_free(replacement);
}

static void republish_state(void) {
    publish_own("$Context", context_as_string());
    publish_own("$ContextPath", context_path_as_list());
}

Expr* context_as_string(void) {
    return expr_new_string(g_current ? g_current : "Global`");
}

Expr* context_path_as_list(void) {
    Expr** args = (Expr**)malloc(g_path_count * sizeof(Expr*));
    for (size_t i = 0; i < g_path_count; i++) args[i] = expr_new_string(g_path[i]);
    Expr* head = expr_new_symbol("List");
    Expr* res = expr_new_function(head, args, g_path_count);
    /* expr_new_function memcpy's the slot pointers into its own array;
     * the caller-owned `args` buffer is no longer needed. */
    free(args);
    return res;
}

/* -------------------- Public accessors -------------------- */

const char* context_current(void) { return g_current ? g_current : "Global`"; }
size_t      context_path_count(void) { return g_path_count; }
const char* context_path_at(size_t i) { return (i < g_path_count) ? g_path[i] : NULL; }

/* -------------------- Resolution -------------------- */

/* Return length of the context prefix of `name`, including the final backtick
 * (e.g. for "foo`bar`x" returns 8). Returns 0 if no backtick present. */
static size_t context_prefix_len(const char* name) {
    size_t last = 0;
    for (size_t i = 0; name[i]; i++) if (name[i] == '`') last = i + 1;
    return last;
}

bool context_is_visible(const char* prefix) {
    if (!prefix || *prefix == '\0') return false;
    if (strcmp(prefix, "System`") == 0) return true; /* System` always visible */
    if (g_current && strcmp(prefix, g_current) == 0) return true;
    return path_contains(prefix);
}

const char* context_display_name(const char* name) {
    if (!name) return name;
    size_t plen = context_prefix_len(name);
    if (plen == 0) return name;
    /* We need a probe string equal to the prefix; do it without allocation. */
    char prefix[128];
    if (plen >= sizeof(prefix)) return name;
    memcpy(prefix, name, plen);
    prefix[plen] = '\0';
    if (context_is_visible(prefix)) return name + plen;
    return name;
}

char* context_resolve_name(const char* raw) {
    if (!raw || *raw == '\0') return ctx_strdup(raw ? raw : "");

    /* System variables ($Foo) bypass context qualification. */
    if (raw[0] == '$') return ctx_strdup(raw);

    /* Leading backtick: relative to current context. */
    if (raw[0] == '`') {
        const char* cur = context_current();
        return ctx_concat(cur, raw + 1);
    }

    if (has_backtick(raw)) {
        /* Absolute qualified name. Canonicalize "System`Foo" -> bare "Foo"
         * because builtins are stored under bare names. */
        const char* sys = "System`";
        size_t sn = 7; /* strlen("System`") */
        if (strncmp(raw, sys, sn) == 0 && raw[sn] != '\0') {
            return ctx_strdup(raw + sn);
        }
        return ctx_strdup(raw);
    }

    /* Bare name. Search {current} U path for an existing qualified match. */
    const char* cur = context_current();
    char* probe = ctx_concat(cur, raw);
    if (symtab_lookup(probe)) return probe;
    free(probe);

    for (size_t i = 0; i < g_path_count; i++) {
        const char* prefix = g_path[i];
        if (strcmp(prefix, "System`") == 0) continue; /* System is checked via bare */
        probe = ctx_concat(prefix, raw);
        if (symtab_lookup(probe)) return probe;
        free(probe);
    }

    /* Fall back: if a bare builtin exists, use it (implicit System`). */
    SymbolDef* bare = symtab_lookup(raw);
    if (bare && bare->builtin_func) return ctx_strdup(raw);

    /* If current context is Global`, keep bare for legacy behavior. Otherwise
     * qualify with the current context. */
    if (strcmp(cur, "Global`") == 0) return ctx_strdup(raw);
    return ctx_concat(cur, raw);
}

/* -------------------- Stack ops -------------------- */

static void frame_push(FrameKind kind) {
    CtxFrame* f = (CtxFrame*)malloc(sizeof(CtxFrame));
    f->kind = kind;
    f->saved_context = ctx_strdup(g_current);
    f->saved_path = path_snapshot(&f->saved_path_count);
    f->next = g_stack;
    g_stack = f;
}

static bool frame_pop(char** out_closed_ctx) {
    if (!g_stack) return false;
    CtxFrame* f = g_stack;
    g_stack = f->next;
    if (out_closed_ctx) *out_closed_ctx = ctx_strdup(g_current);
    free(g_current);
    g_current = f->saved_context;
    path_restore(f->saved_path, f->saved_path_count);
    free(f);
    return true;
}

bool context_begin(const char* ctx) {
    if (!ctx || !ends_with_tick(ctx)) return false;
    frame_push(FRAME_BEGIN);
    char* new_ctx;
    if (ctx[0] == '`') {
        /* Relative context: extend current. */
        new_ctx = ctx_concat(g_current, ctx + 1);
    } else {
        new_ctx = ctx_strdup(ctx);
    }
    free(g_current);
    g_current = new_ctx;
    republish_state();
    return true;
}

bool context_begin_package(const char* ctx, const char** needs, size_t needs_count) {
    if (!ctx || !ends_with_tick(ctx)) return false;
    if (ctx[0] == '`') return false; /* BeginPackage requires absolute context */
    frame_push(FRAME_PACKAGE);
    free(g_current);
    g_current = ctx_strdup(ctx);
    path_clear();
    path_append(g_current);
    path_append("System`");
    for (size_t i = 0; i < needs_count; i++) {
        if (needs[i] && ends_with_tick(needs[i]) && !path_contains(needs[i])) {
            path_append(needs[i]);
        }
    }
    republish_state();
    return true;
}

bool context_end(char** out_ended_ctx) {
    bool ok = frame_pop(out_ended_ctx);
    if (ok) republish_state();
    return ok;
}

bool context_end_package(char** out_ended_ctx) {
    if (!g_stack) return false;
    char* closed = NULL;
    if (!frame_pop(&closed)) return false;
    /* Prepend the closed package context to $ContextPath (if not already on it). */
    if (closed && !path_contains(closed)) path_prepend(closed);
    if (out_ended_ctx) *out_ended_ctx = closed;
    else free(closed);
    republish_state();
    return true;
}

/* -------------------- Builtins -------------------- */

/* Context[]               -> current context
 * Context[sym]            -> context of the symbol
 * Context["name"]         -> context of the symbol named "name" if it exists */
static Expr* builtin_context(Expr* res) {
    size_t argc = res->data.function.arg_count;
    if (argc == 0) {
        Expr* out = expr_new_string(context_current());
        return out;
    }
    if (argc != 1) return NULL;

    Expr* a = res->data.function.args[0];
    const char* name = NULL;
    if (a->type == EXPR_SYMBOL) {
        name = a->data.symbol;
    } else if (a->type == EXPR_STRING) {
        name = a->data.string;
    } else {
        return NULL;
    }

    /* If the name already has a backtick, peel off its prefix as the answer. */
    size_t plen = context_prefix_len(name);
    if (plen > 0) {
        char* prefix = (char*)malloc(plen + 1);
        memcpy(prefix, name, plen);
        prefix[plen] = '\0';
        Expr* out = expr_new_string(prefix);
        free(prefix);
        return out;
    }

    /* Bare name: look it up. If it exists as a builtin, report "System`".
     * If it exists as a plain symbol with rules/docstring, report "Global`".
     * If it doesn't exist at all (string form), emit a notfound-style
     * warning and return the expression unevaluated. */
    SymbolDef* def = symtab_lookup(name);
    if (def) {
        const char* home = (def->builtin_func != NULL) ? "System`" : "Global`";
        Expr* out = expr_new_string(home);
        return out;
    }

    if (a->type == EXPR_STRING) {
        fprintf(stderr, "Context::notfound: Symbol %s not found.\n", name);
        return NULL;
    }

    /* Unassigned bare symbol: treat as future Global`. */
    return expr_new_string("Global`");
}

static const char* string_arg(Expr* e) {
    return (e && e->type == EXPR_STRING) ? e->data.string : NULL;
}

static Expr* builtin_begin(Expr* res) {
    if (res->data.function.arg_count != 1) return NULL;
    const char* ctx = string_arg(res->data.function.args[0]);
    if (!ctx) return NULL;
    if (!context_begin(ctx)) {
        fprintf(stderr, "Begin::cxt: Context specification %s is not valid.\n", ctx);
        return NULL;
    }
    return expr_new_string(context_current());
}

static Expr* builtin_begin_package(Expr* res) {
    size_t argc = res->data.function.arg_count;
    if (argc < 1 || argc > 2) return NULL;
    const char* ctx = string_arg(res->data.function.args[0]);
    if (!ctx) return NULL;

    const char** needs = NULL;
    size_t needs_count = 0;
    Expr* needs_list = (argc == 2) ? res->data.function.args[1] : NULL;
    if (needs_list) {
        if (needs_list->type != EXPR_FUNCTION
            || needs_list->data.function.head->type != EXPR_SYMBOL
            || needs_list->data.function.head->data.symbol != SYM_List) {
            return NULL;
        }
        needs_count = needs_list->data.function.arg_count;
        needs = (const char**)malloc(needs_count * sizeof(char*));
        for (size_t i = 0; i < needs_count; i++) {
            Expr* ne = needs_list->data.function.args[i];
            if (ne->type != EXPR_STRING) { free((void*)needs); return NULL; }
            needs[i] = ne->data.string;
        }
    }

    bool ok = context_begin_package(ctx, needs, needs_count);
    if (needs) free((void*)needs);
    if (!ok) {
        fprintf(stderr, "BeginPackage::cxt: Context specification %s is not valid.\n", ctx);
        return NULL;
    }
    return expr_new_string(context_current());
}

static Expr* builtin_end(Expr* res) {
    if (res->data.function.arg_count != 0) return NULL;
    char* closed = NULL;
    if (!context_end(&closed)) {
        fprintf(stderr, "End::noctx: No previous context to revert to.\n");
        return NULL;
    }
    Expr* out = expr_new_string(closed ? closed : "");
    free(closed);
    return out;
}

static Expr* builtin_end_package(Expr* res) {
    if (res->data.function.arg_count != 0) return NULL;
    char* closed = NULL;
    if (!context_end_package(&closed)) {
        fprintf(stderr, "EndPackage::noctx: No previous package to close.\n");
        return NULL;
    }
    free(closed);
    /* EndPackage returns Null (an unevaluated Symbol named "Null"). */
    return expr_new_symbol("Null");
}

/* -------------------- Init / Shutdown -------------------- */

void context_init(void) {
    if (g_current) return; /* idempotent */
    g_current = ctx_strdup("Global`");
    path_append("Global`");
    path_append("System`");

    symtab_add_builtin("Context",       builtin_context);
    symtab_add_builtin("BeginPackage",  builtin_begin_package);
    symtab_add_builtin("Begin",         builtin_begin);
    symtab_add_builtin("End",           builtin_end);
    symtab_add_builtin("EndPackage",    builtin_end_package);

    /* Attributes. Context must not evaluate its symbol argument. */
    symtab_get_def("Context")->attributes      |= ATTR_HOLDFIRST | ATTR_PROTECTED;
    symtab_get_def("BeginPackage")->attributes |= ATTR_PROTECTED;
    symtab_get_def("Begin")->attributes        |= ATTR_PROTECTED;
    symtab_get_def("End")->attributes          |= ATTR_PROTECTED;
    symtab_get_def("EndPackage")->attributes   |= ATTR_PROTECTED;

    /* Docstrings. */
    symtab_set_docstring("Context",
        "Context[] gives the current context ($Context).\n"
        "Context[sym] gives the context in which sym resides.\n"
        "Context[\"name\"] gives the context of the symbol named \"name\" if it exists.");
    symtab_set_docstring("BeginPackage",
        "BeginPackage[\"ctx`\"] sets the current context to \"ctx`\" and restricts\n"
        "$ContextPath to {\"ctx`\", \"System`\"}, matching Mathematica's package\n"
        "prologue.\nBeginPackage[\"ctx`\", {\"need1`\", ...}] additionally prepends the\n"
        "listed contexts to $ContextPath.");
    symtab_set_docstring("Begin",
        "Begin[\"ctx`\"] sets the current context ($Context) to \"ctx`\", saving\n"
        "the previous value for End[] to restore. If the argument starts with a\n"
        "backtick it is interpreted relative to the current context:\n"
        "Begin[\"`Private`\"] inside \"MyPkg`\" yields \"MyPkg`Private`\".");
    symtab_set_docstring("End",
        "End[] restores the context that was active before the matching Begin[]\n"
        "and returns the closed context as a string.");
    symtab_set_docstring("EndPackage",
        "EndPackage[] restores the state saved by BeginPackage and prepends the\n"
        "just-closed package context to $ContextPath so its exported symbols are\n"
        "visible under short names. Returns Null.");

    /* $Context / $ContextPath are ordinary symbols whose OwnValues track the
     * internal state. Protect them so the user cannot reassign them directly;
     * mutation goes through Begin/BeginPackage/End/EndPackage. */
    republish_state();
    symtab_get_def("$Context")->attributes     |= ATTR_PROTECTED;
    symtab_get_def("$ContextPath")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("$Context",
        "$Context is a string giving the current context. New symbols are created\n"
        "in this context unless otherwise qualified. Modified by Begin[],\n"
        "BeginPackage[], End[], and EndPackage[].");
    symtab_set_docstring("$ContextPath",
        "$ContextPath is a list of contexts used (in order) to resolve bare\n"
        "identifiers to existing qualified symbols. Modified by BeginPackage[]\n"
        "and EndPackage[].");
}

void context_shutdown(void) {
    while (g_stack) {
        CtxFrame* f = g_stack;
        g_stack = f->next;
        free(f->saved_context);
        for (size_t i = 0; i < f->saved_path_count; i++) free(f->saved_path[i]);
        free(f->saved_path);
        free(f);
    }
    free(g_current); g_current = NULL;
    path_clear();
    free(g_path); g_path = NULL; g_path_cap = 0;
}
