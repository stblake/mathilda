#include "context.h"
#include "symtab.h"
#include "attr.h"
#include "expr.h"
#include "sym_names.h"
#include "sym_intern.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "message.h"   /* mth_message: Quiet/Check funnel */
#include "print.h"     /* expr_to_string_outputform for ::cxlist / ::cxset */
#include "assoc.h"     /* $ContextAliases is an Association */
#include "common.h"    /* builtin_arg_error */
#include "readwrite.h" /* readwrite_get: Needs loads through Get */

/* -------------------- State -------------------- */

/* The context state IS the OwnValues of $Context and $ContextPath: every reader
 * below consults them directly, so any write -- Set, AppendTo, Block, Clear,
 * Begin/End -- takes effect at once and the symbols never disagree with the
 * resolver.  Set/SetDelayed are validated first (context_assign, called from
 * apply_assignment); a value written another way is read defensively: a
 * $Context that is not a context string reads as "Global`", and entries of
 * $ContextPath that are not context strings are skipped.  Pointers handed out
 * (context_current, context_path_at) borrow the OwnValue and stay valid until
 * the variable is next assigned.
 *
 * $Packages (the contexts of loaded packages; BeginPackage prepends) and
 * $ContextAliases (an Association "alias`" -> "ctx`", consulted when a
 * qualified name is resolved) follow the same convention. */

typedef enum {
    FRAME_BEGIN,
    FRAME_PACKAGE
} FrameKind;

typedef struct CtxFrame {
    FrameKind kind;
    Expr*     saved_context;   /* owned copies of the values on entry (or NULL) */
    Expr*     saved_path;
    struct CtxFrame* next;
} CtxFrame;

static CtxFrame* g_stack = NULL;
static bool      g_inited = false;

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

static bool is_context_string(const Expr* e) {
    return e && e->type == EXPR_STRING && ends_with_tick(e->data.string);
}

static bool is_list(const Expr* e) {
    return e && e->type == EXPR_FUNCTION &&
           e->data.function.head->type == EXPR_SYMBOL &&
           e->data.function.head->data.symbol.name == SYM_List;
}

/* Borrowed OwnValue of a symbol, or NULL. */
static Expr* own_value(const char* name) {
    Rule* r = symtab_get_own_values(name);
    return r ? r->replacement : NULL;
}

/* $ContextPath's List (borrowed), or NULL when it is not a List. */
static Expr* path_list(void) {
    Expr* v = own_value("$ContextPath");
    return is_list(v) ? v : NULL;
}

static bool path_contains(const char* s) {
    Expr* l = path_list();
    if (!l) return false;
    for (size_t i = 0; i < l->data.function.arg_count; i++) {
        Expr* e = l->data.function.args[i];
        if (is_context_string(e) && strcmp(e->data.string, s) == 0) return true;
    }
    return false;
}

static bool is_rule(const Expr* e) {
    return e && e->type == EXPR_FUNCTION && e->data.function.arg_count == 2 &&
           e->data.function.head->type == EXPR_SYMBOL &&
           (e->data.function.head->data.symbol.name == SYM_Rule ||
            e->data.function.head->data.symbol.name == SYM_RuleDelayed);
}

/* True when the string list `l` holds context string `s`. */
static bool list_has_string(const Expr* l, const char* s) {
    if (!is_list(l)) return false;
    for (size_t i = 0; i < l->data.function.arg_count; i++) {
        Expr* e = l->data.function.args[i];
        if (is_context_string(e) && strcmp(e->data.string, s) == 0) return true;
    }
    return false;
}

/* -------------------- $Context / $ContextPath publishing -------------------- */

/* Replace the OwnValue of the named symbol with `replacement` (consumed). */
static void publish_own(const char* name, Expr* replacement) {
    symtab_clear_symbol(name);
    Expr* pat = expr_new_symbol(name);
    symtab_add_own_value(name, pat, replacement);
    expr_free(pat);
    expr_free(replacement);
}

/* A fresh List of the given context strings. */
static Expr* make_path(const char** ctxs, size_t n) {
    Expr** args = (Expr**)malloc((n ? n : 1) * sizeof(Expr*));
    for (size_t i = 0; i < n; i++) args[i] = expr_new_string(ctxs[i]);
    Expr* res = expr_new_function(expr_new_symbol(SYM_List), args, n);
    free(args);
    return res;
}

/* A fresh List: `first`, then the context strings of `rest` (borrowed; may be
 * NULL or a non-List) other than `first`. */
static Expr* prepend_unique(const char* first, const Expr* rest) {
    size_t n = is_list(rest) ? rest->data.function.arg_count : 0;
    const char** v = (const char**)malloc((n + 1) * sizeof(char*));
    size_t k = 0;
    v[k++] = first;
    for (size_t i = 0; i < n; i++) {
        Expr* e = rest->data.function.args[i];
        if (is_context_string(e) && strcmp(e->data.string, first) != 0)
            v[k++] = e->data.string;
    }
    Expr* out = make_path(v, k);
    free((void*)v);
    return out;
}

/* Put `ctx` at the front of $ContextPath unless it is already on it. */
static void path_prepend(const char* ctx) {
    if (path_contains(ctx)) return;
    publish_own("$ContextPath", prepend_unique(ctx, path_list()));  /* copies first */
}

bool context_packages_contains(const char* ctx) {
    return list_has_string(own_value("$Packages"), ctx);
}

/* Record a loaded package: move/put `ctx` at the front of $Packages. */
static void packages_prepend(const char* ctx) {
    publish_own("$Packages", prepend_unique(ctx, own_value("$Packages")));
}

/* Rewrite a leading context alias of the qualified name `raw`
 * ("al`x" -> "Pkg`x" when $ContextAliases holds "al`" -> "Pkg`"); the longest
 * matching alias wins.  Returns a malloc'd name, or NULL when no alias
 * applies.  $ContextAliases is read defensively: an Association or a List of
 * rules, with non-context-string entries skipped. */
static char* alias_rewrite(const char* raw) {
    Expr* a = own_value("$ContextAliases");
    if (!a || a->type != EXPR_FUNCTION || a->data.function.arg_count == 0 ||
        a->data.function.head->type != EXPR_SYMBOL ||
        (a->data.function.head->data.symbol.name != SYM_Association &&
         a->data.function.head->data.symbol.name != SYM_List))
        return NULL;
    const char* target = NULL;
    size_t best = 0;
    for (size_t i = 0; i < a->data.function.arg_count; i++) {
        Expr* r = a->data.function.args[i];
        if (!is_rule(r)) continue;
        Expr* k = r->data.function.args[0];
        Expr* v = r->data.function.args[1];
        if (!is_context_string(k) || !is_context_string(v)) continue;
        size_t kn = strlen(k->data.string);
        if (kn > best && strncmp(raw, k->data.string, kn) == 0) {
            best = kn;
            target = v->data.string;
        }
    }
    return target ? ctx_concat(target, raw + best) : NULL;
}

Expr* context_as_string(void) {
    return expr_new_string(context_current());
}

Expr* context_path_as_list(void) {
    Expr* l = path_list();
    return l ? expr_copy(l) : make_path(NULL, 0);
}

/* -------------------- Public accessors -------------------- */

const char* context_current(void) {
    Expr* v = own_value("$Context");
    return is_context_string(v) ? v->data.string : "Global`";
}

size_t context_path_count(void) {
    Expr* l = path_list();
    size_t n = 0;
    if (l)
        for (size_t i = 0; i < l->data.function.arg_count; i++)
            n += is_context_string(l->data.function.args[i]);
    return n;
}

const char* context_path_at(size_t i) {
    Expr* l = path_list();
    if (!l) return NULL;
    for (size_t k = 0; k < l->data.function.arg_count; k++) {
        Expr* e = l->data.function.args[k];
        if (is_context_string(e) && i-- == 0) return e->data.string;
    }
    return NULL;
}

const char* context_home(const char* bare, const SymbolDef* def) {
    if (intern_is_mathilda(bare)) return "Mathilda`";
    if (def && def->builtin_func != NULL) return "System`";
    bool spelled = bare[0] == '$' || (bare[0] >= 'A' && bare[0] <= 'Z');
    return (spelled && intern_is_system(bare)) ? "System`" : "Global`";
}

/* -------------------- Assignment (Set / SetDelayed) -------------------- */

bool context_is_state_var(const char* name) {
    return name && name[0] == '$' &&
           (strcmp(name, "$Context") == 0 || strcmp(name, "$ContextPath") == 0);
}

bool context_assign(const char* name, Expr* value) {
    bool ok;
    if (strcmp(name, "$Context") == 0) {
        ok = is_context_string(value);
    } else {
        ok = is_list(value);
        for (size_t i = 0; ok && i < value->data.function.arg_count; i++)
            ok = is_context_string(value->data.function.args[i]);
    }
    if (ok) {
        publish_own(name, expr_copy(value));
        return true;
    }
    char* txt = expr_to_string_outputform(value);
    if (strcmp(name, "$Context") == 0)
        mth_message("$Context", "cxset",
                    "Cannot set $Context to %s; value is not a parsable context name ending in `.",
                    txt ? txt : "?");
    else
        mth_message("$ContextPath", "cxlist",
                    "Cannot set $ContextPath to %s; value must be a list of strings ending in `.",
                    txt ? txt : "?");
    free(txt);
    return false;
}

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
    if (strcmp(prefix, context_current()) == 0) return true;
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
/* Developer` names that are Mathilda builtins under their bare name: code
 * written for the reference kernel spells the packing helpers with their
 * Developer` context, and that spelling must reach the same builtin. */
static bool is_developer_alias(const char* name) {
    static const char* const names[] = { "ToPackedArray", "FromPackedArray", "PackedArrayQ" };
    for (size_t i = 0; i < sizeof names / sizeof names[0]; i++)
        if (strcmp(name, names[i]) == 0) return true;
    return false;
}


char* context_resolve_name(const char* raw) {
    if (!raw || *raw == '\0') return ctx_strdup(raw ? raw : "");

    /* A leading '$' is NOT special: `$x` inside a package is qualified into
     * the package context like any other name (issue #104). System `$`
     * variables ($Version, $Failed, $Assumptions, ...) stay bare via the
     * System fall-back below, since each has a SymbolDef or the System flag. */

    /* Leading backtick: relative to current context. */
    if (raw[0] == '`') {
        const char* cur = context_current();
        return ctx_concat(cur, raw + 1);
    }

    if (has_backtick(raw)) {
        /* Absolute qualified name: first expand a $ContextAliases alias, then
         * canonicalize "System`Foo" -> bare "Foo" because builtins are stored
         * under bare names. */
        char* aliased = alias_rewrite(raw);
        const char* q = aliased ? aliased : raw;
        const char* sys = "System`";
        size_t sn = 7; /* strlen("System`") */
        const char* mth = "Mathilda`";
        size_t mn = 9; /* strlen("Mathilda`") */
        char* out;
        const char* glb = "Global`";
        size_t gn = 7; /* strlen("Global`") */
        if (strncmp(q, sys, sn) == 0 && q[sn] != '\0')
            out = ctx_strdup(q + sn);
        else if (strncmp(q, mth, mn) == 0 && q[mn] != '\0' && intern_is_mathilda(q + mn))
            out = ctx_strdup(q + mn);                 /* Mathilda`Foo is stored bare too */
        else if (strncmp(q, "Developer`", 10) == 0 && is_developer_alias(q + 10))
            out = ctx_strdup(q + 10);                 /* Developer`ToPackedArray IS ToPackedArray */
        else if (strncmp(q, glb, gn) == 0 && q[gn] != '\0' && !has_backtick(q + gn) &&
                 strcmp(context_home(q + gn, symtab_lookup(q + gn)), "Global`") == 0)
            /* Global` symbols are stored under bare keys: Global`x IS x (unless
             * the bare name is a System` / Mathilda` symbol, which Global`Sin
             * must not alias). */
            out = ctx_strdup(q + gn);
        else
            out = ctx_strdup(q);
        free(aliased);
        return out;
    }

    /* Bare name. Search {current} U path for an existing qualified match. */
    const char* cur = context_current();
    char* probe = ctx_concat(cur, raw);
    if (symtab_lookup(probe)) return probe;
    free(probe);

    size_t np = context_path_count();
    for (size_t i = 0; i < np; i++) {
        const char* prefix = context_path_at(i);
        if (strcmp(prefix, "System`") == 0) continue; /* System is checked via bare */
        probe = ctx_concat(prefix, raw);
        if (symtab_lookup(probe)) return probe;
        free(probe);
    }

    /* Fall back: resolve to the bare (implicit System`) name when it denotes a
     * System symbol. System` is always on the conceptual context path, and its
     * symbols are stored under bare names. A name qualifies as System if it
     * either already has a SymbolDef (registered builtins, and option symbols
     * that carry a docstring) OR was flagged System at kernel-init time
     * (intern_mark_all_system). The flag is essential because many System
     * symbols -- List, Symbol, Rule, Integer, Infinity, Automatic, None, the
     * pattern heads, ... -- are pure interned symbols with no SymbolDef until
     * first use. Without it, bare heads like `_List` / `_Symbol` or defaults
     * like Automatic inside a package (a non-Global current context) were
     * mis-qualified to the package's private context, so e.g. `allOpts_List`
     * never matched a real System`List. In the Global case this is a no-op:
     * the legacy branch below already returns the bare name. */
    SymbolDef* bare = symtab_lookup(raw);
    if (bare || intern_is_system(raw) || intern_is_mathilda(raw)) {
        /* Bare keys hold System`, Mathilda` AND Global` symbols alike, so check
         * the bare symbol's home context is visible: System` always is; a
         * Mathilda-only builtin or a Global` symbol only while its context is
         * current or on $ContextPath.  Inside BeginPackage (path {ctx`,
         * System`}) a package's own x or DataType is therefore a new symbol in
         * its context -- not Global`x left behind by earlier code, nor the
         * kernel's DataType. */
        const char* home = context_home(raw, bare);
        if (strcmp(home, "System`") == 0 || context_is_visible(home))
            return ctx_strdup(raw);
    }

    /* If current context is Global`, keep bare for legacy behavior. Otherwise
     * qualify with the current context. */
    if (strcmp(cur, "Global`") == 0) return ctx_strdup(raw);
    return ctx_concat(cur, raw);
}

/* -------------------- Stack ops -------------------- */

static Expr* copy_or_null(Expr* e) { return e ? expr_copy(e) : NULL; }

static void frame_push(FrameKind kind) {
    CtxFrame* f = (CtxFrame*)malloc(sizeof(CtxFrame));
    f->kind = kind;
    f->saved_context = copy_or_null(own_value("$Context"));
    f->saved_path = copy_or_null(own_value("$ContextPath"));
    f->next = g_stack;
    g_stack = f;
}

/* Restore a saved value (consumed); NULL clears the variable. */
static void restore_own(const char* name, Expr* saved) {
    if (saved) publish_own(name, saved);
    else symtab_clear_symbol(name);
}

static bool frame_pop(char** out_closed_ctx) {
    if (!g_stack) return false;
    CtxFrame* f = g_stack;
    g_stack = f->next;
    if (out_closed_ctx) *out_closed_ctx = ctx_strdup(context_current());
    restore_own("$Context", f->saved_context);
    restore_own("$ContextPath", f->saved_path);
    free(f);
    return true;
}

bool context_begin(const char* ctx) {
    if (!ctx || !ends_with_tick(ctx)) return false;
    /* Relative context: extend the current one. */
    char* new_ctx = (ctx[0] == '`') ? ctx_concat(context_current(), ctx + 1)
                                    : ctx_strdup(ctx);
    frame_push(FRAME_BEGIN);
    publish_own("$Context", expr_new_string(new_ctx));
    free(new_ctx);
    return true;
}

bool context_begin_package(const char* ctx, const char** needs, size_t needs_count) {
    if (!ctx || !ends_with_tick(ctx)) return false;
    if (ctx[0] == '`') return false; /* BeginPackage requires absolute context */
    frame_push(FRAME_PACKAGE);
    /* $ContextPath = {ctx, "System`", needs...} (needs deduplicated). */
    const char** path = (const char**)malloc((needs_count + 2) * sizeof(char*));
    size_t n = 0;
    path[n++] = ctx;
    path[n++] = "System`";
    for (size_t i = 0; i < needs_count; i++) {
        if (!needs[i] || !ends_with_tick(needs[i])) continue;
        bool dup = false;
        for (size_t k = 0; k < n && !dup; k++) dup = strcmp(path[k], needs[i]) == 0;
        if (!dup) path[n++] = needs[i];
    }
    publish_own("$Context", expr_new_string(ctx));
    publish_own("$ContextPath", make_path(path, n));
    free((void*)path);
    packages_prepend(ctx);
    return true;
}

bool context_end(char** out_ended_ctx) {
    return frame_pop(out_ended_ctx);
}

bool context_end_package(char** out_ended_ctx) {
    if (!g_stack) return false;
    char* closed = NULL;
    if (!frame_pop(&closed)) return false;
    /* Prepend the closed package context to $ContextPath (if not already on it). */
    if (closed) path_prepend(closed);
    if (out_ended_ctx) *out_ended_ctx = closed;
    else free(closed);
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
        name = a->data.symbol.name;
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

    /* Bare string name: resolve it through {$Context} U $ContextPath exactly
     * as the parser would, so Context["fooF"] after EndPackage reports the
     * package context of an existing Foo`fooF (it used to say notfound). */
    if (a->type == EXPR_STRING) {
        char* q = context_resolve_name(name);
        size_t qlen = q ? context_prefix_len(q) : 0;
        if (qlen > 0 && symtab_lookup(q)) {
            q[qlen] = '\0';
            Expr* out = expr_new_string(q);
            free(q);
            return out;
        }
        free(q);
    }

    /* Bare name: look it up. A builtin, or a name flagged System at kernel
     * init (Pi, Automatic, $Version, $Failed, ...), reports "System`"; any
     * other existing symbol reports "Global`". The init flag also covers
     * scratch variables the C side interns (x, y, t, ...), so it is trusted
     * only for names spelled like System names: a leading '$' or capital.
     * If the name doesn't exist at all (string form), emit a notfound-style
     * warning and return unevaluated. */
    SymbolDef* def = symtab_lookup(name);
    const char* home = context_home(name, def);
    if (def || strcmp(home, "System`") == 0) return expr_new_string(home);

    if (a->type == EXPR_STRING) {
        mth_message("Context", "notfound", "Symbol %s not found.", name);
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
        mth_message("Begin", "cxt", "Context specification %s is not valid.", ctx);
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
            || needs_list->data.function.head->data.symbol.name != SYM_List) {
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
        mth_message("BeginPackage", "cxt", "Context specification %s is not valid.", ctx);
        return NULL;
    }
    return expr_new_string(context_current());
}

static Expr* builtin_end(Expr* res) {
    if (res->data.function.arg_count != 0) return NULL;
    char* closed = NULL;
    if (!context_end(&closed)) {
        mth_message("End", "noctx", "No previous context to revert to.");
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
        mth_message("EndPackage", "noctx", "No previous package to close.");
        return NULL;
    }
    free(closed);
    /* EndPackage returns Null (an unevaluated Symbol named "Null"). */
    return expr_new_symbol(SYM_Null);
}

/* -------------------- Needs -------------------- */

/* A loadable context name: one or more symbol-name segments, each followed by
 * a backtick ("Pkg`", "A`B`").  Non-ASCII bytes count as letters. */
static bool is_valid_context_name(const char* s) {
    if (!s || !*s || !ends_with_tick(s)) return false;
    bool seg_start = true;
    for (const unsigned char* p = (const unsigned char*)s; *p; p++) {
        unsigned char c = *p;
        if (c == '`') {
            if (seg_start) return false;            /* empty segment */
            seg_start = true;
            continue;
        }
        bool letter = (c >= 'A' && c <= 'Z') || (c >= 'a' && c <= 'z') ||
                      c == '$' || c >= 0x80;
        bool digit = c >= '0' && c <= '9';
        if (seg_start ? !letter : !(letter || digit)) return false;
        seg_start = false;
    }
    return true;
}

static void needs_msg(const char* tag, const char* fmt, Expr* res) {
    char* txt = expr_to_string_outputform(res);
    mth_message("Needs", tag, fmt, txt ? txt : "Needs[...]");
    free(txt);
}

/* $ContextAliases["alias`"] = "ctx`", keeping the variable an Association. */
static void add_alias(const char* alias, const char* ctx) {
    Expr* cur = own_value("$ContextAliases");
    Expr* key = expr_new_string(alias);
    Expr* updated = NULL;
    if (cur && is_association(cur)) {
        Expr* val = expr_new_string(ctx);
        updated = assoc_set_key(cur, key, val);
        if (!updated) expr_free(val);
    }
    if (!updated) {
        /* Missing, a List of rules, or garbage: rebuild from its valid rules. */
        size_t n = (cur && cur->type == EXPR_FUNCTION) ? cur->data.function.arg_count : 0;
        Expr** rules = (Expr**)malloc((n + 1) * sizeof(Expr*));
        size_t k = 0;
        for (size_t i = 0; i < n; i++)
            if (is_rule(cur->data.function.args[i])) rules[k++] = cur->data.function.args[i];
        Expr** rargs = (Expr**)malloc(2 * sizeof(Expr*));
        rargs[0] = expr_copy(key);
        rargs[1] = expr_new_string(ctx);
        Expr* rule = expr_new_function(expr_new_symbol(SYM_Rule), rargs, 2);
        free(rargs);
        rules[k++] = rule;
        updated = assoc_from_rules(rules, k);
        expr_free(rule);
        free(rules);
    }
    expr_free(key);
    publish_own("$ContextAliases", updated);
}

/* Needs["ctx`"]                 Get["ctx`"] unless ctx` is in $Packages, then
 *                               put ctx` on $ContextPath.
 * Needs["ctx`", "file"]         the same, loading "file" instead.
 * Needs["ctx`" -> "alias`"]     load, keep ctx` off $ContextPath (undoing the
 *                               package's own EndPackage), and record the
 *                               alias in $ContextAliases.
 * Needs["ctx`" -> None]         load only.
 * Returns Null; when the load does not create ctx`, emits Needs::nocont and
 * returns what Get returned ($Failed if the file could not be opened). */
static Expr* builtin_needs(Expr* res) {
    size_t argc = res->data.function.arg_count;
    if (argc < 1 || argc > 2) return builtin_arg_error("Needs", argc, 1, 2);

    Expr* spec = res->data.function.args[0];
    const char* ctx = NULL;
    const char* alias = NULL;
    bool rule_form = false;
    if (spec->type == EXPR_STRING) {
        ctx = spec->data.string;
    } else if (is_rule(spec) && spec->data.function.args[0]->type == EXPR_STRING) {
        Expr* rhs = spec->data.function.args[1];
        rule_form = true;
        ctx = spec->data.function.args[0]->data.string;
        if (rhs->type == EXPR_STRING && is_valid_context_name(rhs->data.string))
            alias = rhs->data.string;
        else if (!(rhs->type == EXPR_SYMBOL && rhs->data.symbol.name == SYM_None))
            ctx = NULL;
    }
    if (!ctx) {
        needs_msg("cxru", "Context or appropriately structured rule expected at position 1 in %s.", res);
        return NULL;
    }
    if (!is_valid_context_name(ctx)) {
        needs_msg("cxt", "Invalid context specified at position 1 in %s. A context must consist of "
                  "valid symbol names separated by and ending with `.", res);
        return NULL;
    }
    const char* file = NULL;
    if (argc == 2) {
        Expr* f = res->data.function.args[1];
        if (f->type != EXPR_STRING) {
            needs_msg("string", "String expected at position 2 in %s.", res);
            return NULL;
        }
        file = f->data.string;
    }

    if (!context_packages_contains(ctx)) {
        /* The rule form keeps ctx` off $ContextPath, so the path the package's
         * EndPackage extends is put back afterwards. */
        Expr* saved_path = rule_form ? copy_or_null(own_value("$ContextPath")) : NULL;
        Expr* got = readwrite_get(file ? file : ctx);
        if (rule_form) restore_own("$ContextPath", saved_path);
        if (!context_packages_contains(ctx)) {
            mth_message("Needs", "nocont",
                        "Context %s was not created when Needs was evaluated.", ctx);
            return got;
        }
        expr_free(got);
    }

    if (!rule_form) path_prepend(ctx);
    if (alias) add_alias(alias, ctx);
    return expr_new_symbol(SYM_Null);
}

/* -------------------- Init / Shutdown -------------------- */

void context_init(void) {
    if (g_inited) return; /* idempotent */
    g_inited = true;

    symtab_add_builtin("Context",       builtin_context);
    symtab_add_builtin("BeginPackage",  builtin_begin_package);
    symtab_add_builtin("Begin",         builtin_begin);
    symtab_add_builtin("End",           builtin_end);
    symtab_add_builtin("EndPackage",    builtin_end_package);
    symtab_add_builtin("Needs",         builtin_needs);

    /* Attributes. Context must not evaluate its symbol argument. */
    symtab_get_def("Context")->attributes      |= ATTR_HOLDFIRST | ATTR_PROTECTED;
    symtab_get_def("BeginPackage")->attributes |= ATTR_PROTECTED;
    symtab_get_def("Begin")->attributes        |= ATTR_PROTECTED;
    symtab_get_def("End")->attributes          |= ATTR_PROTECTED;
    symtab_get_def("EndPackage")->attributes   |= ATTR_PROTECTED;
    symtab_get_def("Needs")->attributes        |= ATTR_PROTECTED;

    /* Docstrings. */
    symtab_set_docstring("Context",
        "Context[] gives the current context ($Context).\n"
        "Context[sym] gives the context in which sym resides.\n"
        "Context[\"name\"] gives the context of the symbol named \"name\" if it exists.");
    symtab_set_docstring("BeginPackage",
        "BeginPackage[\"ctx`\"] sets the current context to \"ctx`\" and restricts\n"
        "$ContextPath to {\"ctx`\", \"System`\"}, matching the standard package\n"
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

    /* $Context / $ContextPath are ordinary, assignable symbols, and their
     * OwnValues ARE the context state (see the State comment above). */
    publish_own("$Context", expr_new_string("Global`"));
    {
        const char* initial[] = { "Global`", "Mathilda`", "System`" };
        publish_own("$ContextPath", make_path(initial, 3));
        const char* packages[] = { "Mathilda`", "System`", "Global`" };
        publish_own("$Packages", make_path(packages, 3));
    }
    publish_own("$ContextAliases", assoc_from_rules(NULL, 0));
    symtab_set_docstring("$Context",
        "$Context is a string giving the current context. New symbols are created\n"
        "in this context unless otherwise qualified. Set by Begin[], BeginPackage[],\n"
        "End[] and EndPackage[], or assigned directly ($Context = \"ctx`\").");
    symtab_set_docstring("$ContextPath",
        "$ContextPath is a list of contexts used (in order) to resolve bare\n"
        "identifiers to existing qualified symbols, and whose symbols print\n"
        "without their context. Set by BeginPackage[] and EndPackage[], or\n"
        "assigned directly (e.g. AppendTo[$ContextPath, \"ctx`\"]); a Set to\n"
        "anything but a list of context strings emits $ContextPath::cxlist.");
}

void context_shutdown(void) {
    while (g_stack) {
        CtxFrame* f = g_stack;
        g_stack = f->next;
        expr_free(f->saved_context);
        expr_free(f->saved_path);
        free(f);
    }
    g_inited = false;
}
