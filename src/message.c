/*
 * message.c -- global message-suppression depth and the Quiet / Check /
 * Message builtins.  See message.h for rationale.
 */
#include "message.h"
#include "term/term_interrupt.h"
#include "expr.h"
#include "symtab.h"
#include "eval.h"
#include "attr.h"
#include "sym_names.h"
#include "print.h"     /* expr_to_string_outputform: Message argument text */
#include "context.h"   /* context_display_name: the printed message head */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* ------------------------------------------------- Quiet and Check frames */
/* Which messages a Quiet silences, and which a Check sees, follow the
 * reference's nesting rule: a message is silenced by the innermost active
 * Quiet whose spec covers it, and is seen only by the Checks entered INSIDE
 * that Quiet.  So Check[Quiet[1/0]; 5, 7] is 5 (the Quiet stands between the
 * Check and the message) while Quiet[Check[1/0; 5, 7]] is 7, and a message no
 * Quiet covers is seen by every enclosing Check (Check[Check[1/0, a], b] is b).
 *
 * The Quiet frames form a stack; frame i holds its spec (NULL = every message:
 * Quiet[expr] and every internal mth_msg_suppress_push).  Its depth is the
 * suppression depth.  A Check frame records the Quiet depth at its entry, its
 * own spec (NULL = any message) and whether a message it sees has fired. */
typedef struct { Expr* spec; } QuietFrame;          /* spec owned, or NULL */
static void g_msg_fired_bump(void);
typedef struct { int quiet_height; int fired; const Expr* spec; } CheckFrame;

static QuietFrame* g_quiet = NULL;
static int g_quiet_cap = 0;
static int g_msg_suppress_depth = 0;                /* number of Quiet frames */
static CheckFrame* g_checks = NULL;
static int g_checks_cap = 0;
static int g_check_depth = 0;                       /* number of Check frames */

/* Make room for frame index `need` in a growable frame array. */
static int frames_reserve(void** arr, int* cap, int need, size_t elt) {
    if (need < *cap) return 1;
    int nc = *cap ? *cap * 2 : 16;
    while (nc <= need) nc *= 2;
    void* g = realloc(*arr, (size_t)nc * elt);
    if (!g) return 0;
    *arr = g;
    *cap = nc;
    return 1;
}

/* Push a Quiet frame; `spec` is adopted (NULL = all messages).  Without memory
 * for the frame the spec widens to "all", which only silences more. */
static void quiet_push(Expr* spec) {
    if (!frames_reserve((void**)&g_quiet, &g_quiet_cap, g_msg_suppress_depth,
                        sizeof(QuietFrame))) {
        if (spec) expr_free(spec);
        spec = NULL;
        if (g_msg_suppress_depth >= g_quiet_cap) { g_msg_suppress_depth++; return; }
    }
    g_quiet[g_msg_suppress_depth++].spec = spec;
}

/* Drop the Quiet frames above depth d. */
static void quiet_truncate(int d) {
    if (d < 0) d = 0;
    while (g_msg_suppress_depth > d) {
        g_msg_suppress_depth--;
        if (g_msg_suppress_depth < g_quiet_cap && g_quiet[g_msg_suppress_depth].spec) {
            expr_free(g_quiet[g_msg_suppress_depth].spec);
            g_quiet[g_msg_suppress_depth].spec = NULL;
        }
    }
}

void mth_msg_suppress_push(void) { quiet_push(NULL); }

void mth_msg_suppress_pop(void) {
    if (g_msg_suppress_depth > 0) quiet_truncate(g_msg_suppress_depth - 1);
}

/* Does message spec `spec` name head::tag?  A spec is All, None, a message
 * name sym::tag (General::tag names tag under every head) or a list of these;
 * anything else -- including a symbol whose VALUE is a list, since Quiet and
 * Check hold the spec -- names nothing, as in the reference.  head == NULL
 * asks about an anonymous diagnostic, which only All covers. */
static int spec_covers(const Expr* spec, const char* head, const char* tag) {
    if (!spec) return 1;
    if (spec->type == EXPR_SYMBOL) return spec->data.symbol.name == SYM_All;
    if (spec->type != EXPR_FUNCTION || spec->data.function.head->type != EXPR_SYMBOL)
        return 0;
    const char* h = spec->data.function.head->data.symbol.name;
    if (h == SYM_List) {
        for (size_t i = 0; i < spec->data.function.arg_count; i++)
            if (spec_covers(spec->data.function.args[i], head, tag)) return 1;
        return 0;
    }
    if (h == SYM_MessageName && head && spec->data.function.arg_count >= 2) {
        const Expr* s = spec->data.function.args[0];
        const Expr* t = spec->data.function.args[1];
        if (s->type != EXPR_SYMBOL || t->type != EXPR_STRING) return 0;
        return strcmp(t->data.string, tag) == 0 &&
               (strcmp(s->data.symbol.name, head) == 0 ||
                strcmp(s->data.symbol.name, "General") == 0);
    }
    return 0;
}

/* Index of the innermost Quiet frame silencing head::tag, or -1. */
static int quiet_cover_index(const char* head, const char* tag) {
    int top = g_msg_suppress_depth < g_quiet_cap ? g_msg_suppress_depth : g_quiet_cap;
    if (g_msg_suppress_depth > top) return g_msg_suppress_depth - 1;   /* OOM frames */
    for (int i = top - 1; i >= 0; i--)
        if (spec_covers(g_quiet[i].spec, head, tag)) return i;
    return -1;
}

/* A message fired: count it, and mark the Checks that see it -- those entered
 * inside the silencing Quiet frame q (all of them when q == -1). */
static void note_fired_for(int q, const char* head, const char* tag) {
    g_msg_fired_bump();
    for (int i = g_check_depth - 1; i >= 0 && g_checks[i].quiet_height > q; i--)
        if (!g_checks[i].spec || spec_covers(g_checks[i].spec, head, tag))
            g_checks[i].fired = 1;
}

/* Non-zero iff a Quiet that silences every message is active: what a caller
 * can ask without naming the message (mth_message_cont, internal probes). */
int mth_msg_suppressed(void) { return quiet_cover_index(NULL, NULL) >= 0; }

/* Save / restore the Quiet and Check frame depths as one value (Quiet depth in
 * the low 16 bits, Check depth above).  TimeConstrained's timeout siglongjmp
 * can unwind out of a Quiet[] or Check[] region (its pop then never runs,
 * leaving messages silenced, or a stale Check frame marked, for the rest of the
 * session); the guard captures the depths before the body and restores them
 * after the jump. */
int mth_msg_suppress_depth_save(void) {
    return (g_msg_suppress_depth & 0xFFFF) | (g_check_depth << 16);
}
void mth_msg_suppress_depth_load(int d) {
    quiet_truncate(d & 0xFFFF);
    int c = (d >> 16) & 0x7FFF;
    if (g_check_depth > c) g_check_depth = c;
}

/* Depth for the `Solve::ifun` advisory scope (see message.h).  Zero => the
 * message may print; > 0 => a complete-solution caller (Reduce) is on the
 * stack, so the "use Reduce" advice is moot and suppressed. */
static int g_ifun_suppress_depth = 0;

void mth_msg_ifun_suppress_push(void) { g_ifun_suppress_depth++; }

void mth_msg_ifun_suppress_pop(void) {
    if (g_ifun_suppress_depth > 0) g_ifun_suppress_depth--;
}

int mth_msg_ifun_suppressed(void) { return g_ifun_suppress_depth > 0; }

/* Save / restore the ifun-advisory depth across a TimeConstrained unwind, the
 * same way mth_msg_suppress_depth_save/load guard the Quiet depth.  A timeout
 * siglongjmp out of builtin_reduce (which pushes this depth, then grinds in the
 * CAD) would otherwise skip its mth_msg_ifun_suppress_pop and leave Solve::ifun
 * silenced for the rest of the session. */
int  mth_msg_ifun_suppress_depth_save(void)  { return g_ifun_suppress_depth; }
void mth_msg_ifun_suppress_depth_load(int d) { g_ifun_suppress_depth = d; }

/* Message-fired counter (see message.h). */
static unsigned long g_msg_fired = 0;
static void g_msg_fired_bump(void) { g_msg_fired++; }
void          mth_msg_note_fired(void)  { note_fired_for(quiet_cover_index(NULL, NULL), NULL, NULL); }
void mth_msg_note_named(const char* head, const char* tag) {
    if (mth_msg_is_off(head, tag)) return;
    note_fired_for(quiet_cover_index(head, tag), head, tag);
}
unsigned long mth_msg_fired_count(void) { return g_msg_fired; }
void          mth_msg_fired_restore(unsigned long n) { if (n < g_msg_fired) g_msg_fired = n; }

/* ------------------------------------------------------------- Off / On */
/* The set of switched-off messages, as "head::tag" keys.  "General::tag"
 * switches tag off for every head; "head::*" (from Off[head]) every tag of
 * head.  Small and rarely non-empty, so a linear scan of a flat array is
 * cheaper than any table -- and the funnel skips it entirely when empty. */
static char** g_off = NULL;
static size_t g_off_n = 0, g_off_cap = 0;

static long off_find(const char* key) {
    for (size_t i = 0; i < g_off_n; i++)
        if (strcmp(g_off[i], key) == 0) return (long)i;
    return -1;
}

static char* off_key(const char* head, const char* tag) {
    size_t n = strlen(head) + strlen(tag) + 3;
    char* k = malloc(n);
    if (k) snprintf(k, n, "%s::%s", head, tag);
    return k;
}

/* Switch head::tag off (on = 0) or back on (on = 1). */
static void off_set(const char* head, const char* tag, int on) {
    char* k = off_key(head, tag);
    if (!k) return;
    long i = off_find(k);
    if (!on && i < 0) {
        if (g_off_n == g_off_cap) {
            size_t nc = g_off_cap ? 2 * g_off_cap : 8;
            char** ng = realloc(g_off, nc * sizeof(char*));
            if (!ng) { free(k); return; }
            g_off = ng; g_off_cap = nc;
        }
        g_off[g_off_n++] = k;
        return;
    }
    if (on && i >= 0) {
        free(g_off[i]);
        g_off[i] = g_off[--g_off_n];
    }
    free(k);
}

int mth_msg_is_off(const char* head, const char* tag) {
    if (g_off_n == 0) return 0;
    char* k;
    int off = 0;
    if ((k = off_key(head, tag)))      { off = off_find(k) >= 0; free(k); }
    if (!off && (k = off_key("General", tag))) { off = off_find(k) >= 0; free(k); }
    if (!off && (k = off_key(head, "*")))      { off = off_find(k) >= 0; free(k); }
    return off;
}

/* ------------------------------------------------------------ $MessageList */
/* $MessageList: HoldForm[sym::tag] for each message PRINTED during the current
 * input line (a message silenced by Quiet or Off is not listed).  It is an
 * ordinary OwnValue, so Block[{$MessageList = {}}, ...] localises it. */
static void message_list_set(Expr* val) {
    Expr* sym = expr_new_symbol("$MessageList");
    symtab_add_own_value("$MessageList", sym, val);
    expr_free(sym);
}

void mth_message_list_reset(void) {
    Expr* empty = expr_new_function(expr_new_symbol(SYM_List), NULL, 0);
    message_list_set(empty);
    expr_free(empty);
}

static void message_list_append(const char* head, const char* tag) {
    Rule* r = symtab_get_own_values("$MessageList");
    Expr* cur = (r && r->replacement) ? r->replacement : NULL;
    size_t n = (cur && cur->type == EXPR_FUNCTION && cur->data.function.head->type == EXPR_SYMBOL &&
                cur->data.function.head->data.symbol.name == SYM_List) ? cur->data.function.arg_count : 0;
    Expr** items = malloc((n + 1) * sizeof(Expr*));
    if (!items) return;
    for (size_t i = 0; i < n; i++) items[i] = expr_copy(cur->data.function.args[i]);
    Expr* mn_args[2] = { expr_new_symbol(head), expr_new_string(tag) };
    Expr* mn = expr_new_function(expr_new_symbol(SYM_MessageName), mn_args, 2);
    items[n] = expr_new_function(expr_new_symbol(SYM_HoldForm), &mn, 1);
    Expr* l = expr_new_function(expr_new_symbol(SYM_List), items, n + 1);
    free(items);
    message_list_set(l);
    expr_free(l);
}

/* ------------------------------------------------------------- The funnel */
/* The single choke-point for user-facing diagnostics: note (always), then
 * print "Head::tag: <formatted>\n" unless Quiet[] (or an extra local mute) is
 * active.  See message.h. */
/* The interactive REPL's style for the printed Head::tag (borrowed static
 * strings; "" -- byte-identical output -- unless the REPL sets them). */
static const char* g_msg_on = "";
static const char* g_msg_off = "";

void mth_message_set_style(const char* head_on, const char* head_off) {
    g_msg_on = head_on ? head_on : "";
    g_msg_off = head_off ? head_off : "";
}

void mth_message_v(int extra_mute, const char* head, const char* tag,
                   const char* fmt, va_list ap) {
    /* extra_mute is a subsystem-local *probe* mute (g_arith_warnings_muted,
     * g_fm_quiet): an internal probe that pokes at divergent forms and treats
     * the diagnostic as sampling noise, NOT a user-facing event.  So it is
     * suppressed ENTIRELY -- neither printed nor noted -- otherwise an enclosing
     * Check would catch the noise (e.g. Check[Limit[Sin[x]/x, x -> 0], bad]
     * must return 1, not bad, though the 0*ComplexInfinity probe fires
     * Infinity::indet internally).  Quiet is different: it is a request to
     * silence the print of a real message, so the firing is still noted, and a
     * Check[] inside the Quiet still sees it (the Quiet[Check[...]] idiom). */
    if (extra_mute) return;
    /* An interactive interrupt is unwinding (src/term/term_interrupt.h):
     * whatever a declining kernel or a cleanup path says now is noise, so it is
     * neither printed nor noted. Inside AbortProtect messages still print. */
    if (mth_interrupt_should_stop()) return;
    /* Off[head::tag]: the message does not exist for this session -- neither
     * printed nor seen by Check. */
    if (mth_msg_is_off(head, tag)) return;
    /* Seen by the Checks inside the Quiet that silences it (all of them when
     * none does); a silenced message is not printed and not in $MessageList. */
    int q = quiet_cover_index(head, tag);
    note_fired_for(q, head, tag);
    if (q >= 0) return;
    message_list_append(head, tag);
    /* FeynCalc`DiracTrick::failmsg prints as DiracTrick::failmsg while
     * FeynCalc` is on $ContextPath, as the reference prints it. */
    fprintf(stderr, "%s%s::%s%s: ", g_msg_on, context_display_name(head), tag, g_msg_off);
    vfprintf(stderr, fmt, ap);
    fputc('\n', stderr);
}

void mth_message(const char* head, const char* tag, const char* fmt, ...) {
    va_list ap;
    va_start(ap, fmt);
    mth_message_v(0, head, tag, fmt, ap);
    va_end(ap);
}

void mth_message_gated(int extra_mute, const char* head, const char* tag,
                       const char* fmt, ...) {
    va_list ap;
    va_start(ap, fmt);
    mth_message_v(extra_mute != 0, head, tag, fmt, ap);
    va_end(ap);
}

void mth_message_line(int extra_mute, const char* line) {
    /* "Head::tag: text\n" -> mth_message_v(extra_mute, Head, tag, "%s", text). */
    const char* sep = strstr(line, "::");
    const char* colon = sep ? strstr(sep + 2, ": ") : NULL;
    if (!sep || !colon || sep - line >= 128 || colon - (sep + 2) >= 64) {
        mth_message_gated(extra_mute, "General", "msg", "%s", line);
        return;
    }
    char head[128], tag[64];
    memcpy(head, line, (size_t)(sep - line));
    head[sep - line] = '\0';
    memcpy(tag, sep + 2, (size_t)(colon - (sep + 2)));
    tag[colon - (sep + 2)] = '\0';
    const char* text = colon + 2;
    size_t tl = strlen(text);
    while (tl > 0 && text[tl - 1] == '\n') tl--;
    mth_message_gated(extra_mute, head, tag, "%.*s", (int)tl, text);
}

void mth_message_linev(int extra_mute, const char* fmt, va_list ap) {
    char stack[512];
    va_list ap2;
    va_copy(ap2, ap);
    int n = vsnprintf(stack, sizeof stack, fmt, ap);
    if (n < 0) { va_end(ap2); return; }
    if ((size_t)n < sizeof stack) {
        va_end(ap2);
        mth_message_line(extra_mute, stack);
        return;
    }
    char* big = malloc((size_t)n + 1);
    if (!big) { va_end(ap2); mth_message_line(extra_mute, stack); return; }
    vsnprintf(big, (size_t)n + 1, fmt, ap2);
    va_end(ap2);
    mth_message_line(extra_mute, big);
    free(big);
}

void mth_message_cont(int extra_mute, const char* fmt, ...) {
    if (extra_mute || mth_msg_suppressed()) return;    /* no note: already fired */
    if (mth_interrupt_should_stop()) return;           /* as in mth_message_v */
    va_list ap;
    va_start(ap, fmt);
    vfprintf(stderr, fmt, ap);
    va_end(ap);
    fputc('\n', stderr);
}

/* ------------------------------------------------- Quiet / Check / Message */

/* Quiet[expr] / Quiet[expr, spec]: evaluate expr with the messages spec names
 * silenced -- every message, or All, None, sym::tag (General::tag: that tag
 * under any head) or a list of these -- and return its value.  HoldAll: the
 * spec is read as written, so a symbol whose value is a list names nothing.
 * A silenced message is hidden from the Checks outside this Quiet (see the
 * frame model above).  Quiet[expr, moff, mon]'s third argument (messages to
 * switch back on inside an outer Quiet) is ignored. */
Expr* builtin_quiet(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count < 1) return NULL;
    Expr* spec = res->data.function.arg_count >= 2 ? expr_copy(res->data.function.args[1]) : NULL;
    int depth = g_msg_suppress_depth;
    quiet_push(spec);
    Expr* val = evaluate(res->data.function.args[0]);
    quiet_truncate(depth);
    return val;   /* an in-flight Throw sentinel propagates unchanged */
}

/* Check[expr, failexpr] / Check[expr, failexpr, spec]: evaluate expr; if a
 * message it can see fired during that evaluation -- one no Quiet inside the
 * Check silenced, and, with a spec, one the spec names -- return the
 * evaluated failexpr, otherwise expr's value.  HoldAll.  A Throw propagates
 * (it is not a message). */
Expr* builtin_check(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count < 2) return NULL;
    size_t argc = res->data.function.arg_count;
    int frame = g_check_depth;
    if (!frames_reserve((void**)&g_checks, &g_checks_cap, frame, sizeof(CheckFrame)))
        return NULL;
    /* Check[expr, failexpr, s1::t1, s2::t2, ...]: several names (or lists of
     * names) are the List of them, as in the reference; spec_covers walks it.
     * One name keeps borrowing the argument. */
    Expr* owned_spec = NULL;
    if (argc > 3) {
        Expr** names = malloc(sizeof(Expr*) * (argc - 2));
        if (!names) return NULL;
        for (size_t i = 2; i < argc; i++) names[i - 2] = expr_copy(res->data.function.args[i]);
        owned_spec = expr_new_function(expr_new_symbol(SYM_List), names, argc - 2);
        free(names);
    }
    g_checks[frame].quiet_height = g_msg_suppress_depth;
    g_checks[frame].fired = 0;
    g_checks[frame].spec = owned_spec ? owned_spec
                         : argc == 3 ? res->data.function.args[2] : NULL;
    g_check_depth = frame + 1;
    Expr* val = evaluate(res->data.function.args[0]);
    int fired = g_check_depth > frame && g_checks[frame].fired;
    g_check_depth = frame;
    if (owned_spec) expr_free(owned_spec);
    if (val && eval_is_inflight_unwind(val)) return val;
    if (fired) {
        if (val) expr_free(val);
        return evaluate(res->data.function.args[1]);
    }
    return val;
}

/* Is e MessageName[sym, "tag"]?  Sets *head / *tag (borrowed). */
static int message_name_parts(const Expr* e, const char** head, const char** tag) {
    if (!e || e->type != EXPR_FUNCTION || e->data.function.head->type != EXPR_SYMBOL ||
        e->data.function.head->data.symbol.name != SYM_MessageName ||
        e->data.function.arg_count < 2) return 0;
    const Expr* s = e->data.function.args[0];
    const Expr* t = e->data.function.args[1];
    if (s->type != EXPR_SYMBOL || t->type != EXPR_STRING) return 0;
    *head = s->data.symbol.name;
    *tag = t->data.string;
    return 1;
}

/* Texts of the standard General:: messages, used when General::tag has no
 * definition of its own (an assignment General::tag = "..." takes precedence).
 * They let Message[f::argx, f, 3] and friends print the usual text. */
static const struct { const char* tag; const char* text; } GENERAL_TEXTS[] = {
    { "argx", "`1` called with `2` arguments; 1 argument is expected." },
    { "argr", "`1` called with 1 argument; `2` arguments are expected." },
    { "argrx", "`1` called with `2` arguments; `3` arguments are expected." },
    { "argt", "`1` called with `2` arguments; `3` or `4` arguments are expected." },
    { "argtu", "`1` called with 1 argument; `2` or `3` arguments are expected." },
    { "argb", "`1` called with `2` arguments; between `3` and `4` arguments are expected." },
    { "argbu", "`1` called with 1 argument; between `2` and `3` arguments are expected." },
    { "argm", "`1` called with `2` arguments; `3` or more arguments are expected." },
    { "argmu", "`1` called with 1 argument; `2` or more arguments are expected." },
    { "args", "`1` called with invalid parameters." },
    { "string", "String expected at position `1` in `2`." },
    { "strse", "String or list of strings expected at position `1` in `2`." },
    { "int", "Integer expected at position `2` in `1`." },
    { "intp", "Positive integer expected at position `2` in `1`." },
    { "intpm", "Positive machine-sized integer expected at position `2` in `1`." },
    { "intm", "Machine-sized integer expected at position `2` in `1`." },
    { "optx", "Unknown option `1` in `2`." },
    { "sym", "Argument `1` at position `2` is expected to be a symbol." },
    { "normal", "Nonatomic expression expected at position `1` in `2`." },
    { "list", "List expected at position `2` in `1`." },
    { "heads", "Heads `1` and `2` at positions `3` and `4` are expected to be the same." },
    { "partw", "Part `1` of `2` does not exist." },
    { "partd", "Part specification `1` is longer than depth of object." },
    { "pkspec1", "The expression `1` cannot be used as a part specification." },
    { "nonopt", "Options expected (instead of `1`) beyond position `2` in `3`. An option must be a rule or a list of rules." },
    { "noopen", "Cannot open `1`." },
    { "nffil", "File not found during `1`." },
    { "ovfl", "Overflow occurred in computation." },
    { "unfl", "Underflow occurred in computation." },
    { "infy", "Infinite expression `1` encountered." },
    { "indet", "Indeterminate expression `1` encountered." },
    { "wrsym", "Symbol `1` is Protected." },
    { "write", "Tag `1` in `2` is Protected." },
    { "level", "Level specification `1` is not of the form n, {n}, or {m, n}." },
    { "stop", "Further output of `1` will be suppressed during this calculation." },
    { "obspkg", "`1` is now obsolete. The legacy version being loaded may conflict with current functionality. See the Compatibility Guide for updating information." },
    { "ivar", "`1` is not a valid variable." },
    { "fnsym", "First argument in `1` is not a symbol or a string naming a symbol." },
    { "pspec", "Part specification `1` is neither a machine-sized integer nor a list of machine-sized integers." },
    { "rvalue", "`1` is not a variable with a value, so its value cannot be changed." },
    { "shdw", "Symbol `1` appears in multiple contexts `2`; definitions in context `3` may shadow or be shadowed by other definitions." },
    { "nonopt1", "Options expected (instead of `1`) beyond position `2` in `3`." },
    { "seqs", "Sequence specification (+n, -n, {+n}, {-n}, {m, n}, or {m, n, s}) expected at position `2` in `1`." },
};

static const char* general_text(const char* tag) {
    for (size_t i = 0; i < sizeof GENERAL_TEXTS / sizeof GENERAL_TEXTS[0]; i++)
        if (strcmp(GENERAL_TEXTS[i].tag, tag) == 0) return GENERAL_TEXTS[i].text;
    return NULL;
}

/* The text of head::tag, or NULL when none is defined (caller frees). */
static char* message_text(const char* head, const char* tag) {
    Expr* a[2] = { expr_new_symbol(head), expr_new_string(tag) };
    Expr* mn = expr_new_function(expr_new_symbol(SYM_MessageName), a, 2);
    Expr* v = evaluate(mn);
    expr_free(mn);
    char* out = NULL;
    if (v && v->type == EXPR_STRING) {
        size_t n = strlen(v->data.string) + 1;
        out = malloc(n);
        if (out) memcpy(out, v->data.string, n);
    }
    if (v) expr_free(v);
    return out;
}

/* Growable output buffer for the formatted message. */
typedef struct { char* p; size_t n, cap; } MsgBuf;
static void mb_put(MsgBuf* b, const char* s, size_t len) {
    if (b->n + len + 1 > b->cap) {
        size_t nc = b->cap ? b->cap : 64;
        while (b->n + len + 1 > nc) nc *= 2;
        char* np = realloc(b->p, nc);
        if (!np) return;
        b->p = np; b->cap = nc;
    }
    memcpy(b->p + b->n, s, len);
    b->n += len;
    b->p[b->n] = '\0';
}

/* Fill the `n` and `` slots of `tmpl` with args[0..nargs) in OutputForm.
 * A slot past the last argument is left as written, after StringForm::sfr. */
static char* message_fill(const char* tmpl, Expr** args, size_t nargs) {
    MsgBuf b = { NULL, 0, 0 };
    mb_put(&b, "", 0);
    size_t seq = 0;
    for (const char* p = tmpl; *p; ) {
        if (*p == '`') {
            const char* q = p + 1;
            while (*q >= '0' && *q <= '9') q++;
            if (*q == '`') {
                size_t k = (q == p + 1) ? ++seq : (size_t)strtoul(p + 1, NULL, 10);
                if (k >= 1 && k <= nargs) {
                    char* t = expr_to_string_outputform(args[k - 1]);
                    if (t) { mb_put(&b, t, strlen(t)); free(t); }
                } else {
                    mth_message("StringForm", "sfr",
                                "Item %zu requested in \"%s\" out of range; %zu items available.",
                                k, tmpl, nargs);
                    mb_put(&b, p, (size_t)(q - p) + 1);
                }
                p = q + 1;
                continue;
            }
        }
        mb_put(&b, p, 1);
        p++;
    }
    return b.p;
}

/* Message[sym::tag, e1, e2, ...] -- emit the message sym::tag with its `n`
 * slots filled from the ei (in OutputForm).  The text is sym::tag's, else
 * General::tag's, else "-- Message text not found --" followed by each ei in
 * parentheses.  It goes through the mth_message funnel, so it prints as
 * "sym::tag: text", is silenced by Quiet and by Off[sym::tag], and is seen by
 * Check unless switched off.  HoldFirst.  Returns Null. */
Expr* builtin_message(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count < 1) return NULL;
    const char* head;
    const char* tag;
    if (!message_name_parts(res->data.function.args[0], &head, &tag)) {
        /* Not sym::tag: still a fired diagnostic for Check. */
        mth_msg_note_fired();
        return expr_new_symbol(SYM_Null);
    }
    if (mth_msg_is_off(head, tag)) return expr_new_symbol(SYM_Null);
    Expr** args = res->data.function.args + 1;
    size_t nargs = res->data.function.arg_count - 1;
    char* tmpl = message_text(head, tag);
    if (!tmpl) tmpl = message_text("General", tag);
    if (!tmpl) {
        const char* g = general_text(tag);
        if (g) {
            size_t n = strlen(g) + 1;
            tmpl = malloc(n);
            if (tmpl) memcpy(tmpl, g, n);
        }
    }
    char* text;
    if (tmpl) {
        text = message_fill(tmpl, args, nargs);
        free(tmpl);
    } else {
        MsgBuf b = { NULL, 0, 0 };
        const char* nf = "-- Message text not found --";
        mb_put(&b, nf, strlen(nf));
        for (size_t i = 0; i < nargs; i++) {
            char* t = expr_to_string_outputform(args[i]);
            mb_put(&b, " (", 2);
            if (t) { mb_put(&b, t, strlen(t)); free(t); }
            mb_put(&b, ")", 1);
        }
        text = b.p;
    }
    mth_message(head, tag, "%s", text ? text : "");
    free(text);
    return expr_new_symbol(SYM_Null);
}

/* Off[s::t, ...] / On[s::t, ...]: switch messages off or on.  A spec is a
 * message name s::t (General::t covers every head), a symbol s (every message
 * of s), or a list of specs.  HoldAll.  Returns Null. */
static void off_on_spec(const Expr* e, int on) {
    const char* head;
    const char* tag;
    if (message_name_parts(e, &head, &tag)) { off_set(head, tag, on); return; }
    if (e->type == EXPR_SYMBOL) { off_set(e->data.symbol.name, "*", on); return; }
    if (e->type == EXPR_FUNCTION && e->data.function.head->type == EXPR_SYMBOL &&
        e->data.function.head->data.symbol.name == SYM_List) {
        for (size_t i = 0; i < e->data.function.arg_count; i++)
            off_on_spec(e->data.function.args[i], on);
        return;
    }
    /* A message group ("name" strings, ...) is not supported. */
    char* txt = expr_to_string_outputform((Expr*)e);
    mth_message(on ? "On" : "Off", "mgre",
                "Message group %s was not resolved to a list of messages of the form "
                "symbol::name or symbol::name::language.", txt ? txt : "");
    free(txt);
}

static Expr* off_on(Expr* res, int on) {
    if (res->type != EXPR_FUNCTION) return NULL;
    for (size_t i = 0; i < res->data.function.arg_count; i++)
        off_on_spec(res->data.function.args[i], on);
    return expr_new_symbol(SYM_Null);
}

Expr* builtin_off(Expr* res) { return off_on(res, 0); }
Expr* builtin_on(Expr* res)  { return off_on(res, 1); }

void message_init(void) {
    symtab_add_builtin("Quiet", builtin_quiet);
    symtab_get_def("Quiet")->attributes |= ATTR_HOLDALL | ATTR_PROTECTED;
    symtab_set_docstring("Quiet",
        "Quiet[expr] evaluates expr with messages suppressed and returns its value.\n"
        "Quiet[expr, {s1::t1, ...}] suppresses only the named messages (General::t names t\n"
        "\tunder every head; All and None also work). A suppressed message is not seen by a\n"
        "\tCheck outside the Quiet, so Check[Quiet[1/0]; 5, 7] is 5.");

    symtab_add_builtin("Check", builtin_check);
    symtab_get_def("Check")->attributes |= ATTR_HOLDALL | ATTR_PROTECTED;
    symtab_set_docstring("Check",
        "Check[expr, failexpr] returns failexpr if a message is generated while "
        "evaluating expr, otherwise the value of expr.\n"
        "Check[expr, failexpr, {s1::t1, ...}] or Check[expr, failexpr, s1::t1, s2::t2, ...] "
        "checks only for the named messages.\n"
        "\tA message suppressed by a Quiet inside the Check does not count.");

    symtab_add_builtin("Message", builtin_message);
    symtab_get_def("Message")->attributes |= ATTR_HOLDFIRST | ATTR_PROTECTED;
    symtab_set_docstring("Message",
        "Message[sym::tag, e1, e2, ...]\n"
        "\tprints the message sym::tag as \"sym::tag: text\", with the slots `1`, `2`, ... (or ``)\n"
        "\tof its text filled by the ei.  The text is sym::tag's, else General::tag's.  Quiet and\n"
        "\tOff[sym::tag] silence it; Check sees it unless it is switched off.");

    mth_message_list_reset();
    symtab_set_docstring("$MessageList",
        "$MessageList\n\tis a list of the messages printed during the current input line, each as\n"
        "\tHoldForm[sym::tag].  Messages silenced by Quiet or Off are not listed.");

    symtab_add_builtin("Off", builtin_off);
    symtab_get_def("Off")->attributes |= ATTR_HOLDALL | ATTR_PROTECTED;
    symtab_set_docstring("Off",
        "Off[sym::tag]\n"
        "\tswitches off the message sym::tag, so it is neither printed nor seen by Check.\n"
        "Off[General::tag] switches tag off for every symbol; Off[sym] every message of sym;\n"
        "Off[s1, s2, ...] several at once.");
    symtab_add_builtin("On", builtin_on);
    symtab_get_def("On")->attributes |= ATTR_HOLDALL | ATTR_PROTECTED;
    symtab_set_docstring("On",
        "On[sym::tag]\n"
        "\tswitches the message sym::tag back on.  On[General::tag], On[sym] and\n"
        "On[s1, s2, ...] undo the matching Off.");
}
