/* sowreap.c — Sow / Reap: dynamic-scope accumulation of intermediate results.
 *
 * Mathematica semantics
 * ---------------------
 * Reap[expr] evaluates expr and returns {value, second}, where every value
 * passed to Sow during that evaluation is collected. Sow[e] returns e and
 * (as a side effect) records e in the nearest enclosing Reap that matches its
 * tag. Because a single Reap may capture many Sows, this needs genuine
 * dynamic state — a stack of active Reap "frames" — unlike Catch/Throw, whose
 * first-throw-wins semantics let it work statelessly via an in-band sentinel.
 *
 * Tags and patterns
 * -----------------
 *   Sow[e]                 -> tag None
 *   Sow[e, tag]            -> single tag
 *   Sow[e, {t1, t2, ...}]  -> one record per element (repeats allowed, so the
 *                             same value can appear several times)
 * A sown (tag, value) is routed to the innermost frame at least one of whose
 * patterns matches tag; within that frame it is appended to every matching
 * pattern-slot, grouped by structurally-identical tag in first-encounter order.
 *
 * Reap[expr]            == Reap[expr, _]      (single pattern _)
 * Reap[expr, patt]      single pattern; `second` is a flat list of entries
 * Reap[expr, {p1..pk}]  `second` has one slot per pattern (nesting one deeper)
 * Reap[expr, patt, f]   each entry is f[tag, {values}] instead of {values}
 *
 * Data structures
 * ---------------
 * Per the design, values and tag-groups are singly linked lists of Expr, which
 * give O(1) ordered append and preserve Sow order exactly. Each Reap frame
 * lives on builtin_reap's C stack and is linked into a file-global stack via
 * `prev`, so push/pop is a pointer write and the frame struct itself needs no
 * heap allocation. The system is single-threaded, so a plain global head is
 * sufficient. Every return path of builtin_reap pops its frame, keeping the
 * stack strictly balanced (LIFO).
 */

#include <stdbool.h>
#include <stdlib.h>

#include "sowreap.h"
#include "expr.h"
#include "eval.h"
#include "match.h"
#include "common.h"
#include "sym_names.h"

/* ------------------------------------------------------------------ */
/* Collector stack                                                     */
/* ------------------------------------------------------------------ */

typedef struct SowNode {
    Expr* value;             /* owned copy of a sown value, or NULL if moved out */
    struct SowNode* next;
} SowNode;

typedef struct SowGroup {
    Expr* tag;               /* owned copy of the grouping tag */
    SowNode* head;
    SowNode* tail;           /* O(1) append */
    struct SowGroup* next;   /* groups in first-encounter order */
} SowGroup;

typedef struct SowSlot {
    Expr* pattern;           /* borrowed from res (or the frame-owned default) */
    SowGroup* head;
    SowGroup* tail;
} SowSlot;

typedef struct ReapFrame {
    SowSlot* slots;          /* nslots entries */
    size_t   nslots;
    bool     is_list;        /* was the pattern argument a List? controls nesting */
    Expr*    f;              /* borrowed handler, or NULL */
    bool     own_pattern;    /* true when slots[0].pattern is our synthesized _ */
    struct ReapFrame* prev;
} ReapFrame;

/* Innermost active Reap, or NULL. Single-threaded: a plain global suffices. */
static ReapFrame* g_reap_top = NULL;

/* ------------------------------------------------------------------ */
/* Frame teardown                                                      */
/* ------------------------------------------------------------------ */

static void sow_group_free(SowGroup* g) {
    while (g) {
        SowGroup* gn = g->next;
        SowNode* n = g->head;
        while (n) {
            SowNode* nn = n->next;
            expr_free(n->value);   /* NULL-safe; already-moved nodes hold NULL */
            free(n);
            n = nn;
        }
        expr_free(g->tag);
        free(g);
        g = gn;
    }
}

static void reap_frame_free_contents(ReapFrame* fr) {
    for (size_t i = 0; i < fr->nslots; i++)
        sow_group_free(fr->slots[i].head);
    if (fr->own_pattern && fr->nslots > 0)
        expr_free(fr->slots[0].pattern);
    free(fr->slots);
}

/* ------------------------------------------------------------------ */
/* Sow                                                                 */
/* ------------------------------------------------------------------ */

/* Append a copy of `value` under `tag` into `slot` (create the group on first
 * encounter of the tag). */
static void slot_record(SowSlot* slot, Expr* tag, Expr* value) {
    SowGroup* g = slot->head;
    while (g && !expr_eq(g->tag, tag))
        g = g->next;
    if (!g) {
        g = malloc(sizeof(SowGroup));
        g->tag = expr_copy(tag);
        g->head = g->tail = NULL;
        g->next = NULL;
        if (slot->tail) slot->tail->next = g; else slot->head = g;
        slot->tail = g;
    }
    SowNode* n = malloc(sizeof(SowNode));
    n->value = expr_copy(value);
    n->next = NULL;
    if (g->tail) g->tail->next = n; else g->head = n;
    g->tail = n;
}

/* Route one (tag, value): find the innermost frame with a matching pattern and
 * record the value in each matching slot of that frame. If no frame matches,
 * the value is silently discarded (Sow still returns it to the caller). */
static void sow_route(Expr* tag, Expr* value) {
    for (ReapFrame* fr = g_reap_top; fr; fr = fr->prev) {
        bool frame_matched = false;
        for (size_t i = 0; i < fr->nslots; i++) {
            MatchEnv* env = env_new();
            bool m = match(tag, fr->slots[i].pattern, env);
            env_free(env);
            if (m) {
                slot_record(&fr->slots[i], tag, value);
                frame_matched = true;
            }
        }
        if (frame_matched)
            return;   /* innermost matching frame consumes the tag */
    }
}

Expr* builtin_sow(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    if (argc < 1 || argc > 2)
        return builtin_arg_error("Sow", argc, 1, 2);

    Expr* value = res->data.function.args[0];   /* already evaluated (non-held) */

    if (argc == 1) {
        Expr* none = expr_new_symbol(SYM_None);
        sow_route(none, value);
        expr_free(none);
    } else {
        Expr* tagarg = res->data.function.args[1];
        if (tagarg->type == EXPR_FUNCTION &&
            tagarg->data.function.head->type == EXPR_SYMBOL &&
            tagarg->data.function.head->data.symbol.name == SYM_List) {
            /* a list of tags: record once per element (repeats allowed) */
            size_t nt = tagarg->data.function.arg_count;
            for (size_t i = 0; i < nt; i++)
                sow_route(tagarg->data.function.args[i], value);
        } else {
            sow_route(tagarg, value);
        }
    }

    return expr_copy(value);
}

/* ------------------------------------------------------------------ */
/* Reap                                                                */
/* ------------------------------------------------------------------ */

/* Build the {v1, v2, ...} list for one group, MOVING the owned value pointers
 * out of the linked list (so the frame teardown won't double-free them). */
static Expr* group_values_list(SowGroup* g) {
    size_t count = 0;
    for (SowNode* n = g->head; n; n = n->next) count++;
    Expr** args = (count > 0) ? malloc(count * sizeof(Expr*)) : NULL;
    size_t i = 0;
    for (SowNode* n = g->head; n; n = n->next) {
        args[i++] = n->value;
        n->value = NULL;            /* moved out */
    }
    Expr* list = expr_new_function(expr_new_symbol(SYM_List), args, count);
    free(args);
    return list;
}

/* Produce the entry for one group: {values} normally, or f[tag, {values}]. */
static Expr* group_entry(SowGroup* g, Expr* f) {
    Expr* vals = group_values_list(g);
    if (!f) return vals;
    Expr** fa = malloc(2 * sizeof(Expr*));
    fa[0] = expr_copy(g->tag);
    fa[1] = vals;
    Expr* call = expr_new_function(expr_copy(f), fa, 2);
    free(fa);
    return eval_and_free(call);     /* evaluate f[tag, {values}] */
}

/* Build the list of entries collected in one slot, in first-encounter order. */
static Expr* slot_entries_list(SowSlot* slot, Expr* f) {
    size_t count = 0;
    for (SowGroup* g = slot->head; g; g = g->next) count++;
    Expr** args = (count > 0) ? malloc(count * sizeof(Expr*)) : NULL;
    size_t i = 0;
    for (SowGroup* g = slot->head; g; g = g->next)
        args[i++] = group_entry(g, f);
    Expr* list = expr_new_function(expr_new_symbol(SYM_List), args, count);
    free(args);
    return list;
}

Expr* builtin_reap(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    if (argc < 1 || argc > 3)
        return builtin_arg_error("Reap", argc, 1, 3);

    ReapFrame frame;
    frame.is_list = false;
    frame.own_pattern = false;
    frame.f = (argc == 3) ? res->data.function.args[2] : NULL;

    /* Resolve the pattern slots. */
    if (argc >= 2) {
        Expr* patt = res->data.function.args[1];
        if (patt->type == EXPR_FUNCTION &&
            patt->data.function.head->type == EXPR_SYMBOL &&
            patt->data.function.head->data.symbol.name == SYM_List) {
            frame.is_list = true;
            frame.nslots = patt->data.function.arg_count;
            frame.slots = (frame.nslots > 0)
                ? malloc(frame.nslots * sizeof(SowSlot)) : NULL;
            for (size_t i = 0; i < frame.nslots; i++) {
                frame.slots[i].pattern = patt->data.function.args[i]; /* borrowed */
                frame.slots[i].head = frame.slots[i].tail = NULL;
            }
        } else {
            frame.nslots = 1;
            frame.slots = malloc(sizeof(SowSlot));
            frame.slots[0].pattern = patt;                            /* borrowed */
            frame.slots[0].head = frame.slots[0].tail = NULL;
        }
    } else {
        /* Reap[expr] == Reap[expr, _]: synthesize a Blank[] pattern we own. */
        frame.nslots = 1;
        frame.slots = malloc(sizeof(SowSlot));
        frame.slots[0].pattern = expr_new_function(expr_new_symbol(SYM_Blank),
                                                   NULL, 0);
        frame.slots[0].head = frame.slots[0].tail = NULL;
        frame.own_pattern = true;
    }

    /* Push, evaluate the held body, pop. */
    frame.prev = g_reap_top;
    g_reap_top = &frame;
    Expr* result = evaluate(res->data.function.args[0]); /* borrows held body */
    g_reap_top = frame.prev;

    /* A Throw crossing the Reap unwinds: discard sown values, propagate. */
    if (eval_is_inflight_throw(result)) {
        reap_frame_free_contents(&frame);
        return result;
    }

    /* Build `second`. */
    Expr* second;
    if (!frame.is_list) {
        second = slot_entries_list(&frame.slots[0], frame.f);
    } else {
        Expr** slot_lists = (frame.nslots > 0)
            ? malloc(frame.nslots * sizeof(Expr*)) : NULL;
        for (size_t i = 0; i < frame.nslots; i++)
            slot_lists[i] = slot_entries_list(&frame.slots[i], frame.f);
        second = expr_new_function(expr_new_symbol(SYM_List),
                                   slot_lists, frame.nslots);
        free(slot_lists);
    }

    reap_frame_free_contents(&frame);

    /* {result, second} — result owned from evaluate(), placed directly. */
    Expr** out = malloc(2 * sizeof(Expr*));
    out[0] = result;
    out[1] = second;
    Expr* wrapped = expr_new_function(expr_new_symbol(SYM_List), out, 2);
    free(out);
    return wrapped;
}
