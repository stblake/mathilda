/*
 * memacct.c -- kernel memory accounting and MemoryConstrained.  See memacct.h
 * for the model (net live bytes, frames, trips, when a trip is actionable).
 */
#include "memacct.h"

#include <stdlib.h>
#include <string.h>

#include "common.h"
#include "eval.h"
#include "funcprog.h"      /* abort_protect_depth */
#include "message.h"
#include "print.h"
#include "sym_names.h"
#include "symtab.h"

int64_t g_mth_mem_live  = 0;
int64_t g_mth_mem_limit = INT64_MAX;
int     g_cf_tripped    = -1;
int     g_mth_mem_active = 0;

/* One active constraint.  A MEMORY frame's `limit` is absolute (live at entry
 * + budget, saturated) and `saved_eff` the effective limit before it, restored
 * on pop; a TIME frame has limit INT64_MAX and its own `deadline` (the
 * effective deadline is the earliest unspent one).  `protect_base` is the
 * AbortProtect depth at entry. */
typedef struct {
    int64_t limit;
    int64_t saved_eff;
    double  deadline;
    int     protect_base;
    bool    is_time;
    bool    spent;
} McFrame;

/* Nesting deeper than this runs unconstrained rather than failing: 256 nested
 * constraints is far past any real program. */
#define MC_MAX_FRAMES 256
static McFrame g_mc_frames[MC_MAX_FRAMES];
static int     g_mc_depth = 0;
static int     g_mem_frames = 0;     /* MEMORY frames among them */

int cf_depth(void) { return g_mc_depth; }

static McFrame* cf_push_common(void) {
    if (g_mc_depth >= MC_MAX_FRAMES) return NULL;
    McFrame* f = &g_mc_frames[g_mc_depth];
    f->limit = INT64_MAX;
    f->saved_eff = g_mth_mem_limit;
    f->deadline = 0.0;
    f->protect_base = abort_protect_depth();
    f->is_time = false;
    f->spent = false;
    return f;
}

int mc_frame_push(int64_t budget) {
    McFrame* f = cf_push_common();
    if (!f) return -1;
    if (budget < 0) budget = 0;
    int64_t live = g_mth_mem_live;
    f->limit = (live > 0 && budget > INT64_MAX - live) ? INT64_MAX : live + budget;
    if (f->limit < g_mth_mem_limit) g_mth_mem_limit = f->limit;
    g_mem_frames++;
    g_mth_mem_active = 1;
    return g_mc_depth++;
}

int cf_time_push(double deadline) {
    McFrame* f = cf_push_common();
    if (!f) return -1;
    f->is_time = true;
    f->deadline = deadline;
    return g_mc_depth++;
}

void cf_unwind_to(int depth) {
    if (depth < 0) depth = 0;
    while (g_mc_depth > depth) {
        g_mc_depth--;
        g_mth_mem_limit = g_mc_frames[g_mc_depth].saved_eff;
        if (!g_mc_frames[g_mc_depth].is_time) g_mem_frames--;
    }
    if (g_cf_tripped >= g_mc_depth) g_cf_tripped = -1;
    g_mth_mem_active = g_mem_frames > 0;
}

static void cf_trip(int i) {
    if (g_cf_tripped < 0 || i < g_cf_tripped) g_cf_tripped = i;
}

void cf_time_note(double now) {
    for (int i = 0; i < g_mc_depth; i++) {
        const McFrame* f = &g_mc_frames[i];
        if (f->is_time && !f->spent && now >= f->deadline) { cf_trip(i); return; }
    }
}

bool cf_time_deadline(double* out) {
    bool any = false;
    for (int i = 0; i < g_mc_depth; i++) {
        const McFrame* f = &g_mc_frames[i];
        if (!f->is_time || f->spent) continue;
        if (!any || f->deadline < *out) *out = f->deadline;
        any = true;
    }
    return any;
}

int cf_time_protect_base(void) {
    for (int i = 0; i < g_mc_depth; i++)
        if (g_mc_frames[i].is_time && !g_mc_frames[i].spent) return g_mc_frames[i].protect_base;
    return -1;
}

bool cf_is_time(int frame) {
    return frame >= 0 && frame < g_mc_depth && g_mc_frames[frame].is_time;
}

void cf_time_spend(int frame) {
    if (!cf_is_time(frame)) return;
    g_mc_frames[frame].spent = true;
    if (g_cf_tripped == frame) g_cf_tripped = -1;
}

/* Trip the outermost frame whose limit `live + extra` passes: the reference
 * answers with the OUTER failexpr when the outer budget is the one broken. */
static void mc_trip_with(int64_t extra) {
    int64_t need = g_mth_mem_live;
    need = (extra > 0 && need > INT64_MAX - extra) ? INT64_MAX : need + extra;
    for (int i = 0; i < g_mc_depth; i++) {
        if (need > g_mc_frames[i].limit) { cf_trip(i); return; }
    }
}

void mc_note_exceeded(void) { mc_trip_with(0); }

bool cf_abort_actionable(void) {
    if (g_cf_tripped < 0 || g_cf_tripped >= g_mc_depth) return false;
    return abort_protect_depth() <= g_mc_frames[g_cf_tripped].protect_base;
}

bool cf_poll(void) {
    if (g_mth_mem_active && g_mth_mem_live > g_mth_mem_limit) mc_note_exceeded();
    return cf_abort_actionable();
}

void cf_clear_trip(void) { g_cf_tripped = -1; }

void cf_restore_trip(int frame) {
    if (frame < 0 || frame >= g_mc_depth) return;
    cf_trip(frame);
}

bool mth_mem_request(size_t bytes) {
    if (!g_mth_mem_active) return true;
    int64_t req = (bytes > (size_t)INT64_MAX) ? INT64_MAX : (int64_t)bytes;
    int64_t live = g_mth_mem_live;
    if (live <= g_mth_mem_limit && req <= g_mth_mem_limit - live) return true;
    mc_trip_with(req);
    return !cf_abort_actionable();
}

bool mth_mem_request_shape(const Expr* spec, size_t elem_bytes) {
    if (!g_mth_mem_active) return true;
    size_t n = 1;
    const Expr* const* dims = &spec;
    size_t nd = 1;
    if (head_is(spec, SYM_List)) {
        dims = (const Expr* const*)spec->data.function.args;
        nd = spec->data.function.arg_count;
    }
    for (size_t i = 0; i < nd; i++) {
        const Expr* d = dims[i];
        if (d->type != EXPR_INTEGER || d->data.integer < 0) return true;
        size_t k = (size_t)d->data.integer;
        if (k != 0 && n > SIZE_MAX / k) n = SIZE_MAX;
        else n *= k;
    }
    if (elem_bytes != 0 && n > SIZE_MAX / elem_bytes) return mth_mem_request(SIZE_MAX);
    return mth_mem_request(n * elem_bytes);
}

/* ------------------------------------------------------------ the builtin */

/* Budget argument: 1 = a positive machine integer (*out), 0 = Infinity (no
 * constraint), -1 = anything else. */
static int mc_parse_budget(const Expr* b, int64_t* out) {
    if (b->type == EXPR_INTEGER) {
        if (b->data.integer <= 0) return -1;
        *out = b->data.integer;
        return 1;
    }
    if (b->type == EXPR_SYMBOL && b->data.symbol.name == SYM_Infinity) return 0;
    if (b->type == EXPR_FUNCTION && head_is(b, SYM_DirectedInfinity) &&
        b->data.function.arg_count == 1) {
        const Expr* s = b->data.function.args[0];
        if (s->type == EXPR_INTEGER && s->data.integer == 1) return 0;
    }
    return -1;
}

/* MemoryConstrained[expr, b] / MemoryConstrained[expr, b, failexpr].  HoldAll:
 * only b is evaluated up front; failexpr only if the constraint fires, after
 * this frame is gone (so it runs under the enclosing budget, if any). */
Expr* builtin_memory_constrained(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    if (argc < 2 || argc > 3)
        return builtin_arg_error("MemoryConstrained", argc, 2, 3);

    Expr* b = eval_and_free(expr_copy(res->data.function.args[1]));
    int64_t budget = 0;
    int kind = mc_parse_budget(b, &budget);
    expr_free(b);
    if (kind < 0) {
        char* s = expr_to_string_outputform(res);
        mth_message("MemoryConstrained", "ipnfm",
                    "Positive machine-sized integer or Infinity expected at position 2 in %s.",
                    s ? s : "MemoryConstrained[...]");
        free(s);
        return NULL;
    }
    if (kind == 0) return evaluate(res->data.function.args[0]);

    int me = mc_frame_push(budget);
    if (me < 0) return evaluate(res->data.function.args[0]);   /* stack full */

    Expr* r = evaluate(res->data.function.args[0]);

    /* A body that ends on one large allocation (2^(10^7)) may never reach
     * another poll; its result is still held, so the overrun is noted here.
     * Read the trip before popping (the pop forgets a trip of its own frame). */
    if (g_mth_mem_live > g_mth_mem_limit) mc_note_exceeded();
    int tripped = g_cf_tripped;
    cf_unwind_to(me);

    if (tripped == me) {
        /* Our budget broke: consume the abort -- whether it already became the
         * sentinel, or the body ended on a declined request before any poll. */
        if (eval_is_inflight_abort(r)) eval_unwind_consumed();
        expr_free(r);
        if (argc == 3) return evaluate(res->data.function.args[2]);
        return expr_new_symbol(SYM_DollarAborted);
    }
    if (tripped >= 0 && tripped < me) {
        /* An enclosing budget broke: hand the abort outwards. */
        if (!eval_is_inflight_abort(r)) {
            expr_free(r);
            r = expr_new_function(expr_new_symbol(SYM_Abort), NULL, 0);
        }
    }
    return r;   /* a plain Abort[] / interrupt / Throw passes through */
}

void memacct_init(void) {
    symtab_add_builtin(SYM_MemoryConstrained, builtin_memory_constrained);
    symtab_get_def(SYM_MemoryConstrained)->attributes |= ATTR_HOLDALL | ATTR_PROTECTED;
}
