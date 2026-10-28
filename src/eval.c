/*
 * eval.c
 *
 * This file implements the core evaluation engine of Mathilda.
 * It follows the "infinite evaluation" semantics of the Mathematica:
 * expressions are repeatedly transformed until they no longer change.
 *
 * The main entry point is evaluate(), which calls evaluate_step() in a loop.
 */

#include "eval.h"
#include "arithmetic.h"   /* g_arith_warnings_muted: internal-probe marker */
#include "context.h"
#include "symtab.h"
#include "ndarray.h"
#include "pack.h"
#include "core.h"
#include "purefunc.h"
#include "print.h"
#include "deriv.h"
#include "sym_names.h"
#include "sym_intern.h"
#include "assoc.h"                  /* assoc_lookup_value — O(1) <|...|>[key] */
#include "assoc_ops.h"              /* Listable threading over associations, Splice */
#include "interp.h"
#include "interval.h"                /* interval_thread_call — Interval[...] threading */
#include "compile/compiled_function.h"
#include "predict.h"   /* src/ml -- fitted models as callables */
#include "compile/autocompile.h"   /* $AutoCompilation */
#include "numloop.h"                 /* $AutoCompilation also gates numloop */
#include "plot_common.h"             /* $RaylibVerbose backing flag (raylib-free) */
#include "opform.h"                  /* h[o...][x] operator (curried) forms */
#include "numeric.h"                 /* numeric_sanitize_nonfinite -- no inf/nan leaks */
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <limits.h>
#include "message.h"   /* mth_message: Quiet/Check funnel */
#include "sparse/sparse.h" /* EXPR_SPARSE (SparseArray raw objects) */

/*
 * The maximum number of evaluation steps to prevent infinite recursion
 * in cases of circular definitions.
 */
#define MAX_ITERATIONS 4096

/* Small-arity fast path for evaluate_step's per-call scratch array. The vast
 * majority of function calls (arithmetic heads, Set, control flow, ...) have a
 * handful of arguments; sizing a stack buffer to cover them avoids a malloc/
 * free pair per function node per evaluation pass — a dominant cost in tight
 * numeric loops. Calls with more args fall back to a heap allocation. */
#define EVAL_SMALL_ARGS 8

/*
 * $RecursionLimit guard. The REPL is single-threaded so a static counter
 * suffices. Each call to evaluate() bumps eval_recursion_depth on entry and
 * decrements on exit; if the depth would exceed eval_recursion_limit we
 * return the expression wrapped in Hold[] (so it cannot re-enter the
 * evaluator) and emit a $RecursionLimit::reclim message.
 *
 * The default of 1024 matches modern Mathematica and leaves comfortable
 * headroom under the typical 8 MB thread stack while still catching
 * pathological recursion (e.g. a self-referential rule like
 * f[x_] := f[x] + 1).
 *
 * Minimum enforced at 20 (Mathematica-compatible) so users cannot brick
 * the evaluator by setting a value below the depth its own bookkeeping
 * needs.
 */
#define DEFAULT_RECURSION_LIMIT 1024
#define MIN_RECURSION_LIMIT     20
static int  eval_recursion_depth = 0;
static int  eval_recursion_limit = DEFAULT_RECURSION_LIMIT;
/* Once the recursion limit is hit anywhere in the evaluation tree, this
 * sticky flag tells *every* enclosing evaluate() loop to stop iterating
 * and return its current value. Without it, the wrap-in-Hold result
 * causes outer fixed-point loops to perceive endless "progress" (each
 * iteration adds another Hold wrapper) and chew through all 4096
 * outer iterations at every level of the unwind. The flag is cleared
 * at the top of each top-level evaluate() call. */
static bool eval_overflow = false;

/* Bumped every time the packed-list gate (step 2.7) materialises an argument.
 * Monotonic and never reset: evaluate()'s fixed-point test brackets each
 * evaluate_step with it, because a step whose only effect was materialising a
 * buffer is invisible to expr_eq. See the use site for why over-reporting is
 * free and under-reporting would not be. */
static uint64_t g_pack_gate_ticks = 0;

int  eval_get_recursion_limit(void) { return eval_recursion_limit; }
int  eval_get_recursion_depth(void) { return eval_recursion_depth; }
void eval_set_recursion_limit(int n) {
    eval_recursion_limit = (n >= MIN_RECURSION_LIMIT) ? n : DEFAULT_RECURSION_LIMIT;
}

/* See eval.h.  Used after a siglongjmp out of evaluate(): the matching
 * decrements never ran, so we restore the depth counter to the value it
 * had before the timed call entered the evaluator.  We also clear
 * eval_overflow so a future evaluate() is not falsely poisoned by the
 * aborted call. */
void eval_unwind_consumed(void);

void eval_reset_recursion_depth(int n) {
    if (n < 0) n = 0;
    eval_recursion_depth = n;
    eval_overflow = false;
    eval_unwind_consumed();   /* a parked sentinel from the abandoned frames */
}

/* M3 phase-3 evaluation clock. Starts at 1 so a freshly-allocated Expr
 * (last_evaluated_at == 0) is never mistaken for "already evaluated".
 * Bumped by symtab.c (add_rule, symtab_clear_symbol) and attr.c
 * (set_attributes, add/remove_single_attribute). 64 bits is enough to
 * absorb roughly 2^64 mutations, far beyond any practical session. */
static uint64_t g_eval_clock = 1;
/* INERT fixed-point epoch: the clock value of the most recent mutation that
 * could change the evaluation of an INERT node (see node_compute_inert).
 * Every clock bump advances it -- conservatively, since a bump site may stand
 * for any global change -- EXCEPT an OwnValue binding of a symbol no INERT
 * stamp relied on (eval_clock_bump_binding), which is what an iterator or a
 * Module local is. */
static uint64_t g_inert_epoch = 1;
uint64_t eval_clock_get(void) { return g_eval_clock; }
void     eval_clock_bump(void) { g_eval_clock++; g_inert_epoch = g_eval_clock; }
static inline bool ground_head_name(const char* n);
static inline bool inert_builtin_head(const char* n);
/* A binding or a DownValue for `def`: an INERT node can depend on it only if
 * the node relied on `def` (flagged inert_observed) or `def` is one of the
 * builtin heads an INERT node may carry, whose rules are not flagged. */
void     eval_clock_bump_binding(SymbolDef* def) {
    g_eval_clock++;
    if (!def) { g_inert_epoch = g_eval_clock; return; }
    if (def->inert_observed || inert_builtin_head(def->symbol_name) ||
        ground_head_name(def->symbol_name)) {
        g_inert_epoch = g_eval_clock;
        def->inert_observed = 0;
    }
}

/* --- Top-level evaluation id ----------------------------------------------
 * Bumped once at the entry of each OUTERMOST evaluate() call (i.e. once per
 * user command / script statement). Unlike the eval clock, it does NOT change
 * across the nested evaluations, fixed-point iterations, or symbol-table churn
 * that happen while a single top-level expression is being reduced. This gives
 * builtins a stable "am I still inside the same top-level evaluation?" token,
 * used to scope per-command state (e.g. Integrate's fail-memo and the
 * once-per-command Integrate::nonelem diagnostic) so it self-invalidates at the
 * next command without any explicit reset. */
static uint64_t g_toplevel_eval_id = 0;
uint64_t eval_toplevel_id(void) { return g_toplevel_eval_id; }

/* --- Ground fixed-point epoch (loop-invariant re-evaluation) ---------------
 * The eval clock is a single global epoch: ANY symbol-table mutation bumps it
 * and invalidates every cached fixed point. That is correct but coarse -- a
 * Do/Table/Fold loop rebinds its iterator (an OwnValue write) every iteration,
 * so a large loop-invariant value bound to a symbol gets fully re-canonicalised
 * O(size) each step even though nothing it depends on changed.
 *
 * `g_last_rule_change_clock` is a SECOND, finer epoch: the clock value at the
 * most recent mutation that can change how a *head* evaluates -- a DownValue
 * add, an attribute/Protect change, a Clear/Remove. Ordinary OwnValue bindings
 * (iterator variables, `x = 5`, numeric temp-bindings) bump only the eval clock
 * and leave this one alone. A GROUND node (see is_ground_* below: a fixed point
 * built solely from literals under the six pure structural constructors) can be
 * re-validated as a fixed point whenever `stamp >= g_last_rule_change_clock`,
 * regardless of the eval clock -- because it references none of the mutable
 * state an OwnValue binding touches. This is what lifts a loop-invariant
 * association/list from O(size)-per-iteration back to O(1). It only ever
 * advances (monotone), so a ground short-circuit is never taken stale. */
static uint64_t g_last_rule_change_clock = 1;
uint64_t eval_rule_epoch_get(void)  { return g_last_rule_change_clock; }
/* Mark the CURRENT clock as a rule change (caller already bumped the clock). */
void     eval_rule_epoch_mark(void) { g_last_rule_change_clock = g_eval_clock; }
/* Bump the clock AND record it as a rule change (attribute/Protect sites). */
void     eval_rule_epoch_bump(void) {
    g_eval_clock++;
    g_last_rule_change_clock = g_eval_clock;
    g_inert_epoch = g_eval_clock;
}

/* The top bit of `last_evaluated_at` is a benign GROUND flag: the clock is a
 * monotone counter that will never reach 2^63, so this steals no range. The
 * low 63 bits remain the fixed-point stamp; every comparison against the eval
 * clock masks the flag off first. See is_ground_now() / node_compute_ground(). */
#define EVAL_GROUND_BIT  (UINT64_C(1) << 63)
/* Bit 62: the INERT flag (see node_compute_inert).  The clock never reaches
 * 2^62 either. */
#define EVAL_INERT_BIT   (UINT64_C(1) << 62)
#define EVAL_STAMP_MASK  (~(EVAL_GROUND_BIT | EVAL_INERT_BIT))
static inline uint64_t eval_stamp_of(const Expr* e) {
    return e->last_evaluated_at & EVAL_STAMP_MASK;
}
static inline bool eval_ground_of(const Expr* e) {
    return (e->last_evaluated_at & EVAL_GROUND_BIT) != 0;
}
/* A FUNCTION node whose GROUND bit is set is a valid fixed point iff no rule
 * change has occurred since it was stamped -- the straddle-safe predicate used
 * identically at short-circuit time and when a parent consumes a child's bit. */
static inline bool eval_ground_valid(const Expr* e) {
    return eval_ground_of(e) && eval_stamp_of(e) >= g_last_rule_change_clock;
}
/* An INERT node is a valid fixed point while nothing it relied on has changed
 * since it was stamped.  The global rule epoch is NOT consulted: a DownValue
 * advances g_inert_epoch itself when it can matter (eval_clock_bump_binding),
 * so storing a result -- res[14] = <135135-term sum> -- does not invalidate
 * every INERT node, the stored value's own stamp included.  Attribute changes,
 * Clear/Remove and every other conservative bump still advance it. */
static inline bool eval_inert_valid(const Expr* e) {
    if (!(e->last_evaluated_at & EVAL_INERT_BIT)) return false;
    return eval_stamp_of(e) >= g_inert_epoch;
}
/* A FUNCTION node is a re-usable fixed point if it was stamped under the live
 * clock (exact hit) OR it is a still-valid ground node (survives OwnValue churn)
 * OR a still-valid inert one (survives the churn of unrelated bindings). */
static inline bool eval_fixed_point_reusable(const Expr* e) {
    return eval_stamp_of(e) == g_eval_clock || eval_ground_valid(e) ||
           eval_inert_valid(e);
}

/* Public, mask-aware accessors for tests (the raw field now carries the flag). */
uint64_t eval_node_stamp(const Expr* e)   { return e ? eval_stamp_of(e) : 0; }
bool     eval_node_is_ground(const Expr* e) { return e ? eval_ground_of(e) : false; }

/* The pure structural constructors (six generic ones plus the graph value
 * heads). Their canonical form is a total,
 * side-effect-free function of their arguments -- they read no mutable global
 * state -- so a fixed point built only from these heads over literal leaves is
 * immutable until one of the heads is itself redefined (which advances the rule
 * epoch). Heads WITH a value-computing builtin (Plus, Sin, RandomReal, ...) are
 * deliberately excluded: they never reach a stamped fixed point AS themselves
 * when reducible, but even a symbolic residue could in principle depend on
 * global state, so we do not trust them. */
static inline bool ground_head(const Expr* h) {
    if (!h || h->type != EXPR_SYMBOL) return false;
    return ground_head_name(h->data.symbol.name);
}
static inline bool ground_head_name(const char* n) {
    return n == SYM_List || n == SYM_Association || n == SYM_Rule
        || n == SYM_RuleDelayed || n == SYM_Complex || n == SYM_Rational
        /* Graph values. DirectedEdge/UndirectedEdge are inert; Graph's builtin
         * is a pure function of its arguments (it validates, and its memo is a
         * cache of that function), so an evaluated graph stays a fixed point
         * across OwnValue churn. Without these, every assignment made the next
         * use of a stored 5x10^5-edge graph re-walk it (~34 ms). A weighted
         * graph's bare EdgeWeight symbol is not ground, so it is not marked. */
        || n == SYM_Graph || n == SYM_DirectedEdge || n == SYM_UndirectedEdge;
}
/* Is `a` ground *right now*? For a FUNCTION we trust its cached bit only if it
 * is still valid (eval_ground_valid); atoms are decided structurally. This is
 * the recurrence used bottom-up when stamping a parent -- O(arity), not O(size). */
static inline bool is_ground_now(const Expr* a) {
    switch (a->type) {
        case EXPR_INTEGER: case EXPR_REAL: case EXPR_BIGINT: case EXPR_STRING:
            return true;
        case EXPR_FUNCTION:
            return eval_ground_valid(a);
        default:            /* bare SYMBOL, NDARRAY, COMPILED, MPFR: conservative */
            return false;
    }
}
/* Compute the GROUND bit for a FUNCTION node reaching a fixed point: a
 * whitelisted head and every argument ground. Non-FUNCTION nodes are never
 * marked (their stamp is never read by the short-circuits). */
static bool node_compute_ground(const Expr* e) {
    if (e->type != EXPR_FUNCTION) return false;
    if (!ground_head(e->data.function.head)) return false;
    size_t n = e->data.function.arg_count;
    Expr* const* args = e->data.function.args;
    for (size_t i = 0; i < n; i++) {
        if (!args[i] || !is_ground_now(args[i])) return false;
    }
    return true;
}

/* --- INERT fixed points ----------------------------------------------------
 * GROUND covers literal data only.  A symbolic fixed point -- the 135135-term
 * sum of g[x[i], x[j]] products a Dirac trace builds -- was re-canonicalised in
 * full after EVERY iterator or Module-local binding, because any OwnValue write
 * bumps the eval clock.  Its value cannot depend on such a binding unless the
 * bound symbol occurs in it: INERT records exactly that.
 *
 * A FUNCTION fixed point is INERT when its head is a pure arithmetic or
 * structural builtin (Plus, Times, Power, or a ground head) or a symbol with no
 * value, rules or builtin at all, and every argument is a number, a string, a
 * symbol with no value and no builtin, or a still-valid INERT/GROUND node.
 * Plus/Times/Power read no mutable global state (the trig-canon suppression
 * flag bumps the clock when it lifts).  Every symbol relied on -- leaf or head
 * -- is flagged inert_observed; an OwnValue binding of a flagged symbol
 * advances g_inert_epoch (eval_clock_bump_binding), and every other clock bump
 * and every rule change advances it unconditionally, so a stale INERT stamp is
 * never reused. */
static inline bool inert_builtin_head(const char* n) {
    return n == SYM_Plus || n == SYM_Times || n == SYM_Power;
}
static inline SymbolDef* inert_sym_def(Expr* s) {
    SymbolDef* d = s->data.symbol.def;
    if (!d) { d = symtab_get_def(s->data.symbol.name); s->data.symbol.def = d; }
    return d;
}
/* A symbol with no definitions of any kind: evaluates to itself as a leaf and
 * leaves a call it heads to the evaluator's generic argument processing. */
static inline bool inert_bare_symbol(const SymbolDef* d) {
    return d && !d->own_values && !d->down_values && !d->sub_values &&
#if UP_VALUES
           !d->up_values &&
#endif
           !d->builtin_func;
}
static inline bool is_inert_now(Expr* a) {
    switch (a->type) {
        case EXPR_INTEGER: case EXPR_REAL: case EXPR_BIGINT: case EXPR_STRING:
            return true;
        case EXPR_SYMBOL: {
            SymbolDef* d = inert_sym_def(a);
            if (!d || d->own_values || d->builtin_func) return false;
            d->inert_observed = 1;
            return true;
        }
        case EXPR_FUNCTION:
            return eval_inert_valid(a) || eval_ground_valid(a);
        default:
            return false;
    }
}
static bool node_compute_inert(Expr* e) {
    if (e->type != EXPR_FUNCTION) return false;
    Expr* h = e->data.function.head;
    if (!h || h->type != EXPR_SYMBOL) return false;
    const char* hn = h->data.symbol.name;
    if (!inert_builtin_head(hn) && !ground_head(h)) {
        SymbolDef* d = inert_sym_def(h);
        if (!inert_bare_symbol(d)) return false;
        d->inert_observed = 1;
    }
    size_t n = e->data.function.arg_count;
    Expr* const* args = e->data.function.args;
    for (size_t i = 0; i < n; i++) {
        if (!args[i] || !is_inert_now(args[i])) return false;
    }
    return true;
}

/* Stamp a List-of-literals tree a builtin has just built as an evaluated
 * fixed point, so the evaluator's re-walk of the returned value is O(1)
 * instead of O(size). Bottom-up, exactly the stamp the fixed-point exit
 * writes: each List node whose arguments are all ground gets the live clock
 * plus the GROUND bit; anything else is left unstamped and evaluates as usual.
 *
 * Sound because such a node IS a fixed point -- unless List itself carries
 * rules, in which case nothing is stamped. A rule added later advances the
 * rule epoch past these stamps and invalidates them, as for any GROUND node.
 * Returns whether `e` ended up stamped. */
static bool stamp_list_literal_rec(Expr* e) {
    if (e->type != EXPR_FUNCTION || !e->data.function.head ||
        e->data.function.head->type != EXPR_SYMBOL ||
        e->data.function.head->data.symbol.name != SYM_List) return false;
    bool all = true;
    Expr* const* args = e->data.function.args;
    for (size_t i = 0, n = e->data.function.arg_count; i < n; i++) {
        const Expr* a = args[i];
        if (!a) { all = false; continue; }
        switch (a->type) {           /* atoms inline: no call per element */
            case EXPR_INTEGER: case EXPR_REAL: case EXPR_BIGINT: case EXPR_STRING:
                break;
            case EXPR_FUNCTION:
                if (!stamp_list_literal_rec(args[i])) all = false;
                break;
            default:
                all = false;
        }
    }
    if (all) e->last_evaluated_at = g_eval_clock | EVAL_GROUND_BIT;
    return all;
}

bool eval_node_reusable(const Expr* e) {
    if (!e) return false;
    if (e->type == EXPR_FUNCTION) return eval_fixed_point_reusable(e);
    if (e->type == EXPR_SYMBOL) {
        /* A symbol with a value is not final: a builtin may be called from C
         * on raw arguments, not only by the evaluator on evaluated ones. */
        SymbolDef* d = inert_sym_def((Expr*)e);
        return d && !d->own_values;
    }
    return e->type == EXPR_INTEGER || e->type == EXPR_REAL || e->type == EXPR_BIGINT ||
           e->type == EXPR_STRING;
}

void eval_mark_fixed_point(Expr* e) {
    if (!e || e->type != EXPR_FUNCTION) return;
    uint64_t flag = node_compute_ground(e) ? EVAL_GROUND_BIT
                  : node_compute_inert(e)  ? EVAL_INERT_BIT : 0;
    e->last_evaluated_at = g_eval_clock | flag;
}

bool eval_stamp_list_literal(Expr* e) {
    if (!e) return false;
    SymbolDef* d = symtab_get_def(SYM_List);
    if (d && (d->down_values || d->up_values || d->own_values)) return false;
    return stamp_list_literal_rec(e);
}

/* ---- Trace collector (nested) --------------------------------------------
 * Trace[expr] returns a list that mirrors the *structure* of expr's
 * evaluation: each argument sub-evaluation that takes >=1 step appears as a
 * nested sublist, and the reassembled intermediate form (f[evaluated_args],
 * before the head's rule fires) appears as a step. This is NOT a change to
 * evaluation semantics -- the evaluator already performs every sub-evaluation;
 * the collector merely *observes* more of it than the earlier flat v1 did.
 *
 * Mechanism: a stack of frames, one per active evaluate() call. Nesting then
 * follows the evaluator's own recursion automatically -- no depth arithmetic.
 *   - evaluate() pushes a frame on entry (when g_tracing) and pops it on exit.
 *   - On pop, a frame that recorded >=1 entry is turned into a List and spliced
 *     as ONE nested entry into its parent frame; an empty frame (atom / no-op
 *     sub-evaluation) contributes nothing and is discarded. The outermost
 *     frame's List becomes g_trace_root_result, which eval_collect_trace hands
 *     back.
 *   - The fixed-point loop records each real current->next rewrite into the
 *     current (top) frame. The "before" form it records is the reassembled
 *     f[evaluated_args] snapshot (pending_reassembled, set by evaluate_step)
 *     when present, else `current` itself (symbol rewrites, atoms).
 *
 * Reentrancy: eval_collect_trace saves/restores g_trace_top / g_trace_root_result
 * / g_tracing on the C stack, so a nested Trace runs on a fresh stack and
 * appears to the outer trace as a single already-reduced value.
 *
 * Ownership: frame entries are owned by the frame (each a fresh expr_copy or a
 * spliced child List) until the frame is turned into a List via
 * expr_new_function (which takes the elements) -- exactly the earlier v1
 * contract. Cost when not tracing: one predicted-false bool read (g_tracing)
 * per evaluate() call, matching the old g_trace_active check. */
typedef struct TraceFrame {
    Expr**  entries;             /* owned; each a recorded form OR nested List */
    size_t  count, cap;
    Expr*   pending_reassembled; /* f[evaluated_args] snapshot for this step    */
    struct TraceFrame* parent;
} TraceFrame;

static TraceFrame* g_trace_top         = NULL;  /* current frame; NULL = none  */
static Expr*       g_trace_root_result = NULL;  /* outermost frame's List      */
static bool        g_tracing           = false; /* true only inside a Trace[]  */

/* Suppression depth: while > 0, evaluate() neither pushes frames nor records.
 * Bumped around machinery whose internal evaluate() calls are NOT user-visible
 * argument sub-evaluations -- a builtin computing its result (e.g. Range folding
 * i+1) and Listable threading evaluating the threaded elements. Argument and
 * head evaluation happen *before* these bumps, so they are still traced; only
 * the rule's own internal reductions are hidden, matching Mathematica (Range[10]
 * and x^{1..10} each appear as a single rewrite, not a decomposed one). */
static int         g_trace_suppress    = 0;

/* How many non-transparent C builtins are executing (see eval_transparent in
 * symtab.h).  0 at the top level; 1 inside a builtin called from user code. */
static int g_builtin_depth = 0;
/* The head whose builtin the evaluator is running (innermost dispatch). */
static const char* g_dispatch_head = NULL;
int  eval_builtin_depth(void) { return g_builtin_depth; }
const char* eval_dispatch_head(void) { return g_dispatch_head; }
int* eval_builtin_depth_counter(void) { return &g_builtin_depth; }
/* The recursion depth of the evaluate() frame whose current step came to rest
 * without any rewrite (set at evaluate_step_core's rest return, read by the
 * evaluate loop right after the step); -1 otherwise. */
static int         g_step_settled_depth = -1;
#define TRACE_ACTIVE() (g_tracing && g_trace_suppress == 0)

/* Append an already-owned entry to a frame, growing on demand. On OOM the
 * entry is freed rather than leaked (the step is dropped). */
static void frame_append(TraceFrame* f, Expr* owned) {
    if (f->count == f->cap) {
        size_t newcap = f->cap ? f->cap * 2 : 8;
        Expr** grown = realloc(f->entries, newcap * sizeof(Expr*));
        if (!grown) { expr_free(owned); return; }
        f->entries = grown;
        f->cap = newcap;
    }
    f->entries[f->count++] = owned;
}

/* Record a form (borrowed) into the current top frame, with consecutive-dedup:
 * a form structurally equal to the last recorded entry is skipped. Dedup turns
 * a symbol chain a->b->c->1 into {a,b,c,1} and prevents a doubled final form
 * when the reassembled form equals the step result (e.g. f[1+1] -> f[2]). */
static void frame_record(Expr* form) {
    TraceFrame* f = g_trace_top;
    if (!f) return;
    if (f->count > 0 && expr_eq(f->entries[f->count - 1], form)) return;
    frame_append(f, expr_copy(form));
}

/* Push a fresh frame for the evaluate() call now entering. */
static void trace_frame_push(void) {
    TraceFrame* f = calloc(1, sizeof(TraceFrame));
    if (!f) return;                 /* OOM: skip this frame level              */
    f->parent = g_trace_top;
    g_trace_top = f;
}

/* Pop the current frame as the evaluate() call exits: an empty frame is
 * discarded; a non-empty one becomes a List spliced into the parent (or, at the
 * outermost level, stored as g_trace_root_result). */
static void trace_frame_pop(void) {
    TraceFrame* f = g_trace_top;
    if (!f) return;
    g_trace_top = f->parent;
    if (f->pending_reassembled) expr_free(f->pending_reassembled);

    if (f->count == 0) {            /* no steps: contribute nothing            */
        free(f->entries);
        free(f);
        return;
    }
    Expr* list = expr_new_function(expr_new_symbol(SYM_List), f->entries, f->count);
    free(f->entries);
    free(f);
    if (g_trace_top) {
        frame_append(g_trace_top, list);          /* nested sublist            */
    } else {
        if (g_trace_root_result) expr_free(g_trace_root_result);
        g_trace_root_result = list;               /* outermost trace           */
    }
}

/* Discard any reassembled snapshot left on the current frame (a step that did
 * not end up rewriting). Called on the fixed-point / loop-exit paths. */
static void trace_clear_pending(void) {
    if (g_trace_top && g_trace_top->pending_reassembled) {
        expr_free(g_trace_top->pending_reassembled);
        g_trace_top->pending_reassembled = NULL;
    }
}

/*
 * eval_classify_return:
 * See the contract in eval.h. Pointer-equality on interned symbols is
 * the dispatch primitive: every Expr_Symbol's `data.symbol.name` field is the
 * canonical interned pointer (sym_intern.c), so checks like
 * `head->data.symbol.name == SYM_Return` and `target->data.symbol.name ==
 * boundary_head` are O(1) and never strcmp.
 *
 * Care is taken to keep this side-effect free: no eval_clock_bump,
 * no expr_free, no allocation when the answer is NONE/PROPAGATE. The
 * single allocation on CONSUME is either a fresh `Null` symbol (for the
 * 0-arg form) or an expr_copy of args[0] (for 1-arg / 2-arg forms).
 * That copy is necessary because the caller will free `e` after
 * yielding the value.
 */
bool eval_is_inflight_throw(const Expr* e) {
    return e && e->type == EXPR_FUNCTION &&
           e->data.function.head->type == EXPR_SYMBOL &&
           e->data.function.head->data.symbol.name == SYM_Throw &&
           e->data.function.arg_count >= 1 && e->data.function.arg_count <= 3;
}

/* Abort[]: the in-flight abort.  Like a Throw it is its own sentinel and rides
 * the ordinary return path (so every frame's cleanup runs), but no Catch
 * consumes it -- only CheckAbort, or the top level, which yields $Aborted. */
bool eval_is_inflight_abort(const Expr* e) {
    return e && e->type == EXPR_FUNCTION &&
           e->data.function.head->type == EXPR_SYMBOL &&
           e->data.function.head->data.symbol.name == SYM_Abort &&
           e->data.function.arg_count == 0;
}

bool eval_is_inflight_unwind(const Expr* e) {
    return eval_is_inflight_throw(e) || eval_is_inflight_abort(e);
}

/* Set when an abort reaches the top level (see eval_take_top_abort). */
static bool g_top_abort = false;

bool eval_take_top_abort(void) {
    bool a = g_top_abort;
    g_top_abort = false;
    return a;
}

bool eval_is_inflight_goto(const Expr* e) {
    return e && e->type == EXPR_FUNCTION &&
           e->data.function.head->type == EXPR_SYMBOL &&
           e->data.function.head->data.symbol.name == SYM_Goto &&
           e->data.function.arg_count == 1;
}

bool eval_is_inflight_break_continue(const Expr* e) {
    return e && e->type == EXPR_FUNCTION &&
           e->data.function.head->type == EXPR_SYMBOL &&
           (e->data.function.head->data.symbol.name == SYM_Break ||
            e->data.function.head->data.symbol.name == SYM_Continue) &&
           e->data.function.arg_count == 0;
}

EvalReturnAction eval_classify_return(Expr* e,
                                      const char* boundary_head,
                                      Expr** out_value) {
    if (out_value) *out_value = NULL;
    if (!e) return EVAL_RETURN_NONE;
    if (e->type != EXPR_FUNCTION) return EVAL_RETURN_NONE;
    if (e->data.function.head->type != EXPR_SYMBOL) return EVAL_RETURN_NONE;
    if (e->data.function.head->data.symbol.name != SYM_Return) return EVAL_RETURN_NONE;

    size_t argc = e->data.function.arg_count;

    /* Return[]: yield Null at the nearest boundary. */
    if (argc == 0) {
        if (out_value) *out_value = expr_new_symbol(SYM_Null);
        return EVAL_RETURN_CONSUME;
    }

    /* Return[expr]: yield expr at the nearest boundary, regardless of
     * boundary_head. */
    if (argc == 1) {
        if (out_value) *out_value = expr_copy(e->data.function.args[0]);
        return EVAL_RETURN_CONSUME;
    }

    /* Return[expr, h, ...]: target the nearest boundary whose head is h.
     * Extra arguments are ignored (Mathematica accepts up to 2 args; we
     * accept the same and leave further args to a caller-side message
     * if the user supplies them).
     *
     * The target must be an EXPR_SYMBOL for the comparison to be
     * meaningful. If it isn't, fall back to PROPAGATE so that the
     * marker isn't accidentally consumed by an arbitrary boundary. */
    Expr* target = e->data.function.args[1];
    if (boundary_head &&
        target->type == EXPR_SYMBOL &&
        target->data.symbol.name == boundary_head) {
        if (out_value) *out_value = expr_copy(e->data.function.args[0]);
        return EVAL_RETURN_CONSUME;
    }
    return EVAL_RETURN_PROPAGATE;
}

/* ---------------------------------------------------------------------------
 * Boolean system variables that mirror a C-side flag.
 *
 * Both switch off an optimisation that is invisible by construction -- the
 * compiled and packed paths are contracted to give the interpreter's answer --
 * so turning one off must change speed and nothing else. That is exactly what
 * makes them worth having: a differential run flips one and diffs every output,
 * and a user who suspects a fast path has it wrong can confirm in one line.
 *
 * A table rather than a strcmp chain, so the whole set is visible at once and
 * adding one is a line. $RecursionLimit keeps its own hook below: it carries an
 * integer with validation and a roll-back, not a boolean.
 * ------------------------------------------------------------------------- */
typedef struct {
    const char* name;
    void (*set)(bool);
    bool (*get)(void);
} EvalSysFlag;

/* $AutoCompilation covers BOTH compiled paths, because a user turning it off is
 * asking for the interpreter and does not care which of the two internal
 * mechanisms would have compiled their loop: the autocompile adapter (Plot,
 * Table, NIntegrate, ...) and numloop (Do/For/While/Map/Nest bodies). */
static void eval_set_autocompilation(bool on) {
    autocompile_set_enabled(on);
    numloop_set_enabled(on);
}
static bool eval_get_autocompilation(void) { return autocompile_enabled(); }

static const EvalSysFlag EVAL_SYSFLAGS[] = {
    { "$AutoCompilation",  eval_set_autocompilation, eval_get_autocompilation },
    { "$AutoArrayPacking", pack_set_enabled,         pack_enabled },
    { "$RaylibVerbose",    raylib_verbose_set,       raylib_verbose_enabled },
};
#define EVAL_N_SYSFLAGS ((size_t)(sizeof(EVAL_SYSFLAGS) / sizeof(EVAL_SYSFLAGS[0])))

/* True when `name` is one of the boolean system variables above. Separate from
 * the sync below so a caller can decide WITHOUT building a probe value: the
 * assignment hook is reached by every `$`-prefixed symbol, and the REPL hooks
 * ($Pre, $PreRead, $Post, $PrePrint, $Epilog) all live in that namespace with
 * held right-hand sides. Evaluating a probe for those was a regression --
 * repl_hooks_tests caught it -- so the name is checked first and only a real flag
 * ever gets one. */
static bool eval_is_sysflag(const char* name) {
    for (size_t i = 0; i < EVAL_N_SYSFLAGS; i++)
        if (strcmp(name, EVAL_SYSFLAGS[i].name) == 0) return true;
    return false;
}

/* Push a candidate value into the matching C flag. Returns false when `name` is
 * not one of these variables (the caller then does nothing).
 *
 * Only True and False are accepted. Anything else is rejected with a ::flagset
 * message and the OwnValue is rolled back to the live C state, so the symbol
 * never lies about which path is running -- the same discipline
 * $RecursionLimit::limset follows. */
static bool eval_sync_sysflag(const char* name, Expr* value) {
    for (size_t i = 0; i < EVAL_N_SYSFLAGS; i++) {
        if (strcmp(name, EVAL_SYSFLAGS[i].name) != 0) continue;
        bool is_true  = value && value->type == EXPR_SYMBOL &&
                        value->data.symbol.name == SYM_True;
        bool is_false = value && value->type == EXPR_SYMBOL &&
                        value->data.symbol.name == SYM_False;
        if (!is_true && !is_false) {
            mth_message(EVAL_SYSFLAGS[i].name, "flagset",
                    "%s can only be set to True or False.", EVAL_SYSFLAGS[i].name);
            Expr* sym  = expr_new_symbol(EVAL_SYSFLAGS[i].name);
            Expr* curr = expr_new_symbol(EVAL_SYSFLAGS[i].get() ? SYM_True : SYM_False);
            symtab_add_own_value(EVAL_SYSFLAGS[i].name, sym, curr);
            expr_free(sym);
            expr_free(curr);
            return true;
        }
        EVAL_SYSFLAGS[i].set(is_true);
        return true;
    }
    return false;
}

/* Register the boolean system variables with their live defaults as OwnValues, so
 * the user can read them back as well as assign. Called from eval_init, i.e.
 * AFTER the environment overrides (MATHILDA_NO_PACK, MATHILDA_NO_AUTOCOMPILE)
 * have been read -- so `$AutoArrayPacking` reports False in a run started with
 * MATHILDA_NO_PACK=1 rather than claiming True and being wrong. */
static void eval_init_sysflags(void) {
    for (size_t i = 0; i < EVAL_N_SYSFLAGS; i++) {
        Expr* sym = expr_new_symbol(EVAL_SYSFLAGS[i].name);
        Expr* val = expr_new_symbol(EVAL_SYSFLAGS[i].get() ? SYM_True : SYM_False);
        symtab_add_own_value(EVAL_SYSFLAGS[i].name, sym, val);
        expr_free(sym);
        expr_free(val);
    }
}

/*
 * eval_init:
 * Registers the user-visible $RecursionLimit symbol with its default value
 * as an OwnValue so that the user can read or assign to it from the REPL.
 * The C-side state is kept in sync via the hook in apply_assignment.
 *
 * Must be called after symtab_init().
 */
void eval_init(void) {
    eval_init_sysflags();

    Expr* sym = expr_new_symbol(SYM_DollarRecursionLimit);
    Expr* val = expr_new_integer(eval_recursion_limit);
    symtab_add_own_value("$RecursionLimit", sym, val);
    expr_free(sym);
    expr_free(val);

    symtab_set_docstring("$RecursionLimit",
        "$RecursionLimit\n"
        "\tgives the maximum length of the evaluation stack -- the maximum\n"
        "\tnumber of nested invocations of the evaluator that can occur.\n"
        "\n"
        "Assigning a positive integer N (>= 20) updates the limit; smaller\n"
        "values are rejected with a $RecursionLimit::limset message.");
}

/*
 * sync_recursion_limit_from_value:
 * Inspect a candidate value (typically the RHS of $RecursionLimit = ...)
 * and, if it is a positive integer >= MIN_RECURSION_LIMIT, push it into
 * the C-level limit. Otherwise emit a $RecursionLimit::limset message and
 * leave the C state untouched. Bigints are clamped to INT_MAX.
 */
static void sync_recursion_limit_from_value(Expr* value) {
    long n = -1;
    if (value->type == EXPR_INTEGER) {
        n = (long)value->data.integer;
    } else if (value->type == EXPR_BIGINT) {
        /* Anything large enough not to fit in a long is far beyond any
         * useful recursion limit; treat it as "huge and acceptable". */
        if (mpz_fits_slong_p(value->data.bigint)) {
            n = mpz_get_si(value->data.bigint);
        } else if (mpz_sgn(value->data.bigint) > 0) {
            n = (long)INT_MAX;
        }
    }

    if (n < MIN_RECURSION_LIMIT) {
        mth_message("$RecursionLimit", "limset", "Cannot set $RecursionLimit to a value below %d.", MIN_RECURSION_LIMIT);
        /* Restore the OwnValue to the current C-side limit so the symbol
         * does not lie about the active value. */
        Expr* sym  = expr_new_symbol(SYM_DollarRecursionLimit);
        Expr* curr = expr_new_integer(eval_recursion_limit);
        symtab_add_own_value("$RecursionLimit", sym, curr);
        expr_free(sym);
        expr_free(curr);
        return;
    }
    if (n > INT_MAX) n = INT_MAX;
    eval_set_recursion_limit((int)n);
}

/*
 * eval_compare_expr_ptrs:
 * Helper for sorting expression arguments when a head has the Orderless attribute.
 * Uses the canonical expr_compare to ensure a stable, deterministic order.
 */
int eval_compare_expr_ptrs(const void* a, const void* b) {
    Expr* ea = *(Expr**)a;
    Expr* eb = *(Expr**)b;
    return expr_compare(ea, eb);
}

/*
 * flatten_args:
 * Implements the Flat (associative) attribute.
 * If a function has the same head as some of its arguments, those arguments
 * are "unwrapped" and their elements are promoted to be direct arguments 
 * of the parent function.
 * Example: f[a, f[b, c], d] -> f[a, b, c, d]
 */
/* Returns true iff the call actually flattened nested same-head children
 * (i.e. produced a structurally different argument list). When false,
 * `e` is byte-for-byte unchanged and the §3.4 fixed-point detector can
 * count this step as a no-op. */
/* Core flatten, assuming `head_name` is ALREADY the interned canonical
 * pointer. The hot evaluator path (which holds an EXPR_SYMBOL's name, always
 * interned) calls this directly to skip a per-call hash on every Flat head. */
static bool eval_flatten_args_interned(Expr* e, const char* head_name) {
    size_t new_count = 0;
    bool needs_flattening = false;

    for (size_t i = 0; i < e->data.function.arg_count; i++) {
        Expr* arg = e->data.function.args[i];
        if (arg->type == EXPR_FUNCTION && arg->data.function.head->type == EXPR_SYMBOL &&
            arg->data.function.head->data.symbol.name == head_name) {
            new_count += arg->data.function.arg_count;
            needs_flattening = true;
        } else {
            new_count++;
        }
    }

    /* If no nested occurrences of the head were found, we are done */
    if (!needs_flattening) return false;

    /* Second pass: allocate new argument array and copy elements */
    Expr** new_args = malloc(sizeof(Expr*) * new_count);
    size_t idx = 0;
    for (size_t i = 0; i < e->data.function.arg_count; i++) {
        Expr* arg = e->data.function.args[i];
        if (arg->type == EXPR_FUNCTION && arg->data.function.head->type == EXPR_SYMBOL &&
            arg->data.function.head->data.symbol.name == head_name) {
            /* Splat nested arguments into the new array */
            for (size_t j = 0; j < arg->data.function.arg_count; j++) {
                new_args[idx++] = expr_copy(arg->data.function.args[j]);
            }
            /* Free the intermediate nested function node */
            expr_free(arg); 
        } else {
            new_args[idx++] = arg;
        }
    }
    
    /* Replace old argument array with the flattened one */
    free(e->data.function.args);
    e->data.function.args = new_args;
    e->data.function.arg_count = new_count;
    expr_invalidate_hash(e);   /* args rewritten in place: drop memoized hash */
    return true;
}

/* Public entry: callers may hand us a C-string literal (e.g. internal_call_impl
 * passes "Plus") rather than the interned canonical pointer, so funnel through
 * the interner before the pointer-compare core. */
bool eval_flatten_args(Expr* e, const char* head_name) {
    return eval_flatten_args_interned(e, intern_symbol(head_name));
}

/*
 * has_list_arg:
 * Helper for the Listable attribute.
 * Checks if any argument of the function is an explicit List[...].
 */
/*
 * Does every DownValue of `def` bind its arguments OPAQUELY -- i.e. is each
 * top-level argument of every rule's LHS a bare `Pattern[sym, Blank[]]`, with no
 * head restriction, no PatternTest, no Condition and no nested structure?
 *
 * WHY THE GATE NEEDS THIS. A user symbol has no packed_aware bit, so
 * `f[packedArray]` materialised the buffer -- correctly in general, because
 * src/match.c cannot descend a buffer to match `f[{a_, b_}]`. But the shape that
 * actually appears in numerical code is `jac[u_] := ...`: one bare pattern
 * variable, bound and substituted whole. Nothing looks inside the value, so
 * nothing needs it to be a List.
 *
 * The cost of not having this was the largest single gap left in the packed
 * surface. A 512x512 Jacobi sweep written the obvious way --
 *     jac[u_] := (RotateLeft[u,{1,0}] + ... )/4.;  Nest[jac, u0, 100]
 * -- materialised 262144 Expr nodes on every one of the 100 iterations: 21.6 s,
 * against 107 ms for Mathematica and 0.19 s for the identical body written as a
 * pure function, which the gate already exempts. Defining a helper function
 * should not be what costs 100x.
 *
 * A bare Blank[] matches a packed list already: matching it binds without
 * inspecting, and Head[] answers List. Anything that could look INSIDE -- a
 * head-restricted Blank[h], a sequence pattern, a literal list of sub-patterns,
 * a test or condition that might itself pattern-match -- is refused here and
 * still materialises. Conservative by construction: a shape this does not
 * recognise keeps today's behaviour exactly.
 */
static bool dv_pattern_is_opaque(const Expr* lhs) {
    if (!lhs || lhs->type != EXPR_FUNCTION) return false;
    /* The LHS head must be the plain symbol, not f[a][b] or a pattern itself. */
    if (lhs->data.function.head->type != EXPR_SYMBOL) return false;
    for (size_t i = 0; i < lhs->data.function.arg_count; i++) {
        const Expr* a = lhs->data.function.args[i];
        /* Exactly Pattern[sym, Blank[]] -- two args, second a 0-arg Blank. */
        if (!a || a->type != EXPR_FUNCTION) return false;
        if (a->data.function.head->type != EXPR_SYMBOL) return false;
        if (a->data.function.head->data.symbol.name != SYM_Pattern) return false;
        if (a->data.function.arg_count != 2) return false;
        const Expr* b = a->data.function.args[1];
        if (!b || b->type != EXPR_FUNCTION) return false;
        if (b->data.function.head->type != EXPR_SYMBOL) return false;
        if (b->data.function.head->data.symbol.name != SYM_Blank) return false;
        if (b->data.function.arg_count != 0) return false;   /* _h looks at the head */
    }
    return true;
}

static bool dv_binds_opaquely(const SymbolDef* def) {
    if (!def || !def->down_values) return false;
    /* A builtin with DownValues layered on top is not covered: the builtin still
     * reads the args itself, and packed_aware is the claim for that. */
    if (def->builtin_func) return false;
    for (const Rule* r = def->down_values; r; r = r->next)
        if (!dv_pattern_is_opaque(r->pattern)) return false;
    return true;
}

/* Classify the arguments of a Listable call in one pass: LA_LIST if any is an
 * explicit List, LA_ASSOC if any is an Association, LA_SPARSE if any is a
 * SparseArray (any combination). */
enum { LA_LIST = 1, LA_ASSOC = 2, LA_SPARSE = 4 };
static int listable_arg_kinds(const Expr* e) {
    int kinds = 0;
    size_t n = e->data.function.arg_count;
    for (size_t i = 0; i < n; i++) {
        const Expr* arg = e->data.function.args[i];
        if (arg->type == EXPR_SPARSE) { kinds |= LA_SPARSE; continue; }
        if (arg->type != EXPR_FUNCTION || arg->data.function.head->type != EXPR_SYMBOL)
            continue;
        const char* h = arg->data.function.head->data.symbol.name;
        if (h == SYM_List) {
            /* List threading wins; only a SparseArray further on still matters
             * (it densifies first), and that is one type test per argument. */
            for (size_t j = i + 1; j < n; j++)
                if (e->data.function.args[j]->type == EXPR_SPARSE) kinds |= LA_SPARSE;
            return kinds | LA_LIST;
        }
        if (h == SYM_Association) kinds |= LA_ASSOC;
    }
    return kinds;
}

/* SparseArray[pos -> vals, ...] -- the constructor's standard spelling --
 * reads a packed position matrix and value vector straight from their buffers
 * (sparse/sparse_build.c). Rule is deliberately not packed-aware (the gate's
 * no-nesting invariant, see pack.c), so evaluating the Rule as an ordinary
 * argument would materialise 3 x 10^6 Exprs for a 10^6-entry array before
 * SparseArray ever ran. So the Rule's two sides are evaluated directly and the
 * Rule rebuilt around them: Rule has no rules of its own to apply, so this is
 * the same value. If SparseArray then declines, it unpacks them in place
 * (builtin_sparsearray), so no packed node is left inside the inert call. */
static Expr* eval_sparse_rule_arg(const Expr* rule) {
    Expr* lhs = evaluate(rule->data.function.args[0]);
    Expr* rhs = evaluate(rule->data.function.args[1]);
    Expr* args[2] = { lhs, rhs };
    Expr* r = expr_new_function(expr_new_symbol(SYM_Rule), args, 2);
    if (is_packed_list(lhs) || is_packed_list(rhs)) return r;
    Expr* e = evaluate(r);    /* nothing packed: the ordinary evaluation */
    expr_free(r);
    return e;
}

/* A Listable call with SparseArray arguments. Alongside a List (or a packed
 * list) the sparse arrays densify -- a sparse array plus a dense one is a
 * dense array -- and the rewritten call threads as Lists on the next step.
 * Otherwise the call threads over the sparse arrays themselves
 * (sparse/sparse_thread.c). NULL leaves the call to the steps below; the call
 * itself back (expr_copy) settles it as it stands, after a message
 * (Thread::tdlen for unequal lengths, SparseArray::dense over the implicit
 * densify cap) -- the List threading below must not then run, since it would
 * treat the SparseArray as a scalar. */
static Expr* listable_sparse(const Expr* res, int la_kind) {
    /* N numericalises a SparseArray itself (numeric.c), keeping it sparse. */
    const Expr* h = res->data.function.head;
    if (h->type == EXPR_SYMBOL && h->data.symbol.name == SYM_N) return NULL;
    size_t n = res->data.function.arg_count;
    bool dense = (la_kind & LA_LIST) != 0;
    for (size_t i = 0; i < n && !dense; i++) dense = is_packed_list(res->data.function.args[i]);
    if (!dense) return sp_thread_listable(res);
    /* Every threaded operand must share the leading length. */
    int64_t len = -1;
    for (size_t i = 0; i < n; i++) {
        const Expr* x = res->data.function.args[i];
        int64_t li = x->type == EXPR_SPARSE ? x->data.sparse->dims[0]
                   : is_packed_list(x) ? x->data.ndarray.dims[0]
                   : (x->type == EXPR_FUNCTION && x->data.function.head->type == EXPR_SYMBOL &&
                      x->data.function.head->data.symbol.name == SYM_List)
                         ? (int64_t)x->data.function.arg_count : -1;
        if (li < 0) continue;
        if (len >= 0 && li != len) {
            char* shown = expr_to_string_outputform((Expr*)res);
            mth_message("Thread", "tdlen", "Objects of unequal length in %s cannot be combined.",
                        shown ? shown : "");
            free(shown);
            return expr_copy((Expr*)res);
        }
        len = li;
    }
    /* Times keeps a zero-background factor's sparsity: v Range[n] is sparse. */
    if (h->type == EXPR_SYMBOL && h->data.symbol.name == SYM_Times) {
        Expr* t = sp_thread_times_dense(res);
        if (t) return t;
    }
    Expr** a = malloc(sizeof(Expr*) * (n ? n : 1));
    for (size_t i = 0; i < n; i++) {
        Expr* x = res->data.function.args[i];
        bool refused = false;
        Expr* d = x->type == EXPR_SPARSE ? sp_gate_dense(x, &refused) : NULL;
        if (x->type == EXPR_SPARSE && !d) {
            for (size_t j = 0; j < i; j++) expr_free(a[j]);
            free(a);
            return refused ? expr_copy((Expr*)res) : NULL;
        }
        a[i] = d ? d : expr_copy(x);
    }
    Expr* out = expr_new_function(expr_copy(res->data.function.head), a, n);
    free(a);
    return out;
}

static bool has_list_arg(Expr* e) {
    for (size_t i = 0; i < e->data.function.arg_count; i++) {
        Expr* arg = e->data.function.args[i];
        if (arg->type == EXPR_FUNCTION &&
            arg->data.function.head->type == EXPR_SYMBOL &&
            arg->data.function.head->data.symbol.name == SYM_List) {
            return true;
        }
    }
    return false;
}

/*
 * apply_listable:
 * Implements automatic threading of functions over lists.
 * If a function is Listable and contains list arguments, it maps itself over them.
 * Example: f[{a, b}, c] -> {f[a, c], f[b, c]}
 */
static Expr* apply_listable(Expr* e) {
    /* Determine the required result length from the list arguments, and verify
     * that every List argument shares that length (Mathematica threads over
     * equal-length lists only). The length may legitimately be 0: threading a
     * Listable function over an empty list yields an empty list, e.g.
     * BernoulliB[{}] -> {} and f[{}, c] -> {}. */
    bool have_list = false;
    size_t list_len = 0;
    for (size_t i = 0; i < e->data.function.arg_count; i++) {
        Expr* arg = e->data.function.args[i];
        if (arg->type == EXPR_FUNCTION &&
            arg->data.function.head->type == EXPR_SYMBOL &&
            arg->data.function.head->data.symbol.name == SYM_List) {
            size_t len = arg->data.function.arg_count;
            if (!have_list) { have_list = true; list_len = len; }
            else if (len != list_len) {
                char* s = expr_to_string(e);
                mth_message("Thread", "tdlen",
                            "Objects of unequal length in %s cannot be combined.", s);
                free(s);
                return NULL;
            }
        }
    }

    if (!have_list) return NULL;

    /* Empty list: thread to an empty list without per-element work. */
    if (list_len == 0) {
        return expr_new_function(expr_new_symbol(SYM_List), NULL, 0);
    }

    /* Construct a new List containing the threaded evaluations */
    Expr** new_list_args = malloc(sizeof(Expr*) * list_len);
    for (size_t j = 0; j < list_len; j++) {
        Expr** new_func_args = malloc(sizeof(Expr*) * e->data.function.arg_count);
        for (size_t i = 0; i < e->data.function.arg_count; i++) {
            Expr* arg = e->data.function.args[i];
            if (arg->type == EXPR_FUNCTION &&
            arg->data.function.head->type == EXPR_SYMBOL &&
            arg->data.function.head->data.symbol.name == SYM_List) {
                /* All list arguments must have identical lengths */
                if (arg->data.function.arg_count != list_len) {
                    char* s = expr_to_string(e);
                    mth_message("Thread", "tdlen",
                                "Objects of unequal length in %s cannot be combined.", s);
                    free(s);
                    for (size_t k = 0; k < i; k++) expr_free(new_func_args[k]);
                    free(new_func_args);
                    for (size_t k = 0; k < j; k++) expr_free(new_list_args[k]);
                    free(new_list_args);
                    return NULL;
                }
                new_func_args[i] = expr_copy(arg->data.function.args[j]);
            } else {
                /* Non-list arguments are repeated for every element */
                new_func_args[i] = expr_copy(arg);
            }
        }
        /* Recursively evaluate each threaded call */
        Expr* tmp = expr_new_function(expr_copy(e->data.function.head), new_func_args, e->data.function.arg_count);
        free(new_func_args);
        new_list_args[j] = evaluate(tmp);
        expr_free(tmp);
    }
    
    Expr* final_res = expr_new_function(expr_new_symbol(SYM_List), new_list_args, list_len);
    free(new_list_args);
    return final_res;
}

/*
 * lhs_arg_contains_pattern:
 * Walks an expression looking for any pattern construct (Blank,
 * BlankSequence, BlankNullSequence, Pattern, Optional, Repeated,
 * RepeatedNull, PatternTest, HoldPattern, Condition).  Returns true
 * if any node is found.
 *
 * Used by the Set/SetDelayed dispatcher to decide whether the LHS arg
 * should be evaluated.  Mathematica's rule semantics scope pattern
 * variables inside the LHS — they must NOT be rewritten by existing
 * DownValues during definition.  Evaluating a pattern-bearing LHS
 * arg via the generic evaluator inadvertently fires earlier DownValues
 * on the held pattern (e.g. a CRC table entry whose LHS reshapes to
 * match an earlier rule's pattern, returning that rule's RHS shape),
 * destroying the rule being installed.  Holding the arg whenever it
 * contains a pattern construct sidesteps that.
 *
 * The `f[x] = 1` (x has an OwnValue c → defines f[c] = 1) case still
 * works because non-pattern LHSes contain no Blank/Pattern node and
 * continue down the evaluation path.
 */
static bool lhs_arg_contains_pattern(Expr* e) {
    if (!e) return false;
    if (e->type != EXPR_FUNCTION) return false;
    if (e->data.function.head && e->data.function.head->type == EXPR_SYMBOL) {
        const char* h = e->data.function.head->data.symbol.name;
        if (h == SYM_Blank || h == SYM_BlankSequence || h == SYM_BlankNullSequence
            || h == SYM_Pattern || h == SYM_Optional || h == SYM_Repeated
            || h == SYM_RepeatedNull || h == SYM_PatternTest
            || h == SYM_HoldPattern || h == SYM_Condition) {
            return true;
        }
        if (lhs_arg_contains_pattern(e->data.function.head)) return true;
    }
    for (size_t i = 0; i < e->data.function.arg_count; i++) {
        if (lhs_arg_contains_pattern(e->data.function.args[i])) return true;
    }
    return false;
}

/* Does the pattern `e` mention a symbol with an OwnValue where evaluation can
 * reach it? An assignment evaluates its LHS arguments, patterns
 * included, so PD = PropagatorDenominator; f[a__PD] := ... is a rule for
 * __PropagatorDenominator. Pattern arguments are otherwise held here (they
 * are inert, and holding keeps every existing definition byte-identical), so
 * only such an argument is evaluated. A pattern name (Pattern's first slot)
 * and a held test or guard are not reachable. */
static bool lhs_pattern_mentions_valued_symbol(Expr* e) {
    if (!e) return false;
    if (e->type == EXPR_SYMBOL) return symtab_get_own_values(e->data.symbol.name) != NULL;
    if (e->type != EXPR_FUNCTION) return false;
    Expr* h = e->data.function.head;
    size_t n = e->data.function.arg_count;
    if (h->type == EXPR_SYMBOL) {
        const char* hn = h->data.symbol.name;
        if (hn == SYM_HoldPattern || hn == SYM_Verbatim) return false;
        if (hn == SYM_Pattern)
            return n == 2 && lhs_pattern_mentions_valued_symbol(e->data.function.args[1]);
        if (hn == SYM_Condition || hn == SYM_PatternTest)
            return n >= 1 && lhs_pattern_mentions_valued_symbol(e->data.function.args[0]);
    } else if (lhs_pattern_mentions_valued_symbol(h)) {
        return true;
    }
    for (size_t i = 0; i < n; i++)
        if (lhs_pattern_mentions_valued_symbol(e->data.function.args[i])) return true;
    return false;
}

/* The parts of a pattern-bearing LHS argument that cannot be patterns,
 * evaluated: an assignment evaluates its LHS arguments, so the reference stores
 * h[y_ - 4] := ... as h[-4 + y_] and the rule fires on h[x - 4]. Held whole,
 * the parse's Plus[y_, Times[-1, 4]] matched nothing, so FeynCalc's Pair rules
 * on dim_Symbol - 4 never fired. Pattern nodes themselves stay unevaluated (see
 * lhs_arg_contains_pattern for why), so only a pattern-free subexpression -- a
 * constant, a symbol with a value, a call on concrete arguments -- is
 * evaluated, exactly as it would be outside a pattern, and so is a head symbol
 * with a value (DOT = Dot; f[DOT[a_, b_]] := ...). A pattern name, a held
 * guard or test, a HoldPattern / Verbatim body and the held arguments of a
 * Hold-attributed head are left as written. Returns a new expression. */
static Expr* lhs_fold_pattern_free(Expr* e) {
    if (!lhs_arg_contains_pattern(e)) return evaluate(e);
    if (e->type != EXPR_FUNCTION) return expr_copy(e);
    Expr* h0 = e->data.function.head;
    size_t n = e->data.function.arg_count;
    /* A head symbol with a value is evaluated too: DOT = Dot;
     * f[DOT[a_, b_]] := ... is a rule on Dot[a_, b_], as in the reference. */
    Expr* nh = (h0->type == EXPR_SYMBOL && symtab_get_own_values(h0->data.symbol.name))
                   ? evaluate(h0) : expr_copy(h0);
    Expr* h = nh;
    uint32_t attrs = ATTR_NONE;
    if (h->type == EXPR_SYMBOL) {
        const char* hn = h->data.symbol.name;
        if (hn == SYM_HoldPattern || hn == SYM_Verbatim || hn == SYM_Blank ||
            hn == SYM_BlankSequence || hn == SYM_BlankNullSequence) {
            expr_free(nh);
            return expr_copy(e);
        }
        attrs = get_attributes(hn);
        /* Pattern[name, p]: the name is held; Condition / PatternTest hold
         * their guard. Their pattern slot is folded. */
        if (hn == SYM_Pattern) attrs = ATTR_HOLDFIRST;
        else if (hn == SYM_Condition || hn == SYM_PatternTest) attrs = ATTR_HOLDREST;
    }
    Expr** args = malloc(sizeof(Expr*) * (n ? n : 1));
    for (size_t i = 0; i < n; i++) {
        bool hold = (attrs & ATTR_HOLDALLCOMPLETE) == ATTR_HOLDALLCOMPLETE ||
                    (attrs & ATTR_HOLDALL) == ATTR_HOLDALL ||
                    (i == 0 && (attrs & ATTR_HOLDFIRST)) ||
                    (i > 0 && (attrs & ATTR_HOLDREST));
        args[i] = hold ? expr_copy(e->data.function.args[i])
                       : lhs_fold_pattern_free(e->data.function.args[i]);
    }
    Expr* out = expr_new_function(nh, args, n);
    free(args);
    return out;
}

/*
 * assignment_target_symbol:
 * Returns the name of the symbol whose OwnValue/DownValue/Part would be
 * written if lhs were used as the LHS of a Set or SetDelayed, or NULL if
 * lhs does not name a specific symbol. Wrappers that do not change the
 * ultimate target (Condition[pat, test], HoldPattern[pat], Part[x, ...])
 * are unwrapped so that Protected can be detected on the underlying head.
 */
static const char* assignment_target_symbol(Expr* lhs) {
    if (!lhs) return NULL;
    if (lhs->type == EXPR_SYMBOL) return lhs->data.symbol.name;
    if (lhs->type == EXPR_FUNCTION &&
        lhs->data.function.head->type == EXPR_SYMBOL &&
        lhs->data.function.arg_count >= 1) {
        const char* h = lhs->data.function.head->data.symbol.name;
        if (h == SYM_Condition || h == SYM_HoldPattern || h == SYM_Part ||
            h == SYM_MessageName) {
            /* f::tag = ... targets f, not MessageName: a usage message may be
             * attached to an unprotected symbol even though MessageName itself
             * is Protected. */
            return assignment_target_symbol(lhs->data.function.args[0]);
        }
        return h;
    }
    return NULL;
}

/*
 * lhs_matches_nd_shape:
 * Validate a (possibly nested) List LHS against the rectangular shape of a
 * packed array, starting at axis `axis`. A packed array is always a rectangular
 * block of numbers, so the only thing that can be malformed in `{...} = <array>`
 * is the LHS -- which this checks from the shape alone, materialising nothing.
 * A symbol binds the whole remaining sub-array (or scalar); a non-List function
 * head is a DownValue target and is accepted; a List LHS must line up
 * element-for-element with `axis`'s extent and recurse one axis deeper; anything
 * else (a literal on the LHS, or a List deeper than the array's rank) is not an
 * assignable target.
 */
static bool lhs_matches_nd_shape(const Expr* lhs, const Expr* nd, int axis) {
    if (!lhs) return false;
    if (lhs->type == EXPR_SYMBOL) return true;
    if (lhs->type != EXPR_FUNCTION) return false;
    if (lhs->data.function.head->type != EXPR_SYMBOL) return false;
    if (lhs->data.function.head->data.symbol.name != SYM_List) return true;
    if (axis >= nd->data.ndarray.rank) return false;
    if ((int64_t)lhs->data.function.arg_count != nd->data.ndarray.dims[axis]) return false;
    for (size_t i = 0; i < lhs->data.function.arg_count; i++)
        if (!lhs_matches_nd_shape(lhs->data.function.args[i], nd, axis + 1)) return false;
    return true;
}

/*
 * is_assignable_lhs:
 * Validate that `lhs` shaped against `rhs` is a structurally legal target
 * for a Set/SetDelayed, including all destructured sub-elements. Used as a
 * pre-flight check on List destructuring so a malformed element (e.g. a
 * literal integer on the LHS) cannot land partial assignments before being
 * detected by the in-loop failure path.
 */
static bool is_assignable_lhs(Expr* lhs, Expr* rhs) {
    if (!lhs) return false;
    if (lhs->type == EXPR_SYMBOL) return true;
    if (lhs->type != EXPR_FUNCTION) return false;
    if (lhs->data.function.head->type != EXPR_SYMBOL) return false;
    const char* h = lhs->data.function.head->data.symbol.name;
    if (h != SYM_List) return true; /* downvalue / part / etc. handled in apply_assignment */

    /* A packed-list RHS (an EXPR_NDARRAY presenting as List) destructures like
     * the List it stands for -- see apply_assignment. Validate the LHS against
     * the array's rectangular shape without materialising any element. */
    if (is_packed_list(rhs)) return lhs_matches_nd_shape(lhs, rhs, 0);

    if (rhs->type == EXPR_FUNCTION &&
        rhs->data.function.head->type == EXPR_SYMBOL &&
        rhs->data.function.head->data.symbol.name == SYM_List) {
        /* List RHS: destructure element-wise -- lengths must match. */
        if (lhs->data.function.arg_count != rhs->data.function.arg_count) return false;
        for (size_t i = 0; i < lhs->data.function.arg_count; i++) {
            if (!is_assignable_lhs(lhs->data.function.args[i], rhs->data.function.args[i])) {
                return false;
            }
        }
        return true;
    }

    /* Non-List RHS: Set threads it over the targets ({a, b} = c binds a = c,
     * b = c), so every element must be assignable against the WHOLE rhs. A
     * nested List element threads recursively. */
    for (size_t i = 0; i < lhs->data.function.arg_count; i++) {
        if (!is_assignable_lhs(lhs->data.function.args[i], rhs)) return false;
    }
    return true;
}

/* `lhs := body /; t1 /; t2` -- the delayed RHS Condition[Condition[body, t1], t2]
 * -- is stored as the pattern Condition[Condition[lhs, t1], t2] with RHS body,
 * so the guards run at match time, innermost first, exactly as the same chain
 * written on the LHS (`lhs /; t1 /; t2 := body`).  Returns the pattern: lhs
 * itself when the RHS has no guard (borrowed), else a new tree the caller
 * frees; *body is set to the RHS to store (borrowed from rhs). */
static Expr* hoist_rhs_conditions(Expr* lhs, Expr* rhs, Expr** body) {
    size_t n = 0;
    Expr* g = rhs;
    while (g->type == EXPR_FUNCTION && g->data.function.head->type == EXPR_SYMBOL &&
           g->data.function.head->data.symbol.name == SYM_Condition && g->data.function.arg_count == 2) {
        n++;
        g = g->data.function.args[0];
    }
    *body = g;
    if (n == 0) return lhs;
    Expr** guards = malloc(n * sizeof(Expr*));
    g = rhs;
    for (size_t k = 0; k < n; k++) {             /* outermost first */
        guards[k] = g->data.function.args[1];
        g = g->data.function.args[0];
    }
    Expr* pat = expr_copy(lhs);
    for (size_t k = n; k-- > 0; ) {              /* innermost guard nearest lhs */
        Expr* a[2] = { pat, expr_copy(guards[k]) };
        pat = expr_new_function(expr_new_symbol(SYM_Condition), a, 2);
    }
    free(guards);
    return pat;
}

/* The symbol a Format target is keyed on: the symbol itself, or the innermost
 * symbol head of a call (k[x_] -> k, k[a][b_] -> k). NULL otherwise. */
static const char* format_target_tag(const Expr* e) {
    while (e && e->type == EXPR_FUNCTION) e = e->data.function.head;
    return (e && e->type == EXPR_SYMBOL) ? e->data.symbol.name : NULL;
}

/* MakeBoxes[lhs, form] or Format[lhs] / Format[lhs, form]: an assignment to
 * either installs FormatValues, not a DownValue of the Protected head. */
static bool is_format_target(const Expr* e) {
    if (e->type != EXPR_FUNCTION || e->data.function.head->type != EXPR_SYMBOL) return false;
    const char* h = e->data.function.head->data.symbol.name;
    size_t n = e->data.function.arg_count;
    return (h == SYM_Format && (n == 1 || n == 2)) || (h == SYM_MakeBoxes && n == 2);
}

/* Install the FormatValues the reference stores for an assignment to a format
 * target, keyed on `tag`:
 *   MakeBoxes[lhs, form] := rhs    ->  MakeBoxes[lhs, form] :> rhs
 *   Format[lhs] := rhs             ->  MakeBoxes[lhs, FormatType_] :> Format[rhs, FormatType]
 *                                      lhs :> rhs
 *   Format[lhs, OutputForm] := rhs ->  lhs :> rhs
 *   Format[lhs, form] := rhs       ->  MakeBoxes[lhs, form] :> Format[rhs, form]
 * The plain lhs :> rhs rules are the ones the printer applies (print.c). */
static void install_format_rules(const char* tag, Expr* target, Expr* rhs, bool delayed) {
    if (target->data.function.head->data.symbol.name == SYM_MakeBoxes) {
        symtab_add_format_value(tag, target, rhs, delayed);
        return;
    }
    Expr* lhs = target->data.function.args[0];
    Expr* form = target->data.function.arg_count == 2 ? target->data.function.args[1] : NULL;
    if (form && form->type == EXPR_SYMBOL && form->data.symbol.name == SYM_OutputForm) {
        symtab_add_format_value(tag, lhs, rhs, delayed);
        return;
    }
    Expr* form_pat = form ? expr_copy(form)
        : expr_new_function(expr_new_symbol(SYM_Pattern),
              (Expr*[]){ expr_new_symbol(SYM_FormatType),
                         expr_new_function(expr_new_symbol(SYM_Blank), NULL, 0) }, 2);
    Expr* mb = expr_new_function(expr_new_symbol(SYM_MakeBoxes),
                                 (Expr*[]){ expr_copy(lhs), form_pat }, 2);
    Expr* fm = expr_new_function(expr_new_symbol(SYM_Format),
                                 (Expr*[]){ expr_copy(rhs),
                                            form ? expr_copy(form) : expr_new_symbol(SYM_FormatType) }, 2);
    symtab_add_format_value(tag, mb, fm, delayed);
    expr_free(mb);
    expr_free(fm);
    if (!form) symtab_add_format_value(tag, lhs, rhs, delayed);
}

#if UP_VALUES
/* Store an N-value N[x] = v / N[x, p] := v on `tag`'s NValues. A
 * one-argument N[x] is stored as N[x, {MachinePrecision, MachinePrecision}],
 * the form the reference stores and N asks for at machine precision. */
static void install_n_value(const char* tag, Expr* lhs, Expr* rhs, bool delayed) {
    if (lhs->data.function.arg_count != 1) {
        symtab_add_n_value(tag, lhs, rhs, delayed);
        return;
    }
    Expr* mp = expr_new_function(expr_new_symbol(SYM_List),
                   (Expr*[]){ expr_new_symbol(SYM_MachinePrecision),
                              expr_new_symbol(SYM_MachinePrecision) }, 2);
    Expr* pat = expr_new_function(expr_new_symbol(SYM_N),
                   (Expr*[]){ expr_copy(lhs->data.function.args[0]), mp }, 2);
    symtab_add_n_value(tag, pat, rhs, delayed);
    expr_free(pat);
}

static bool is_n_target(const Expr* e) {
    return e->type == EXPR_FUNCTION && e->data.function.head->type == EXPR_SYMBOL &&
           e->data.function.head->data.symbol.name == SYM_N &&
           (e->data.function.arg_count == 1 || e->data.function.arg_count == 2);
}
#endif

/*
 * apply_assignment:
 * Helper to handle the 'Set' (=) and 'SetDelayed' (:=) primitives.
 * Supports recursive list destructuring.
 * Example: {x, y} = {1, 2}
 *
 * Returns true if the caller should respond as if the assignment
 * succeeded (return the RHS for Set, Null for SetDelayed). Attempts to
 * assign to a Protected symbol emit a Set::wrsym message, leave state
 * unchanged, and still return true so the caller yields the RHS, matching
 * Mathematica semantics.
 */
static bool apply_assignment(Expr* lhs, Expr* rhs, bool is_delayed) {
    /* Options[sym] = {name -> value, ...} redefines the symbol's default
     * option settings. Intercepted before the Protected guard below because
     * Options itself is Protected; this writes the dedicated options store
     * rather than a DownValue (mirrors SetOptions). The RHS must be a List of
     * Rule/RuleDelayed with symbol/string names, else the assignment is left
     * unevaluated. */
    if (lhs->type == EXPR_FUNCTION
        && lhs->data.function.head->type == EXPR_SYMBOL
        && lhs->data.function.head->data.symbol.name == SYM_Options
        && lhs->data.function.arg_count == 1
        && lhs->data.function.args[0]->type == EXPR_SYMBOL) {
        if (rhs->type == EXPR_FUNCTION
            && rhs->data.function.head->type == EXPR_SYMBOL
            && rhs->data.function.head->data.symbol.name == SYM_List) {
            bool all_rules = true;
            for (size_t i = 0; i < rhs->data.function.arg_count; i++) {
                Expr* r = rhs->data.function.args[i];
                if (!(r->type == EXPR_FUNCTION
                      && r->data.function.head->type == EXPR_SYMBOL
                      && (r->data.function.head->data.symbol.name == SYM_Rule
                          || r->data.function.head->data.symbol.name == SYM_RuleDelayed)
                      && r->data.function.arg_count == 2
                      && (r->data.function.args[0]->type == EXPR_SYMBOL
                          || r->data.function.args[0]->type == EXPR_STRING))) {
                    all_rules = false;
                    break;
                }
            }
            if (all_rules) {
                symtab_set_options(lhs->data.function.args[0]->data.symbol.name,
                                   expr_copy(rhs));
                eval_clock_bump();
                return true;
            }
        }
        return false;
    }

    /* Format[lhs] = rhs / Format[lhs, form] := rhs: FormatValues of the symbol
     * lhs is keyed on, as the reference stores them. Format is Protected, so
     * this is intercepted before the Protected guard. (An untagged
     * MakeBoxes[...] := rhs stays a DownValue of MakeBoxes, as in the
     * reference; tag /: MakeBoxes[...] := rhs is handled by TagSet.) */
    if (is_format_target(lhs) && lhs->data.function.head->data.symbol.name == SYM_Format) {
        const char* tag = format_target_tag(lhs->data.function.args[0]);
        if (tag) {
            if (get_attributes(tag) & ATTR_PROTECTED) {
                char* lt = expr_to_string_outputform(lhs->data.function.args[0]);
                mth_message(is_delayed ? "SetDelayed" : "Set", "write",
                            "Tag %s in %s is Protected.", tag, lt ? lt : "");
                free(lt);
                return true;
            }
            install_format_rules(tag, lhs, rhs, is_delayed);
            return true;
        }
    }

#if UP_VALUES
    /* N[x] = v, N[f[...]] := v (and N[x, p] = v): an N-value of the symbol x or
     * the head f, kept in NValues as in the reference; N consults them
     * (numeric.c n_value). N itself is Protected, so this is intercepted
     * before the Protected guard. */
    if (is_n_target(lhs)) {
        Expr* x = lhs->data.function.args[0];
        const char* tag = NULL;
        if (x->type == EXPR_SYMBOL) tag = x->data.symbol.name;
        else if (x->type == EXPR_FUNCTION && x->data.function.head->type == EXPR_SYMBOL)
            tag = x->data.function.head->data.symbol.name;
        if (tag) {
            if (get_attributes(tag) & ATTR_PROTECTED) {
                char* lt = expr_to_string_outputform(lhs);
                mth_message(is_delayed ? "SetDelayed" : "Set", "write",
                            "Tag %s in %s is Protected.", tag, lt ? lt : "");
                free(lt);
                return true;
            }
            install_n_value(tag, lhs, rhs, is_delayed);
            return true;
        }
    }

    /* UpValues[sym] = list / DownValues[sym] = list / OwnValues[sym] = list:
     * replace the whole rule list from a List of (Rule|RuleDelayed) rules.
     * Intercepted before the Protected guard because these heads are themselves
     * Protected (mirrors the Options[sym]=... intercept above). */
    if (lhs->type == EXPR_FUNCTION
        && lhs->data.function.head->type == EXPR_SYMBOL
        && lhs->data.function.arg_count == 1
        && lhs->data.function.args[0]->type == EXPR_SYMBOL
        && (lhs->data.function.head->data.symbol.name == SYM_UpValues
            || lhs->data.function.head->data.symbol.name == SYM_DownValues
            || lhs->data.function.head->data.symbol.name == SYM_OwnValues
            || lhs->data.function.head->data.symbol.name == SYM_FormatValues
            || lhs->data.function.head->data.symbol.name == SYM_NValues)) {
        const char* h = lhs->data.function.head->data.symbol.name;
        const char* target = lhs->data.function.args[0]->data.symbol.name;
        /* A list-assignment mutates `target`'s definitions, so it must honour
         * Protected just like `target[...] = rhs` would (otherwise
         * `DownValues[Sin] = {...}` would silently overwrite a builtin). */
        if (get_attributes(target) & ATTR_PROTECTED) {
            /* Mathematica reports this as Set::write ("Tag f in DownValues[f] is
             * Protected."), not Set::wrsym -- match it so `::write` greps hit. */
            mth_message(is_delayed ? "SetDelayed" : "Set", "write",
                        "Tag %s in %s[%s] is Protected.", target, h, target);
            return true;
        }
        if (rhs->type == EXPR_FUNCTION
            && rhs->data.function.head->type == EXPR_SYMBOL
            && rhs->data.function.head->data.symbol.name == SYM_List) {
            if (h == SYM_UpValues)          symtab_set_up_values(target, rhs);
            else if (h == SYM_FormatValues) symtab_set_format_values(target, rhs);
            else if (h == SYM_NValues)      symtab_set_n_values(target, rhs);
            else if (h == SYM_DownValues)   symtab_set_down_values(target, rhs);
            else                          symtab_set_own_values(target, rhs);
            return true;
        }
        return false;   /* RHS is not a List: leave unevaluated */
    }
#endif

    /* Block writes to Protected symbols. List destructuring is recursed
     * into below and each child runs through apply_assignment again, so
     * per-element protection checks happen naturally -- we only skip the
     * outer check when lhs itself is a List. */
    bool lhs_is_list = (lhs->type == EXPR_FUNCTION &&
                        lhs->data.function.head->type == EXPR_SYMBOL &&
                        lhs->data.function.head->data.symbol.name == SYM_List);
    if (!lhs_is_list) {
        const char* target = assignment_target_symbol(lhs);
        if (target && (get_attributes(target) & ATTR_PROTECTED)) {
            /* A symbol LHS is ::wrsym; a rule for a protected tag is ::write,
             * naming the LHS, as the reference reports them. */
            if (lhs->type == EXPR_SYMBOL) {
                mth_message(is_delayed ? "SetDelayed" : "Set", "wrsym",
                            "Symbol %s is Protected.", target);
            } else {
                char* lt = expr_to_string_outputform(lhs);
                mth_message(is_delayed ? "SetDelayed" : "Set", "write",
                            "Tag %s in %s is Protected.", target, lt ? lt : "");
                free(lt);
            }
            return true;
        }
    }

    if (lhs->type == EXPR_SYMBOL && context_is_state_var(lhs->data.symbol.name)) {
        /* $Context / $ContextPath ARE the resolver's state, so they are
         * validated before the write and always stored evaluated (a
         * SetDelayed is evaluated here, once). An invalid value emits
         * ::cxset / ::cxlist and keeps the old one; Set still returns rhs. */
        Expr* probe = is_delayed ? evaluate(expr_copy(rhs)) : expr_copy(rhs);
        context_assign(lhs->data.symbol.name, probe);
        expr_free(probe);
        return true;
    }

    if (lhs->type == EXPR_SYMBOL) {
        /* Standard symbol assignment */
        symtab_add_own_value_as(lhs->data.symbol.name, lhs, rhs, is_delayed);

        /* Special system variables: keep their C-side mirror state in sync.
         * Set has HoldFirst, so for `$RecursionLimit = expr` the rhs is
         * already evaluated; for SetDelayed, we evaluate a copy here so the
         * limit reflects the value the user expects to see when they read
         * the symbol back. If validation fails, the OwnValue is rolled back
         * to the current C-side limit. */
        if (strcmp(lhs->data.symbol.name, "$RecursionLimit") == 0) {
            Expr* probe = is_delayed ? evaluate(expr_copy(rhs)) : expr_copy(rhs);
            sync_recursion_limit_from_value(probe);
            expr_free(probe);
        } else if (lhs->data.symbol.name[0] == '$' &&
                   eval_is_sysflag(lhs->data.symbol.name)) {
            /* Sigil first so an ordinary symbol assignment pays one byte compare,
             * then the NAME -- never a probe for a symbol that is not a flag. */
            Expr* probe = is_delayed ? evaluate(expr_copy(rhs)) : expr_copy(rhs);
            eval_sync_sysflag(lhs->data.symbol.name, probe);
            expr_free(probe);
        }
        return true;
    } else if (lhs->type == EXPR_FUNCTION &&
               lhs->data.function.head->type == EXPR_SYMBOL &&
               lhs->data.function.head->data.symbol.name == SYM_List) {
        /* A List LHS is either DESTRUCTURED (a List RHS of matching length,
         * element for element) or THREADED (any other RHS is broadcast to each
         * target: {a, b} = c binds a = c, b = c -- Wolfram Set semantics). It
         * NEVER installs a DownValue on List, so this branch always returns
         * here rather than falling through below.
         *
         * A packed-array RHS is an EXPR_NDARRAY (present_as List), not a List
         * node. Set is a packed-aware head so a whole-value binding
         * (x = Range[10^6]) keeps its argument packed; but destructuring is the
         * one assignment path that reads the RHS *structure*, and the
         * transparency gate leaves a packed argument intact for an aware head.
         * So normalise a packed RHS to a List of its top-level slices first --
         * slices stay packed, so {xc, yc} = {Range[m], Range[n]} binds packed
         * vectors rather than materialising one Expr per element. */
        Expr* rhs_ds = rhs;
        bool own_rhs_ds = false;
        if (is_packed_list(rhs)) {
            Expr* sliced = ndarray_unpack_top_level(rhs);
            if (sliced) { rhs_ds = sliced; own_rhs_ds = true; }
        }

        bool rhs_is_list = (rhs_ds->type == EXPR_FUNCTION &&
                            rhs_ds->data.function.head->type == EXPR_SYMBOL &&
                            rhs_ds->data.function.head->data.symbol.name == SYM_List);

        /* Pre-flight every element so a malformed target (e.g. a literal
         * integer on the LHS) fails the whole assignment before any sibling is
         * bound -- partial assignments would otherwise leak past the failure.
         * In the threaded case each element pairs with the whole rhs; in the
         * destructured case with the matching rhs element. */
        bool all_ok = true;
        size_t n = lhs->data.function.arg_count;
        if (rhs_is_list && n != rhs_ds->data.function.arg_count) {
            all_ok = false;                         /* length mismatch: leave unevaluated */
        } else {
            for (size_t i = 0; i < n; i++) {
                Expr* rhs_i = rhs_is_list ? rhs_ds->data.function.args[i] : rhs_ds;
                if (!is_assignable_lhs(lhs->data.function.args[i], rhs_i)) {
                    all_ok = false;
                    break;
                }
            }
            if (all_ok) {
                for (size_t i = 0; i < n; i++) {
                    Expr* rhs_i = rhs_is_list ? rhs_ds->data.function.args[i] : rhs_ds;
                    if (!apply_assignment(lhs->data.function.args[i], rhs_i, is_delayed)) {
                        all_ok = false;
                    }
                }
            }
        }
        if (own_rhs_ds) expr_free(rhs_ds);
        return all_ok;
    } else if (lhs->type == EXPR_FUNCTION) {
        if (lhs->data.function.head->type == EXPR_SYMBOL && lhs->data.function.head->data.symbol.name == SYM_Part) {
            Expr* expr_part_assign(Expr* lhs, Expr* rhs); // Forward declare or include part.h
            Expr* assigned = expr_part_assign(lhs, rhs);
            if (assigned) {
                expr_free(assigned);
                return true;
            }
            return false;
        } else if (lhs->data.function.head->type == EXPR_SYMBOL) {
            /* Pattern-based assignment (DownValues) */
            /* We use the entire lhs as the pattern, and its head as the key */
            const char* symbol_name = lhs->data.function.head->data.symbol.name;

            /* Association element assignment: `assoc[k1, ..., kn] = val` when
             * the head symbol currently holds an Association -- Wolfram's
             * assoc[key] = val sugar (a = <||>; a["x"] = 5 mutates a, it does
             * NOT install a DownValue). Reroute to the Part machinery with each
             * key wrapped in Key[...] (assoc[k] === assoc[[Key[k]]]); the keys
             * are evaluated first (a = <||>; k = "x"; a[k] = 5 sets a["x"]).
             * Guard on the symbol's stored OwnValue being an Association so an
             * ordinary DownValue definition (f[x_] := ..., g[1] = 5) is
             * untouched -- expr_part_assign re-reads and writes the symbol. */
            {
                Expr* expr_part_assign(Expr* lhs_p, Expr* rhs_p);
                Rule* ov = symtab_get_own_values(symbol_name);
                if (ov && ov->replacement &&
                    ov->replacement->type == EXPR_FUNCTION &&
                    ov->replacement->data.function.head->type == EXPR_SYMBOL &&
                    ov->replacement->data.function.head->data.symbol.name == SYM_Association) {
                    size_t nk = lhs->data.function.arg_count;
                    Expr** pargs = malloc(sizeof(Expr*) * (nk + 1));
                    pargs[0] = expr_copy(lhs->data.function.head);
                    for (size_t i = 0; i < nk; i++) {
                        Expr* kev = evaluate(lhs->data.function.args[i]);   /* evaluate borrows */
                        Expr** karg = malloc(sizeof(Expr*));
                        karg[0] = kev;
                        pargs[i + 1] = expr_new_function(expr_new_symbol(SYM_Key), karg, 1);
                        free(karg);
                    }
                    Expr* part_lhs = expr_new_function(expr_new_symbol(SYM_Part), pargs, nk + 1);
                    free(pargs);
                    Expr* assigned = expr_part_assign(part_lhs, rhs);
                    expr_free(part_lhs);
                    if (assigned) { expr_free(assigned); eval_clock_bump(); return true; }
                    /* Could not place the key (e.g. a nested path through a
                     * missing intermediate association): leave the Set
                     * unevaluated rather than installing a spurious DownValue. */
                    return false;
                }
            }

            /* f::usage = "..." additionally registers the string as f's
             * docstring so ?f and Information[f] surface it. The message is
             * still installed as a DownValue on MessageName below, which makes
             * MessageName[f, "usage"] (i.e. f::usage) retrievable. */
            if (symbol_name == SYM_MessageName &&
                lhs->data.function.arg_count == 2 &&
                lhs->data.function.args[0]->type == EXPR_SYMBOL &&
                lhs->data.function.args[1]->type == EXPR_STRING &&
                strcmp(lhs->data.function.args[1]->data.string, "usage") == 0 &&
                rhs->type == EXPR_STRING) {
                symtab_set_docstring(lhs->data.function.args[0]->data.symbol.name,
                                     rhs->data.string);
            }

            Expr* actual_pattern = lhs;
            Expr* actual_rhs = rhs;

            /* If the RHS is Condition[body, test], move the condition to the LHS.
             * This makes f[x_] := body /; test equivalent to f[x_] /; test := body.
             * This is standard Mathematica semantics. */
            if (is_delayed)
                actual_pattern = hoist_rhs_conditions(lhs, rhs, &actual_rhs);

            if (symbol_name == SYM_Condition && actual_pattern->data.function.arg_count == 2) {
                /* `lhs /; c1 /; c2` nests to the LEFT: peel every guard. */
                Expr* inner_lhs = actual_pattern->data.function.args[0];
                while (inner_lhs->type == EXPR_FUNCTION && inner_lhs->data.function.head->type == EXPR_SYMBOL &&
                       inner_lhs->data.function.head->data.symbol.name == SYM_Condition &&
                       inner_lhs->data.function.arg_count == 2)
                    inner_lhs = inner_lhs->data.function.args[0];
                if (inner_lhs->type == EXPR_FUNCTION && inner_lhs->data.function.head->type == EXPR_SYMBOL) {
                    symbol_name = inner_lhs->data.function.head->data.symbol.name;
                } else if (inner_lhs->type == EXPR_SYMBOL) {
                    symbol_name = inner_lhs->data.symbol.name;
                }
            }
            /* HoldPattern[f[...]] := rhs (also under a guard, and with a guard
             * inside): a rule for f, stored with its HoldPattern, as the
             * reference stores it. HoldPattern[s] = v is an OwnValue of s.
             * The rule used to land on HoldPattern itself, silently. */
            if (symbol_name == SYM_HoldPattern) {
                Expr* inner = lhs;
                while (inner->type == EXPR_FUNCTION && inner->data.function.head->type == EXPR_SYMBOL &&
                       (inner->data.function.head->data.symbol.name == SYM_HoldPattern ||
                        inner->data.function.head->data.symbol.name == SYM_Condition) &&
                       inner->data.function.arg_count >= 1)
                    inner = inner->data.function.args[0];
                if (inner->type == EXPR_SYMBOL && actual_pattern == lhs) {
                    return apply_assignment(inner, rhs, is_delayed);
                }
                if (inner->type == EXPR_FUNCTION && inner->data.function.head->type == EXPR_SYMBOL &&
                    inner->data.function.head->data.symbol.name != SYM_List &&
                    inner->data.function.head->data.symbol.name != SYM_Part)
                    symbol_name = inner->data.function.head->data.symbol.name;
            }
            /* `f[x_][y_] /; c := ...`: the guarded LHS has a compound head --
             * a SubValue of its innermost head symbol. */
            if (symbol_name == SYM_Condition) {
                Expr* inner = actual_pattern->data.function.args[0];
                while (inner->type == EXPR_FUNCTION && inner->data.function.head->type == EXPR_SYMBOL &&
                       inner->data.function.head->data.symbol.name == SYM_Condition &&
                       inner->data.function.arg_count == 2)
                    inner = inner->data.function.args[0];
                const char* sub_sym = symtab_innermost_head(inner);
                if (sub_sym) {
                    if (get_attributes(sub_sym) & ATTR_PROTECTED) {
                        mth_message(is_delayed ? "SetDelayed" : "Set", "write",
                                    "Tag %s in the assignment is Protected.", sub_sym);
                        if (actual_pattern != lhs) expr_free(actual_pattern);
                        return true;
                    }
                    symtab_add_sub_value(sub_sym, actual_pattern, actual_rhs, is_delayed);
                    if (actual_pattern != lhs) expr_free(actual_pattern);
                    return true;
                }
            }
            symtab_add_down_value_as(symbol_name, actual_pattern, actual_rhs, is_delayed);
            if (actual_pattern != lhs) expr_free(actual_pattern);
            return true;
        } else if (lhs->data.function.head->type == EXPR_FUNCTION) {
            /* f[a][b] = v / f[x_][y_] := rhs: a SubValue of the innermost head
             * symbol f (any depth: g[a][b][c]). */
            const char* sub_sym = symtab_innermost_head(lhs);
            if (!sub_sym) return false;
            if (get_attributes(sub_sym) & ATTR_PROTECTED) {
                mth_message(is_delayed ? "SetDelayed" : "Set", "write",
                            "Tag %s in the assignment is Protected.", sub_sym);
                return true;
            }
            Expr* actual_rhs = rhs;
            Expr* actual_pattern = is_delayed ? hoist_rhs_conditions(lhs, rhs, &actual_rhs) : lhs;
            symtab_add_sub_value(sub_sym, actual_pattern, actual_rhs, is_delayed);
            if (actual_pattern != lhs) expr_free(actual_pattern);
            return true;
        }
    }
    return false;
}

/* Evaluate the ARGUMENTS of an assignment LHS to canonical form, holding pattern
 * constructs and Hold-attributed positions, exactly as the Set primitive does
 * (so `f[x] = v` with x=c installs f[c]=v, while a pattern LHS is preserved).
 * Returns a new node (sets *freed=true, caller frees it) or `lhs` unchanged
 * (*freed=false) when lhs is not a function. Shared by Set and, under UP_VALUES,
 * by UpSet / TagSet / TagUnset. */
static Expr* build_assignment_target(Expr* lhs, bool* freed) {
    *freed = false;
    if (lhs->type != EXPR_FUNCTION) return lhs;

    Expr** eval_args = malloc(sizeof(Expr*) * lhs->data.function.arg_count);
    bool is_part = (lhs->data.function.head->type == EXPR_SYMBOL &&
                    lhs->data.function.head->data.symbol.name == SYM_Part);
    bool is_list = (lhs->data.function.head->type == EXPR_SYMBOL &&
                    lhs->data.function.head->data.symbol.name == SYM_List);
    uint32_t lhs_attrs = ATTR_NONE;
    if (lhs->data.function.head->type == EXPR_SYMBOL)
        lhs_attrs = get_attributes(lhs->data.function.head->data.symbol.name);

    for (size_t i = 0; i < lhs->data.function.arg_count; i++) {
        bool hold = false;
        if ((lhs_attrs & ATTR_HOLDALLCOMPLETE) == ATTR_HOLDALLCOMPLETE) hold = true;
        else if ((lhs_attrs & ATTR_HOLDALL) == ATTR_HOLDALL) hold = true;
        else if (i == 0 && (lhs_attrs & ATTR_HOLDFIRST)) hold = true;
        else if (i > 0 && (lhs_attrs & ATTR_HOLDREST)) hold = true;
        if (is_part && i == 0) hold = true;
        bool fold = !hold && lhs_arg_contains_pattern(lhs->data.function.args[i]) &&
                    !lhs_pattern_mentions_valued_symbol(lhs->data.function.args[i]);
        if (is_list) {
            Expr* child = lhs->data.function.args[i];
            if (child->type == EXPR_SYMBOL) hold = true;
            else if (child->type == EXPR_FUNCTION &&
                     child->data.function.head->type == EXPR_SYMBOL &&
                     child->data.function.head->data.symbol.name == SYM_List) hold = true;
        }
        eval_args[i] = hold ? expr_copy(lhs->data.function.args[i])
                     : fold ? lhs_fold_pattern_free(lhs->data.function.args[i])
                            : evaluate(lhs->data.function.args[i]);
    }
    /* A compound head f[a][b] (a SubValue target) has its own arguments
     * evaluated the same way: x = 5; f[x][1] = 3 defines f[5][1]. */
    Expr* head;
    if (lhs->data.function.head->type == EXPR_FUNCTION) {
        bool hfreed = false;
        head = build_assignment_target(lhs->data.function.head, &hfreed);
        if (!hfreed) head = expr_copy(head);
    } else {
        head = expr_copy(lhs->data.function.head);
        /* The head is evaluated too: h = g; h[x_] := ... defines g, and so
         * does a Block/Module-local holding a computed symbol
         * (ch = ToExpression["cache" <> name]; ch[a_] := ...). Only a head
         * with an OwnValue can change, so the common case pays one lookup.
         * Only a symbol value retargets: a = <||>; a[k] = v is an element
         * assignment into the Association (apply_assignment's own path), and
         * a number or pure function is no assignable head. */
        if (head->type == EXPR_SYMBOL && symtab_get_own_values(head->data.symbol.name)) {
            Expr* hv = evaluate(head);
            if (hv->type == EXPR_SYMBOL) {
                expr_free(head);
                head = hv;
            } else {
                expr_free(hv);
            }
        }
    }
    Expr* target = expr_new_function(head, eval_args, lhs->data.function.arg_count);
    free(eval_args);
    *freed = true;
    return target;
}

Expr* eval_hoist_rhs_conditions(Expr* lhs, Expr* rhs, Expr** body) {
    return hoist_rhs_conditions(lhs, rhs, body);
}

/* Public face of build_assignment_target for Unset, which lives in core.c:
 * `f[i] =.` with i = 2 removes the f[2] rule, exactly as `f[i] = v` set it. */
Expr* eval_assignment_target(Expr* lhs, bool* freed) {
    return build_assignment_target(lhs, freed);
}

#if UP_VALUES
#define UP_MAX_LEVEL1 64

/* Heads whose held LHS must not be scanned for UpValues during evaluation --
 * otherwise an all-heads upvalue (_[g[x_]] ^= ...) could hijack the held LHS of
 * an assignment, e.g. Set[g[y], 1]. */
/* Set / SetDelayed / Unset are NOT on this list: an upvalue defined for them
 * (`dt /: Set[dt[a_, b_], v_] := ...`, `g /: Set[g, v_] := ...`) intercepts the
 * assignment, as the reference does -- packages use it to keep bookkeeping
 * (FeynCalc's DataType, $NonComm, ScalarProduct[p, q] = ...).  A rule only
 * fires when its whole pattern matches the Set call, so an all-heads upvalue
 * like _[g[x_]] ^= v (one argument) cannot hijack Set[g[y], 1]. */
static inline bool is_assignment_primitive(const char* head_name) {
    return head_name == SYM_UpSet || head_name == SYM_UpSetDelayed ||
           head_name == SYM_TagSet || head_name == SYM_TagSetDelayed ||
           head_name == SYM_TagUnset;
}

/* The interned "dispatch head" of a level-one LHS element -- the symbol whose
 * up_values an occurrence of this element would consult. A bare symbol is itself;
 * an ordinary call g[...] gives g; and a pattern that constrains an argument's
 * head -- a_h (Pattern[a, Blank[h]]) or _h (Blank[h] / __h / ___h) -- gives h, so
 * `mod /: a_mod + b_mod := ...` keys on mod. Returns NULL for a head-free blank
 * (a_ / _) or a non-symbol-headed element. */
static const char* element_dispatch_head(Expr* a) {
    if (!a) return NULL;
    if (a->type == EXPR_SYMBOL) return a->data.symbol.name;
    if (a->type != EXPR_FUNCTION || !a->data.function.head) return NULL;
    /* A compound head g[x][y] dispatches on its innermost symbol head, g --
     * `clc /: MakeBoxes[clc[x___][y___], fmt] := ...` keys on clc. */
    if (a->data.function.head->type == EXPR_FUNCTION) {
        Expr* hh = a->data.function.head;
        while (hh->type == EXPR_FUNCTION && hh->data.function.head) hh = hh->data.function.head;
        return hh->type == EXPR_SYMBOL ? hh->data.symbol.name : NULL;
    }
    if (a->data.function.head->type != EXPR_SYMBOL) return NULL;
    const char* h = a->data.function.head->data.symbol.name;
    if (h == SYM_Pattern && a->data.function.arg_count == 2)
        return element_dispatch_head(a->data.function.args[1]);   /* a_h -> head of the blank; a_ -> NULL */
    if (h == SYM_Blank || h == SYM_BlankSequence || h == SYM_BlankNullSequence)
        return (a->data.function.arg_count == 1 &&                /* _h -> h ; bare _ -> NULL */
                a->data.function.args[0]->type == EXPR_SYMBOL)
               ? a->data.function.args[0]->data.symbol.name : NULL;
    /* Pattern-construct wrappers dispatch on the WRAPPED pattern's head, never the
     * wrapper itself (x_Integer?EvenQ -> Integer, x_:0 -> head of x_, p/;c -> head
     * of p). Alternatives has no single head. Otherwise g[...] -> g. */
    if (h == SYM_PatternTest || h == SYM_Condition || h == SYM_Optional)
        return (a->data.function.arg_count >= 1)
               ? element_dispatch_head(a->data.function.args[0]) : NULL;
    if (h == SYM_Alternatives) return NULL;
    return h;                                                     /* g[...] -> g */
}

/* Distinct interned dispatch heads of the level-one elements of an UpSet LHS (its
 * arguments -- NOT the head of lhs, since UpSet keys on operands, not operator).
 * A top-level Condition (lhs /; test) is unwrapped first, since the test is not an
 * operand. Returns the count. */
static int collect_level_one_symbols(Expr* lhs, const char** out, int cap) {
    int n = 0;
    if (!lhs) return 0;
    while (lhs->type == EXPR_FUNCTION && lhs->data.function.head->type == EXPR_SYMBOL &&
        (lhs->data.function.head->data.symbol.name == SYM_Condition ||
         lhs->data.function.head->data.symbol.name == SYM_HoldPattern) &&
        lhs->data.function.arg_count >= 1)
        lhs = lhs->data.function.args[0];
    if (lhs->type != EXPR_FUNCTION) return 0;
    for (size_t i = 0; i < lhs->data.function.arg_count && n < cap; i++) {
        const char* name = element_dispatch_head(lhs->data.function.args[i]);
        if (!name) continue;
        int seen = 0;
        for (int k = 0; k < n; k++) if (out[k] == name) { seen = 1; break; }
        if (!seen) out[n++] = name;
    }
    return n;
}

/* Would `lhs ^= rhs` have placed its rule on `tag` ALONE? Then Definition can
 * render the upvalue in UpSet form (lhs ^= rhs); otherwise it needs the TagSet
 * form (tag /: lhs = rhs), since UpSet would also have tagged the other
 * operands. A literal atom operand (2, "s", 1.5) counts as another tag -- its
 * head (Integer, String, Real) is what UpSet would try to tag. */
bool upvalue_renders_as_upset(const char* tag, Expr* lhs) {
    const char* syms[UP_MAX_LEVEL1];
    int n = collect_level_one_symbols(lhs, syms, UP_MAX_LEVEL1);
    if (n != 1 || syms[0] != tag) return false;
    while (lhs->type == EXPR_FUNCTION && lhs->data.function.head->type == EXPR_SYMBOL &&
        lhs->data.function.head->data.symbol.name == SYM_Condition &&
        lhs->data.function.arg_count >= 1)
        lhs = lhs->data.function.args[0];
    if (lhs->type != EXPR_FUNCTION) return false;
    for (size_t i = 0; i < lhs->data.function.arg_count; i++) {
        ExprType t = lhs->data.function.args[i]->type;
        if (t != EXPR_SYMBOL && t != EXPR_FUNCTION) return false;
    }
    return true;
}

typedef enum { TAG_NOTFOUND, TAG_OWN, TAG_DOWN, TAG_UP, TAG_SUB } TagPos;

/* Classify where interned name `tag` occurs in `lhs`, per the Mathematica rule
 * for TagSet / TagUnset ("defines upvalues, downvalues or subvalues as
 * appropriate"). Names are interned, so identity compares. */
static TagPos classify_tag_position(const char* tag, Expr* lhs) {
    /* A Condition wrapper (lhs /; test) is part of the pattern, not a position
     * for the tag -- classify against the underlying lhs so e.g.
     * `xr /: xr + y_ /; y > -2 := ...` finds xr inside the held Plus. */
    while (lhs->type == EXPR_FUNCTION &&
        lhs->data.function.head->type == EXPR_SYMBOL &&
        (lhs->data.function.head->data.symbol.name == SYM_Condition ||
         lhs->data.function.head->data.symbol.name == SYM_HoldPattern) &&
        lhs->data.function.arg_count >= 1)
        lhs = lhs->data.function.args[0];
    if (lhs->type == EXPR_SYMBOL)
        return (lhs->data.symbol.name == tag) ? TAG_OWN : TAG_NOTFOUND;
    if (lhs->type != EXPR_FUNCTION) return TAG_NOTFOUND;
    Expr* head = lhs->data.function.head;
    if (head->type == EXPR_SYMBOL && head->data.symbol.name == tag) return TAG_DOWN;
    if (head->type == EXPR_FUNCTION && symtab_innermost_head(lhs) == tag) return TAG_SUB;
    for (size_t i = 0; i < lhs->data.function.arg_count; i++) {
        if (element_dispatch_head(lhs->data.function.args[i]) == tag) return TAG_UP;
    }
    return TAG_NOTFOUND;
}

/* A delayed RHS of the form `body /; test` has its Condition moved onto the LHS
 * pattern, so `lhs ^:= body /; test` behaves as `lhs /; test ^:= body` (identical
 * to how SetDelayed is handled in apply_assignment). Returns the pattern to store:
 * `lhs` unchanged (*owned=false) or a new Condition[lhs, test] (*owned=true, the
 * caller frees it); *body is set to the RHS to store (borrowed either way). */
static Expr* move_rhs_condition(Expr* lhs, Expr* rhs, bool is_delayed,
                                bool* owned, Expr** body) {
    *body = rhs;
    if (!is_delayed) { *owned = false; return lhs; }
    Expr* pat = hoist_rhs_conditions(lhs, rhs, body);
    *owned = (pat != lhs);
    return pat;
}

/* lhs ^= rhs / lhs ^:= rhs : install the rule on the up_values of every distinct
 * level-one symbol of lhs. rhs is already evaluated (UpSet is HoldFirst) or held
 * (UpSetDelayed is HoldAll); either way symtab_add_up_value copies it in.
 * Returns true if there was at least one symbol to attach to. */
static bool apply_up_assignment(Expr* lhs, Expr* rhs, bool is_delayed) {
    const char* syms[UP_MAX_LEVEL1];
    int n = collect_level_one_symbols(lhs, syms, UP_MAX_LEVEL1);
    if (n == 0) {
        /* No argument symbol to attach to (e.g. `f[x_] ^= rhs`): message + leave
         * the UpSet unevaluated, matching Mathematica's UpSet::nosym. */
        mth_message(is_delayed ? "UpSetDelayed" : "UpSet", "nosym",
                    "No new symbol was found on which to place the up-value.");
        return false;
    }
    bool owned; Expr* body;
    Expr* pat = move_rhs_condition(lhs, rhs, is_delayed, &owned, &body);
    for (int i = 0; i < n; i++) {
        if (get_attributes(syms[i]) & ATTR_PROTECTED) {
            mth_message(is_delayed ? "UpSetDelayed" : "UpSet", "write",
                        "Tag %s in the assignment is Protected.", syms[i]);
            continue;
        }
        if (is_n_target(pat)) install_n_value(syms[i], pat, body, is_delayed);   /* N[z] ^= v */
        else symtab_add_up_value_as(syms[i], pat, body, is_delayed);
    }
    if (owned) expr_free(pat);
    return true;
}

/* Does symbol `tag` (interned) occur anywhere in e, heads included? */
static bool tag_occurs_in(const Expr* e, const char* tag) {
    if (e->type == EXPR_SYMBOL) return e->data.symbol.name == tag;
    if (e->type != EXPR_FUNCTION) return false;
    if (tag_occurs_in(e->data.function.head, tag)) return true;
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (tag_occurs_in(e->data.function.args[i], tag)) return true;
    return false;
}

/* tag /: lhs = rhs  (and := ): install on tag as own/down/up value per the tag's
 * position in lhs. Returns true when a rule was installed, false when the
 * assignment was refused with a message (tag not a symbol, Protected, absent or
 * too deep). The caller (branch 6c) yields the RHS for TagSet either way, as
 * the Set family does on an invalid left-hand side (`3 = 4` emits Set::setraw
 * and returns 4), and for TagSetDelayed Null on success, $Failed on refusal. */
static bool apply_tag_assignment(Expr* tag, Expr* lhs, Expr* rhs, bool is_delayed) {
    const char* head = is_delayed ? "TagSetDelayed" : "TagSet";
    if (tag->type != EXPR_SYMBOL) {
        mth_message(head, "sym", "Argument at position 1 is expected to be a symbol.");
        return false;
    }
    const char* tagn = tag->data.symbol.name;
    if (get_attributes(tagn) & ATTR_PROTECTED) {
        mth_message(head, "write", "Tag %s is Protected.", tagn);
        return false;
    }
    /* tag /: MakeBoxes[lhs, form] := rhs and tag /: Format[lhs] := rhs define
     * FormatValues of tag, not UpValues: the reference never consults them in
     * evaluation, only in output. */
    if (is_format_target(lhs) && tag_occurs_in(lhs, tagn)) {
        install_format_rules(tagn, lhs, rhs, is_delayed);
        return true;
    }
    /* tag /: N[x] = v: an N-value of tag (the reference's NValues). */
    if (is_n_target(lhs) && tag_occurs_in(lhs->data.function.args[0], tagn)) {
        install_n_value(tagn, lhs, rhs, is_delayed);
        return true;
    }
    TagPos pos = classify_tag_position(tagn, lhs);
    if (pos == TAG_NOTFOUND) {
        /* Not at a position a rule can be keyed on. If the tag still occurs in
         * lhs it is merely too deep (h[sq[k[s_]]]: k sits two levels down, so
         * no evaluation of h[...] would ever consult k); otherwise it is absent. */
        char* shown = expr_to_string_outputform(lhs);
        if (tag_occurs_in(lhs, tagn))
            mth_message(head, "tagpos",
                        "Tag %s in %s is too deep for an assigned rule to be found.",
                        tagn, shown ? shown : "");
        else
            mth_message(head, "tagnf", "Tag %s not found in %s.",
                        tagn, shown ? shown : "");
        free(shown);
        return false;
    }
    /* Lift a delayed `body /; test` onto the LHS before storing, like SetDelayed. */
    bool owned; Expr* body;
    Expr* pat = move_rhs_condition(lhs, rhs, is_delayed, &owned, &body);
    if (pos == TAG_OWN)       symtab_add_own_value_as(tagn, pat, body, is_delayed);
    else if (pos == TAG_DOWN) symtab_add_down_value_as(tagn, pat, body, is_delayed);
    else if (pos == TAG_SUB)  symtab_add_sub_value(tagn, pat, body, is_delayed);
    else                      symtab_add_up_value_as(tagn, pat, body, is_delayed); /* TAG_UP */
    if (owned) expr_free(pat);
    return true;
}

/* tag /: lhs =. : remove the matching own/down/up rule on tag. */
static void apply_tag_unset(Expr* tag, Expr* lhs) {
    if (tag->type != EXPR_SYMBOL) {
        mth_message("TagUnset", "sym",
                    "Argument at position 1 is expected to be a symbol.");
        return;
    }
    const char* tagn = tag->data.symbol.name;
    switch (classify_tag_position(tagn, lhs)) {
        case TAG_OWN:  symtab_remove_matching_rule(tagn, lhs, true);  break;
        case TAG_DOWN: symtab_remove_matching_rule(tagn, lhs, false); break;
        case TAG_UP:   symtab_remove_matching_up_value(tagn, lhs);    break;
        default: break;  /* SubValue / not found: quietly nothing to remove */
    }
}
#endif /* UP_VALUES */

/*
 * flatten_sequences:
 * Flattens any Sequence[...] heads found in the arguments of e.
 * Returns true iff the args list was actually rewritten.
 */
static bool flatten_sequences(Expr* e) {
    if (e->type != EXPR_FUNCTION) return false;

    size_t new_count = 0;
    bool found_sequence = false;
    bool found_splice = false;   /* Splice[list, h]: see eval_splice_args */
    for (size_t i = 0; i < e->data.function.arg_count; i++) {
        Expr* arg = e->data.function.args[i];
        if (arg->type == EXPR_FUNCTION && arg->data.function.head->type == EXPR_SYMBOL &&
            arg->data.function.head->data.symbol.name == SYM_Sequence) {
            new_count += arg->data.function.arg_count;
            found_sequence = true;
        } else {
            if (arg->type == EXPR_FUNCTION && arg->data.function.head->type == EXPR_SYMBOL &&
                arg->data.function.head->data.symbol.name == SYM_Splice)
                found_splice = true;
            new_count++;
        }
    }

    /* A lone Sequence[x] changes structure without changing arg_count, so we
     * cannot gate on (new_count == arg_count) -- test for any Sequence head. */
    if (!found_sequence) return found_splice && eval_splice_args(e);
    
    Expr** new_args = malloc(sizeof(Expr*) * new_count);
    size_t k = 0;
    for (size_t i = 0; i < e->data.function.arg_count; i++) {
        Expr* arg = e->data.function.args[i];
        if (arg->type == EXPR_FUNCTION && arg->data.function.head->type == EXPR_SYMBOL &&
            arg->data.function.head->data.symbol.name == SYM_Sequence) {
            for (size_t j = 0; j < arg->data.function.arg_count; j++) {
                new_args[k++] = expr_copy(arg->data.function.args[j]);
            }
            expr_free(arg);
        } else {
            new_args[k++] = arg;
        }
    }
    
    free(e->data.function.args);
    e->data.function.args = new_args;
    e->data.function.arg_count = new_count;
    expr_invalidate_hash(e);   /* Sequence splice rewrote args in place */
    if (found_splice) eval_splice_args(e);
    return true;
}

/*
 * strip_nothing:
 * Removes the special symbol `Nothing` (and any `Nothing[...]` form) from the
 * arguments of a List. Per WL, `Nothing` is the identity element of list
 * construction and vanishes from any list it appears in — the idiom behind
 * `Table[If[cond, x, Nothing], ...]` to conditionally build a list. List-specific:
 * for a non-List head `Nothing` is an ordinary symbol. Returns true iff the args
 * were rewritten.
 */
static bool strip_nothing(Expr* e) {
    if (e->type != EXPR_FUNCTION) return false;
    Expr* head = e->data.function.head;
    if (head->type != EXPR_SYMBOL || head->data.symbol.name != SYM_List) return false;
    size_t n = e->data.function.arg_count;
    size_t keep = 0;
    for (size_t i = 0; i < n; i++) {
        Expr* a = e->data.function.args[i];
        bool is_nothing =
            (a->type == EXPR_SYMBOL && a->data.symbol.name == SYM_Nothing) ||
            (a->type == EXPR_FUNCTION && a->data.function.head->type == EXPR_SYMBOL &&
             a->data.function.head->data.symbol.name == SYM_Nothing);
        if (!is_nothing) keep++;
    }
    if (keep == n) return false;   /* no Nothing present — unchanged */
    Expr** na = (keep > 0) ? (Expr**)malloc(sizeof(Expr*) * keep) : NULL;
    size_t k = 0;
    for (size_t i = 0; i < n; i++) {
        Expr* a = e->data.function.args[i];
        bool is_nothing =
            (a->type == EXPR_SYMBOL && a->data.symbol.name == SYM_Nothing) ||
            (a->type == EXPR_FUNCTION && a->data.function.head->type == EXPR_SYMBOL &&
             a->data.function.head->data.symbol.name == SYM_Nothing);
        if (is_nothing) expr_free(a);
        else            na[k++] = a;
    }
    free(e->data.function.args);
    e->data.function.args = na;
    e->data.function.arg_count = keep;
    expr_invalidate_hash(e);
    return true;
}

/*
 * evaluate_step:
 * Performs exactly one level of evaluation transformation.
 *
 * `changed` is an out-parameter set to true iff a real rewrite fired
 * during this step (M3 §3.4 — eager early-exit fixed-point loop). The
 * outer loop in evaluate() uses it to skip the O(tree) expr_eq compare.
 * See the contract in eval.h. NULL is permitted for callers that do not
 * care about the signal.
 */
/* Does any argument of the call `res` hold a dense buffer (a packed List or a
 * visible NDArray)? The packed-array gate in evaluate_step has work to do only
 * then. */
/* Collect the distinct symbols of `e` (heads included) into names[]. Returns
 * false when there are more than `cap`. */
static bool eval_collect_symbols(const Expr* e, const char** names, size_t* n, size_t cap) {
    if (e->type == EXPR_SYMBOL) {
        for (size_t i = 0; i < *n; i++) if (names[i] == e->data.symbol.name) return true;
        if (*n == cap) return false;
        names[(*n)++] = e->data.symbol.name;
        return true;
    }
    if (e->type != EXPR_FUNCTION) return true;
    if (!eval_collect_symbols(e->data.function.head, names, n, cap)) return false;
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (!eval_collect_symbols(e->data.function.args[i], names, n, cap)) return false;
    return true;
}

static bool eval_mentions_any(const Expr* e, const char* const* names, size_t n) {
    if (e->type == EXPR_SYMBOL) {
        for (size_t i = 0; i < n; i++) if (names[i] == e->data.symbol.name) return true;
        return false;
    }
    if (e->type != EXPR_FUNCTION) return false;
    if (eval_mentions_any(e->data.function.head, names, n)) return true;
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (eval_mentions_any(e->data.function.args[i], names, n)) return true;
    return false;
}

/* Does `e` mention any symbol that occurs in `lhs` (the assignment target)?
 * Conservatively true for a left side with very many symbols. */
static bool eval_mentions_symbols_of(const Expr* e, const Expr* lhs) {
    const char* names[32];
    size_t n = 0;
    if (!eval_collect_symbols(lhs, names, &n, 32)) return true;
    return eval_mentions_any(e, names, n);
}

/* Are `res`'s arguments still exactly the `n` pointers in `snap`? */
static inline bool eval_args_untouched(const Expr* res, Expr* const* snap, size_t n) {
    return res->type == EXPR_FUNCTION && res->data.function.arg_count == n &&
           (n == 0 || memcmp(snap, res->data.function.args, sizeof(Expr*) * n) == 0);
}

/* Does the sparse gate have anything to convert in this call: a SparseArray
 * argument, or (SPG_DEEP) one inside a List argument, two levels down? */
static bool eval_sparse_in(const Expr* e, int depth) {
    if (e->type == EXPR_SPARSE) return true;
    if (depth == 0 || e->type != EXPR_FUNCTION || e->data.function.head->type != EXPR_SYMBOL ||
        e->data.function.head->data.symbol.name != SYM_List)
        return false;
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (eval_sparse_in(e->data.function.args[i], depth - 1)) return true;
    return false;
}
static bool eval_sparse_gate_applies(const Expr* res, int mode_bits) {
    int depth = (mode_bits & SPG_DEEP) ? 2 : 0;
    for (size_t i = 0; i < res->data.function.arg_count; i++)
        if (eval_sparse_in(res->data.function.args[i], depth)) return true;
    return false;
}

static inline bool eval_any_ndarray_arg(const Expr* res) {
    for (size_t i = 0; i < res->data.function.arg_count; i++)
        if (res->data.function.args[i]->type == EXPR_NDARRAY) return true;
    return false;
}

/* Is `v` a value that evaluate_step provably leaves unchanged, so that the
 * fixed-point loop may stop on it without running a confirming step?
 *   - Numeric and string atoms: evaluate_step returns expr_copy for them.
 *     Symbols are excluded (an OwnValue may rewrite them), and so are NDArray
 *     and compiled objects (rare, so not worth reasoning about the pack gate).
 *   - Complex[re, im] with both parts machine Reals and im != 0: builtin_complex
 *     returns NULL for it (only a zero imaginary part collapses). Its parts are
 *     atoms, so no UpValue can attach. Complex is Protected, so no DownValue can
 *     fire either. This is the shape every machine-complex special-function
 *     result takes.
 * Cheap: a type switch, and for Complex two more loads. */
static inline bool eval_is_terminal_value(const Expr* v) {
    switch (v->type) {
        case EXPR_INTEGER: case EXPR_REAL: case EXPR_BIGINT: case EXPR_STRING:
#ifdef USE_MPFR
        case EXPR_MPFR:
#endif
            return true;
        case EXPR_FUNCTION: {
            const Expr* h = v->data.function.head;
            if (h->type != EXPR_SYMBOL || h->data.symbol.name != SYM_Complex ||
                v->data.function.arg_count != 2)
                return false;
            const Expr* re = v->data.function.args[0];
            const Expr* im = v->data.function.args[1];
            return re->type == EXPR_REAL && im->type == EXPR_REAL &&
                   im->data.real != 0.0;
        }
        default:
            return false;
    }
}

/* ---- Swallowed unwinds -----------------------------------------------------
 *
 * A Throw or Abort is a sentinel riding the ordinary return path.  The argument
 * loop passes one upward, but a builtin (or the matcher) that evaluates user
 * code and only inspects the result -- a Select criterion, a Sort comparator,
 * a _?test or /; guard -- would otherwise drop it: Catch[Select[{1,2,3},
 * If[# == 2, Throw[9], True] &]] gave {1, 3}.  So where a USER CALLBACK is
 * evaluated (eval_callback / eval_park_unwind: the matcher's guards and tests,
 * the predicate / comparator / key-function sites of the list builtins), a
 * sentinel it yields is parked in g_pending_unwind.  Parking is deliberately
 * not done for every evaluate(): Mathilda's symbolic engines evaluate user
 * expressions speculatively (Sum probes its summand symbolically even for
 * numeric bounds) and rightly ignore what such a probe throws.  Each
 * evaluate_step frame saves and clears the slot on entry, and on exit, if the
 * slot is set but the frame's own result is not a sentinel (it was swallowed),
 * the parked sentinel becomes the result.  A genuine consumer -- Catch taking a
 * Throw, CheckAbort taking an Abort -- clears the slot (eval_unwind_consumed). */
static Expr* g_pending_unwind = NULL;

void eval_unwind_consumed(void) {
    if (g_pending_unwind) { expr_free(g_pending_unwind); g_pending_unwind = NULL; }
}

static Expr* evaluate_step_core(Expr* e, bool* changed);

Expr* evaluate_step(Expr* e, bool* changed) {
    Expr* saved = g_pending_unwind;
    g_pending_unwind = NULL;
    Expr* r = evaluate_step_core(e, changed);
    if (g_pending_unwind) {
        if (r && eval_is_inflight_unwind(r)) {
            expr_free(g_pending_unwind);
        } else {
            if (r) expr_free(r);
            r = g_pending_unwind;
            if (changed) *changed = true;
        }
        g_pending_unwind = NULL;
    }
    g_pending_unwind = saved;
    return r;
}

static Expr* evaluate_step_core(Expr* e, bool* changed) {
    /* Local sink so we can write through `*changed` unconditionally. */
    bool sink = false;
    if (!changed) changed = &sink;
    *changed = false;

    if (!e) return NULL;

    switch (e->type) {        /* Atomics evaluate to themselves */
        case EXPR_INTEGER:
        case EXPR_REAL:
        case EXPR_STRING:
        case EXPR_BIGINT:
        case EXPR_NDARRAY:        /* dense ndarray: an atomic value */
        case EXPR_COMPILED:       /* compiled function object: an atomic value */
        case EXPR_SPARSE:         /* SparseArray raw object: never re-evaluated */
#ifdef USE_MPFR
        case EXPR_MPFR:
#endif
            return expr_copy(e);

        case EXPR_SYMBOL: {
            /* Check for immediate assignments (OwnValues) like x = 5 */
            Expr* own = apply_own_values(e);
            if (own) { *changed = true; return own; }
            return expr_copy(e);
        }
            
        case EXPR_FUNCTION: {
            /* 1. Evaluate the head recursively (e.g. f[x][y]).
             * Refcount sharing means evaluate() returns the same pointer
             * when nothing rewrote — pointer-inequality is a sound
             * "head changed" signal, and we only fall through to the
             * structural compare if needed. */
            Expr* orig_head = e->data.function.head;
            Expr* head = evaluate(orig_head);
            if (head != orig_head) *changed = true;

            uint32_t attrs = ATTR_NONE;
            /* Phase 3a: resolve the head's definition ONCE, then thread it
             * through attribute lookup, DownValue dispatch, and builtin dispatch
             * below -- instead of re-resolving (re-hashing) the same head up to
             * three times per node per evaluation pass. The def node is stable
             * (Phase 2: never freed/reallocated), so the pointer stays valid
             * even if the symbol is redefined mid-step; each use reads its
             * fields fresh. */
            SymbolDef* hdef = NULL;
            if (head->type == EXPR_SYMBOL) {
                /* Phase 3b: the symbol node caches its resolved def. First touch
                 * looks it up (one hash on the unified table); every later pass
                 * over this same node is a pointer load. Writing the cache on a
                 * possibly-shared node is safe -- it is benign metadata, and all
                 * sharers are the same symbol with the same def. */
                hdef = head->data.symbol.def;
                if (!hdef) {
                    hdef = symtab_get_def(head->data.symbol.name);
                    head->data.symbol.def = hdef;
                }
                attrs = get_attributes_def(hdef);
            } else if (head->type == EXPR_FUNCTION) {
                attrs = pure_function_attributes(head);
            }
            
            /* 2. Handle 'Hold' attributes.
             *
             * HoldAllComplete is like HoldAll but additionally suppresses
             * Sequence flattening, Unevaluated stripping, Flat flattening,
             * and (eventually) UpValues lookup. Inside a HoldAllComplete
             * head, Evaluate[expr] does NOT force evaluation. Built-ins and
             * DownValues attached to the head still apply -- this is what
             * lets Hold-style heads do useful work (e.g. Length[Hold[a,b,c]]
             * via the down code on Hold). */
            bool hold_all_complete = (attrs & ATTR_HOLDALLCOMPLETE) != 0;

            /* Evaluate arguments unless suppressed by HoldFirst, HoldRest,
             * HoldAll, or HoldAllComplete.
             *
             * Pointer-identity check against the original arg signals
             * "did sub-evaluation change anything"; with refcount sharing
             * (M3 phase-2) evaluate() returns the same pointer for
             * already-stable inputs. Stripping an Evaluate[] wrapper is
             * itself a rewrite even if the wrapped expression evaluates
             * to itself — flag it explicitly. */
            size_t argc = e->data.function.arg_count;
            Expr* new_args_stack[EVAL_SMALL_ARGS];
            Expr** new_args = (argc <= EVAL_SMALL_ARGS)
                                  ? new_args_stack
                                  : malloc(sizeof(Expr*) * argc);
            /* Pointers of Unevaluated[...] wrappers that landed in a HELD slot.
             * Per WMA, such wrappers are NOT stripped below (the argument was
             * never going to be evaluated, so Unevaluated has nothing to do).
             * We track by pointer identity, not index, because flatten_sequences
             * can shift argument positions before the strip pass runs; held
             * Unevaluated nodes are never Sequence, so their pointers survive it.
             * Held Unevaluated wrappers are vanishingly rare (essentially never
             * in numeric code), so this tracking array is allocated lazily on
             * first sighting rather than on every function node. */
            Expr** held_uneval = NULL;
            size_t held_uneval_count = 0;
            for (size_t i = 0; i < e->data.function.arg_count; i++) {
                bool hold = hold_all_complete;
                if (i == 0 && (attrs & ATTR_HOLDFIRST)) hold = true;
                if (i > 0 && (attrs & ATTR_HOLDREST)) hold = true;

                /* Track whether this slot was produced by *evaluation* (vs a
                 * held plain copy). Only an evaluated argument can be an
                 * in-flight Throw that must short-circuit the call; a held
                 * copy of a literal Throw[...] (e.g. Hold[Throw[1]]) is inert
                 * and must NOT short-circuit. */
                bool arg_evaluated = false;
                Expr* orig_arg = e->data.function.args[i];
                if (hold) {
                    /* Check for Evaluate[...] - overrides HoldFirst/HoldRest/HoldAll
                     * but NOT HoldAllComplete. Arity 1 forces the single argument;
                     * any other arity forces each argument and splices them via a
                     * Sequence, which flatten_sequences resolves below (so
                     * Evaluate[] vanishes and Evaluate[a, b] -> a, b). */
                    if (!hold_all_complete &&
                        orig_arg->type == EXPR_FUNCTION &&
                        orig_arg->data.function.head->type == EXPR_SYMBOL &&
                        orig_arg->data.function.head->data.symbol.name == SYM_Evaluate) {
                        size_t ne = orig_arg->data.function.arg_count;
                        if (ne == 1) {
                            new_args[i] = evaluate(orig_arg->data.function.args[0]);
                        } else {
                            Expr** seq = malloc(sizeof(Expr*) * ne);
                            for (size_t k = 0; k < ne; k++)
                                seq[k] = evaluate(orig_arg->data.function.args[k]);
                            new_args[i] = expr_new_function(
                                expr_new_symbol(SYM_Sequence), seq, ne);
                            free(seq);
                        }
                        arg_evaluated = true;
                        *changed = true; /* Evaluate wrapper stripped */
                    } else {
                        new_args[i] = expr_copy(orig_arg);
                    }
                    /* Record a held Unevaluated[...] wrapper so the strip pass
                     * below leaves it intact (WMA: held wrappers are not removed).
                     * Covers both the plain held copy and the
                     * Hold[Evaluate[Unevaluated[x]]] case where the override
                     * evaluates to Unevaluated[x] in a held slot. */
                    if (new_args[i]->type == EXPR_FUNCTION &&
                        new_args[i]->data.function.head->type == EXPR_SYMBOL &&
                        new_args[i]->data.function.head->data.symbol.name == SYM_Unevaluated &&
                        new_args[i]->data.function.arg_count == 1) {
                        /* Lazily allocate: capacity argc is always enough since
                         * there is at most one held wrapper per argument slot. */
                        if (!held_uneval)
                            held_uneval = malloc(sizeof(Expr*) * (argc ? argc : 1));
                        if (held_uneval) held_uneval[held_uneval_count++] = new_args[i];
                    }
                } else {
                    /* Atom fast path. A raw atom (number/string/ndarray/
                     * compiled) always evaluates to itself -- evaluate() would
                     * just expr_copy it after paying the recursion-depth, trace,
                     * deadline-check and fixed-point-loop overhead per call.
                     * Skipping straight to expr_copy is exactly equivalent (same
                     * pointer, so `*changed` stays false, and an atom is never an
                     * in-flight Throw/Goto so the sentinel check is moot). This
                     * is the hot path when a large List of numbers is
                     * re-evaluated -- e.g. every pass of `list //. rule` re-walks
                     * the whole (mostly unchanged) list. SYMBOL is excluded: it
                     * may carry an OwnValue and must go through evaluate(). */
                    switch (orig_arg->type) {
                        case EXPR_INTEGER: case EXPR_REAL: case EXPR_STRING:
                        case EXPR_BIGINT: case EXPR_NDARRAY: case EXPR_COMPILED:
                        case EXPR_SPARSE:
#ifdef USE_MPFR
                        case EXPR_MPFR:
#endif
                            new_args[i] = expr_copy(orig_arg);
                            continue;   /* atoms are never Throw/Goto sentinels */
                        default: break;
                    }
                    if (i == 0 && head->type == EXPR_SYMBOL &&
                        head->data.symbol.name == SYM_SparseArray &&
                        orig_arg->type == EXPR_FUNCTION && orig_arg->data.function.arg_count == 2 &&
                        orig_arg->data.function.head->type == EXPR_SYMBOL &&
                        orig_arg->data.function.head->data.symbol.name == SYM_Rule)
                        new_args[i] = eval_sparse_rule_arg(orig_arg);
                    else
                        new_args[i] = evaluate(orig_arg);
                    arg_evaluated = true;
                    if (new_args[i] != orig_arg) *changed = true;
                }

                /* Catch/Throw and Goto/Label: an evaluated argument that is an
                 * in-flight Throw or Goto short-circuits the entire call. Free
                 * the siblings produced so far and the evaluated head, then hand
                 * the sentinel up through the normal return path -- this frame
                 * and every enclosing frame still run their own cleanup (this is
                 * why longjmp is not used). `res` has not been built yet, so
                 * nothing else is owned here. Runs before `res` is built and
                 * before Orderless sorting, so first-throw-wins holds for
                 * Plus/Times arguments too. A Goto bubbles the same way until an
                 * enclosing CompoundExpression consumes it (see
                 * builtin_compoundexpression). */
                if (arg_evaluated && (eval_is_inflight_unwind(new_args[i]) ||
                                      eval_is_inflight_goto(new_args[i]))) {
                    Expr* sentinel = new_args[i];
                    for (size_t j = 0; j < i; j++) expr_free(new_args[j]);
                    if (new_args != new_args_stack) free(new_args);
                    free(held_uneval);   /* free(NULL) is a no-op */
                    expr_free(head);
                    *changed = true;
                    return sentinel;
                }
            }
            
            Expr* res = expr_new_function(head, new_args, argc);
            if (new_args != new_args_stack) free(new_args);

            /* The node as argument evaluation left it. If nothing after this
             * point touches it -- no Sequence splice, Unevaluated strip, Flat
             * (which flattens one level per step), Orderless, gate or builtin
             * rewriting arguments in place (And stores its evaluated held
             * arguments) -- and it comes to rest, a further step would only
             * replay this one: see g_step_settled_depth. */
            enum { SETTLE_SNAP_MAX = 64 };
            Expr* settle_snap[SETTLE_SNAP_MAX];
            const size_t settle_argc = res->data.function.arg_count;
            const bool settle_ok = settle_argc <= SETTLE_SNAP_MAX;
            if (settle_ok && settle_argc > 0)   /* args is NULL for f[] */
                memcpy(settle_snap, res->data.function.args, sizeof(Expr*) * settle_argc);
            /* Messages fired so far: a head that reports a message while
             * declining has rejected the call (see the post-gate). */
            const unsigned long step_msgs_before = mth_msg_fired_count();

            /* Trace: snapshot the reassembled f[evaluated_args] form, before
             * flatten/Flat/Listable/Orderless and before any rule fires, so the
             * trace can show e.g. `8 + 16 + 1` rather than the pre-arg-eval
             * `2^3 + 4^2 + 1`. Stored on the current frame (not a global) so a
             * builtin that itself re-enters evaluate() cannot clobber it. */
            if (TRACE_ACTIVE() && g_trace_top) {
                if (g_trace_top->pending_reassembled)
                    expr_free(g_trace_top->pending_reassembled);
                g_trace_top->pending_reassembled = expr_copy(res);
            }
            
    /* 2.5 Flatten Sequences - must happen before attributes.
     * Suppressed for heads carrying SequenceHold (e.g. Set/SetDelayed/Rule/
     * RuleDelayed, so assignments and rules can return Sequence objects) and for
     * HoldAllComplete heads (HoldAllComplete implies SequenceHold). Any
     * user-defined head with SequenceHold is honored automatically. */
    if (hold_all_complete || (attrs & ATTR_SEQUENCEHOLD)) {
        /* SequenceHold / HoldAllComplete leaves Sequence intact */
    } else {
        if (flatten_sequences(res)) *changed = true;
    }

    /* Strip `Nothing` from Lists (the list-construction identity element). Runs
     * for every evaluated List; genuinely held Lists are not evaluated, so their
     * `Nothing`s survive until released, matching WL. */
    if (strip_nothing(res)) *changed = true;

            /* 2.6 Strip Unevaluated wrappers.
             * f[Unevaluated[expr]] passes expr (unevaluated) to f, with the
             * wrapper removed. This runs AFTER flatten_sequences so that any
             * Sequence[...] directly inside Unevaluated is preserved (e.g.
             * Length[Unevaluated[Sequence[a,b]]] gives 2 because Sequence is
             * not flattened into Length's argument list).
             * Per WMA semantics, the wrapper is stripped ONLY in positions that
             * would otherwise be evaluated (non-held slots): the wrapper's job is
             * to temporarily hold an argument that a non-Hold head would evaluate.
             * In genuinely held slots -- HoldFirst/HoldRest/HoldAll on that
             * position, or content forced there by Evaluate -- the wrapper is
             * left intact (so f[Unevaluated[1+2]] with f HoldAll stays
             * f[Unevaluated[1+2]], and Hold[Evaluate[Unevaluated[1+2]]] stays
             * Hold[Unevaluated[1+2]]). Those held wrappers were recorded by
             * pointer in held_uneval above. HoldAllComplete heads keep every
             * wrapper and skip this pass entirely. Stripping merely removes the
             * wrapper; the exposed content is not itself evaluated, so the head's
             * hold attributes are respected on the next pass. */
            if (!hold_all_complete) {
                for (size_t i = 0; i < res->data.function.arg_count; i++) {
                    Expr* arg = res->data.function.args[i];
                    if (arg->type == EXPR_FUNCTION &&
                        arg->data.function.head->type == EXPR_SYMBOL &&
                        arg->data.function.head->data.symbol.name == SYM_Unevaluated &&
                        arg->data.function.arg_count == 1) {
                        bool held_here = false;
                        for (size_t k = 0; k < held_uneval_count; k++) {
                            if (held_uneval[k] == arg) { held_here = true; break; }
                        }
                        if (held_here) continue; /* held slot: keep wrapper */
                        Expr* stripped = expr_copy(arg->data.function.args[0]);
                        expr_free(arg);
                        res->data.function.args[i] = stripped;
                        *changed = true; /* Unevaluated wrapper removed */
                    }
                }
            }
            free(held_uneval);

            /* 2.7 PACKED-LIST TRANSPARENCY GATE.
             *
             * A packed list is an ordinary List stored as a dense buffer -- one
             * EXPR_NDARRAY node with no args[] array (see pack.h). Roughly 7100
             * sites in this tree read data.function.args directly, so a head
             * that has not been taught about packing either crashes on one
             * (arg_count aliases the payload's data pointer) or, worse, takes
             * its "not a function" branch and answers confidently wrongly --
             * Count[packed, _Real] would say 0.
             *
             * So a packed list never reaches such a head: materialise here, and
             * every unaware consumer is correct by construction. That covers
             * AtomQ, Count, Level, Position, Cases, ReplaceAll, ListQ, Insert,
             * and any user DownValue with a structural pattern, none of which
             * need a line of packing-specific code.
             *
             * THE POSITION IS FORCED, in both directions:
             *   - after flatten_sequences and the Unevaluated strip, because
             *     f[Sequence[packed]] has no packed argument until the splice
             *     happens (both passes key on EXPR_FUNCTION, so a packed node
             *     travels through them untouched);
             *   - before Flat, Listable, Orderless, DownValues and the builtin.
             *     Listable is the subtle one: has_list_arg tests for an
             *     EXPR_FUNCTION headed List and returns FALSE for a packed
             *     node, so an unaware Listable head would silently skip
             *     threading altogether. Gating first fixes that with no change
             *     to has_list_arg -- while an AWARE head still skips threading,
             *     which is what lets its ndarray kernel fire below.
             *
             * Sits outside the `head is a symbol` block, so f[x][y] is covered
             * too. Exactly ONE non-symbol head is exempt: a pure Function/&.
             * There, substitution drops the packed value into the body and
             * every head in the body is gated on the next pass, so nothing can
             * come to rest nested -- and materialising here would cost
             * (#^2 &)[packed] its fast path.
             *
             * The exemption must not be widened to "head is any EXPR_FUNCTION".
             * An unevaluated application like gg[xx][packed] has no builtin to
             * run, so it comes to REST with a packed node inside a plain
             * function node, and the shallow gate at the enclosing level never
             * looks that deep. Measured, with the exemption too broad:
             * Count[gg[xx][p], _Real, 2] gave 0 against 4, LeafCount 3 against
             * 7, Cases {} against the elements, and a ReplaceAll did nothing.
             *
             * The `!down_values` term is load-bearing, not caution: Protected
             * is opt-in per symbol, so even an aware head can carry a user
             * DownValue, and the matcher cannot descend a buffer. */
            /* A List being BUILT out of packed rows: absorb them into one array
             * instead of materialising every one.
             *
             * `List` is not, and must not be, a packed-aware head -- a plain List
             * node holding EXPR_NDARRAY elements is precisely the malformed shape
             * the gate below exists to prevent, and ~7100 sites would read it
             * wrongly. But the answer is not to unpack the rows either. It is to
             * pack the WHOLE thing: pack_sniff already absorbs already-packed rows
             * (that is how Table[i j, {i,300},{j,300}] becomes rank 2), so a list
             * of n packed vectors is exactly a rank-2 buffer, and the result is a
             * genuine packed list that every consumer handles correctly.
             *
             * Without this, a function that returns several arrays destroys them
             * all on the way out, and the caller pays on the NEXT operation, not
             * this one. Measured on the N-body step, which returns six 1024-vectors:
             * the first call took 42 ms with a packed argument and the second
             * 5.75 s, because its own result came back as six plain Lists and
             * Outer, the elementwise arithmetic and the reductions all fell off the
             * buffer path together. 137x, from the return statement.
             *
             * pack_offer declines (ragged rows, mixed exact/inexact, under the
             * threshold) by returning the node unchanged, and the gate below then
             * materialises the rows exactly as before -- so this can only turn a
             * materialise into a pack, never into anything new. */
            if (pack_any_created() && head->type == EXPR_SYMBOL &&
                head->data.symbol.name == SYM_List) {
                bool any_packed = false;
                for (size_t i = 0; i < res->data.function.arg_count; i++)
                    if (is_packed_list(res->data.function.args[i])) { any_packed = true; break; }
                if (any_packed) {
                    Expr* packed = pack_offer(res);
                    if (packed != res) { *changed = true; return packed; }
                }
            }

            /* The whole gate below is a no-op unless some argument IS a buffer:
             * pack_lift_listable_args needs an NDArray to lift toward, and the
             * materialise loop only touches packed arguments. Once any array has
             * been packed in the session (pack_any_created), every call to every
             * Listable head would otherwise pay the scans: has_nd_kernel,
             * has_list_arg, dv_binds_opaquely, and one is_packed_list call per
             * argument. A scalar F[x] is the hot case that pays for nothing here,
             * so one type scan decides it. */
            if (pack_any_created() && eval_any_ndarray_arg(res)) {
                bool pure_fn = head->type == EXPR_FUNCTION &&
                               head->data.function.head->type == EXPR_SYMBOL &&
                               head->data.function.head->data.symbol.name == SYM_Function;
                /* A Listable head holding BOTH a plain List argument and a packed
                 * one. Threading fires (has_list_arg sees the plain List) and
                 * treats the buffer as a scalar, so it broadcasts instead of
                 * threading elementwise: NDArray[{1.,2.,3.}] * {10,20,30} gave
                 * the 3x3 outer product where {1.,2.,3.} * {10,20,30} gives
                 * {10., 40., 90.}. Materialising makes every list argument thread
                 * the same way, at the cost of the buffer path for a mixed call.
                 *
                 * The bug is older than packing -- it was always there for an
                 * explicit NDArray[...] -- but automatic packing is what makes it
                 * reachable from code that never mentions an array, so it belongs
                 * to the gate. It surfaced as ListConvolve's reference
                 * computation being off by 103 once RandomReal[{-1,1},{40,40}]
                 * started packing. */
                /* A CompiledFunction object is aware, and had to be told so
                 * explicitly: its head is an EXPR_COMPILED with no SymbolDef, so
                 * it read as unaware and the gate materialised the very buffer
                 * the compiled boundary exists to borrow. That cost the whole
                 * Phase-5 win invisibly -- f[Range[1., 200000.]] measured 75x
                 * slower than f[NDArray[Range[1., 200000.]]], the two differing
                 * only in `present_as`. compiled_function_apply reads the buffer
                 * directly and returns one with the argument's presentation; if
                 * it declines, cf_fallback re-runs the body through the
                 * evaluator, where every head inside is gated on the next pass. */
                bool compiled_head = head->type == EXPR_COMPILED;
                /* An InterpolatingFunction object applied to a packed array of
                 * query points is aware for the same reason a CompiledFunction
                 * is: interp_apply's vectorised 1-D path reads the buffer
                 * directly (a batch ifn[array] is scipy's cs(array)), so
                 * materialising the points into 10^5 boxed Exprs first is pure
                 * loss.  Integer query points are read exactly (ndt_get to 2^53)
                 * and give the same real result as materialising, so int64 is
                 * fine too. */
                bool interp_head = head->type == EXPR_FUNCTION &&
                                   head->data.function.head->type == EXPR_SYMBOL &&
                                   head->data.function.head->data.symbol.name ==
                                       SYM_InterpolatingFunction;
                /* MIXED PACKED/PLAIN: pack the List UP, do not unpack the buffer
                 * DOWN.
                 *
                 * `listable_mixed` below exists because apply_listable walks its
                 * arguments as ordinary Lists, so a buffer alongside a plain List
                 * has to be materialised for threading to work. That is correct
                 * and it is also ruinous: the cost is set by the LARGEST operand,
                 * and one plain operand forfeits the buffer path for the whole
                 * call. On 10^6 elements `a + b` measured 418 ms with `b` plain
                 * against 2.1 ms with both packed -- the same 200x as the fully
                 * unpacked `b + b`, from a single unpacked argument. Plain
                 * numeric Lists are easy to come by: anything under
                 * PACK_MIN_ELEMENTS, a literal, or any producer without a packed
                 * path.
                 *
                 * Packing instead is value-preserving by contract (pack.h), and
                 * it hands the call to the head's own buffer kernel rather than
                 * to threading. Restricted to shapes those kernels actually
                 * cover: equal shapes, which every ND kernel handles, or -- for a
                 * head claiming packed_broadcast_ok -- a lower-rank operand whose
                 * dims prefix the higher one's. Any other combination is left
                 * exactly as it was, to thread. */
                /* A head with a matching-arity NDArray kernel is aware of a
                 * packed buffer even when it ALSO carries DownValues: the kernel
                 * runs at step 4b, AFTER DownValues (step 4), so a guarded
                 * symbolic rule still gets first crack and a numeric array falls
                 * through to the kernel. Without this, any rule-bearing head --
                 * the whole Bessel family -- materialised its buffer and threaded
                 * scalar by scalar (BesselI/BesselK over 10^6 elements cost tens
                 * of seconds). The .m rules are guarded with !NDArrayQ[z] so they
                 * decline on a buffer rather than fire element-wise. */
                size_t nd_argc = res->data.function.arg_count;
                bool has_nd_kernel = hdef &&
                    ((nd_argc == 1 && hdef->ndarray_unary_kernel) ||
                     (nd_argc == 2 && hdef->ndarray_binary_kernel));
                if ((attrs & ATTR_LISTABLE) && hdef && hdef->packed_aware &&
                    (!hdef->down_values || has_nd_kernel) && pack_any_created())
                    pack_lift_listable_args(res, hdef->packed_broadcast_ok != 0,
                                            hdef->packed_int64_ok != 0, changed);
                bool listable_mixed = (attrs & ATTR_LISTABLE) && has_list_arg(res);
                bool aware = ((hdef && hdef->packed_aware
                                    && (!hdef->down_values || has_nd_kernel))
                              || pure_fn || compiled_head || interp_head
                              || dv_binds_opaquely(hdef)) && !listable_mixed;
                /* A head that binds opaquely is exact on an int64 buffer for
                 * exactly the reason it is aware at all: it reads no element.
                 * `f[x_] := body` binds the whole value and substitutes it, and
                 * whatever the body then does is gated at the next head.
                 *
                 * Without this the exemption covered only float64, because a
                 * user symbol's packed_int64_ok is false -- so an INTEGER grid
                 * still materialised at every helper call. That is the whole of
                 * the vectorised Game of Life benchmark, whose grid is integer:
                 * `probe[q_] := NDArrayQ[q]` answered False for a packed integer
                 * argument while answering True for a real one. */
                bool int64_ok = pure_fn || compiled_head || interp_head ||
                                (hdef && hdef->packed_int64_ok) ||
                                dv_binds_opaquely(hdef);
                for (size_t i = 0; i < res->data.function.arg_count; i++) {
                    Expr* pa = res->data.function.args[i];
                    if (!is_packed_list(pa)) continue;
                    /* An int64 buffer needs the stronger claim. Most of the ND
                     * layer reads through ndt_get, which is exact only to 2^53
                     * and was written when only Compile[] could make an integer
                     * buffer -- so without this, Total[{1,2,3}] came back as
                     * 6. and Sin[{1,2,3}] as {0,0,0}. Materialising is always
                     * an option; being wrong is not. */
                    if (aware &&
                        (pa->data.ndarray.dtype != NDT_INT64 || int64_ok))
                        continue;
                    pack_gate_report(head, pa);          /* MATHILDA_PACK_DIAG=gate */
                    res->data.function.args[i] = ndarray_to_nested_list(pa);
                    expr_free(pa);
                    *changed = true;
                    g_pack_gate_ticks++;
                }
            }

            /* A VISIBLE NDArray beside a plain List that the lift above could
             * not absorb -- unequal length, or a symbolic element. The gate
             * materialises a packed argument in that position but by design
             * never touches a visible one, so threading below treated the array
             * as a scalar and broadcast it: NDArray[{1., 2., 3.}] + {1, 2} gave
             * two whole arrays instead of Thread::tdlen, and
             * NDArray[{1., 2.}] + {x, y} gave {nd + x, nd + y} instead of
             * {1. + x, 2. + y}. Delist it, so the call threads exactly as the
             * List would. A same-shape numeric List was already lifted to a
             * buffer above and never reaches here. */
            if ((attrs & ATTR_LISTABLE) && eval_any_ndarray_arg(res) && has_list_arg(res)) {
                for (size_t i = 0; i < res->data.function.arg_count; i++) {
                    Expr* pa = res->data.function.args[i];
                    if (!is_ndarray(pa) || is_packed_list(pa)) continue;
                    res->data.function.args[i] = ndarray_to_nested_list(pa);
                    expr_free(pa);
                    *changed = true;
                }
            }

            /* The SparseArray transparency gate (sparse/sparse_gate.c): a head
             * that opted in by sparse_mode sees an array argument densified,
             * densified and re-sparsified, or split into its rows -- after its
             * registered native kernel (sparse_native), if any, declines. A
             * builtin that reads sparse arguments itself (SPG_NATIVE) gets
             * the gate only if it declines -- see the post-gate below. */
            if (hdef && hdef->sparse_native && sp_any_created() &&
                eval_sparse_gate_applies(res, hdef->sparse_mode)) {
                /* A native sparse kernel (sparse/sparse_*.c) gets first try. */
                Expr* nat = ((SpNativeFn)hdef->sparse_native)(res);
                if (nat == res) {
                    /* The call itself back (expr_copy): the kernel settled it
                     * as it stands, after its own message (Dot::dotsh). It
                     * rests here -- neither the gate nor a further step
                     * re-asks and repeats the message. */
                    expr_free(nat);
                    g_step_settled_depth = eval_recursion_depth;
                    return res;
                }
                if (nat) { expr_free(res); *changed = true; return nat; }
            }
            if (hdef && (hdef->sparse_mode & SPG_MODE_MASK) &&
                !(hdef->sparse_mode & SPG_NATIVE) && sp_any_created() &&
                eval_sparse_gate_applies(res, hdef->sparse_mode)) {
                bool rest = false;
                g_trace_suppress++;
                Expr* g = sp_gate(res, hdef->sparse_mode, &rest);
                g_trace_suppress--;
                if (g) { expr_free(res); *changed = true; return g; }
                if (rest) {
                    if (settle_ok && eval_args_untouched(res, settle_snap, settle_argc))
                        g_step_settled_depth = eval_recursion_depth;
                    return res;
                }
            }

            /* 3. Apply structural and semantic attributes.
             * Order follows Withoff §3.1: Flat → Sequence (already done above) →
             * Listable → Orderless. Flat must run before Listable so that lists
             * exposed by flattening get threaded, e.g. Plus[Plus[a,{1,2}],3]
             * → Plus[a,{1,2},3] → {Plus[a,1,3], Plus[a,2,3]}. */

            /* Flat: associative flattening (requires symbolic head, suppressed by HoldAllComplete) */
            if (head->type == EXPR_SYMBOL && (attrs & ATTR_FLAT) && !hold_all_complete) {
                /* head is an EXPR_SYMBOL, so its name is already interned —
                 * call the core directly and skip the redundant hash. */
                if (eval_flatten_args_interned(res, head->data.symbol.name)) *changed = true;
            }

            /* Listable: automatic threading over Lists, else over the values
             * of Association arguments (assoc_thread_listable; the List case
             * wins when both are present, so a List stays the outer level).
             * One pass classifies the arguments, so the List path pays only a
             * pointer compare per argument for the association check. */
            int la_kind = (attrs & ATTR_LISTABLE) ? listable_arg_kinds(res) : 0;
            if (la_kind & LA_SPARSE) {
                g_trace_suppress++;
                Expr* sp_res = listable_sparse(res, la_kind);
                g_trace_suppress--;
                if (sp_res == res) {
                    /* settled after its message: no List threading, no re-ask */
                    expr_free(sp_res);
                    g_step_settled_depth = eval_recursion_depth;
                    return res;
                }
                if (sp_res) {
                    expr_free(res);
                    *changed = true; /* threaded over (or densified) SparseArrays */
                    return sp_res;
                }
            }
            if (la_kind == LA_ASSOC) {
                Expr* assoc_res = assoc_thread_listable(res);
                if (assoc_res) {
                    expr_free(res);
                    *changed = true;
                    return assoc_res;
                }
            }
            if (la_kind & LA_LIST) {
                /* Trace: threading evaluates the threaded elements internally;
                 * hide those sub-evaluations so x^{1..10} shows as a single
                 * rewrite to {x,x^2,...} rather than a decomposed one. */
                g_trace_suppress++;
                Expr* list_res = apply_listable(res);
                g_trace_suppress--;
                if (list_res) {
                    expr_free(res);
                    *changed = true; /* List threading reshaped the call */
                    return list_res;
                }
            }

            /* SubValues: a call with a compound head f[a][b] (the head already
             * evaluated) tries UpValues of its arguments, then the SubValue
             * rules of its innermost head symbol f. */
            if (head->type == EXPR_FUNCTION && !hold_all_complete) {
                /* The innermost head symbol's def via the node's cached pointer
                 * (no hash): this branch runs for every pure-function call. */
                Expr* inner = head;
                while (inner->type == EXPR_FUNCTION && inner->data.function.head)
                    inner = inner->data.function.head;
                SymbolDef* sdef = NULL;
                if (inner->type == EXPR_SYMBOL) {
                    sdef = inner->data.symbol.def;
                    if (!sdef) sdef = symtab_lookup(inner->data.symbol.name);
                }
#if UP_VALUES
                if (symtab_up_value_count) {
                    Expr* up = apply_up_values(res);
                    if (up) {
                        expr_free(res);
                        *changed = true;
                        return up;
                    }
                }
#endif
                if (sdef && sdef->sub_values) {
                    Expr* sub = apply_sub_values(sdef, res);
                    if (sub) {
                        expr_free(res);
                        *changed = true;
                        return sub;
                    }
                }
            }
            if (head->type == EXPR_SYMBOL) {
                const char* head_name = head->data.symbol.name;

                /* Orderless: commutative sorting. Pre-check whether the
                 * args are already in canonical order so the §3.4 detector
                 * can skip a no-op qsort on stable expressions
                 * (Plus[a,b,c] re-evaluating, etc.).
                 *
                 * Plus and Times canonicalize their OWN argument order inside
                 * builtin_plus/builtin_times (which sort the collapsed group
                 * set before returning), so the generic ORDERLESS sort here is
                 * redundant for them. Worse, on a large sum/product it is the
                 * dominant cost: sorting N raw terms with the heavy polynomial
                 * expr_compare is O(N log N) with a per-compare constant orders
                 * of magnitude above an atom compare, even though the collapsed
                 * result has only a handful of distinct terms (e.g. Total of
                 * 10^5 monomials over ~11 powers). Skip the whole block — the
                 * builtin always runs for arg_count >= 2 (arity 0/1 is handled
                 * earlier) and sorts its small output, and a held/unevaluated
                 * Plus was never sorted anyway. */
                if ((attrs & ATTR_ORDERLESS) &&
                    head_name != SYM_Plus && head_name != SYM_Times) {
                    bool already_sorted = true;
                    for (size_t i = 1; i < res->data.function.arg_count; i++) {
                        if (eval_compare_expr_ptrs(&res->data.function.args[i - 1],
                                                   &res->data.function.args[i]) > 0) {
                            already_sorted = false;
                            break;
                        }
                    }
                    if (!already_sorted) {
                        expr_orderless_sort(res->data.function.args, res->data.function.arg_count);
                        expr_invalidate_hash(res);   /* args reordered in place */
                        *changed = true;
                    }
                }

                /* 4. Apply user-defined DownValues FIRST (Withoff §3.1).
                 * In Mathematica's evaluation pipeline, user-defined
                 * DownValues take precedence over internal "down code"
                 * (built-in implementations). This lets a user override
                 * a built-in for non-Protected symbols, while Protected
                 * symbols (which is most builtin-bearing heads) are
                 * unaffected because apply_assignment refuses to install
                 * DownValues on a Protected target. */
#if UP_VALUES
                /* 3b. UpValues fire BEFORE the enclosing head's DownValues and
                 * builtin (Withoff). HoldAllComplete suppresses them; HoldAll does
                 * not (its held args are still present to be scanned). Skipped for
                 * the assignment primitives so an all-heads upvalue cannot hijack a
                 * held assignment LHS. symtab_up_value_count makes the common
                 * no-upvalue case a single load+branch. */
                if (!hold_all_complete && symtab_up_value_count &&
                    !is_assignment_primitive(head_name)) {
                    /* NOTE: on an Orderless head with upvalue candidates,
                     * apply_up_values may reorder res's args IN PLACE even when it
                     * returns NULL (no rule fired). That is intentional and safe --
                     * the Plus/Times builtin re-sorts below -- but it means res's
                     * argument order must not be relied on after this call. */
                    Expr* up = apply_up_values(res);
                    if (up) {
                        expr_free(res);
                        *changed = true; /* UpValue rule fired */
                        return up;
                    }
                }
#endif
                Expr* down = apply_down_values_def(hdef, res);
                if (down) {
                    expr_free(res);
                    *changed = true; /* DownValue rule fired */
                    return down;
                }

                /* 4b. NDArray element-wise fast path. Listable functions with a
                 * registered machine kernel (elementary/special functions) map
                 * directly over an NDArray argument's flat buffer at C speed
                 * instead of falling through to the slow List-threading path
                 * (which NDArrays don't even trigger). A kernel failure on any
                 * element (pole/overflow) degrades faithfully to the List path,
                 * so results always match f[{...}]. Runs after user DownValues
                 * (so overrides still win) and before the builtin (whose scalar
                 * numeric paths ignore NDArrays anyway). */
                if (hdef) {
                    size_t na = res->data.function.arg_count;
                    Expr** aa = res->data.function.args;
                    Expr* nd = NULL;
                    if (na == 1 && hdef->ndarray_unary_kernel && is_ndarray(aa[0])) {
                        nd = ndarray_map_unary(aa[0], hdef->ndarray_unary_kernel);
                        if (!nd) nd = ndarray_delist_and_reeval(res);
                    } else if (na == 2 && hdef->ndarray_binary_kernel &&
                               (is_ndarray(aa[0]) || is_ndarray(aa[1]))) {
                        /* BOTH operands may be arrays. The test used to be an
                         * XOR -- one array, one broadcast scalar -- so a kernel
                         * registered for a genuinely two-argument function was
                         * unreachable whenever both arguments were arrays, which
                         * is the ordinary way to call it: ArcTan[v, w] is
                         * numpy's arctan2 of two vectors and Beta[p, q] takes
                         * two. Both stayed unevaluated on a visible NDArray and
                         * boxed every element on a packed List, with NDKB_ArcTan
                         * and NDKB_Beta sitting registered and unused. */
                        nd = (is_ndarray(aa[0]) && is_ndarray(aa[1]))
                                 ? ndarray_map_binary2(aa[0], aa[1], hdef->ndarray_binary_kernel)
                                 : ndarray_map_binary(aa[0], aa[1], hdef->ndarray_binary_kernel);
                        if (!nd) nd = ndarray_delist_and_reeval(res);
                    } else if (na >= 3 && hdef->ndarray_nary_kernel && (attrs & ATTR_LISTABLE)) {
                        /* An n-ary Listable head (LerchPhi, Hypergeometric2F1,
                         * PolyLog[n, p, z], ...) over one or more array operands
                         * threads exactly as the List path does, in one pass
                         * over the buffers. */
                        bool any_nd = false;
                        for (size_t i = 0; i < na && !any_nd; i++) any_nd = is_ndarray(aa[i]);
                        if (any_nd) {
                            nd = ndarray_map_nary(aa, na, (const NDNaryKernel*)hdef->ndarray_nary_kernel);
                            if (!nd) nd = ndarray_delist_and_reeval(res);
                        }
                    }
                    if (nd) {
                        expr_free(res);
                        *changed = true; /* NDArray fast path rewrote the call */
                        return nd;
                    }
                }

                /* 5. Call C-level Built-in Functions (internal "down code") */
                if (hdef && hdef->builtin_func) {
                    /* Trace: a builtin's internal evaluate() calls (e.g. Range
                     * folding i+1) are the rule's own computation, not user
                     * argument sub-evaluations -- hide them so the builtin shows
                     * as a single rewrite. Arguments were already evaluated (and
                     * traced) above, before this bump. */
                    g_trace_suppress++;
                    unsigned long msg_before = mth_msg_fired_count();
                    int opaque = !hdef->eval_transparent;
                    g_builtin_depth += opaque;
                    const char* outer_head = g_dispatch_head;
                    g_dispatch_head = (res->type == EXPR_FUNCTION &&
                                       res->data.function.head->type == EXPR_SYMBOL)
                                    ? res->data.function.head->data.symbol.name : NULL;
                    Expr* ret = hdef->builtin_func(res);
                    /* 5a. Mixed exact/inexact numeric arguments. A numeric
                     * builtin that declined f[2., Sqrt[2]] is re-called once
                     * with its exact numeric arguments numericalized to the
                     * inexact ones' precision (numeric_mixed_args_call), so
                     * Zeta[2., Sqrt[2]] is a number as Sin[2.] is. Flat heads
                     * (Plus, Times) run their own contagion pass, and a decline
                     * that already emitted a message is a rejection of the
                     * call, not of its exactness, so neither is retried. */
                    if (!ret && (attrs & ATTR_NUMERICFUNCTION) && !(attrs & ATTR_FLAT) &&
                        res->type == EXPR_FUNCTION && res->data.function.arg_count >= 2 &&
                        mth_msg_fired_count() == msg_before) {
                        Expr* mixed = numeric_mixed_args_call(res, attrs);
                        if (mixed) {
                            ret = hdef->builtin_func(mixed);
                            expr_free(mixed);
                        }
                    }
                    g_builtin_depth -= opaque;
                    g_dispatch_head = outer_head;
                    g_trace_suppress--;
                    if (ret) {
                        /* No inf/nan may leak: a machine fast path that
                         * overflowed a double past even MPFR's exponent range
                         * would otherwise hand back `inf.0`/`nan`, which are
                         * libm artifacts, not Mathilda values. Rewrite to
                         * Overflow[]/-Overflow[]/Indeterminate. Not gated on
                         * NumericFunction -- BernoulliB / EulerE and friends are
                         * only Listable, yet must never leak either -- the cheap
                         * type check IS the gate, so the finite path is
                         * untouched. (Representable overflows are promoted to an
                         * extended-exponent MPFR real inside each head, upstream
                         * of here.) */
                        if (numeric_result_has_nonfinite(ret)) {
                            ret = numeric_sanitize_nonfinite(ret);
                        }
                        expr_free(res);
                        *changed = true; /* Built-in produced a rewrite */
                        return ret;
                    }
                }

                /* 5b. Interval[...] threading. A NumericFunction applied to an
                 * Interval argument that no builtin handled (Erfi, ExpIntegralEi,
                 * PolyLog, the piecewise heads, ...) threads by certified
                 * monotonicity — interval_thread_call() tries the tight bespoke
                 * handlers, then interval-evaluates the symbolic derivative to
                 * prove the sign (src/interval.c). Runs after DownValues and the
                 * builtin, so user overrides and the tight elementary handlers
                 * already win. The common no-interval path costs only a bit test
                 * plus a cheap argument scan on this already-cold branch. */
                if (hdef && (attrs & ATTR_NUMERICFUNCTION) &&
                    res->type == EXPR_FUNCTION) {
                    for (size_t i = 0; i < res->data.function.arg_count; i++) {
                        if (is_interval(res->data.function.args[i])) {
                            Expr* ivr = interval_thread_call(res);
                            if (ivr) {
                                expr_free(res);
                                *changed = true; /* Interval threading rewrote */
                                return ivr;
                            }
                            break;
                        }
                    }
                }

                /* 6. Special primitives (Set, SetDelayed) */
                if ((head_name == SYM_Set || head_name == SYM_SetDelayed) && res->data.function.arg_count == 2) {
                    Expr* lhs = res->data.function.args[0];
                    Expr* rhs = res->data.function.args[1];
                    int is_delayed = (head_name == SYM_SetDelayed);

                    /* Evaluate the LHS arguments to canonical form (so f[x]=1 with
                     * x=c defines f[c]=1), holding patterns and Hold positions. */
                    bool free_target = false;
                    Expr* target_lhs = build_assignment_target(lhs, &free_target);

                    if (apply_assignment(target_lhs, rhs, is_delayed)) {
                        /* The value of x = rhs is rhs, already evaluated as this
                         * step's argument. It needs another evaluation only if
                         * the assignment can have changed it -- if it mentions a
                         * symbol of the left side (y = y + 1 re-evaluates, and
                         * recurses, as it should). Otherwise it is final: the
                         * re-evaluation repeated every message of a declining
                         * call (z = Det[{{1, 2}}] printed Det::matsq twice). */
                        bool final_value = !is_delayed && settle_ok && settle_argc == 2 &&
                                           settle_snap[1] == rhs &&
                                           !eval_mentions_symbols_of(rhs, target_lhs);
                        Expr* ret = is_delayed ? expr_new_symbol(SYM_Null)
                                  : final_value ? expr_copy(rhs) : evaluate(rhs);
                        if (free_target) expr_free(target_lhs);
                        expr_free(res);
                        *changed = true; /* Set/SetDelayed installed a rule */
                        if (final_value) g_step_settled_depth = eval_recursion_depth;
                        return ret;
                    }
                    if (free_target) expr_free(target_lhs);
                }
#if UP_VALUES
                /* 6b. UpSet / UpSetDelayed (lhs ^= / ^:= rhs). Same shape as Set
                 * (HoldFirst / HoldAll); installs on the operands' up_values. */
                else if ((head_name == SYM_UpSet || head_name == SYM_UpSetDelayed) &&
                         res->data.function.arg_count == 2) {
                    Expr* lhs = res->data.function.args[0];
                    Expr* rhs = res->data.function.args[1];
                    int is_delayed = (head_name == SYM_UpSetDelayed);
                    bool free_target = false;
                    Expr* target_lhs = build_assignment_target(lhs, &free_target);
                    if (apply_up_assignment(target_lhs, rhs, is_delayed)) {
                        Expr* ret = is_delayed ? expr_new_symbol(SYM_Null) : evaluate(rhs);
                        if (free_target) expr_free(target_lhs);
                        expr_free(res);
                        *changed = true;
                        return ret;
                    }
                    if (free_target) expr_free(target_lhs);
                }
                /* 6c. TagSet / TagSetDelayed (f /: lhs = / := rhs). Both HoldAll,
                 * so for the immediate case evaluate rhs here before storing. */
                else if ((head_name == SYM_TagSet || head_name == SYM_TagSetDelayed) &&
                         res->data.function.arg_count == 3) {
                    Expr* tag = res->data.function.args[0];
                    Expr* lhs = res->data.function.args[1];
                    Expr* rhs = res->data.function.args[2];
                    int is_delayed = (head_name == SYM_TagSetDelayed);
                    bool free_target = false;
                    Expr* target_lhs = build_assignment_target(lhs, &free_target);
                    Expr* store_rhs = is_delayed ? rhs : evaluate(rhs); /* owned iff immediate */
                    bool installed = apply_tag_assignment(tag, target_lhs, store_rhs, is_delayed);
                    Expr* ret = !is_delayed ? store_rhs
                              : expr_new_symbol(installed ? SYM_Null : SYM_DollarFailed);
                    if (free_target) expr_free(target_lhs);
                    expr_free(res);
                    *changed = true;
                    return ret;
                }
                /* 6d. TagUnset (f /: lhs =.). HoldAll. */
                else if (head_name == SYM_TagUnset && res->data.function.arg_count == 2) {
                    Expr* tag = res->data.function.args[0];
                    Expr* lhs = res->data.function.args[1];
                    bool free_target = false;
                    Expr* target_lhs = build_assignment_target(lhs, &free_target);
                    apply_tag_unset(tag, target_lhs);
                    if (free_target) expr_free(target_lhs);
                    expr_free(res);
                    *changed = true;
                    return expr_new_symbol(SYM_Null);
                }
#endif

                /* OneIdentity is intentionally NOT rewritten at evaluation
                 * time. In Mathematica it is purely a pattern-matching
                 * attribute: it lets f[x_, y_:def] match a literal `a`.
                 * The 1-arg collapse f[x] -> x is the responsibility of
                 * each head's builtin (Plus, Times, Power, GCD, LCM, And,
                 * Or, Dot all handle the n==1 case explicitly), so a
                 * user-defined OneIdentity head like
                 *   SetAttributes[g, OneIdentity]
                 * leaves g[x] as g[x] rather than rewriting to x.
                 * The pattern-matching half lives in src/match.c (search
                 * for ATTR_ONEIDENTITY). */
            } else if (head->type == EXPR_FUNCTION && head->data.function.head->type == EXPR_SYMBOL &&
                       head->data.function.head->data.symbol.name == SYM_Function) {

                /* 7. Apply Pure Function */
                Expr* applied = apply_pure_function(head, res->data.function.args, res->data.function.arg_count);
                if (applied) {
                    expr_free(res);
                    *changed = true; /* Pure Function applied */
                    return applied;
                }
            } else if (head->type == EXPR_FUNCTION && head->data.function.head->type == EXPR_SYMBOL &&
                       res->data.function.arg_count >= 1 &&
                       ml_model_apply_probe(head)) {
                /* 7a-pre. A fitted machine-learning model as a callable:
                 * PredictorFunction[...][features] predicts, and
                 * PredictorFunction[...]["Coefficients"] reads a property. This sits in
                 * the same composite-head chain as Function[...][args] just above and
                 * Association[...][key] just below, which is the point -- a trained
                 * model is a callable object, and this codebase already has an idiom
                 * for that, so no new evaluation concept is introduced. See
                 * src/ml/predict.h for why the model is a plain EXPR_FUNCTION rather
                 * than a new node type. */
                Expr* applied = ml_model_apply(head, res->data.function.args,
                                               res->data.function.arg_count);
                if (applied) {
                    expr_free(res);
                    *changed = true;
                    return applied;
                }
            } else if (head->type == EXPR_FUNCTION && head->data.function.head->type == EXPR_SYMBOL &&
                       head->data.function.head->data.symbol.name == SYM_Association &&
                       res->data.function.arg_count >= 1) {
                /* 7a. Association as accessor: <|...|>[key] looks the key up,
                 * giving the value or Missing["KeyAbsent", key] -- the idiomatic
                 * Wolfram accessor, complementing Lookup and Part. The argument
                 * is always a LITERAL key: unlike Part and Lookup, the accessor
                 * does not unwrap Key[...], so <|"a" -> 1|>[Key["a"]] is
                 * Missing["KeyAbsent", Key["a"]] (Mathematica 15).
                 * Multi-key <|...|>[k1, k2, ...] is nested lookup: the value for
                 * k1 is then applied to the remaining keys. */
                Expr* lookup_key = res->data.function.args[0];
                Expr* found = assoc_lookup_value(head, lookup_key);  /* O(1) via key index */
                Expr* out;
                if (found) {
                    out = expr_copy(found);
                } else {
                    Expr* margs[2] = { expr_new_string("KeyAbsent"), expr_copy(lookup_key) };
                    out = expr_new_function(expr_new_symbol(SYM_Missing), margs, 2);
                }
                size_t nkeys = res->data.function.arg_count;
                if (nkeys > 1 && found) {
                    /* Apply the retrieved value to the remaining keys and let the
                     * evaluator continue (nested associations recurse here). */
                    size_t rest = nkeys - 1;
                    Expr** rest_args = malloc(sizeof(Expr*) * rest);
                    for (size_t i = 0; i < rest; i++)
                        rest_args[i] = expr_copy(res->data.function.args[i + 1]);
                    Expr* nested = expr_new_function(out, rest_args, rest);
                    free(rest_args);
                    out = nested;
                }
                expr_free(res);
                *changed = true;
                return out;
            } else if (head->type == EXPR_FUNCTION && head->data.function.head->type == EXPR_SYMBOL &&
                       head->data.function.head->data.symbol.name == SYM_Key &&
                       head->data.function.arg_count == 1 &&
                       res->data.function.arg_count == 1 &&
                       res->data.function.args[0]->type == EXPR_FUNCTION &&
                       res->data.function.args[0]->data.function.head->type == EXPR_SYMBOL &&
                       res->data.function.args[0]->data.function.head->data.symbol.name == SYM_Association) {
                /* 7a'. Key[k][assoc] operator form: extract the value at key k,
                 * giving the value or Missing["KeyAbsent", k]. This is the
                 * curried complement of assoc[Key[k]], and it is what makes
                 * record pipelines like GroupBy[records, Key["field"]] and
                 * SortBy[records, Key["field"]] work. */
                Expr* key   = head->data.function.args[0];
                Expr* assoc = res->data.function.args[0];
                Expr* found = assoc_lookup_value(assoc, key);   /* O(1) via key index */
                Expr* out;
                if (found) {
                    out = expr_copy(found);
                } else {
                    Expr* margs[2] = { expr_new_string("KeyAbsent"), expr_copy(key) };
                    out = expr_new_function(expr_new_symbol(SYM_Missing), margs, 2);
                }
                expr_free(res);
                *changed = true;
                return out;
            } else if (head->type == EXPR_FUNCTION && head->data.function.head->type == EXPR_SYMBOL &&
                       head->data.function.head->data.symbol.name == SYM_Derivative &&
                       res->data.function.arg_count == 1 &&
                       res->data.function.args[0]->type == EXPR_FUNCTION &&
                       res->data.function.args[0]->data.function.head->type == EXPR_SYMBOL &&
                       res->data.function.args[0]->data.function.head->data.symbol.name == SYM_Function) {
                /* 7b. Derivative[n1,...,nm][Function[...]] reduces to a new
                 * Function whose body has been differentiated. Without this
                 * step, f'[x] for a pure-function f would remain stuck as
                 * Derivative[1][Function[...]][x]. */
                Expr* reduced = derivative_of_pure_function(head, res->data.function.args[0]);
                if (reduced) {
                    expr_free(res);
                    *changed = true; /* Derivative-of-Function reduced */
                    return reduced;
                }
            } else if (head->type == EXPR_FUNCTION && head->data.function.head->type == EXPR_SYMBOL &&
                       head->data.function.head->data.symbol.name == SYM_Derivative &&
                       res->data.function.arg_count == 1 &&
                       res->data.function.args[0]->type == EXPR_FUNCTION &&
                       res->data.function.args[0]->data.function.head->type == EXPR_SYMBOL &&
                       res->data.function.args[0]->data.function.head->data.symbol.name == SYM_InterpolatingFunction) {
                /* 7b''. Derivative[d1,...,dm][InterpolatingFunction[...]] reduces
                 * to a fresh InterpolatingFunction carrying the accumulated
                 * derivative orders, which evaluates the mixed partial when
                 * applied. This makes ifun'[x] and D[ifun[x],x] work. */
                Expr* reduced = interp_make_derivative(head, res->data.function.args[0]);
                if (reduced) {
                    expr_free(res);
                    *changed = true; /* Derivative-of-InterpolatingFunction reduced */
                    return reduced;
                }
            } else if (head->type == EXPR_FUNCTION && head->data.function.head->type == EXPR_SYMBOL &&
                       head->data.function.head->data.symbol.name == SYM_Derivative &&
                       res->data.function.arg_count == 1 &&
                       res->data.function.args[0]->type == EXPR_SYMBOL) {
                /* 7b'. Derivative[n1,...,nm][f] where f is a symbol with
                 * DownValues. Reduce by synthesising
                 *     Function[{t1,...,tm}, f[t1,...,tm]]
                 * after the DownValue rewrite, then differentiate via the
                 * pure-function pipeline. This is what makes f'[x] work
                 * for user-defined f[x_] := body. */
                Expr* reduced = derivative_of_symbol(head, res->data.function.args[0]);
                if (reduced) {
                    expr_free(res);
                    *changed = true; /* Derivative-of-symbol reduced */
                    return reduced;
                }
            } else if (head->type == EXPR_FUNCTION && head->data.function.head->type == EXPR_SYMBOL &&
                       head->data.function.head->data.symbol.name == SYM_Composition &&
                       head->data.function.arg_count >= 1) {
                /* 7c. Composition[f1, ..., fn][args...] -> f1[f2[...[fn[args...]]]].
                 * The innermost call carries all the user-supplied arguments;
                 * each outer fk wraps the previous result as a single argument. */
                size_t nf = head->data.function.arg_count;
                size_t na = res->data.function.arg_count;
                Expr** call_args = malloc(sizeof(Expr*) * na);
                for (size_t i = 0; i < na; i++) {
                    call_args[i] = expr_copy(res->data.function.args[i]);
                }
                Expr* inner = expr_new_function(
                    expr_copy(head->data.function.args[nf - 1]),
                    call_args, na);
                free(call_args);
                for (size_t k = nf - 1; k > 0; k--) {
                    Expr* one[1] = { inner };
                    inner = expr_new_function(
                        expr_copy(head->data.function.args[k - 1]),
                        one, 1);
                }
                expr_free(res);
                *changed = true; /* Composition unrolled */
                return inner;
            } else if (head->type == EXPR_FUNCTION && head->data.function.head->type == EXPR_SYMBOL &&
                       head->data.function.head->data.symbol.name == SYM_InterpolatingFunction) {
                /* 7d. InterpolatingFunction[domain, table][x] -> interpolated
                 * value. The object itself is a normal form; only its
                 * application is reduced here. interp_apply returns NULL for a
                 * symbolic / out-of-form argument, leaving the call intact. */
                Expr* applied = interp_apply(head, res->data.function.args, res->data.function.arg_count);
                if (applied) {
                    expr_free(res);
                    *changed = true; /* InterpolatingFunction evaluated */
                    return applied;
                }
            } else if (head->type == EXPR_FUNCTION && res->data.function.arg_count == 1 &&
                       opform_matches(head)) {
                /* 7f. Operator (curried) form of a registered builtin:
                 * Select[crit][x] -> Select[x, crit], Map[f][x] -> Map[f, x],
                 * Insert[e, n][x] -> Insert[x, e, n], ... (table in opform.c).
                 * NULL (head declined the rewritten call) keeps h[o][x]. */
                Expr* applied = opform_apply(head, res->data.function.args[0]);
                if (applied) {
                    expr_free(res);
                    *changed = true;
                    return applied;
                }
            } else if (head->type == EXPR_SPARSE) {
                /* 7g. SparseArray[...]["prop"] -> property (sparse_builtins.c);
                 * NULL (after SparseArray::nomthd for an unknown name) leaves
                 * it unevaluated. */
                Expr* applied = sp_apply(head, res->data.function.args,
                                         res->data.function.arg_count);
                if (applied) {
                    expr_free(res);
                    *changed = true;
                    return applied;
                }
            } else if (head->type == EXPR_COMPILED) {
                /* 7e. CompiledFunction[...][args] -> numeric result (or the
                 * interpreter fallback for symbolic args).  Returns NULL only on
                 * an arity mismatch, leaving the application unevaluated. */
                Expr* applied = compiled_function_apply(head->data.compiled,
                                                        res->data.function.args,
                                                        res->data.function.arg_count);
                if (applied) {
                    expr_free(res);
                    *changed = true; /* CompiledFunction applied */
                    return applied;
                }
            }

            /* A head with a native sparse path declined this call: the
             * transparency gate is its fallback (see the pre-gate above). */
            if (hdef && (hdef->sparse_mode & SPG_NATIVE) && sp_any_created() &&
                eval_sparse_gate_applies(res, hdef->sparse_mode)) {
                bool rest = false;
                g_trace_suppress++;
                Expr* g = sp_gate(res, hdef->sparse_mode, &rest);
                g_trace_suppress--;
                if (g) { expr_free(res); *changed = true; return g; }
            }

            /* THE POST-GATE. Nothing above rewrote this node, so it is coming to
             * REST -- and a packed argument that is still here was kept for a
             * fast path that did not fire. Whatever the head is, it has now
             * declined the buffer, so leaving one inside an inert expression can
             * only differ from the plain list.
             *
             * Measured: Mod[buffer] came to rest as Mod[{1., 2., 3.}] where the
             * plain Mod[{1., 2., 3.}] threads to {Mod[1.], Mod[2.], Mod[3.]} --
             * an aware Listable head skips threading so its kernel can fire, and
             * no kernel matches that arity. Materialising here and reporting the
             * change lets the next pass thread it exactly as the list does.
             *
             * Converges in one extra pass: the retry finds no packed argument
             * left, so nothing changes and the node rests for real. */
            bool post_gated = false;
            if (pack_any_created()) {
                for (size_t i = 0; i < res->data.function.arg_count; i++) {
                    Expr* pa = res->data.function.args[i];
                    if (!is_packed_list(pa)) continue;
                    pack_gate_report(res->data.function.head, pa);
                    res->data.function.args[i] = ndarray_to_nested_list(pa);
                    expr_free(pa);
                    *changed = true;
                    post_gated = true;
                    g_pack_gate_ticks++;
                }
            }

            /* Nothing rewrote the node: it comes to rest here, and a further
             * step would only replay this one (see g_step_settled_depth). A
             * non-Listable head that declined a buffer WITH a message
             * (Dot::dotsh on a packed vector) rejected the call: its
             * materialised form rests too, rather than being re-asked and
             * printing the message again. (A Listable head's materialised
             * call still threads on the next step.) */
            if ((!post_gated && settle_ok && eval_args_untouched(res, settle_snap, settle_argc)) ||
                (post_gated && !(attrs & ATTR_LISTABLE) &&
                 mth_msg_fired_count() != step_msgs_before))
                g_step_settled_depth = eval_recursion_depth;
            return res;
        }
    }
    return expr_copy(e);
}

/*
 * evaluate:
 * The primary evaluator loop. Repeatedly applies evaluate_step until the 
 * expression reaches a fixed point (no further changes) or the iteration 
 * limit is reached.
 */
/* Convert an uncaught in-flight Throw that reached the top level into its WL
 * result: Throw[v,t,f] -> f[v,t] (evaluated); Throw[v] / Throw[v,t] ->
 * Hold[Throw[...]] (inert, so feeding it back does not re-throw). Emits
 * Throw::nocatch on stderr -- the channel eval already uses for
 * $RecursionLimit/$IterationLimit (there is no Message[] builtin). Takes
 * ownership of `thr`. */
static Expr* eval_report_uncaught_throw(Expr* thr) {
    char* s = expr_to_string(thr);
    mth_message("Throw", "nocatch", "Uncaught %s returned to top level.", s ? s : "Throw[...]");
    free(s);
    if (thr->data.function.arg_count == 3) {
        Expr* fa[2] = { expr_copy(thr->data.function.args[0]),
                        expr_copy(thr->data.function.args[1]) };
        Expr* call = expr_new_function(expr_copy(thr->data.function.args[2]), fa, 2);
        expr_free(thr);
        Expr* out = evaluate(call);
        expr_free(call);
        return out;
    }
    Expr* one[1] = { thr };
    return expr_new_function(expr_new_symbol(SYM_Hold), one, 1);
}

/* An abort that reaches the top level ends the input: the result is $Aborted,
 * and the front end learns of it through eval_take_top_abort (a -file script
 * stops there).  Takes ownership of the sentinel. */
static Expr* eval_report_uncaught_abort(Expr* ab) {
    expr_free(ab);
    g_top_abort = true;
    return expr_new_symbol(SYM_DollarAborted);
}

/* An in-flight Goto that survives to the top level found no matching Label in
 * any enclosing CompoundExpression. Emit Goto::nolabel on stderr (same channel
 * as Throw::nocatch) and return the Goto[tag] node unchanged -- unlike Throw,
 * the inert node itself is the WL result, so no rewrite is needed. Takes
 * ownership of and returns `g`. */
static Expr* eval_report_uncaught_goto(Expr* g) {
    char* s = expr_to_string(g);
    mth_message("Goto", "nolabel", "%s found no matching Label.", s ? s : "Goto[...]");
    free(s);
    return g;
}

/* A Break[]/Continue[] that survives to the top level had no enclosing Do/For/
 * While to consume it. Emit <head>::nofwd (same stderr channel as the other
 * flow-control reporters) and rewrite to Hold[Break[]] / Hold[Continue[]] so
 * feeding the result back does not re-trigger the marker. Takes ownership of
 * and consumes `e`, returning the Hold[] wrapper. */
static Expr* eval_report_uncaught_break_continue(Expr* e) {
    const char* h = e->data.function.head->data.symbol.name;  /* Break | Continue */
    char* s = expr_to_string(e);
    mth_message(h, "nofwd",
            "No enclosing For, While, Until or Do found for %s.",
            s ? s : (h == SYM_Break ? "Break[]" : "Continue[]"));
    free(s);
    Expr* one[1] = { e };
    return expr_new_function(expr_new_symbol(SYM_Hold), one, 1);
}

/* See eval.h. */
Expr* eval_report_uncaught(Expr* e) {
    if (eval_is_inflight_throw(e))          return eval_report_uncaught_throw(e);
    if (eval_is_inflight_abort(e))          return eval_report_uncaught_abort(e);
    if (eval_is_inflight_goto(e))           return eval_report_uncaught_goto(e);
    if (eval_is_inflight_break_continue(e)) return eval_report_uncaught_break_continue(e);
    return e;
}

static Expr* evaluate_core(Expr* e);

Expr* evaluate(Expr* e) {
    Expr* r = evaluate_core(e);
    /* Module locals that escaped and have since been dropped are reclaimed
     * when the outermost evaluation ends, or here in a nested evaluation once
     * enough have queued up (symtab.c, "Timing"; a Module reclaims its own
     * locals itself). */
    if (eval_recursion_depth == 0) symtab_temp_sweep_if_queued();
    else if (g_temp_pending_n >= g_temp_sweep_at) symtab_temp_sweep();
    return r;
}

void eval_park_unwind(const Expr* r) {
    if (eval_recursion_depth > 0 && r && eval_is_inflight_unwind(r) && !g_arith_warnings_muted) {
        if (g_pending_unwind) expr_free(g_pending_unwind);
        g_pending_unwind = expr_copy((Expr*)r);
    }
}

Expr* eval_callback(Expr* call) {
    Expr* r = evaluate(call);
    eval_park_unwind(r);
    return r;
}

static Expr* evaluate_core(Expr* e) {
    if (!e) return NULL;

    /* M3 phase-3 timestamp early-exit. If this expression has been fully
     * evaluated under the current symbol-table state (clock unchanged
     * since), return an inc-ref'd view immediately and skip the entire
     * fixed-point loop and all of evaluate_step. This lifts the cost of
     * a re-evaluation from O(tree size) to O(1). Atoms and bare symbols
     * already short-circuit cheaply inside evaluate_step (atom returns
     * expr_copy, symbol-no-OwnValue returns expr_copy), so we limit the
     * pre-check to FUNCTION nodes -- both to avoid an extra branch on
     * the common atom path and because atoms are never expensive to
     * "re-evaluate" anyway. */
    if (e->type == EXPR_FUNCTION && eval_fixed_point_reusable(e)) {
        return expr_copy(e);
    }

    bool is_top_level = (eval_recursion_depth == 0);
    if (is_top_level) { eval_overflow = false; g_toplevel_eval_id++; }

    /* Guard the C stack: when nested evaluate() calls would exceed the
     * recursion limit, wrap the input in Hold[] so it stops re-entering
     * the evaluator, set the sticky overflow flag so all enclosing
     * fixed-point loops bail, and emit a message exactly once per
     * top-level evaluation. */
    if (eval_recursion_depth >= eval_recursion_limit) {
        if (!eval_overflow) {
            mth_message("$RecursionLimit", "reclim", "Recursion depth of %d exceeded.", eval_recursion_limit);
        }
        eval_overflow = true;
        Expr** wrap = malloc(sizeof(Expr*));
        wrap[0] = expr_copy(e);
        Expr* held = expr_new_function(expr_new_symbol(SYM_Hold), wrap, 1);
        free(wrap);
        return held;
    }

    eval_recursion_depth++;

    /* Trace: this evaluate() call gets its own frame; nested sub-evaluations
     * push child frames and splice back on exit (see trace_frame_pop). While
     * suppressed (inside a builtin / Listable threading) nothing is pushed, so
     * that whole subtree stays out of the trace. */
    bool trace_here = TRACE_ACTIVE();
    if (trace_here) trace_frame_push();

    Expr* current = expr_copy(e);
    Expr* next = NULL;
    int iterations = 0;

    while (iterations < MAX_ITERATIONS) {
        /* TimeConstrained's cooperative wall-clock backstop.  No-op
         * unless we're inside an active TimeConstrained call; on
         * hosts where ITIMER_PROF/SIGPROF are reliable the signal
         * normally fires first and this check stays a cheap inactive
         * read.  On hosts where they aren't (WSL 1), this is what
         * actually enforces the deadline.  When tripped, the call
         * siglongjmp's out of this loop straight to TimeConstrained's
         * sigsetjmp; no further cleanup runs here, exactly matching
         * the signal-handler path.
         *
         * Throttled: clock_gettime on every loop iteration was ~10% of
         * DSolve's field-arithmetic profile.  The SIGPROF async timer is the
         * PRIMARY deadline mechanism (and inside async-deferred regions, where
         * SIGPROF is muted, allocation safepoints via tc_alloc_safepoint are),
         * so this cooperative wall-clock read only needs to be sampled
         * periodically.  Once every 128 iterations is far finer than any
         * user-visible deadline yet reads the clock 128x less often; the very
         * first iteration (tick 0) still checks so a pre-elapsed deadline trips
         * at once.  A single global tick (evaluate is single-threaded) gives a
         * uniform sampling rate across all nested evaluate() frames. */
        {
            static unsigned tc_check_tick = 0;
            if ((tc_check_tick++ & 127u) == 0) tc_check_deadline();
        }

        /* In-loop timestamp fixed-point exit.  The entry short-circuit at the top
         * of evaluate() only catches an ALREADY-stamped INPUT; a stamped FUNCTION
         * can also become `current` MID-loop — most importantly a canonical
         * value (e.g. an Association) reached through an OwnValue symbol
         * (`a = <|...|>; ... a ...`).  There evaluate(a) rewrites the symbol to
         * its stored value and then, without this check, evaluate_step would
         * rebuild that value O(tree size) every time — re-canonicalising the
         * association and discarding its cached key index — even though it is
         * already at a fixed point.  If `current` is a FUNCTION fully evaluated
         * under the current clock, stop here.  Same invariant as the entry check
         * (a stamp is set only on a clean fixed-point exit and is invalidated by
         * any symbol-table mutation via the clock), so this is a pure speedup. */
        if (current->type == EXPR_FUNCTION && eval_fixed_point_reusable(current)) {
            eval_recursion_depth--;
            if (trace_here) { trace_clear_pending(); trace_frame_pop(); }
            if (is_top_level) current = eval_report_uncaught(current);
            return current;
        }

        bool step_changed = false;
        uint64_t gate_mark = g_pack_gate_ticks;
        g_step_settled_depth = -1;
        next = evaluate_step(current, &step_changed);
        bool gate_fired = (g_pack_gate_ticks != gate_mark);
        /* The step changed the node (its head or arguments evaluated, Flat /
         * Orderless reshaped it) and then came to REST: no builtin, rule or
         * gate rewrote it. Another step would re-evaluate the same stamped
         * arguments and re-ask the same declining builtin -- re-emitting its
         * message (Det[m] printed Det::matsq twice) for nothing. So the step's
         * result is the fixed point. Only the rest return of THIS frame's step
         * sets this depth; a nested evaluate runs deeper, and a step that
         * returns by any other path leaves -1 or a deeper depth here. */
        bool step_settled = (g_step_settled_depth == eval_recursion_depth);

        /* M3 phase-4 (§3.4): eager early exit. evaluate_step signals via
         * the `step_changed` out-parameter whether any rewrite fired
         * during the step (head re-evaluation, arg evaluation, Sequence
         * flatten, Unevaluated strip, Flat, Listable, Orderless,
         * DownValue, built-in, special primitive, pure Function,
         * Derivative-of-Function, Composition-unfold). When nothing
         * fired, the result is structurally identical to `current` and
         * we are at a fixed point — skip the O(tree) expr_eq compare.
         *
         * Some built-ins (Plus, Times, ...) unconditionally rebuild
         * their output even when no terms combined; those trip the
         * change flag without producing a structural difference. Use
         * expr_eq as a fallback in the changed-true branch so those
         * "false positives" still converge in one iteration, matching
         * the old semantics. Cost: identical to the pre-§3.4 path on
         * the slow case; the win is the cheap boolean fast-path on the
         * common case where nothing fires (atoms, bare symbols, fully
         * reduced functions). */
        bool is_fixed_point;
        if (step_changed && !trace_here && eval_is_terminal_value(next)) {
            /* The step rewrote `current` into a value that is a fixed point by
             * construction -- a numeric/string atom, or a canonical machine
             * Complex (see eval_is_terminal_value). The loop would otherwise run
             * one more evaluate_step on it purely to discover "nothing changed",
             * and pay an expr_eq to compare a FUNCTION against an atom. That
             * confirming pass is the dominant per-call overhead of a numeric
             * builtin -- F[x_machine] -> Real -- so take `next` as the fixed
             * point directly. Same result, same stamp, one step fewer. Not
             * taken while tracing: the rewrite branch below records the step. */
            expr_free(current);
            current = next;
            next = NULL;
            is_fixed_point = true;
        } else if (step_changed && step_settled && !trace_here) {
            /* (A pack-gate materialisation in this step is kept too: `next`
             * is the newer representation, as the gate_fired branch below
             * would also choose.) */
            expr_free(current);
            current = next;
            next = NULL;
            is_fixed_point = true;
        } else {
            is_fixed_point = !step_changed || expr_eq(current, next);
        }
        if (is_fixed_point && gate_fired && next) {
            /* The packed-list gate materialised a buffer somewhere in this step,
             * and expr_eq is DELIBERATELY blind to that -- a packed list and its
             * plain form are the same value, which is the whole point. So a step
             * whose only effect was materialising looks like no progress, and
             * keeping `current` would throw the work away and hand a packed node
             * to a head that has just declared it cannot read one. Measured
             * before this: Table[i j, {i,300}, {j,300}] came back as a List of
             * 300 packed rows (Dimensions {300}, and Total[m, 2] a list instead
             * of a number) because the outer List's materialisation was
             * discarded exactly here.
             *
             * `current` and `next` are equal, so taking `next` is value-neutral;
             * it just keeps the newer representation. The tick counter only ever
             * over-reports (a nested evaluate's gate bumps it too), and an
             * over-report costs one pointer choice between two equal values. */
            expr_free(current);
            current = next;
            next = NULL;
        }
        if (is_fixed_point) {
            expr_free(next);
            /* M3 phase-3: stamp the fully-evaluated result with the
             * current clock so a subsequent evaluate(current) hits the
             * early-exit above. We deliberately stamp ONLY on a clean
             * fixed-point exit; the iteration-cap and recursion-overflow
             * paths below leave the timestamp untouched, so a later
             * evaluator gets a fresh chance to make progress. The
             * write is benign metadata, so it is safe even when
             * `current` is shared (refcount > 1). */
            if (!eval_overflow) {
                /* Stamp with the live clock, plus the GROUND flag when this is a
                 * whitelisted-constructor node over ground args -- so a later
                 * evaluate() can re-validate it as a fixed point even after the
                 * eval clock has churned (loop-invariant O(1) re-check). */
                uint64_t flag = node_compute_ground(current) ? EVAL_GROUND_BIT
                              : node_compute_inert(current)  ? EVAL_INERT_BIT : 0;
                current->last_evaluated_at = g_eval_clock | flag;
            }
            eval_recursion_depth--;
            /* Trace: the last step didn't rewrite; drop its reassembled
             * snapshot and finalize this frame (splice into parent or discard). */
            if (trace_here) { trace_clear_pending(); trace_frame_pop(); }
            /* An in-flight Throw that survives to the top level is uncaught
             * (any enclosing Catch would have consumed it at depth >= 1). */
            if (is_top_level) current = eval_report_uncaught(current);
            return current;
        }

        /* Trace: a real top-level rewrite (current -> next). Record the
         * "before" form -- the reassembled f[evaluated_args] snapshot when
         * evaluate_step produced one, else `current` (symbol/atom rewrites) --
         * then the result. frame_record's consecutive-dedup collapses chains
         * and avoids doubling a result equal to its reassembled form. */
        if (trace_here && g_trace_top) {
            Expr* before = g_trace_top->pending_reassembled
                               ? g_trace_top->pending_reassembled : current;
            frame_record(before);
            frame_record(next);
            trace_clear_pending();
        }

        /* Prepare for the next iteration */
        expr_free(current);
        current = next;
        iterations++;

        /* If a deeper call hit the recursion limit, the rewrites
         * above are no longer making real progress -- bail out so the
         * unwind doesn't burn $IterationLimit at every level. */
        if (eval_overflow) break;
    }

    if (iterations >= MAX_ITERATIONS) {
        fprintf(stderr, "$IterationLimit exceeded\n");
    }
    eval_recursion_depth--;
    /* Trace: iteration-cap / overflow exit -- finalize this frame too. */
    if (trace_here) { trace_clear_pending(); trace_frame_pop(); }
    if (is_top_level && eval_is_inflight_throw(current))
        current = eval_report_uncaught_throw(current);
    else if (is_top_level && eval_is_inflight_abort(current))
        current = eval_report_uncaught_abort(current);
    else if (is_top_level && eval_is_inflight_goto(current))
        current = eval_report_uncaught_goto(current);
    else if (is_top_level && eval_is_inflight_break_continue(current))
        current = eval_report_uncaught_break_continue(current);
    return current;
}

/* See eval.h. Evaluate `held_expr` to a fixed point while recording the nested
 * structure of its evaluation, returning a fresh List (with nested sublists for
 * each sub-evaluation that took a step). `held_expr` is BORROWED -- evaluate()
 * copies it and the caller retains ownership. Reentrant: the previous trace
 * state (frame stack, root result, tracing flag) is saved on the C stack and
 * restored on return, so a nested Trace runs on a fresh stack and appears to the
 * outer trace as a single already-reduced value. */
Expr* eval_collect_trace(Expr* held_expr) {
    if (!held_expr) return NULL;

    TraceFrame* saved_top      = g_trace_top;
    Expr*       saved_result   = g_trace_root_result;
    bool        saved_on       = g_tracing;
    int         saved_suppress = g_trace_suppress;
    g_trace_top         = NULL;    /* fresh stack for this Trace */
    g_trace_root_result = NULL;
    g_tracing           = true;
    g_trace_suppress    = 0;       /* Trace itself runs under the caller's
                                    * builtin suppression; clear it so the
                                    * traced expression is actually recorded. */

    /* Bump the eval clock so the timestamp early-exit (see top of evaluate())
     * cannot short-circuit an already-stamped root -- a full re-evaluation is
     * required to observe its steps. Trace is thereby non-memoized (it
     * invalidates the eval cache once per call). */
    eval_clock_bump();

    Expr* final = evaluate(held_expr);   /* borrowed in; fresh copy out */
    expr_free(final);

    Expr* list = g_trace_root_result;    /* outermost frame's collected List */
    if (!list)                           /* no step taken -> {} */
        list = expr_new_function(expr_new_symbol(SYM_List), NULL, 0);

    g_trace_top         = saved_top;     /* restore the outer trace state */
    g_trace_root_result = saved_result;
    g_tracing           = saved_on;
    g_trace_suppress    = saved_suppress;
    return list;
}
