#include "match.h"
#include "options.h"   /* optionvalue_inject_context: OptionValue in /; guards */
#include <gmp.h>
#include "part.h"
#include "modular.h"
#include "eval.h"
#include "print.h"
#include "symtab.h"
#include "attr.h"
#include "default_helper.h"
#include "sym_names.h"
#include "message.h"
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <limits.h>

/* An option argument: Rule[name,val] / RuleDelayed[name,val] with a
 * symbol/string name. Used by OptionsPattern matching. */
static bool mo_is_option_rule(const Expr* e) {
    if (!e || e->type != EXPR_FUNCTION) return false;
    const Expr* h = e->data.function.head;
    if (!h || h->type != EXPR_SYMBOL) return false;
    if (h->data.symbol.name != SYM_Rule && h->data.symbol.name != SYM_RuleDelayed) return false;
    if (e->data.function.arg_count != 2) return false;
    const Expr* lhs = e->data.function.args[0];
    return lhs && (lhs->type == EXPR_SYMBOL || lhs->type == EXPR_STRING);
}
static bool mo_is_list(const Expr* e) {
    return e && e->type == EXPR_FUNCTION && e->data.function.head->type == EXPR_SYMBOL
        && e->data.function.head->data.symbol.name == SYM_List;
}

/* Bounded free-list pool of MatchEnv objects. The DownValue dispatch path
 * (apply_down_values_def) allocates and frees one env per candidate rule, so a
 * recursive program like fib[n_]:=fib[n-1]+fib[n-2] churns millions of them.
 * env_new's three mallocs dominated that path; recycling the whole struct
 * INCLUDING its symbols/values arrays across matches makes reuse allocation-free
 * after warm-up. Single-threaded REPL, so a plain static free-list is safe and
 * reentrant (nested matches just draw more from the pool). Mirrors the bounded
 * Expr node pool in expr.c; pooled envs are program-lifetime "still reachable",
 * not leaked. env_free preserves its per-binding cleanup (the strdup'd keys and
 * the value refs are still released); only the arrays and the struct are kept. */
#define ENV_POOL_CAP 256
static MatchEnv* g_env_pool[ENV_POOL_CAP];
static size_t    g_env_pool_len = 0;

MatchEnv* env_new(void) {
    if (g_env_pool_len > 0) {
        MatchEnv* env = g_env_pool[--g_env_pool_len];
        env->count = 0;                 /* arrays + capacity retained */
        env->callback = NULL;
        env->callback_data = NULL;
        env->top_pattern = NULL;
        return env;
    }
    MatchEnv* env = malloc(sizeof(MatchEnv));
    env->count = 0;
    env->capacity = 8;
    env->callback = NULL;
    env->callback_data = NULL;
    env->top_pattern = NULL;
    env->symbols = malloc(sizeof(char*) * env->capacity);
    env->values = malloc(sizeof(Expr*) * env->capacity);
    return env;
}

/* KEY OWNERSHIP: env BORROWS the symbol-name pointer, it does not own it.
 * Every env_set caller passes either an interned symbol name (permanent, owned
 * by the interner) or a string literal ("$OptionsPattern$") -- both outlive any
 * MatchEnv, which is per-match-attempt. So the key is stored directly with no
 * strdup, and env_free / env_rollback never free it. That removes one malloc +
 * one free per binding (millions in recursive rewriting), on top of the struct
 * pool above. Comparison is pointer-first (interned keys hit immediately -- the
 * hot path, since the matcher and replace_bindings both pass
 * sym->data.symbol.name) with a strcmp fallback so the literal callers, whose
 * pointer differs per translation unit, still match correctly. */
void env_free(MatchEnv* env) {
    if (!env) return;
    for (size_t i = 0; i < env->count; i++)
        expr_free(env->values[i]);            /* keys are borrowed -- never freed */
    env->count = 0;
    if (g_env_pool_len < ENV_POOL_CAP) {
        g_env_pool[g_env_pool_len++] = env;   /* recycle struct + arrays */
        return;
    }
    free(env->symbols);
    free(env->values);
    free(env);
}

void env_set(MatchEnv* env, const char* symbol, Expr* value) {
    for (size_t i = 0; i < env->count; i++) {
        if (env->symbols[i] == symbol || strcmp(env->symbols[i], symbol) == 0) {
            expr_free(env->values[i]);
            env->values[i] = expr_copy(value);
            return;
        }
    }
    if (env->count >= env->capacity) {
        env->capacity *= 2;
        env->symbols = realloc(env->symbols, sizeof(char*) * env->capacity);
        env->values = realloc(env->values, sizeof(Expr*) * env->capacity);
    }
    env->symbols[env->count] = (char*)symbol;   /* borrowed; see note above */
    env->values[env->count] = expr_copy(value);
    env->count++;
}

static void env_rollback(MatchEnv* env, size_t saved_count) {
    while (env->count > saved_count) {
        env->count--;
        expr_free(env->values[env->count]);   /* keys are borrowed -- never freed */
    }
}

Expr* env_get(MatchEnv* env, const char* symbol) {
    for (size_t i = 0; i < env->count; i++) {
        if (env->symbols[i] == symbol || strcmp(env->symbols[i], symbol) == 0) {
            return env->values[i];
        }
    }
    return NULL;
}

static bool is_pattern(Expr* e, Expr** sym_out, Expr** pat_out) {
    if (e->type == EXPR_FUNCTION && e->data.function.head->type == EXPR_SYMBOL) {
        if (e->data.function.head->data.symbol.name == SYM_Pattern) {
            if (e->data.function.arg_count == 2) {
                if (sym_out) *sym_out = e->data.function.args[0];
                if (pat_out) *pat_out = e->data.function.args[1];
                return true;
            }
        }
    }
    return false;
}

static bool is_blank(Expr* e, Expr** head_out) {
    if (e->type == EXPR_FUNCTION && e->data.function.head->type == EXPR_SYMBOL) {
        if (e->data.function.head->data.symbol.name == SYM_Blank) {
            if (e->data.function.arg_count == 0) {
                if (head_out) *head_out = NULL;
                return true;
            } else if (e->data.function.arg_count == 1) {
                if (head_out) *head_out = e->data.function.args[0];
                return true;
            }
        }
    }
    return false;
}

static bool is_sequence_blank(Expr* e, Expr** head_out, int* min_len) {
    if (e->type == EXPR_FUNCTION && e->data.function.head->type == EXPR_SYMBOL) {
        const char* head = e->data.function.head->data.symbol.name;
        if (head == SYM_BlankSequence) {
            if (min_len) *min_len = 1;
            if (head_out) *head_out = (e->data.function.arg_count == 1) ? e->data.function.args[0] : NULL;
            return true;
        }
        if (head == SYM_BlankNullSequence) {
            if (min_len) *min_len = 0;
            if (head_out) *head_out = (e->data.function.arg_count == 1) ? e->data.function.args[0] : NULL;
            return true;
        }
    }
    return false;
}

/* Read one repeat bound of a Repeated / RepeatedNull spec. A machine integer
 * is taken as is; a BigInt (necessarily huge -- normalised integers that fit
 * stay EXPR_INTEGER) or one beyond INT_MAX is as good as unbounded, since no
 * argument list can be that long, and becomes -1; `Infinity` is -1 when
 * `allow_inf`. Anything else -- a symbol, a negative count, a real -- is not
 * a valid bound and yields false. */
static bool repeat_bound(Expr* b, bool allow_inf, int* out) {
    if (b->type == EXPR_INTEGER) {
        if (b->data.integer < 0) return false;
        *out = (b->data.integer > INT_MAX) ? -1 : (int)b->data.integer;
        return true;
    }
    if (b->type == EXPR_BIGINT) {
        if (mpz_sgn(b->data.bigint) < 0) return false;
        *out = -1;
        return true;
    }
    if (allow_inf && b->type == EXPR_SYMBOL && b->data.symbol.name == SYM_Infinity) {
        *out = -1;
        return true;
    }
    return false;
}

/* Repeated::range / Repeated::order for a malformed repeat spec, as WL. */
static void repeat_spec_message(Expr* e, bool order) {
    char* txt = expr_to_string(e);
    if (order)
        mth_message(e->data.function.head->data.symbol.name, "order",
                    "Range out of order in %s.", txt ? txt : "");
    else
        mth_message(e->data.function.head->data.symbol.name, "range",
                    "Range specification in integers (max or {min, max}) expected "
                    "at position 2 in %s.", txt ? txt : "");
    free(txt);
}

/* Recognise Repeated[p] / RepeatedNull[p] with an optional repeat spec:
 *   p..  / p...                    1.. / 0.. repeats (unbounded)
 *   Repeated[p, max]               min .. max   (min = 1, or 0 for RepeatedNull)
 *   Repeated[p, {n}]               exactly n
 *   Repeated[p, {min, max}]        min .. max   (max may be Infinity)
 * On success fills rep_pat / min_len / max_len (-1 = no upper bound). An
 * invalid spec (negative or symbolic count, min > max, wrong shape) returns
 * false: the expression is then not a repeat pattern at all and is matched
 * structurally, so it matches nothing -- it must never silently behave as an
 * unbounded repeat. With `warn` (the matcher's dispatch sites) the malformed
 * spec is reported as Repeated::range / Repeated::order, as in WL. */
static bool is_repeated_ex(Expr* e, Expr** rep_pat, int* min_len, int* max_len, bool warn) {
    if (e->type != EXPR_FUNCTION || e->data.function.head->type != EXPR_SYMBOL) return false;
    const char* head = e->data.function.head->data.symbol.name;
    bool is_rep = (head == SYM_Repeated);
    bool is_rep_null = (head == SYM_RepeatedNull);
    if (!is_rep && !is_rep_null) return false;
    size_t argc = e->data.function.arg_count;
    if (argc < 1 || argc > 2) return false;

    int lo = is_rep ? 1 : 0;
    int hi = -1;   /* -1 means infinity */
    if (argc == 2) {
        Expr* spec = e->data.function.args[1];
        if (spec->type == EXPR_FUNCTION && spec->data.function.head->type == EXPR_SYMBOL
            && spec->data.function.head->data.symbol.name == SYM_List) {
            size_t sc = spec->data.function.arg_count;
            if (sc == 1) {
                if (!repeat_bound(spec->data.function.args[0], false, &lo)) goto bad_range;
                if (lo == -1) return false;   /* exactly 10^30 args: valid, unsatisfiable */
                hi = lo;
            } else if (sc == 2) {
                if (!repeat_bound(spec->data.function.args[0], false, &lo)) goto bad_range;
                if (lo == -1) return false;   /* likewise at least 10^30 */
                if (!repeat_bound(spec->data.function.args[1], true, &hi)) goto bad_range;
            } else {
                goto bad_range;
            }
        } else if (!repeat_bound(spec, true, &hi)) {
            goto bad_range;
        }
        if (hi != -1 && lo > hi) {
            if (warn) repeat_spec_message(e, true);
            return false;
        }
    }
    *rep_pat = e->data.function.args[0];
    *min_len = lo;
    *max_len = hi;
    return true;

bad_range:
    if (warn) repeat_spec_message(e, false);
    return false;
}

static bool is_repeated(Expr* e, Expr** rep_pat, int* min_len, int* max_len) {
    return is_repeated_ex(e, rep_pat, min_len, max_len, false);
}

/* True when pattern element `p` matches a *variable* number of subject
 * elements -- a BlankSequence (__) or BlankNullSequence (___), optionally
 * wrapped in Pattern[name, ...]. These are the elements whose Orderless
 * subset enumeration is exponential, so for an Orderless head they must be
 * matched only after the fixed-arity elements have pruned. */
static bool pat_is_variable_seq(Expr* p) {
    Expr* sym = NULL;
    Expr* inner = NULL;
    if (!is_pattern(p, &sym, &inner)) inner = p;
    Expr* h = NULL;
    int ml = 0;
    return is_sequence_blank(inner, &h, &ml);
}

/* True when pattern element `p` is (or wraps) an Optional -- x_. or
 * Optional[...]. Optional carries a *positional* default resolved via
 * get_default_value(pat_head, original_index, ...), so a pattern list that
 * contains one must NOT be reordered: the index would shift. */
static bool pat_has_optional(Expr* p) {
    Expr* cur = p;
    while (cur && cur->type == EXPR_FUNCTION
           && cur->data.function.head
           && cur->data.function.head->type == EXPR_SYMBOL
           && cur->data.function.arg_count >= 1) {
        const char* h = cur->data.function.head->data.symbol.name;
        if (h == SYM_Optional) return true;
        if (h == SYM_Pattern && cur->data.function.arg_count == 2) {
            cur = cur->data.function.args[1];
        } else if (h == SYM_Shortest || h == SYM_Longest) {
            cur = cur->data.function.args[0];
        } else {
            break;
        }
    }
    return false;
}

/* True when `p` can never match a group `flat_head[x1, ..., xk]` (k >= 2): a
 * literal atom, or a call on a fixed symbol head other than `flat_head` that is
 * neither a pattern construct nor OneIdentity (so no default can collapse it). */
static bool flat_elem_rejects_group(const Expr* p, const char* flat_head) {
    if (p->type != EXPR_FUNCTION) return true;
    const Expr* hd = p->data.function.head;
    if (hd->type != EXPR_SYMBOL) return false;
    const char* h = hd->data.symbol.name;
    if (h == flat_head
        || h == SYM_Blank || h == SYM_BlankSequence || h == SYM_BlankNullSequence
        || h == SYM_Pattern || h == SYM_HoldPattern || h == SYM_Alternatives
        || h == SYM_Condition || h == SYM_PatternTest || h == SYM_Optional
        || h == SYM_Shortest || h == SYM_Longest || h == SYM_Repeated
        || h == SYM_RepeatedNull || h == SYM_Except || h == SYM_Verbatim
        || h == SYM_PatternSequence || h == SYM_OptionsPattern)
        return false;
    SymbolDef* d = symtab_get_def(h);
    return !(d && (d->attributes & ATTR_ONEIDENTITY));
}

/* True when `p` contains no pattern construct at any depth (heads included),
 * so it can only match an expression structurally equal to itself. */
static bool pat_is_literal(const Expr* p) {
    if (p->type != EXPR_FUNCTION) return true;
    const Expr* hd = p->data.function.head;
    if (hd->type == EXPR_SYMBOL) {
        const char* h = hd->data.symbol.name;
        if (h == SYM_Blank || h == SYM_BlankSequence || h == SYM_BlankNullSequence ||
            h == SYM_Pattern || h == SYM_Alternatives || h == SYM_Except ||
            h == SYM_PatternTest || h == SYM_Condition || h == SYM_Optional ||
            h == SYM_Repeated || h == SYM_RepeatedNull || h == SYM_HoldPattern ||
            h == SYM_Verbatim || h == SYM_Longest || h == SYM_Shortest ||
            h == SYM_OptionsPattern || h == SYM_PatternSequence ||
            h == SYM_OrderlessPatternSequence || h == SYM_KeyValuePattern)
            return false;
    } else if (!pat_is_literal(hd)) {
        return false;
    }
    for (size_t i = 0; i < p->data.function.arg_count; i++)
        if (!pat_is_literal(p->data.function.args[i])) return false;
    return true;
}

/* ------------------------------------------------------------------
 *  Span analysis for elements of an ordered Flat call pattern
 *  (Dot, **, And, Or, StringJoin, a user head with only Flat).
 *
 *  In a Flat head a non-sequence element normally absorbs a GROUP of 1..n
 *  operands (shortest first), so `ff[x_, y_]` against ff[a1..an] tries every
 *  split.  These predicates let the matcher -- and the rule-application
 *  leftover search (match_flat_leftover) -- skip splits that provably cannot
 *  succeed WITHOUT changing which split is found first:
 *
 *   - flat_elem_single_span: the element can only ever consume exactly ONE
 *     operand (a literal, or a call like op[j_] whose head is not the Flat
 *     head nor a pattern construct).  It is an "anchor": in an ordered head it
 *     sits at one known operand position.
 *   - flat_elem_screenable: an anchor whose standalone match (fresh, empty
 *     env) is a NECESSARY condition for its match under any earlier bindings.
 *     Pattern-variable reuse only adds equality constraints and a PatternTest
 *     sees only the element, so this holds unless the element contains a
 *     construct that is not monotone in the bindings: a Condition (its guard
 *     can read an earlier binding) or an Except (more bindings make its inner
 *     pattern STRICTER, i.e. the Except looser).
 *   - flat_elem_min_span: a lower bound on the operands an element consumes.
 * ------------------------------------------------------------------ */

/* True when pattern element `p` (an argument of a `flat_head` call pattern) can
 * only ever match exactly ONE operand -- never a group, never zero.  Mirrors
 * the matcher's own Flat-element span pruning in match_args_internal: a literal
 * atom, or a call whose head is neither `flat_head` nor a pattern construct and
 * which cannot OneIdentity-absorb a group.  Conservative: false when unsure.
 * PatternSequence / OrderlessPatternSequence are spliced or reordered by the
 * matcher rather than matched against one operand, so they are excluded. */
static bool flat_elem_single_span(Expr* p, const char* flat_head) {
    if (!p) return false;
    if (p->type == EXPR_FUNCTION && p->data.function.head &&
        p->data.function.head->type == EXPR_SYMBOL &&
        p->data.function.head->data.symbol.name == SYM_Pattern &&
        p->data.function.arg_count == 2)
        p = p->data.function.args[1];              /* name:inner -> inner */
    if (p->type != EXPR_FUNCTION) return true;     /* literal atom */
    Expr* hd = p->data.function.head;
    if (!hd || hd->type != EXPR_SYMBOL) return false;
    const char* h = hd->data.symbol.name;
    if (h == flat_head || h == SYM_Blank || h == SYM_BlankSequence ||
        h == SYM_BlankNullSequence || h == SYM_Pattern || h == SYM_HoldPattern ||
        h == SYM_Alternatives || h == SYM_Condition || h == SYM_PatternTest ||
        h == SYM_Optional || h == SYM_Shortest || h == SYM_Longest ||
        h == SYM_Repeated || h == SYM_RepeatedNull || h == SYM_Except ||
        h == SYM_Verbatim || h == SYM_OptionsPattern ||
        h == SYM_PatternSequence || h == SYM_OrderlessPatternSequence)
        return false;
    if ((get_attributes(h) & ATTR_ONEIDENTITY) && p->data.function.arg_count == 1)
        return false;
    return true;
}

/* True when `p` contains (at any depth, heads included) a call whose head is
 * the symbol `h1` or `h2` (h2 may be NULL). */
static bool pat_contains_head(const Expr* p, const char* h1, const char* h2) {
    if (!p || p->type != EXPR_FUNCTION) return false;
    const Expr* hd = p->data.function.head;
    if (hd && hd->type == EXPR_SYMBOL &&
        (hd->data.symbol.name == h1 || (h2 && hd->data.symbol.name == h2)))
        return true;
    if (pat_contains_head(hd, h1, h2)) return true;
    for (size_t i = 0; i < p->data.function.arg_count; i++)
        if (pat_contains_head(p->data.function.args[i], h1, h2)) return true;
    return false;
}

/* Condition and Except: the constructs whose match is not monotone in the
 * earlier bindings. */
static bool pat_has_nonmonotone(const Expr* p) {
    return pat_contains_head(p, SYM_Condition, SYM_Except);
}

static bool flat_elem_screenable(Expr* p, const char* flat_head) {
    return flat_elem_single_span(p, flat_head) && !pat_has_nonmonotone(p);
}

/* An anchor that is a call op[...] (after its name) can only match a call
 * with head op -- unless op is OneIdentity or the pattern carries an Optional
 * default (x^n_., op[x_:d, y_] may match a non-op operand).  Returns op when
 * that head shortcut is sound, else NULL.  Computed once per anchor, not per
 * operand, since it costs an attribute lookup and a pattern walk. */
static const char* flat_anchor_head(Expr* anchor) {
    Expr* a = anchor;
    if (a->type == EXPR_FUNCTION && a->data.function.head->type == EXPR_SYMBOL &&
        a->data.function.head->data.symbol.name == SYM_Pattern &&
        a->data.function.arg_count == 2)
        a = a->data.function.args[1];
    if (a->type != EXPR_FUNCTION || !a->data.function.head ||
        a->data.function.head->type != EXPR_SYMBOL) return NULL;
    const char* h = a->data.function.head->data.symbol.name;
    if ((get_attributes(h) & ATTR_ONEIDENTITY) || pat_contains_head(a, SYM_Optional, NULL))
        return NULL;
    return h;
}

/* Standalone (fresh-env) match of one operand against an anchor element, with
 * the head shortcut `fast_head` (from flat_anchor_head, or NULL) tried first. */
static bool flat_anchor_accepts(Expr* operand, Expr* anchor, const char* fast_head) {
    if (fast_head) {
        switch (operand->type) {
            case EXPR_INTEGER: case EXPR_BIGINT: case EXPR_REAL:
            case EXPR_SYMBOL:  case EXPR_STRING:
#ifdef USE_MPFR
            case EXPR_MPFR:
#endif
                return false;      /* an atom is never a fast_head[...] call */
            case EXPR_FUNCTION:
                if (!operand->data.function.head ||
                    operand->data.function.head->type != EXPR_SYMBOL ||
                    operand->data.function.head->data.symbol.name != fast_head)
                    return false;
                break;
            default:               /* packed arrays, associations, ...: ask match() */
                break;
        }
    }
    MatchEnv* probe = env_new();
    bool ok = match(operand, anchor, probe);
    env_free(probe);
    return ok;
}

/* Lower bound on the number of operands element `p` consumes in a Flat call:
 * 1 for an anchor, a (named / tested / guarded) Blank or BlankSequence; 0 for
 * anything else (BlankNullSequence, Optional, Repeated, PatternSequence, ...),
 * which is always a sound bound. */
static size_t flat_elem_min_span(Expr* p, const char* flat_head) {
    if (flat_elem_single_span(p, flat_head)) return 1;
    for (;;) {
        if (!p || p->type != EXPR_FUNCTION || !p->data.function.head ||
            p->data.function.head->type != EXPR_SYMBOL) return 0;
        const char* h = p->data.function.head->data.symbol.name;
        if (h == SYM_Pattern && p->data.function.arg_count == 2)
            p = p->data.function.args[1];
        else if ((h == SYM_PatternTest || h == SYM_Condition) &&
                 p->data.function.arg_count == 2)
            p = p->data.function.args[0];
        else
            return (h == SYM_Blank || h == SYM_BlankSequence) ? 1 : 0;
    }
}

/* Could the k-operand group  pat_head[subset[0..k)]  (k > 1) be SameQ to the
 * existing binding `existing`?  Exactly expr_eq(existing, <that group>) except
 * that it never builds the group, and it answers "maybe" (true) for a packed /
 * NDArray binding, which expr_eq compares by its own rules.  NULL (unbound)
 * is "maybe" too. */
static bool flat_group_may_equal(const Expr* existing, const Expr* pat_head,
                                 Expr** subset, size_t k) {
    if (!existing || existing->type == EXPR_NDARRAY) return true;
    if (existing->type != EXPR_FUNCTION) return false;   /* group is a call */
    if (existing->data.function.arg_count != k) return false;
    if (!expr_eq(existing->data.function.head, pat_head)) return false;
    for (size_t i = 0; i < k; i++)
        if (!expr_eq(existing->data.function.args[i], subset[i])) return false;
    return true;
}

/* Cap on the top-level pattern-element count for which the Orderless
 * constraint-reordering buffer is stack-allocated. Real patterns have a
 * handful of top-level elements; above this we skip the (purely
 * performance) reordering and match in source order. */
#define MATCH_REORDER_CAP 64

typedef struct ParentMatch {
    Expr** exprs;
    size_t n_exprs;
    Expr** pats;
    size_t n_pats;
    Expr* condition;
    Expr* pat_head;
    size_t total_pats;
    struct ParentMatch* parent;
} ParentMatch;

static bool match_internal(Expr* expr, Expr* pattern, MatchEnv* env, ParentMatch* parent);
/* Is `h` a OneIdentity symbol? Out of line so the recursive matcher's frame
 * does not grow (deep-nesting cases sit close to the C stack limit). */
static bool head_is_oneidentity(const Expr* h) {
    return h->type == EXPR_SYMBOL &&
           (get_attributes(h->data.symbol.name) & ATTR_ONEIDENTITY);
}

static bool match_internal_impl(Expr* expr, Expr* pattern, MatchEnv* env, ParentMatch* parent);
static bool match_args_internal(Expr** exprs, size_t n_exprs, Expr** pats, size_t n_pats, MatchEnv* env, Expr* condition, Expr* pat_head, size_t total_pats, ParentMatch* parent);
/* Success continuation for rep_units_match: called once the units cover the
 * whole span; returning false makes the search try the next way to split it. */
typedef bool (*RepCont)(MatchEnv* env, void* data);
static bool rep_units_match(Expr** xs, size_t n, Expr** unit, size_t m,
                            int min_reps, int max_reps,
                            MatchEnv* env, Expr* pat_head,
                            RepCont k, void* kdata);
static bool rep_unit_split(Expr* rep_pat, Expr** slot, Expr*** unit, size_t* m);
static bool pat_is_single_elem(Expr* p);
static void pat_width(Expr* p, size_t* lo, size_t* hi);

/* C-stack-overflow guard for the mutually-recursive matcher.  match_internal
 * is the per-expression-level entry (the descent that overflows on deeply
 * nested subjects/patterns), so counting live depth here bounds the stack.
 * The bound is the SAME user-visible $RecursionLimit the evaluator honours, so
 * it is adjustable and WL-consistent rather than a hidden magic cap; exceeding
 * it degrades to a graceful non-match instead of a SIGSEGV.  Raising
 * $RecursionLimit far past the C-stack capacity can still overflow, exactly as
 * for the evaluator itself. */
static int match_recursion_depth = 0;
static bool match_internal(Expr* expr, Expr* pattern, MatchEnv* env, ParentMatch* parent) {
    if (match_recursion_depth >= eval_get_recursion_limit()) return false;
    match_recursion_depth++;
    bool r = match_internal_impl(expr, pattern, env, parent);
    match_recursion_depth--;
    return r;
}

/* Evaluate a Condition guard with the current bindings substituted.
 *
 * Mathematica parses nested `/;` right-associatively, so
 *   f[a, b] /; c1 /; c2
 * becomes
 *   Condition[f[a, b], Condition[c1, c2]]
 * with the guard being `Condition[c1, c2]`.  Mathilda's `Condition`
 * head has no evaluator rule that collapses `Condition[True, True]`
 * to `True` (it's primarily a pattern wrapper), so a naive "evaluate
 * and check `result == True`" gives False here.  We treat any chain
 * of Condition[...] guards as conjunction: every leaf must evaluate
 * to True.  Equivalent to substituting `And` for `Condition` at the
 * guard level, without changing the global Condition semantics. */
/* Recursively AND the leaves of a (possibly nested) Condition guard.
 *
 * `Condition` is HoldAll, so a nested guard like Condition[c1, c2] is NOT
 * reduced by the top-level evaluate() -- its arguments arrive unevaluated.
 * Each leaf that is neither literally True nor another Condition must
 * therefore be evaluated on its own before we can judge it.  Recursion (vs a
 * fixed-size stack) means an arbitrarily deep `/;` chain is handled without a
 * magic cap. */
static bool guard_all_true(Expr* g) {
    if (!g) return false;
    if (g->type == EXPR_SYMBOL && g->data.symbol.name == SYM_True) return true;
    if (g->type == EXPR_FUNCTION && g->data.function.head
        && g->data.function.head->type == EXPR_SYMBOL
        && g->data.function.head->data.symbol.name == SYM_Condition
        && g->data.function.arg_count == 2) {
        return guard_all_true(g->data.function.args[0])
            && guard_all_true(g->data.function.args[1]);
    }
    /* evaluate() borrows its argument (it copies internally), so the temporary
     * we hand it is ours to free -- keep it in a variable rather than inlining
     * expr_copy(g), which would orphan one node per non-True guard leaf. */
    Expr* gc = expr_copy(g);
    Expr* leaf = eval_callback(gc);
    expr_free(gc);
    bool t = leaf && leaf->type == EXPR_SYMBOL && leaf->data.symbol.name == SYM_True;
    if (leaf) expr_free(leaf);
    return t;
}

/* Head of the definition whose rule is being matched (see match_set_rule_head),
 * so a `/;` guard can resolve OptionValue like the rule's RHS does. */
static const char* g_guard_rule_head = NULL;

const char* match_set_rule_head(const char* head) {
    const char* prev = g_guard_rule_head;
    g_guard_rule_head = head;
    return prev;
}

static bool eval_guard_true(Expr* guard, MatchEnv* env) {
    Expr* expanded = replace_bindings(guard, env);
    /* A guard of a rule whose LHS carries OptionsPattern[] sees OptionValue
     * exactly as the RHS does: `f[x_, OptionsPattern[]] := x /; OptionValue[o] == 1`
     * stores the condition in the pattern, so it is evaluated here, before
     * the RHS is ever built. */
    Expr* opts = g_guard_rule_head ? env_get(env, "$OptionsPattern$") : NULL;
    if (opts) {
        const char* head = g_guard_rule_head;
        Expr* oh = env_get(env, "$OptionsPatternHead$");
        if (oh && oh->type == EXPR_SYMBOL) head = oh->data.symbol.name;
        Expr* resolved = optionvalue_inject_context(expanded, head, opts);
        expr_free(expanded);
        expanded = resolved;
    }
    Expr* result = eval_callback(expanded);
    expr_free(expanded);
    if (!result) return false;
    bool ok = guard_all_true(result);
    expr_free(result);
    return ok;
}

/* Evaluate the guard chain of a Condition-wrapped SEQUENCE pattern element such
 * as `block___ /; c1 /; c2`, i.e. Condition[Condition[base, c1], c2] where the
 * innermost `base` is the (possibly named) sequence blank / Repeated -- NOT a
 * guard.  eval_guard_true cannot be handed the whole node (its guard_all_true
 * would try to evaluate `base` as a boolean leaf and fail); instead we descend
 * args[0] to skip the base and evaluate only the args[1] guards.  Mathematica's
 * `p /; c1 /; c2` checks the inner guard first, so we recurse before evaluating
 * this layer's guard, short-circuiting on the first failure (this is what makes
 * `block___ /; Length[{block}] > 0 /; Total[{block}] == k` test the cheap length
 * guard before the sum).  Bindings for the sequence variable must already be in
 * `env`.  Returns true once the base (a non-Condition) is reached. */
static bool eval_seq_conditions(Expr* p, MatchEnv* env) {
    if (p && p->type == EXPR_FUNCTION && p->data.function.head
        && p->data.function.head->type == EXPR_SYMBOL
        && p->data.function.head->data.symbol.name == SYM_Condition
        && p->data.function.arg_count == 2) {
        return eval_seq_conditions(p->data.function.args[0], env)   /* inner /; first */
            && eval_guard_true(p->data.function.args[1], env);
    }
    return true;   /* reached the base pattern: no guards left */
}

#include "part.h" // for expr_head

/* Head-typed Blank check: does `expr` have head `b_head`?  b_head == NULL is a
 * plain Blank[] and matches anything.  Covers every atom type -- Integer and
 * BigInt both report Integer, Real and MPFR both report Real, plus Symbol and
 * String -- and the Rational/Complex "structured atoms", which are EXPR_FUNCTION
 * nodes with a Rational/Complex head (is_atom() returns true for them, which is
 * why the old get_expr_head_borrowed()-based check missed `_Rational`/`_Complex`
 * entirely). */
static bool blank_head_matches(Expr* expr, Expr* b_head) {
    if (!b_head) return true;
    if (!expr) return false;
    if (expr->type == EXPR_FUNCTION && expr->data.function.head) {
        return expr_eq(expr->data.function.head, b_head);
    }
    if (b_head->type != EXPR_SYMBOL) return false;
    const char* hn = b_head->data.symbol.name;
    switch (expr->type) {
        case EXPR_INTEGER:
        case EXPR_BIGINT:  return hn == SYM_Integer;
        case EXPR_REAL:    return hn == SYM_Real;
#ifdef USE_MPFR
        case EXPR_MPFR:    return hn == SYM_Real;
#endif
        case EXPR_SYMBOL:  return hn == SYM_Symbol;
        case EXPR_STRING:  return hn == SYM_String;
        case EXPR_SPARSE:  return hn == SYM_SparseArray;
        default:           return false;
    }
}

/* PatternTest predicate: is test_func[elem] === True?  Shared by the
 * single-element PatternTest handler and the __?test / ___?test per-element
 * sequence check. */
static bool apply_pattern_test(Expr* elem, Expr* test_func) {
    Expr* call_args[1] = { expr_copy(elem) };
    Expr* test_call = expr_new_function(expr_copy(test_func), call_args, 1);
    Expr* result = eval_callback(test_call);
    expr_free(test_call);
    bool ok = result && result->type == EXPR_SYMBOL
              && result->data.symbol.name == SYM_True;
    if (result) expr_free(result);
    return ok;
}

static bool next_combination(int* comb, int n, int k) {
    if (k == 0) return false;
    int i = k - 1;
    while (i >= 0 && comb[i] == n - k + i) {
        i--;
    }
    if (i < 0) return false;
    comb[i]++;
    for (int j = i + 1; j < k; j++) {
        comb[j] = comb[j - 1] + 1;
    }
    return true;
}

/* Split `exprs` into the `comb`-selected `subset` and the `remainder`. `subset`
 * may be NULL: callers that only consume the remainder (an unnamed, untyped
 * sequence blank -- a plain `__`/`___` whose matched elements are never
 * inspected) pass NULL to skip the copy while the comb index still advances so
 * the remainder is partitioned correctly. */
static void extract_subset(Expr** exprs, size_t n_exprs, int* comb, int k, Expr** subset, Expr** remainder) {
    int c_idx = 0;
    int r_idx = 0;
    for (int i = 0; i < (int)n_exprs; i++) {
        if (c_idx < k && comb[c_idx] == i) {
            if (subset) subset[c_idx] = exprs[i];
            c_idx++;
        } else {
            remainder[r_idx++] = exprs[i];
        }
    }
}

// Call parent if exists, otherwise return true (or check top condition)

static ParentMatch top_level_sentinel = {0};

/* Bind every pattern name in `p` that the match left unbound to Sequence[].
 * A name can stay unbound only when it sits in an Alternatives branch that was
 * not taken (or under Except); it then stands for
 * nothing on the right-hand side: 1 /. (x_ | y_) :> {x, y} is {1}. Done once,
 * at the point a whole match is accepted -- not when the branch is chosen --
 * so a later occurrence can still bind the name: f[1, 2] /. f[x_ | y_, y_]
 * gives x -> 1, y -> 2. Only pattern positions are walked: Verbatim content
 * is literal data, and the test / condition / default operand of PatternTest,
 * Condition and Optional is an expression whose own patterns (a nested
 * MatchQ, say) belong to an inner scope. Runs only on success, so failed
 * attempts -- the DownValue hot path -- pay nothing. */
static void bind_leftover_names(Expr* p, MatchEnv* env) {
    if (p->type != EXPR_FUNCTION) return;
    Expr* h = p->data.function.head;
    size_t n = p->data.function.arg_count;
    if (h->type == EXPR_SYMBOL) {
        const char* hn = h->data.symbol.name;
        if (hn == SYM_Verbatim) return;
        if (hn == SYM_Pattern && n == 2 && p->data.function.args[0]->type == EXPR_SYMBOL) {
            const char* name = p->data.function.args[0]->data.symbol.name;
            if (!env_get(env, name)) {
                Expr* empty = expr_new_function(expr_new_symbol(SYM_Sequence), NULL, 0);
                env_set(env, name, empty);
                expr_free(empty);
            }
            bind_leftover_names(p->data.function.args[1], env);
            return;
        }
        if ((hn == SYM_PatternTest || hn == SYM_Condition || hn == SYM_Optional) && n >= 1) {
            bind_leftover_names(p->data.function.args[0], env);
            return;
        }
    } else {
        bind_leftover_names(h, env);
    }
    for (size_t i = 0; i < n; i++)
        bind_leftover_names(p->data.function.args[i], env);
}

static bool call_parent(MatchEnv* env, ParentMatch* parent) {
    if (parent == &top_level_sentinel) {
        size_t saved = env->count;
        if (env->top_pattern) bind_leftover_names(env->top_pattern, env);
        if (env->callback) {
            if (env->callback(env, env->callback_data)) return true;
            env_rollback(env, saved);
            return false;
        }
        return true;
    }
    if (parent) {
        return match_args_internal(parent->exprs, parent->n_exprs, parent->pats, parent->n_pats, env, parent->condition, parent->pat_head, parent->total_pats, parent->parent);
    }
    return true;
}

/* RepCont adapter: continue with the enclosing match. */
static bool rep_cont_parent(MatchEnv* env, void* data) {
    return call_parent(env, (ParentMatch*)data);
}

/* True iff `e` is a two-argument Rule/RuleDelayed. */
static bool kvp_is_rule2(Expr* e) {
    return e->type == EXPR_FUNCTION && e->data.function.head->type == EXPR_SYMBOL &&
           (e->data.function.head->data.symbol.name == SYM_Rule ||
            e->data.function.head->data.symbol.name == SYM_RuleDelayed) &&
           e->data.function.arg_count == 2;
}

/* Backtracking matcher for KeyValuePattern requirements: assign requirement
 * `ri` to some entry of `subj` whose key matches ki and value matches pi, then
 * recurse to the next requirement; on failure, undo the bindings and try
 * another entry. When all requirements are satisfied (ri == nreq), call_parent
 * finalizes — and if that fails, the search backtracks to try other
 * assignments. This makes e.g. KeyValuePattern[{k_ -> _, _ -> k_}] solvable. */
static bool kvp_match_reqs(Expr** reqs, size_t nreq, size_t ri, Expr* subj,
                           MatchEnv* env, ParentMatch* parent) {
    if (ri == nreq) return call_parent(env, parent);
    Expr* rule = reqs[ri];
    if (!kvp_is_rule2(rule)) return false;
    Expr* kpat = rule->data.function.args[0];
    Expr* vpat = rule->data.function.args[1];
    for (size_t i = 0; i < subj->data.function.arg_count; i++) {
        Expr* entry = subj->data.function.args[i];
        if (!kvp_is_rule2(entry)) continue;
        size_t es = env->count;
        if (match_internal(entry->data.function.args[0], kpat, env, NULL) &&
            match_internal(entry->data.function.args[1], vpat, env, NULL) &&
            kvp_match_reqs(reqs, nreq, ri + 1, subj, env, parent)) {
            return true;
        }
        env_rollback(env, es);
    }
    return false;
}

/* If `e` is an OrderlessPatternSequence[...] (optionally wrapped in a
 * Pattern[name, ...]), return the OPS node and set *bind_out to the binding
 * name symbol (or NULL). Otherwise return NULL. */
static Expr* ops_unwrap(Expr* e, Expr** bind_out) {
    *bind_out = NULL;
    Expr* cur = e;
    if (cur && cur->type == EXPR_FUNCTION && cur->data.function.head
        && cur->data.function.head->type == EXPR_SYMBOL
        && cur->data.function.head->data.symbol.name == SYM_Pattern
        && cur->data.function.arg_count == 2) {
        *bind_out = cur->data.function.args[0];
        cur = cur->data.function.args[1];
    }
    if (cur && cur->type == EXPR_FUNCTION && cur->data.function.head
        && cur->data.function.head->type == EXPR_SYMBOL
        && cur->data.function.head->data.symbol.name == SYM_OrderlessPatternSequence) {
        return cur;
    }
    *bind_out = NULL;
    return NULL;
}

/* Assign the elements qs[qi..nq) of an OrderlessPatternSequence to the
 * still-unused elements of `run` (the contiguous block of arguments the OPS
 * occupies), in any order. A single-element pattern takes one element; a
 * literal one (no pattern construct) takes only the FIRST equal element, as
 * every equal choice leaves the same state; a variable-width one (x__, x___,
 * p..) takes an order-preserving subset, and the last pattern takes whatever
 * is left. Once every pattern is placed and the whole run consumed, the OPS
 * name (if any) binds the run in its ORIGINAL order and matching continues
 * with the arguments after the run against pats[1..]. */
static bool ops_assign(Expr** run, size_t L, bool* used,
                       Expr** qs, size_t qi, size_t nq,
                       Expr** tail, size_t n_tail,
                       Expr** pats, size_t n_pats, MatchEnv* env,
                       Expr* condition, Expr* pat_head, size_t total_pats,
                       ParentMatch* parent, Expr* bind_sym) {
    if (qi == nq) {
        for (size_t i = 0; i < L; i++) if (!used[i]) return false;
        size_t saved = env->count;
        if (bind_sym && bind_sym->type == EXPR_SYMBOL) {
            Expr* seqv = expr_new_function(expr_new_symbol(SYM_Sequence), NULL, L);
            for (size_t i = 0; i < L; i++) seqv->data.function.args[i] = expr_copy(run[i]);
            Expr* existing = env_get(env, bind_sym->data.symbol.name);
            bool okbind = true;
            if (existing) okbind = expr_eq(seqv, existing);
            else env_set(env, bind_sym->data.symbol.name, seqv);
            expr_free(seqv);
            if (!okbind) return false;
        }
        if (match_args_internal(tail, n_tail, pats + 1, n_pats - 1, env,
                                condition, pat_head, total_pats, parent))
            return true;
        env_rollback(env, saved);
        return false;
    }
    Expr* q = qs[qi];
    if (pat_is_single_elem(q)) {
        bool lit = pat_is_literal(q);
        for (size_t p = 0; p < L; p++) {
            if (used[p]) continue;
            if (lit && !expr_eq(run[p], q)) continue;
            size_t saved = env->count;
            used[p] = true;
            if (match_internal(run[p], q, env, NULL)
                && ops_assign(run, L, used, qs, qi + 1, nq, tail, n_tail, pats, n_pats,
                              env, condition, pat_head, total_pats, parent, bind_sym))
                return true;
            used[p] = false;
            env_rollback(env, saved);
            if (lit) break;
        }
        return false;
    }
    /* Variable width: choose `w` of the free elements (kept in run order). */
    size_t free_n = 0;
    for (size_t i = 0; i < L; i++) if (!used[i]) free_n++;
    size_t lo, hi;
    pat_width(q, &lo, &hi);
    size_t wmin = lo, wmax = hi < free_n ? hi : free_n;
    if (qi + 1 == nq) { if (free_n < lo || free_n > hi) return false; wmin = wmax = free_n; }
    int* freeidx = malloc((free_n ? free_n : 1) * sizeof(int));
    size_t fw = 0;
    for (size_t i = 0; i < L; i++) if (!used[i]) freeidx[fw++] = (int)i;
    Expr** sub = malloc((free_n ? free_n : 1) * sizeof(Expr*));
    int* comb = malloc((free_n ? free_n : 1) * sizeof(int));
    bool r = false;
    for (size_t w = wmin; w <= wmax && !r; w++) {
        for (size_t i = 0; i < w; i++) comb[i] = (int)i;
        do {
            for (size_t i = 0; i < w; i++) sub[i] = run[freeidx[comb[i]]];
            size_t saved = env->count;
            if (match_args_internal(sub, w, &q, 1, env, NULL, pat_head, 1, NULL)) {
                for (size_t i = 0; i < w; i++) used[freeidx[comb[i]]] = true;
                r = ops_assign(run, L, used, qs, qi + 1, nq, tail, n_tail, pats, n_pats,
                               env, condition, pat_head, total_pats, parent, bind_sym);
                for (size_t i = 0; i < w; i++) used[freeidx[comb[i]]] = false;
            }
            if (r) break;
            env_rollback(env, saved);
        } while (next_combination(comb, (int)free_n, (int)w));
    }
    free(comb); free(sub); free(freeidx);
    return r;
}


static bool match_internal_impl(Expr* expr, Expr* pattern, MatchEnv* env, ParentMatch* parent) {
    if (!pattern) return false;
    if (!expr) return false;

    /* HoldPattern[p] is equivalent to p for matching purposes. The
     * HoldAll attribute on HoldPattern prevents the inner expression from
     * being evaluated when the rule LHS is formed, so patterns that would
     * otherwise simplify (e.g. _+_ becoming 2 Blank[]) are preserved. */
    if (pattern->type == EXPR_FUNCTION && pattern->data.function.head->type == EXPR_SYMBOL &&
        pattern->data.function.head->data.symbol.name == SYM_HoldPattern &&
        pattern->data.function.arg_count == 1) {
        return match_internal(expr, pattern->data.function.args[0], env, parent);
    }

    /* Verbatim[p] matches expr iff expr is structurally identical to p, with p
     * taken LITERALLY -- its pattern constructs (Blank, Pattern, ...) are not
     * interpreted.  So Verbatim[x_] matches only the literal Pattern[x,Blank[]],
     * not an arbitrary expression. */
    if (pattern->type == EXPR_FUNCTION && pattern->data.function.head->type == EXPR_SYMBOL &&
        pattern->data.function.head->data.symbol.name == SYM_Verbatim &&
        pattern->data.function.arg_count == 1) {
        if (expr_eq(expr, pattern->data.function.args[0])) return call_parent(env, parent);
        return false;
    }

    /* Longest[p] / Shortest[p] at a non-sequence position are transparent: the
     * min/max preference only bites inside a sequence slot (handled by the
     * wrapper peeling in match_args_internal).  Unwrapping here lets a bare
     * top-level MatchQ[e, Longest[p]] / Shortest[p] work. */
    if (pattern->type == EXPR_FUNCTION && pattern->data.function.head->type == EXPR_SYMBOL &&
        (pattern->data.function.head->data.symbol.name == SYM_Longest ||
         pattern->data.function.head->data.symbol.name == SYM_Shortest) &&
        pattern->data.function.arg_count == 1) {
        return match_internal(expr, pattern->data.function.args[0], env, parent);
    }

    /* KeyValuePattern[{k1 -> p1, ...}] (or a single k -> p) matches an
     * association — or a list of rules — that contains, for each required
     * ki -> pi, an entry whose key matches ki and whose value matches pi.
     * Requirements are matched with backtracking (kvp_match_reqs), so a
     * consistent assignment is found whenever one exists even when requirements
     * share bound variables. An empty spec matches any association. This is a
     * new pattern head, so it cannot affect existing matching. */
    if (pattern->type == EXPR_FUNCTION && pattern->data.function.head->type == EXPR_SYMBOL &&
        pattern->data.function.head->data.symbol.name == SYM_KeyValuePattern &&
        pattern->data.function.arg_count == 1) {
        Expr* subj = expr;
        bool subj_ok = subj->type == EXPR_FUNCTION && subj->data.function.head->type == EXPR_SYMBOL &&
            (subj->data.function.head->data.symbol.name == SYM_Association ||
             subj->data.function.head->data.symbol.name == SYM_List);
        if (!subj_ok) return false;

        Expr* spec = pattern->data.function.args[0];
        Expr** reqs; size_t nreq;
        if (spec->type == EXPR_FUNCTION && spec->data.function.head->type == EXPR_SYMBOL &&
            spec->data.function.head->data.symbol.name == SYM_List) {
            reqs = spec->data.function.args; nreq = spec->data.function.arg_count;
        } else {
            reqs = &spec; nreq = 1;   /* a single rule */
        }

        size_t saved = env->count;
        if (kvp_match_reqs(reqs, nreq, 0, subj, env, parent)) return true;
        env_rollback(env, saved);
        return false;
    }

    // Handle Except
    if (pattern->type == EXPR_FUNCTION && pattern->data.function.head->type == EXPR_SYMBOL &&
        pattern->data.function.head->data.symbol.name == SYM_Except) {
        if (pattern->data.function.arg_count == 1) {
            size_t saved_env_count = env->count;
            if (match_internal(expr, pattern->data.function.args[0], env, NULL)) {
                env_rollback(env, saved_env_count);
                return false;
            }
            return call_parent(env, parent);
        }
        if (pattern->data.function.arg_count == 2) {
            size_t saved_env_count = env->count;
            if (!match_internal(expr, pattern->data.function.args[1], env, NULL)) {
                env_rollback(env, saved_env_count);
                return false;
            }
            size_t inner_saved = env->count;
            if (match_internal(expr, pattern->data.function.args[0], env, NULL)) {
                env_rollback(env, saved_env_count);
                return false;
            }
            env_rollback(env, inner_saved);
            return call_parent(env, parent);
        }
    }

    // Handle Condition
    if (pattern->type == EXPR_FUNCTION && pattern->data.function.head->type == EXPR_SYMBOL &&
        pattern->data.function.head->data.symbol.name == SYM_Condition) {
        if (pattern->data.function.arg_count != 2) return false;

        Expr* inner_pat = pattern->data.function.args[0];
        Expr* cond = pattern->data.function.args[1];

        /* The "structured" branch threads cond through match_args_internal
         * so the guard runs once all bindings are established. It only
         * applies when inner_pat is a literal-headed function call (so we
         * can route argument matching). When inner_pat is itself a
         * pattern wrapper -- a nested Condition/PatternTest/HoldPattern/
         * Verbatim guard OR a pattern object whose head only *looks* like
         * a function call (Pattern[x, Blank[]] for `x_`, Blank[], etc.) --
         * we MUST fall through to the simpler recursive path below, which
         * re-enters this handler for the inner wrapper. Without this gate,
         * e.g. `f[a_,b_] /; c1 /; c2` would try to structurally match
         * input's head against `Condition`, and `e_ /; cond` would match
         * input's head against `Pattern`, both wrong dispatches. */
        bool inner_is_wrapper = false;
        if (inner_pat->type == EXPR_FUNCTION && inner_pat->data.function.head
            && inner_pat->data.function.head->type == EXPR_SYMBOL) {
            const char* ih = inner_pat->data.function.head->data.symbol.name;
            if (ih == SYM_Condition || ih == SYM_PatternTest
                || ih == SYM_HoldPattern || ih == SYM_Verbatim
                || ih == SYM_Pattern || ih == SYM_Blank
                || ih == SYM_BlankSequence || ih == SYM_BlankNullSequence
                || ih == SYM_Alternatives || ih == SYM_Optional
                || ih == SYM_Repeated || ih == SYM_RepeatedNull
                || ih == SYM_Shortest || ih == SYM_Longest
                || ih == SYM_Except || ih == SYM_KeyValuePattern) {
                inner_is_wrapper = true;
            }
        }

        if (!inner_is_wrapper
            && inner_pat->type == EXPR_FUNCTION && expr->type == EXPR_FUNCTION) {
            size_t saved_env_count = env->count;
            if (match_internal(expr->data.function.head, inner_pat->data.function.head, env, NULL)) {
                bool res = match_args_internal(expr->data.function.args, expr->data.function.arg_count,
                                     inner_pat->data.function.args, inner_pat->data.function.arg_count, env, cond, inner_pat->data.function.head, inner_pat->data.function.arg_count, parent);
                if (!res) env_rollback(env, saved_env_count);
                return res;
            }
            env_rollback(env, saved_env_count);
            /* A head mismatch is final -- unless the inner head is OneIdentity,
             * where f[p1, ..., pn] can still match a lone non-f subject through
             * an Optional's default (h[g] against Times[a_., h[b_]] as a -> 1).
             * That reading lives at the end of match_internal, which the
             * generic path below reaches with the guard as continuation. */
            if (!head_is_oneidentity(inner_pat->data.function.head)) return false;
        }

        /* The guard is the inner match's continuation (an empty argument level
         * carrying `cond`, then the parent), so a guard that fails backtracks into
         * the inner pattern's other assignments instead of ending the match. */
        size_t saved_env_count = env->count;
        ParentMatch guard_cont = { NULL, 0, NULL, 0, cond, NULL, 0, parent };
        if (!match_internal(expr, inner_pat, env, &guard_cont)) {
            env_rollback(env, saved_env_count);
            return false;
        }
        return true;
    }

    Expr* b_head = NULL;
    if (is_blank(pattern, &b_head)) {
        if (blank_head_matches(expr, b_head)) return call_parent(env, parent);
        return false;
    }

    int min_len = 0;
    if (is_sequence_blank(pattern, &b_head, &min_len)) {
        if (min_len > 1) return false;
        if (blank_head_matches(expr, b_head)) return call_parent(env, parent);
        return false;
    }

    /* Top-level Repeated / RepeatedNull. A `p..` / `p...` pattern matches a
     * sequence; at the top level the single subject `expr` is a length-1
     * sequence, so it matches iff 1 is within the [min, max] repeat bounds and
     * `expr` matches the inner pattern. Without this, `MatchQ[a, a..]` was
     * False (should be True), and nested repeats routed through here --
     * `Repeated[RepeatedNull[a]]` checks each element with a top-level
     * `match_internal(elem, RepeatedNull[a])` -- also failed. */
    {
        Expr* rep_pat = NULL;
        int rmin = 0, rmax = -1;
        /* No warning here: a malformed spec inside a pattern list is
         * reported once at the sequence site, which then falls back to
         * matching it as a single element -- and so lands here again. */
        if (is_repeated(pattern, &rep_pat, &rmin, &rmax)) {
            size_t saved = env->count;
            Expr* slot = NULL; Expr** unit = NULL; size_t um = 0;
            if (rep_unit_split(rep_pat, &slot, &unit, &um)) {
                /* A variable-width body ((PatternSequence[...])..): the lone
                 * subject is a one-argument sequence for the units to cover. */
                if (rep_units_match(&expr, 1, unit, um, rmin, rmax, env, NULL,
                                    rep_cont_parent, parent)) return true;
                env_rollback(env, saved);
            } else if (rmin <= 1 && (rmax == -1 || rmax >= 1)) {
                if (match_internal(expr, rep_pat, env, NULL)) {
                    if (call_parent(env, parent)) return true;
                }
                env_rollback(env, saved);
            }
            return false;
        }
    }

    Expr* p_sym = NULL;
    Expr* p_pat = NULL;
    if (is_pattern(pattern, &p_sym, &p_pat)) {
        if (p_sym->type == EXPR_SYMBOL) {
            /* Bind (or check) the name FIRST, then match the inner pattern with
             * the real continuation. The name's value is `expr` whatever the
             * inner pattern binds, so binding early changes no result; but
             * matching the inner pattern with a NULL continuation (as this
             * once did) committed to its first assignment, so a parent that
             * failed could not backtrack into the others: ReplaceList[2 x,
             * p:x_.*y_ :> ...] found one of the three matches. */
            size_t saved_env_count = env->count;
            Expr* existing = env_get(env, p_sym->data.symbol.name);
            if (existing) {
                if (!expr_eq(expr, existing)) return false;
            } else {
                env_set(env, p_sym->data.symbol.name, expr);
            }
            if (match_internal(expr, p_pat, env, parent)) return true;
            env_rollback(env, saved_env_count);
            return false;
        }
    }

    // Handle Alternatives
    if (pattern->type == EXPR_FUNCTION && pattern->data.function.head->type == EXPR_SYMBOL &&
        pattern->data.function.head->data.symbol.name == SYM_Alternatives) {
        size_t saved_env_count = env->count;
        for (size_t i = 0; i < pattern->data.function.arg_count; i++) {
            if (match_internal(expr, pattern->data.function.args[i], env, parent)) return true;
            env_rollback(env, saved_env_count);
        }
        return false;
    }

    // Handle PatternTest
    if (pattern->type == EXPR_FUNCTION && pattern->data.function.head->type == EXPR_SYMBOL &&
        pattern->data.function.head->data.symbol.name == SYM_PatternTest) {
        if (pattern->data.function.arg_count != 2) return false;
        size_t saved_env_count = env->count;
        if (!match_internal(expr, pattern->data.function.args[0], env, NULL)) {
            return false;
        }

        if (!apply_pattern_test(expr, pattern->data.function.args[1])) {
            env_rollback(env, saved_env_count);
            return false;
        }
        if (call_parent(env, parent)) return true;
        env_rollback(env, saved_env_count);
        return false;
    }

    if (expr->type == pattern->type) {
        switch (expr->type) {
            case EXPR_INTEGER:
                if (expr->data.integer == pattern->data.integer) return call_parent(env, parent);
                return false;
            case EXPR_REAL:
                if (expr->data.real == pattern->data.real) return call_parent(env, parent);
                return false;
            case EXPR_SYMBOL: {
                Expr* bound = env_get(env, pattern->data.symbol.name);
                if (bound) {
                    if (expr_eq(expr, bound)) return call_parent(env, parent);
                    return false;
                }
                if (strcmp(expr->data.symbol.name, pattern->data.symbol.name) == 0) return call_parent(env, parent);
                return false;
            }
            case EXPR_STRING:
                if (strcmp(expr->data.string, pattern->data.string) == 0) return call_parent(env, parent);
                return false;
            case EXPR_BIGINT:
                if (mpz_cmp(expr->data.bigint, pattern->data.bigint) == 0) return call_parent(env, parent);
                return false;
#ifdef USE_MPFR
            case EXPR_MPFR:
                /* Match identical-precision MPFR values bit-for-bit. */
                if (mpfr_get_prec(expr->data.mpfr) == mpfr_get_prec(pattern->data.mpfr)
                    && mpfr_equal_p(expr->data.mpfr, pattern->data.mpfr)) {
                    return call_parent(env, parent);
                }
                return false;
#endif
            case EXPR_NDARRAY:
                /* A literal NDArray pattern matches an identical dense array. */
                if (expr_eq(expr, pattern)) return call_parent(env, parent);
                return false;
            case EXPR_COMPILED:
                /* CompiledFunction objects match by identity. */
                if (expr_eq(expr, pattern)) return call_parent(env, parent);
                return false;
            case EXPR_SPARSE:
                /* A SparseArray is a raw object: a literal one matches only a
                 * SameQ array, never by destructuring. */
                if (expr_eq(expr, pattern)) return call_parent(env, parent);
                return false;
            case EXPR_FUNCTION: {
                size_t saved_env_count = env->count;
                if (match_internal(expr->data.function.head, pattern->data.function.head, env, NULL)) {
                    if (match_args_internal(expr->data.function.args, expr->data.function.arg_count,
                                   pattern->data.function.args, pattern->data.function.arg_count, env, NULL, pattern->data.function.head, pattern->data.function.arg_count, parent)) {
                        return true;
                    }
                }
                env_rollback(env, saved_env_count);
                break;
            }
        }
    }

    /* OneIdentity: f[p1, ..., pn] may match a lone x as if it were f[x] --
     * but only through DEFAULTS: some pi must be
     * Optional, so the others can match x while it takes its default
     * (x matches Times[n_., x_] as n -> 1). Without an Optional, f[x_] does
     * not match x, nor f[a, b] (which would otherwise bind x -> f[a, b]), nor
     * f[] under Flat; and f[a] gets one match from the ordinary path above,
     * not a second one through f[f[a]]. */
    if (pattern->type == EXPR_FUNCTION && pattern->data.function.head->type == EXPR_SYMBOL) {
        SymbolDef* def = symtab_get_def(pattern->data.function.head->data.symbol.name);
        bool any_optional = false;
        if (def && (def->attributes & ATTR_ONEIDENTITY)) {
            for (size_t i = 0; i < pattern->data.function.arg_count && !any_optional; i++)
                any_optional = pat_has_optional(pattern->data.function.args[i]);
        }
        /* Under Flat, a subject that already has this head is covered by the
         * ordinary path (the defaulted element leaves the rest to group into
         * one f[...] operand), so the f[f[...]] reading would only repeat that
         * match: ReplaceList[a + x, x_. + y_ :> ...] lists {0, a + x} once. */
        if (any_optional && (def->attributes & ATTR_FLAT) && expr->type == EXPR_FUNCTION
            && expr_eq(expr->data.function.head, pattern->data.function.head))
            any_optional = false;
        if (any_optional) {
            size_t saved_env_count = env->count;
            Expr* args[1] = { expr };
            if (match_args_internal(args, 1, pattern->data.function.args, pattern->data.function.arg_count, env, NULL, pattern->data.function.head, pattern->data.function.arg_count, parent)) {
                return true;
            }
            env_rollback(env, saved_env_count);
        }
    }

    return false;
}

bool match(Expr* expr, Expr* pattern, MatchEnv* env) {
    Expr* outer = env->top_pattern;          /* match() may nest via conditions */
    env->top_pattern = pattern;
    bool r = match_internal(expr, pattern, env, &top_level_sentinel);
    env->top_pattern = outer;
    return r;
}

/* True when pattern `p` always consumes exactly ONE argument of a sequence:
 * a blank, a literal, or a structural pattern, possibly named, tested,
 * conditioned or alternated. Sequence blanks, repeats, Optional (which may
 * consume zero), Longest/Shortest and (Orderless)PatternSequence have a
 * variable width and return false. Used to give a Repeated[PatternSequence]
 * unit a fixed width so it can be matched chunk by chunk. */
static bool pat_is_single_elem(Expr* p) {
    while (p->type == EXPR_FUNCTION && p->data.function.head->type == EXPR_SYMBOL) {
        const char* h = p->data.function.head->data.symbol.name;
        size_t ac = p->data.function.arg_count;
        if (h == SYM_BlankSequence || h == SYM_BlankNullSequence
            || h == SYM_Repeated || h == SYM_RepeatedNull || h == SYM_Optional
            || h == SYM_Longest || h == SYM_Shortest
            || h == SYM_PatternSequence || h == SYM_OrderlessPatternSequence)
            return false;
        if (h == SYM_Pattern && ac == 2) { p = p->data.function.args[1]; continue; }
        if ((h == SYM_Condition || h == SYM_PatternTest) && ac == 2) {
            p = p->data.function.args[0]; continue;
        }
        if (h == SYM_HoldPattern && ac == 1) { p = p->data.function.args[0]; continue; }
        if (h == SYM_Alternatives) {
            for (size_t i = 0; i < ac; i++)
                if (!pat_is_single_elem(p->data.function.args[i])) return false;
            return true;
        }
        return true;
    }
    return true;
}

/* Static bounds on how many arguments pattern `p` can consume in a sequence:
 * [*lo, *hi], *hi == SIZE_MAX when unbounded. Conservative -- the true widths
 * always lie inside -- so it only prunes spans that could never match. Lets
 * rep_units_match try just the feasible unit widths instead of every span. */
static void pat_width(Expr* p, size_t* lo, size_t* hi) {
    *lo = 1; *hi = 1;
    if (p->type != EXPR_FUNCTION || p->data.function.head->type != EXPR_SYMBOL) return;
    const char* h = p->data.function.head->data.symbol.name;
    size_t ac = p->data.function.arg_count;
    if (h == SYM_Pattern && ac == 2) { pat_width(p->data.function.args[1], lo, hi); return; }
    if ((h == SYM_Condition || h == SYM_PatternTest) && ac == 2) {
        pat_width(p->data.function.args[0], lo, hi); return;
    }
    if ((h == SYM_HoldPattern || h == SYM_Longest || h == SYM_Shortest) && ac >= 1) {
        pat_width(p->data.function.args[0], lo, hi); return;
    }
    if (h == SYM_Optional && ac >= 1) {
        pat_width(p->data.function.args[0], lo, hi); *lo = 0; return;
    }
    if (h == SYM_BlankSequence) { *lo = 1; *hi = SIZE_MAX; return; }
    if (h == SYM_BlankNullSequence) { *lo = 0; *hi = SIZE_MAX; return; }
    if (h == SYM_Repeated || h == SYM_RepeatedNull) {
        Expr* rp = NULL; int rmin = 0, rmax = -1;
        if (!is_repeated(p, &rp, &rmin, &rmax)) return;
        size_t ilo, ihi;
        pat_width(rp, &ilo, &ihi);
        *lo = ilo * (size_t)rmin;
        *hi = (rmax == -1 || ihi == SIZE_MAX) ? SIZE_MAX
              : (ihi > SIZE_MAX / ((size_t)rmax + 1) ? SIZE_MAX : ihi * (size_t)rmax);
        if (rmax != -1 && rmax == 0) *hi = 0;
        return;
    }
    if (h == SYM_PatternSequence || h == SYM_OrderlessPatternSequence) {
        size_t tlo = 0, thi = 0;
        for (size_t i = 0; i < ac; i++) {
            size_t a, b;
            pat_width(p->data.function.args[i], &a, &b);
            tlo += a;
            thi = (thi == SIZE_MAX || b == SIZE_MAX) ? SIZE_MAX : thi + b;
        }
        *lo = tlo; *hi = thi;
        return;
    }
    if (h == SYM_Alternatives && ac > 0) {
        size_t tlo = SIZE_MAX, thi = 0;
        for (size_t i = 0; i < ac; i++) {
            size_t a, b;
            pat_width(p->data.function.args[i], &a, &b);
            if (a < tlo) tlo = a;
            if (b > thi) thi = b;
        }
        *lo = tlo; *hi = thi;
        return;
    }
}

/* Does xs[0..n) split into r consecutive units, min_reps <= r <= max_reps
 * (max_reps -1 = unbounded), each unit matching the pattern sequence
 * unit[0..m)? This is the sequence semantics of a Repeated whose body can
 * span a variable number of arguments: {a,b,a,b} matches
 * {PatternSequence[a,b]..} (unit = a, b), {c,a,b} matches
 * {(PatternSequence[a,b] | c)..} and {1,2,1,2} matches {x__..} with x bound
 * to Sequence[1,2] (unit = the one body pattern). Bindings made by the units
 * accumulate in `env`, so a named sub-pattern must bind the same value in
 * every repetition, exactly as for single-argument repeats; on failure
 * everything this call bound is rolled back.
 *
 * When the unit is a PatternSequence of single-argument patterns (the common
 * PatternSequence[a, b] / PatternSequence[x_, _Integer]) its width is fixed
 * and the split is unique, so the units are checked in one linear pass.
 * Otherwise the split is searched depth-first over the width of each unit,
 * with an explicit stack (no recursion, so a long subject cannot exhaust the
 * C stack); each unit takes its first inner match at a given width. */
static bool rep_cont_call(RepCont k, void* kdata, MatchEnv* env) {
    if (!k) return true;
    size_t mark = env->count;
    if (k(env, kdata)) return true;
    env_rollback(env, mark);
    return false;
}

static bool rep_units_match(Expr** xs, size_t n, Expr** unit, size_t m,
                            int min_reps, int max_reps,
                            MatchEnv* env, Expr* pat_head,
                            RepCont k, void* kdata) {
    size_t saved = env->count;
    if (n == 0) {
        /* Zero units suffice when the bounds allow it; otherwise min_reps
         * units must each match nothing -- one empty match stands for all of
         * them, since repeating it binds nothing new. */
        if (min_reps == 0) return rep_cont_call(k, kdata, env);
        if (max_reps == 0) return false;
        if (match_args_internal(xs, 0, unit, m, env, NULL, pat_head, m, NULL)
            && rep_cont_call(k, kdata, env)) return true;
        env_rollback(env, saved);
        return false;
    }

    bool fixed = (m > 0);
    for (size_t i = 0; i < m && fixed; i++) fixed = pat_is_single_elem(unit[i]);
    if (fixed) {
        if (n % m != 0) return false;
        size_t reps = n / m;
        if (reps < (size_t)min_reps) return false;
        if (max_reps != -1 && reps > (size_t)max_reps) return false;
        for (size_t off = 0; off < n; off += m) {
            if (!match_args_internal(xs + off, m, unit, m, env, NULL, pat_head, m, NULL)) {
                env_rollback(env, saved);
                return false;
            }
        }
        if (rep_cont_call(k, kdata, env)) return true;
        env_rollback(env, saved);
        return false;
    }

    /* Depth-first search over unit widths. Frame d records where unit d
     * started, the width it took, and the env mark to roll back to when that
     * width is abandoned. A unit always consumes at least one argument here
     * (n > 0; a zero-width unit would repeat forever without progress), and
     * only widths inside the unit's static bounds are tried, which keeps a
     * bounded unit such as (PatternSequence[a, b] | c) linear in n. */
    size_t wlo = 0, whi = 0;
    for (size_t i = 0; i < m; i++) {
        size_t a, b;
        pat_width(unit[i], &a, &b);
        wlo += a;
        whi = (whi == SIZE_MAX || b == SIZE_MAX) ? SIZE_MAX : whi + b;
    }
    if (wlo < 1) wlo = 1;
    if (whi == 0) return false;   /* only zero-width units, and n > 0 */
    size_t cap = 16, depth = 0;
    size_t* f_off = malloc(cap * sizeof(size_t));
    size_t* f_w = malloc(cap * sizeof(size_t));
    size_t* f_env = malloc(cap * sizeof(size_t));
    size_t off = 0, next_w = wlo;
    bool ok = false;
    for (;;) {
        if (off == n) {
            /* A complete cover: accept it only if the continuation does too;
             * otherwise fall through to the dead-end step and re-split, so
             * f[a, a] /. f[x__..] can give x -> a (two units) and, on
             * backtracking, x -> Sequence[a, a] (one unit). */
            if (depth >= (size_t)min_reps && rep_cont_call(k, kdata, env)) { ok = true; break; }
        } else if (max_reps == -1 || depth < (size_t)max_reps) {
            bool pushed = false;
            size_t wmax = (whi < n - off) ? whi : n - off;
            for (size_t w = next_w; w <= wmax; w++) {
                size_t mark = env->count;
                if (match_args_internal(xs + off, w, unit, m, env, NULL, pat_head, m, NULL)) {
                    if (depth == cap) {
                        cap *= 2;
                        f_off = realloc(f_off, cap * sizeof(size_t));
                        f_w = realloc(f_w, cap * sizeof(size_t));
                        f_env = realloc(f_env, cap * sizeof(size_t));
                    }
                    f_off[depth] = off; f_w[depth] = w; f_env[depth] = mark;
                    depth++;
                    off += w;
                    next_w = wlo;
                    pushed = true;
                    break;
                }
                env_rollback(env, mark);
            }
            if (pushed) continue;
        }
        /* Dead end: abandon the last unit's width and try the next one. */
        if (depth == 0) break;
        depth--;
        env_rollback(env, f_env[depth]);
        off = f_off[depth];
        next_w = f_w[depth] + 1;
    }
    free(f_off); free(f_w); free(f_env);
    if (!ok) env_rollback(env, saved);
    return ok;
}

/* Split the body of a repeat into the unit pattern sequence rep_units_match
 * walks: the arguments of a PatternSequence, otherwise the body itself as a
 * one-pattern sequence (*slot holds it, so `unit` may point at it). Returns
 * false for a body that always consumes exactly one argument -- the caller's
 * per-argument loop handles that case faster. */
static bool rep_unit_split(Expr* rep_pat, Expr** slot, Expr*** unit, size_t* m) {
    if (pat_is_single_elem(rep_pat)) return false;
    if (rep_pat->type == EXPR_FUNCTION && rep_pat->data.function.head->type == EXPR_SYMBOL
        && rep_pat->data.function.head->data.symbol.name == SYM_PatternSequence) {
        *unit = rep_pat->data.function.args;
        *m = rep_pat->data.function.arg_count;
    } else {
        *slot = rep_pat;
        *unit = slot;
        *m = 1;
    }
    return true;
}

/* The commit step of a sequence element (x__, x___, p.., ...) that has taken
 * the k arguments `subset`: bind its name(s) to Sequence[subset] -- a name
 * already bound (nonlinear reuse, or the outer/inner names of x:(b___/;c))
 * must agree -- evaluate its element-level Condition guard, then match the
 * rest of the argument list. Bindings it adds are left for the caller to roll
 * back on failure. Also used as a RepCont, so a Repeated unit search can
 * resume when the rest fails. */
typedef struct {
    Expr** subset; size_t k;
    Expr* p_sym; Expr* extra_bind_sym; Expr* seq_cond;
    Expr** remainder; size_t n_rem;
    Expr** pats; size_t n_pats;
    Expr* condition; Expr* pat_head; size_t total_pats; ParentMatch* parent;
} SeqCommit;

static bool seq_commit(MatchEnv* env, void* data) {
    SeqCommit* c = (SeqCommit*)data;
    if (c->p_sym) {
        Expr* seq_val = expr_new_function(expr_new_symbol(SYM_Sequence), NULL, c->k);
        for (size_t i = 0; i < c->k; i++) seq_val->data.function.args[i] = expr_copy(c->subset[i]);
        bool bound_ok = true;
        Expr* existing = env_get(env, c->p_sym->data.symbol.name);
        if (existing) { if (!expr_eq(seq_val, existing)) bound_ok = false; }
        else env_set(env, c->p_sym->data.symbol.name, seq_val);
        if (bound_ok && c->extra_bind_sym) {
            Expr* ex2 = env_get(env, c->extra_bind_sym->data.symbol.name);
            if (ex2) { if (!expr_eq(seq_val, ex2)) bound_ok = false; }
            else env_set(env, c->extra_bind_sym->data.symbol.name, seq_val);
        }
        expr_free(seq_val);
        if (!bound_ok) return false;
    }
    /* Element-level Condition guard (block___ /; ...): evaluated once, now
     * that the sequence variable is bound. */
    if (c->seq_cond && !eval_seq_conditions(c->seq_cond, env)) return false;
    return match_args_internal(c->remainder, c->n_rem, c->pats, c->n_pats, env,
                               c->condition, c->pat_head, c->total_pats, c->parent);
}

/* OptionsPattern slot `opt_pat` taking exprs[0..take) (each an option rule
 * or a list of them) and the later pattern elements matching the rest; see
 * the OptionsPattern branch of match_args_internal. */
static bool opt_match_prefix(Expr** exprs, size_t take, size_t n_exprs, Expr** pats, size_t n_pats,
                             Expr* opt_pat, Expr* bind_sym, Expr* extra_bind_sym, Expr* seq_cond,
                             MatchEnv* env, Expr* condition, Expr* pat_head, size_t total_pats,
                             ParentMatch* parent) {
        size_t cap = (take ? take : 1), cnt = 0;
        Expr** collected = malloc(sizeof(Expr*) * cap);
        bool ok = true;
        /* Each argument is a rule or a list of rules, lists nested to any
         * depth ({{a -> 5}}), flattened in order, as the reference accepts. */
        size_t sp = 0, scap = 8;
        Expr** stack = malloc(sizeof(Expr*) * scap);
        for (size_t i = take; i-- > 0; ) {
            if (sp == scap) { scap *= 2; stack = realloc(stack, sizeof(Expr*) * scap); }
            stack[sp++] = exprs[i];
        }
        while (sp > 0 && ok) {
            Expr* a = stack[--sp];
            if (mo_is_option_rule(a)) {
                if (cnt == cap) { cap *= 2; collected = realloc(collected, sizeof(Expr*) * cap); }
                collected[cnt++] = a;
            } else if (mo_is_list(a)) {
                for (size_t j = a->data.function.arg_count; j-- > 0; ) {
                    if (sp == scap) { scap *= 2; stack = realloc(stack, sizeof(Expr*) * scap); }
                    stack[sp++] = a->data.function.args[j];
                }
            } else {
                ok = false;
            }
        }
        free(stack);
        if (ok) {
            /* OptionsPattern[{o -> v, ...}] (or a single rule) carries its own
             * defaults: append them AFTER the passed options, so OptionValue's
             * first-match lookup sees a passed value first and the default
             * otherwise. A symbol argument is handled below. */
            Expr* dflt = (opt_pat->data.function.arg_count >= 1)
                         ? opt_pat->data.function.args[0] : NULL;
            size_t nd = 0;
            if (dflt && mo_is_list(dflt)) nd = dflt->data.function.arg_count;
            else if (dflt && mo_is_option_rule(dflt)) nd = 1;
            Expr** items = malloc(sizeof(Expr*) * (cnt + nd ? cnt + nd : 1));
            for (size_t i = 0; i < cnt; i++) items[i] = expr_copy(collected[i]);
            for (size_t i = 0; i < nd; i++)
                items[cnt + i] = expr_copy(mo_is_list(dflt) ? dflt->data.function.args[i] : dflt);
            Expr* optlist = expr_new_function(expr_new_symbol(SYM_List), items, cnt + nd);
            free(items);
            size_t saved = env->count;
            env_set(env, "$OptionsPattern$", optlist);  /* env_set copies value */
            expr_free(optlist);
            /* If the pattern is OptionsPattern[s] (a target symbol), record s so
             * OptionValue resolves an unpassed option's DEFAULT against s's
             * Options, not the enclosing definition's own (A9). A bare
             * OptionsPattern[] leaves this unset -> the enclosing head is used. */
            if (opt_pat->data.function.arg_count >= 1
                && opt_pat->data.function.args[0]->type == EXPR_SYMBOL) {
                Expr* tgt = expr_copy(opt_pat->data.function.args[0]);
                env_set(env, "$OptionsPatternHead$", tgt);
                expr_free(tgt);
            }
            /* If the OptionsPattern was named (opts : OptionsPattern[]), also
             * bind that symbol to a Sequence of the matched options so the RHS
             * can splice them -- Sequence @@ ..., {opts}, or passing opts down
             * to another call. WL binds the name to the option sequence; an
             * empty match binds Sequence[] which vanishes on splice. Without
             * this only OptionValue (via $OptionsPattern$) works, not the
             * common opts-forwarding idiom. */
            /* The name binds the arguments as written -- {a -> 1}, b -> 2 stays
             * a list and a rule, not the flattened rules -- as in the reference
             * (f[x_, o:OptionsPattern[]] := {o}; f[1, {a->1}, b->2] is
             * {{a -> 1}, b -> 2}). */
            Expr* named[2] = { bind_sym, extra_bind_sym };
            for (int nb = 0; nb < 2; nb++) {
                if (!named[nb] || named[nb]->type != EXPR_SYMBOL) continue;
                Expr** sitems = malloc(sizeof(Expr*) * (take ? take : 1));
                for (size_t i = 0; i < take; i++) sitems[i] = expr_copy(exprs[i]);
                Expr* seq = expr_new_function(expr_new_symbol(SYM_Sequence), sitems, take);
                free(sitems);
                env_set(env, named[nb]->data.symbol.name, seq);
                expr_free(seq);
            }
            /* The rest of the pattern list matches the arguments after the prefix. */
            if ((!seq_cond || eval_seq_conditions(seq_cond, env)) &&
                match_args_internal(exprs + take, n_exprs - take, pats + 1, n_pats - 1, env,
                                    condition, pat_head, total_pats, parent)) {
                free(collected);
                return true;
            }
            env_rollback(env, saved);
        }
        free(collected);
        return false;
}

static bool match_args_internal(Expr** exprs, size_t n_exprs, Expr** pats, size_t n_pats, MatchEnv* env, Expr* condition, Expr* pat_head, size_t total_pats, ParentMatch* parent) {
    
    if (n_pats == 0) {
        if (n_exprs != 0) return false;

        if (condition && !eval_guard_true(condition, env)) return false;

        return call_parent(env, parent);
    }

    /* PatternSequence[q1, ..., qm] as a pattern element expands in place to the
     * m sub-patterns q1..qm; PatternSequence[] matches zero arguments.  Splice
     * it into the pattern list and recurse.  total_pats grows by the net number
     * of spliced elements so Optional positional defaults downstream keep a
     * consistent arity.  (A named x:PatternSequence[...] is not spliced here --
     * it needs the sequence bound to a name, handled separately.) */
    if (pats[0]->type == EXPR_FUNCTION && pats[0]->data.function.head->type == EXPR_SYMBOL
        && pats[0]->data.function.head->data.symbol.name == SYM_PatternSequence) {
        Expr* ps = pats[0];
        size_t m = ps->data.function.arg_count;
        size_t new_n = m + (n_pats - 1);
        Expr** newpats = malloc(sizeof(Expr*) * (new_n ? new_n : 1));
        for (size_t i = 0; i < m; i++) newpats[i] = ps->data.function.args[i];
        for (size_t i = 1; i < n_pats; i++) newpats[m + i - 1] = pats[i];
        bool r = match_args_internal(exprs, n_exprs, newpats, new_n, env, condition,
                                     pat_head, total_pats + m - 1, parent);
        free(newpats);
        return r;
    }

    /* Sequence alternatives: an unnamed Alternatives element with a branch
     * that can span other than exactly one argument (PatternSequence[a, b] | c,
     * __ | x_) is not a single-argument test -- each branch is spliced into
     * the pattern list in turn, like PatternSequence above, so the branch takes
     * as many arguments as it needs and the rest of the list continues after
     * it. {a, b, d} matches {PatternSequence[a, b] | c, d}. Plain alternatives
     * (a | b | _h) keep the cheaper single-argument path in match_internal. */
    if (pats[0]->type == EXPR_FUNCTION && pats[0]->data.function.head->type == EXPR_SYMBOL
        && pats[0]->data.function.head->data.symbol.name == SYM_Alternatives
        && !pat_is_single_elem(pats[0])) {
        Expr* alts = pats[0];
        Expr** newpats = malloc(sizeof(Expr*) * n_pats);
        for (size_t i = 1; i < n_pats; i++) newpats[i] = pats[i];
        size_t saved = env->count;
        bool r = false;
        for (size_t a = 0; a < alts->data.function.arg_count && !r; a++) {
            newpats[0] = alts->data.function.args[a];
            r = match_args_internal(exprs, n_exprs, newpats, n_pats, env, condition,
                                    pat_head, total_pats, parent);
            if (!r) env_rollback(env, saved);
        }
        free(newpats);
        return r;
    }

    /* OrderlessPatternSequence[q1, ..., qm] in first position: it occupies a
     * CONTIGUOUS run of the arguments starting here -- as any sequence
     * element does -- and the qi match the run's elements in any order
     * (f[2, 1, b] matches f[OrderlessPatternSequence[1, 2], ___], f[1, b, 2]
     * does not). Run lengths are tried shortest first; when the OPS is the
     * last element it must take everything that is left. Literal qi are
     * placed first, so each settles on one element before the variable ones
     * enumerate (one match, not two, for f[a, a] against OPS[x_, a]). */
    {
        Expr* ops_bind = NULL;
        Expr* ops = ops_unwrap(pats[0], &ops_bind);
        if (ops) {
            size_t nq = ops->data.function.arg_count;
            Expr** qs = malloc((nq ? nq : 1) * sizeof(Expr*));
            size_t w = 0;
            for (size_t i = 0; i < nq; i++)
                if (pat_is_literal(ops->data.function.args[i])) qs[w++] = ops->data.function.args[i];
            for (size_t i = 0; i < nq; i++)
                if (!pat_is_literal(ops->data.function.args[i])) qs[w++] = ops->data.function.args[i];
            size_t lo, hi;
            pat_width(ops, &lo, &hi);
            size_t lmin = lo, lmax = hi < n_exprs ? hi : n_exprs;
            if (n_pats == 1) { lmin = n_exprs; if (n_exprs < lo || n_exprs > hi) lmax = 0, lmin = 1; }
            bool r = false;
            bool* used = calloc(n_exprs ? n_exprs : 1, sizeof(bool));
            for (size_t L = lmin; L <= lmax && !r; L++) {
                r = ops_assign(exprs, L, used, qs, 0, nq, exprs + L, n_exprs - L,
                               pats, n_pats, env, condition, pat_head, total_pats,
                               parent, ops_bind);
            }
            free(used);
            free(qs);
            return r;
        }
    }

    /* Constraint ordering for Orderless AC-matching.  In an Orderless head a
     * leading variable-length sequence blank (__/___) forces the matcher to
     * enumerate every subset of the subject for that blank *before* a later,
     * highly-selective fixed element gets a chance to fail-fast -- making the
     * match exponential on large sums/products.  (Concretely, the trig
     * Pythagorean simplification rules `Cosh[x_]^2 - Sinh[x_]^2 + r___`
     * applied to a ~20-term expansion cost ~2 s per rule, and the enclosing
     * Simplify hung.)  Pattern-element order is semantically irrelevant for an
     * Orderless head, so stably move the variable-length sequence blanks to
     * the end; the fixed elements then prune first (O(n) instead of O(2^n)).
     * Skipped when any element carries an Optional positional default, whose
     * index would shift under reordering. */
    Expr* reorder_storage[MATCH_REORDER_CAP];
    if (n_pats >= 2 && pat_head && pat_head->type == EXPR_SYMBOL) {
        SymbolDef* hd = symtab_get_def(pat_head->data.symbol.name);
        if (hd && (hd->attributes & ATTR_ORDERLESS)) {
            bool any_opt = false, need = false, seen_seq = false;
            for (size_t i = 0; i < n_pats; i++) {
                if (pat_has_optional(pats[i])) { any_opt = true; break; }
                if (pat_is_variable_seq(pats[i])) seen_seq = true;
                else if (seen_seq) need = true;  /* fixed element after a blank */
            }
            if (!any_opt && need) {
                /* Fast path: reorder into the stack buffer for the common case.
                 * For a very wide pattern (> MATCH_REORDER_CAP elements) fall
                 * back to a heap buffer and recurse once so the fail-fast
                 * ordering still applies -- the recursive call sees fixed
                 * elements first, so `need` is false and it does not re-reorder.
                 * This keeps the O(2^n) blowup fix working at any width without
                 * an arbitrary element cap or a per-call malloc in the common
                 * case. */
                Expr** dst = reorder_storage;
                Expr** heap = NULL;
                if (n_pats > MATCH_REORDER_CAP) {
                    heap = malloc(n_pats * sizeof(Expr*));
                    dst = heap;
                }
                size_t w = 0;
                for (size_t i = 0; i < n_pats; i++)
                    if (!pat_is_variable_seq(pats[i])) dst[w++] = pats[i];
                for (size_t i = 0; i < n_pats; i++)
                    if (pat_is_variable_seq(pats[i])) dst[w++] = pats[i];
                if (heap) {
                    bool r = match_args_internal(exprs, n_exprs, heap, n_pats, env,
                                                 condition, pat_head, total_pats, parent);
                    free(heap);
                    return r;
                }
                pats = dst;
            }
        }
    }

    Expr* p = pats[0];

    bool is_optional = false;
    bool is_shortest = false;
    bool is_longest = false;
    Expr* opt_container = p;
    Expr* bind_sym = NULL;
    Expr* current_p = p;

    while (current_p->type == EXPR_FUNCTION && current_p->data.function.head->type == EXPR_SYMBOL && current_p->data.function.arg_count >= 1) {
        const char* head = current_p->data.function.head->data.symbol.name;
        if (head == SYM_Pattern && current_p->data.function.arg_count == 2) {
            bind_sym = current_p->data.function.args[0];
            current_p = current_p->data.function.args[1];
        } else if (head == SYM_Shortest) {
            is_shortest = true;
            current_p = current_p->data.function.args[0];
        } else if (head == SYM_Longest) {
            is_longest = true;
            is_shortest = false;
            current_p = current_p->data.function.args[0];
        } else if (head == SYM_Optional) {
            is_optional = true;
            opt_container = current_p;
            current_p = current_p->data.function.args[0];
        } else {
            break;
        }
    }
    
    Expr* opt_pat = current_p;

    /* __?test / ___?test: a PatternTest wrapping a sequence blank applies the
     * predicate to EVERY matched element.  Strip the PatternTest here (leaving
     * the bare sequence blank / named sequence blank for the machinery below)
     * and remember the test; the subset loop re-applies it per element.  A
     * PatternTest on a non-sequence pattern is deliberately left intact so
     * match_internal's single-element PatternTest handler runs instead. */
    Expr* seq_test = NULL;
    if (opt_pat->type == EXPR_FUNCTION && opt_pat->data.function.head->type == EXPR_SYMBOL
        && opt_pat->data.function.head->data.symbol.name == SYM_PatternTest
        && opt_pat->data.function.arg_count == 2) {
        Expr* pt_inner = opt_pat->data.function.args[0];
        Expr* tsym = NULL, *tpat = NULL;
        Expr* seqcheck = pt_inner;
        if (is_pattern(pt_inner, &tsym, &tpat)) seqcheck = tpat;
        Expr* dummy_h = NULL; int dummy_ml = 0;
        if (is_sequence_blank(seqcheck, &dummy_h, &dummy_ml)) {
            seq_test = opt_pat->data.function.args[1];
            opt_pat = pt_inner;
        }
    }

    /* block___ /; c   (also nested  block___ /; c1 /; c2,  and named
     * x:(block___ /; c)):  a Condition wrapping a SEQUENCE blank / Repeated.
     * match_internal's Condition handler can only gate a single matched element,
     * so a multi-element `x__ /; test` silently failed to match (only the
     * width-1 case slipped through, because a BlankSequence also matches one
     * element).  Strip the Condition layer(s) here -- leaving the bare (possibly
     * named) sequence blank / Repeated for the machinery below -- and remember
     * the guard chain in `seq_cond`; it is evaluated once, after the sequence is
     * bound, at each commit point.  Left intact when the inner pattern is NOT a
     * sequence blank / Repeated, so the well-tested single-element `x_ /; test`
     * path is unchanged.  A surviving inner name under an already-peeled outer
     * name (x:(block___ /; c) binds both x and block to the same sequence) is
     * remembered in `extra_bind_sym`. */
    Expr* seq_cond = NULL;
    Expr* extra_bind_sym = NULL;
    if (opt_pat->type == EXPR_FUNCTION && opt_pat->data.function.head->type == EXPR_SYMBOL
        && opt_pat->data.function.head->data.symbol.name == SYM_Condition
        && opt_pat->data.function.arg_count == 2) {
        Expr* base = opt_pat;
        while (base->type == EXPR_FUNCTION && base->data.function.head->type == EXPR_SYMBOL
               && base->data.function.head->data.symbol.name == SYM_Condition
               && base->data.function.arg_count == 2) {
            base = base->data.function.args[0];
        }
        Expr* iname = NULL, *ipat = NULL;
        Expr* seqcheck = base;
        if (is_pattern(base, &iname, &ipat)) seqcheck = ipat;
        Expr* dh = NULL; int dml = 0;
        Expr* drp = NULL; int drmin = 0, drmax = -1;
        bool opts_pat = seqcheck->type == EXPR_FUNCTION
            && seqcheck->data.function.head->type == EXPR_SYMBOL
            && seqcheck->data.function.head->data.symbol.name == SYM_OptionsPattern;
        if (opts_pat) {
            /* (opts:OptionsPattern[]) /; test: the guard is checked once the
             * options are bound (the OptionsPattern block below), as for a
             * sequence blank. It used to never match: FeynCalc's
             * Isolate[ex_, (opts:OptionsPattern[]) /; opts =!= {}] was dead. */
            seq_cond = opt_pat;
            if (iname && !bind_sym) bind_sym = iname;
            else if (iname) extra_bind_sym = iname;
            opt_pat = seqcheck;
        } else if (is_sequence_blank(seqcheck, &dh, &dml)
            || is_repeated(seqcheck, &drp, &drmin, &drmax)) {
            seq_cond = opt_pat;              /* guard chain, evaluated at commit */
            if (iname && bind_sym) {         /* x:(block___ /; c): two names, one seq */
                extra_bind_sym = iname;
                opt_pat = seqcheck;          /* bare blank/Repeated; outer name is p_sym */
            } else {
                opt_pat = base;              /* bare, or single inner-named Pattern[b,..] */
            }
        }
    }

    /* Named PatternSequence: x:PatternSequence[q1..qm] binds x to the Sequence
     * of the arguments consumed by q1..qm.  (The unnamed form is spliced away
     * before the reorder block; this handles the bound form, which the peeling
     * loop leaves as bind_sym=x, opt_pat=PatternSequence[...].)  We enumerate
     * the prefix length k, require exprs[0..k) to fully match the sub-patterns,
     * bind x to that prefix, and recurse on the remainder.  Note: the prefix
     * sub-match takes the first internal assignment of q1..qm; it does not
     * backtrack the internal split against the outer remainder, which only
     * matters for PatternSequence bodies that themselves contain sequence
     * blanks. */
    /* The same applies to a named sequence alternative, x:(PatternSequence[a,
     * b] | c): the body is matched as a one-pattern sequence (the unnamed
     * splice above then picks the branch) and x binds what it consumed. */
    Expr* named_alt_slot = opt_pat;
    bool named_seq_alt = opt_pat->type == EXPR_FUNCTION
        && opt_pat->data.function.head->type == EXPR_SYMBOL
        && opt_pat->data.function.head->data.symbol.name == SYM_Alternatives
        && !pat_is_single_elem(opt_pat);
    if (bind_sym && bind_sym->type == EXPR_SYMBOL
        && opt_pat->type == EXPR_FUNCTION && opt_pat->data.function.head->type == EXPR_SYMBOL
        && (opt_pat->data.function.head->data.symbol.name == SYM_PatternSequence
            || named_seq_alt)) {
        Expr** qs = named_seq_alt ? &named_alt_slot : opt_pat->data.function.args;
        size_t m = named_seq_alt ? 1 : opt_pat->data.function.arg_count;
        for (size_t k = 0; k <= n_exprs; k++) {
            size_t saved = env->count;
            /* Full boolean sub-match of exprs[0..k) against q1..qm (parent NULL
             * => success when both sub-lists are exhausted, without firing the
             * outer callback). */
            if (match_args_internal(exprs, k, qs, m, env, NULL, pat_head, m, NULL)) {
                Expr* seqv = expr_new_function(expr_new_symbol(SYM_Sequence), NULL, k);
                for (size_t i = 0; i < k; i++)
                    seqv->data.function.args[i] = expr_copy(exprs[i]);
                Expr* existing = env_get(env, bind_sym->data.symbol.name);
                bool okbind = true;
                if (existing) { if (!expr_eq(seqv, existing)) okbind = false; }
                else env_set(env, bind_sym->data.symbol.name, seqv);
                expr_free(seqv);
                if (okbind && match_args_internal(exprs + k, n_exprs - k, pats + 1,
                        n_pats - 1, env, condition, pat_head, total_pats, parent))
                    return true;
            }
            env_rollback(env, saved);
        }
        return false;
    }

    /* OptionsPattern[]: consume the remaining option arguments (each a
     * Rule/RuleDelayed with a symbol/string name, or a List of such rules),
     * flatten them, and bind the collected options under the reserved key
     * "$OptionsPattern$" so OptionValue can resolve them when the rule fires.
     * It is the least-specific, variable-arity slot, so it greedily takes the
     * whole trailing option run; any remaining non-option argument makes this
     * rule fail to match (e.g. f[OptionsPattern[]] does not match f[x]). */
    if (opt_pat->type == EXPR_FUNCTION
        && opt_pat->data.function.head->type == EXPR_SYMBOL
        && opt_pat->data.function.head->data.symbol.name == SYM_OptionsPattern) {
        /* Not necessarily the last pattern element: s[a__, OptionsPattern[]]
         * on an Orderless head is sorted to s[OptionsPattern[], a__], and
         * t[OptionsPattern[], a_] matches t[x] in the reference. So the slot
         * takes a PREFIX of option arguments -- the longest first, down to
         * none -- and the rest goes to the later pattern elements. As the last
         * element only the whole run can succeed, as before. */
        size_t run = 0;                    /* leading arguments that are options */
        while (run < n_exprs && (mo_is_option_rule(exprs[run]) || mo_is_list(exprs[run]))) run++;
        for (size_t take = run + 1; take-- > 0; ) {
            if (n_pats == 1 && take != n_exprs) continue;   /* last slot: everything */
            if (opt_match_prefix(exprs, take, n_exprs, pats, n_pats, opt_pat, bind_sym,
                                 extra_bind_sym, seq_cond, env, condition, pat_head,
                                 total_pats, parent))
                return true;
        }
        return false;
    }

    if (is_optional && is_shortest) {
        size_t saved_env_count = env->count;
        Expr* def_val = NULL;
        if (opt_container->data.function.arg_count == 2) {
            def_val = expr_copy(opt_container->data.function.args[1]);
        } else {
            def_val = get_default_value(pat_head, total_pats - n_pats + 1, total_pats);
        }

        if (def_val) {
            Expr* p_sym = bind_sym;
            Expr* inner_p = opt_pat;
            if (!p_sym) is_pattern(opt_pat, &p_sym, &inner_p);
            bool consistent = true;
            if (p_sym && p_sym->type == EXPR_SYMBOL) {
                Expr* existing = env_get(env, p_sym->data.symbol.name);
                if (existing) {
                    if (!expr_eq(existing, def_val)) consistent = false;
                } else {
                    env_set(env, p_sym->data.symbol.name, def_val);
                }
            }
            expr_free(def_val);


            if (consistent && match_args_internal(exprs, n_exprs, pats + 1, n_pats - 1, env, condition, pat_head, total_pats, parent)) {
                return true;
            }
            env_rollback(env, saved_env_count);
        }
    }

    Expr* inner_p = opt_pat;
    Expr* p_sym = bind_sym;
    int min_len = 0;
    int max_len = -1;
    Expr* b_head = NULL;
    Expr* rep_pat = NULL;

    bool is_seq = false;
    bool is_rep = false;
    if (!p_sym && is_pattern(opt_pat, &p_sym, &inner_p)) {
        is_seq = is_sequence_blank(inner_p, &b_head, &min_len);
        if (!is_seq) is_rep = is_repeated_ex(inner_p, &rep_pat, &min_len, &max_len, true);
    } else {
        is_seq = is_sequence_blank(inner_p, &b_head, &min_len);
        if (!is_seq) is_rep = is_repeated_ex(inner_p, &rep_pat, &min_len, &max_len, true);
    }


    bool is_orderless = false;
    bool is_flat = false;
    if (pat_head && pat_head->type == EXPR_SYMBOL) {
        SymbolDef* def = symtab_get_def(pat_head->data.symbol.name);
        if (def) {
            if (def->attributes & ATTR_FLAT) is_flat = true;
            if (def->attributes & ATTR_ORDERLESS) is_orderless = true;
        }
    }

    size_t min_k, max_k;
    if (is_seq) {
        min_k = min_len;
        max_k = n_exprs;
        /* The LAST pattern element must consume every remaining argument --
         * nothing follows it to take the rest -- so only k = n_exprs can
         * succeed. Fixing it here avoids enumerating all C(n_exprs, k) subsets
         * for each k < n_exprs under an Orderless head: that is the 2^n blowup
         * that made e.g. Plus[x__, y__] over a wide sum hang, because the
         * trailing y__ tried every subset before taking all of them. */
        if (n_pats == 1) {
            if ((size_t)min_len <= n_exprs) { min_k = max_k = n_exprs; }
            else { min_k = 1; max_k = 0; }   /* not enough args for this blank */
        }
    } else if (is_rep && !pat_is_single_elem(rep_pat)) {
        /* A variable-width body (PatternSequence[...].., (p | PatternSequence[
         * ...]).., x__..): the repeat bounds count UNITS, not arguments, and a
         * unit's width is only known while matching, so offer every span and
         * let rep_units_match decide. The last pattern element must consume
         * everything that is left. */
        min_k = 0;
        max_k = n_exprs;
        if (n_pats == 1) min_k = n_exprs;
    } else if (is_rep) {
        min_k = min_len;
        max_k = (max_len == -1 || max_len > (int)n_exprs) ? n_exprs : (size_t)max_len;
        if (n_pats == 1) {
            /* Same "last element consumes all" constraint, within the repeat
             * bounds; k = n_exprs is the only candidate. */
            if ((size_t)min_len <= n_exprs
                && (max_len == -1 || (size_t)max_len >= n_exprs)) {
                min_k = max_k = n_exprs;
            } else { min_k = 1; max_k = 0; }
        }
    } else if (is_flat) {
        min_k = 1;
        max_k = n_exprs;
        /* Pruning: if inner_p is a specific Function pattern whose head
         * is neither pat_head's symbol nor a pattern-construct head that
         * may match anything, then no constructed
         * Function[pat_head, subset...] (k > 1) can ever match it -- it
         * would fail at the head check. The only exception is when
         * inner_p's head has OneIdentity AND has exactly one pattern arg
         * (e.g., Plus[a_]), in which case ONEIDENTITY lets the
         * constructed multi-factor value be absorbed into that single
         * arg. Without this guard the orderless+Flat enumeration is
         * exponential in the factor count even though k > 1 cannot
         * succeed -- dominates Simplify on multi-factor trig products. */
        if (inner_p && inner_p->type == EXPR_FUNCTION
            && inner_p->data.function.head
            && inner_p->data.function.head->type == EXPR_SYMBOL
            && pat_head && pat_head->type == EXPR_SYMBOL) {
            const char* inner_hn = inner_p->data.function.head->data.symbol.name;
            if (inner_hn != pat_head->data.symbol.name
                && inner_hn != SYM_Blank && inner_hn != SYM_BlankSequence
                && inner_hn != SYM_BlankNullSequence && inner_hn != SYM_Pattern
                && inner_hn != SYM_HoldPattern
                && inner_hn != SYM_Alternatives && inner_hn != SYM_Condition
                && inner_hn != SYM_PatternTest && inner_hn != SYM_Optional
                && inner_hn != SYM_Shortest && inner_hn != SYM_Longest
                && inner_hn != SYM_Repeated && inner_hn != SYM_RepeatedNull) {
                SymbolDef* idef = symtab_get_def(inner_hn);
                bool has_oid = idef && (idef->attributes & ATTR_ONEIDENTITY);
                /* OneIdentity lets a group be absorbed by a one-argument
                 * pattern, or through defaults -- Power[y_, n_.] takes the
                 * group 2 x as y with n -> 1 -- so keep the group spans then. */
                bool oid_can_absorb = has_oid
                    && inner_p->data.function.arg_count == 1;
                for (size_t ai = 0; has_oid && !oid_can_absorb
                                    && ai < inner_p->data.function.arg_count; ai++)
                    oid_can_absorb = pat_has_optional(inner_p->data.function.args[ai]);
                /* Through defaults, the group is absorbed by what is left once
                 * every optional argument takes its default: nothing, or one
                 * argument that may itself match pat_head[...].  Two or more
                 * required arguments keep inner_p's own head, and a required
                 * call on a fixed non-OneIdentity head (Cosh[a_] in
                 * c_. Cosh[a_]) cannot equal the group -- without this a sum of
                 * coefficient-times-Cosh terms enumerated every subset. */
                if (oid_can_absorb && inner_p->data.function.arg_count > 1) {
                    size_t n_req = 0;
                    Expr* req = NULL;
                    for (size_t ai = 0; ai < inner_p->data.function.arg_count; ai++)
                        if (!pat_has_optional(inner_p->data.function.args[ai])) {
                            n_req++;
                            req = inner_p->data.function.args[ai];
                        }
                    if (n_req >= 2 || (req && flat_elem_rejects_group(req, pat_head->data.symbol.name)))
                        oid_can_absorb = false;
                }
                if (!oid_can_absorb) {
                    max_k = (n_exprs > 0) ? 1 : 0;
                }
            }
        } else if (inner_p && inner_p->type == EXPR_FUNCTION &&
                   inner_p->data.function.head &&
                   inner_p->data.function.head->type != EXPR_SYMBOL) {
            /* A call with a pattern (or compound) head, g_[y__] or _[x_],
             * matches one operand, never a group pat_head[subset]: in the
             * reference ReplaceList[ff[a[1], b[2], c[3]], ff[g_[y__], c[3]]
             * :> ...] (ff Flat, Orderless) is {}, and an UpValue
             * q[z_, x_] f_[y__] never tries f -> Times over the other
             * factors (FeynCalc's DiracTrace PairContract rules). */
            max_k = (n_exprs > 0) ? 1 : 0;
        } else if (inner_p && inner_p->type != EXPR_FUNCTION) {
            /* A literal atom (integer/real/symbol/string/bigint) as a
             * pattern element can only match a single subject element that
             * is structurally equal to it: a constructed pat_head[subset]
             * with k > 1 is a Function and can never equal an atom, so the
             * Flat subset enumeration (max_k = n_exprs) is pure wasted work,
             * exponential in the argument count. Cap it at 1. This is the
             * atomic analogue of the Function-pattern pruning above; e.g.
             * the Pythagorean rule `1 - Cos[x_]^2 + r___`, whose literal `1`
             * element otherwise drove a 2^n subset scan over a large sum. */
            max_k = (n_exprs > 0) ? 1 : 0;
        }
    } else {
        min_k = 1;
        max_k = (n_exprs > 0) ? 1 : 0;
    }

    /* Anchor screen (ordered Flat head only).  When this element can take a
     * variable span and the NEXT element is a screenable anchor (it consumes
     * exactly the one operand that follows this element's span, and its
     * standalone match is necessary -- see flat_elem_screenable), a span k is
     * viable only if exprs[k] exists and passes the anchor on its own.  Skipping
     * the other spans removes the group build + downstream match for each, so
     * `And[x_, Not[x_]]` against a long conjunction costs one cheap operand
     * test per span instead of an And[...] group per span.  Only spans that
     * would fail anyway are skipped, in the same order, so the first match --
     * and therefore the result -- is unchanged. */
    Expr* anchor_next = NULL;
    const char* anchor_head = NULL;
    if (is_flat && !is_orderless && n_pats >= 2 && min_k < max_k &&
        pat_head && pat_head->type == EXPR_SYMBOL &&
        flat_elem_screenable(pats[1], pat_head->data.symbol.name)) {
        anchor_next = pats[1];
        anchor_head = flat_anchor_head(anchor_next);
        if (max_k >= n_exprs) {
            if (n_exprs == 0) { min_k = 1; max_k = 0; }  /* no operand for the anchor */
            else max_k = n_exprs - 1;
        }
    }

    if (min_k <= max_k) {
        for (size_t step = 0; step <= max_k - min_k; step++) {
            size_t k = is_longest ? (max_k - step) : (min_k + step);
            if (k > n_exprs) continue;
            if (anchor_next && !flat_anchor_accepts(exprs[k], anchor_next, anchor_head)) continue;

            /* comb (the index combination) is only consumed by extract_subset
             * and advanced by next_combination -- both Orderless-only below.
             * For an ordered head the consumed run is always the prefix, so skip
             * allocating and filling it entirely. */
            int* comb = NULL;
            if (k > 0 && is_orderless) {
                comb = malloc(k * sizeof(int));
                for (int i = 0; i < (int)k; i++) comb[i] = i;
            }

        /* The matched elements of a plain, unnamed, untyped sequence blank
         * (`__`/`___` with no bound name, head restriction, or per-element
         * test, and not a Repeated) are never read -- only the remainder is
         * recursed into. Skip allocating and filling `subset` in that case; for
         * a leading `___` scanned across a long list this halves the per-attempt
         * memory traffic (e.g. the duplicate-search `{___, x_, ___, x_, ___}`). */
        bool need_subset = !(is_seq && !b_head && !seq_test && !p_sym);
        /* A literal element (no pattern construct, no name) of an Orderless
         * head can only take an argument equal to itself, and every such
         * argument leaves the same bindings and an equal remainder. So it is
         * matched against the FIRST equal argument only:
         * ReplaceList[o[a, a], o[a, x_] :> x] is {a}, not {a, a}.
         * Unequal candidates are skipped without building the continuation. */
        bool lit_elem = is_orderless && k == 1 && !is_seq && !is_rep && !p_sym
                        && pat_is_literal(inner_p);
        bool lit_done = false;
        do {
            /* For an ordered (non-Orderless) head the do-while runs exactly once
             * with comb = [0..k-1] (next_combination is gated on is_orderless
             * below): the consumed elements are the contiguous prefix
             * exprs[0..k-1] and the remainder is the contiguous suffix exprs[k..].
             * Borrow both as slices of the caller's array instead of
             * malloc+extract_subset -- this removes the O(n) remainder copy paid
             * at EVERY pattern element, which (once atomic re-eval became cheap)
             * dominated `list //. {a_, b_, r___}`. The borrowed arrays are
             * read-only here (type checks, and expr_copy into seq_val/matched_val)
             * and in the recursive callee, so sharing is safe. `_owned` holds the
             * heap copies to free; they stay NULL in the borrow case, so the
             * free sites below are no-ops. Only the Orderless path still copies. */
            Expr** subset = NULL;
            Expr** remainder = NULL;
            Expr** subset_owned = NULL;
            Expr** remainder_owned = NULL;
            if (!is_orderless) {
                subset = exprs;                 /* only [0..k-1] are read */
                remainder = exprs + k;          /* contiguous suffix */
            } else {
                if (k > 0 && need_subset) subset_owned = malloc(k * sizeof(Expr*));
                if (n_exprs - k > 0) remainder_owned = malloc((n_exprs - k) * sizeof(Expr*));
                extract_subset(exprs, n_exprs, comb, (int)k, subset_owned, remainder_owned);
                subset = subset_owned;
                remainder = remainder_owned;
            }

            size_t saved_env = env->count;

            if (is_seq || is_rep) {
                bool type_ok = true;
                bool rep_variable = false;
                Expr* rep_slot = NULL; Expr** rep_unit = NULL; size_t rep_um = 0;
                if (p_sym) {
                    /* Already-bound sequence name of a different length: the
                     * Sequence built below can never equal it -- reject before
                     * the per-element tests and the O(k) build. */
                    Expr* ex = env_get(env, p_sym->data.symbol.name);
                    if (ex && ex->type == EXPR_FUNCTION && ex->data.function.head &&
                        ex->data.function.head->type == EXPR_SYMBOL &&
                        ex->data.function.head->data.symbol.name == SYM_Sequence &&
                        ex->data.function.arg_count != k)
                        type_ok = false;
                }
                if (!type_ok) {
                    /* fall through to the rollback below */
                } else if (is_seq) {
                    if (b_head) {
                        for (size_t i = 0; i < k; i++) {
                            if (!blank_head_matches(subset[i], b_head)) { type_ok = false; break; }
                        }
                    }
                    /* __?test / ___?test: the PatternTest applies to every
                     * matched element of the sequence. */
                    if (type_ok && seq_test) {
                        for (size_t i = 0; i < k; i++) {
                            if (!apply_pattern_test(subset[i], seq_test)) { type_ok = false; break; }
                        }
                    }
                } else if (is_rep) {
                    if (rep_unit_split(rep_pat, &rep_slot, &rep_unit, &rep_um)) {
                        rep_variable = true;    /* matched below, with a continuation */
                    } else {
                        for (size_t i = 0; i < k; i++) {
                            if (!match_internal(subset[i], rep_pat, env, NULL)) { type_ok = false; break; }
                        }
                    }
                }

                if (type_ok) {
                    SeqCommit sc = { subset, k, p_sym, extra_bind_sym, seq_cond,
                                     remainder, n_exprs - k, pats + 1, n_pats - 1,
                                     condition, pat_head, total_pats, parent };
                    /* A variable-width repeat hands the commit step to the unit
                     * search as its continuation, so a failure after it (or
                     * ReplaceList asking for more) re-splits the units. */
                    bool done = rep_variable
                        ? rep_units_match(subset, k, rep_unit, rep_um, min_len, max_len, env,
                                          is_orderless ? NULL : pat_head, seq_commit, &sc)
                        : seq_commit(env, &sc);
                    if (done) {
                        if (subset_owned) free(subset_owned);
                        if (remainder_owned) free(remainder_owned);
                        if (comb) free(comb);
                        return true;
                    }
                }
            } else {
                Expr* matched_val = NULL;
                if (k > 1 && p_sym && !flat_group_may_equal(
                        env_get(env, p_sym->data.symbol.name), pat_head, subset, k)) {
                    /* Already-bound name (x_ ... x_): the group head[subset]
                     * cannot equal the binding, so bind_ok below would fail.
                     * Skip it without building the k-operand group -- for a
                     * repeated name in a Flat head this turns the per-split
                     * O(k) build into an O(1) reject. */
                } else if (k == 1) {
                    matched_val = expr_copy(subset[0]);
                } else if (k > 1) {
                    matched_val = expr_new_function(expr_copy(pat_head), NULL, k);
                    for (size_t i = 0; i < k; i++) {
                        matched_val->data.function.args[i] = expr_copy(subset[i]);
                    }
                }

                if (matched_val && lit_elem && !expr_eq(matched_val, inner_p)) {
                    expr_free(matched_val);     /* cannot equal the literal */
                    matched_val = NULL;
                } else if (lit_elem) {
                    lit_done = true;            /* the one candidate worth trying */
                }

                if (matched_val) {
                    size_t saved_env_inner = env->count;
                    /* Bind the element's name first, then match the inner pattern
                     * with the REST OF THIS LEVEL as its continuation.  Matching
                     * the inner pattern with no continuation committed to its first
                     * internal assignment: a nested Orderless Times could never
                     * re-split when the remaining arguments or the rule's /;
                     * condition then failed, so IntegrateTable[a_ f_, x_] /;
                     * FreeQ[a, x] missed z x^k whenever x^k sorted first (and
                     * MatchQ[g[z x^k], g[a_ f_] /; FreeQ[a, x]] was False).  With
                     * the continuation, a failure downstream backtracks into the
                     * nested enumeration, as in Mathematica; a pattern with a single
                     * inner assignment behaves exactly as before. */
                    bool bind_ok = true;
                    if (p_sym) {
                        Expr* existing = env_get(env, p_sym->data.symbol.name);
                        if (existing) bind_ok = expr_eq(matched_val, existing);
                        else env_set(env, p_sym->data.symbol.name, matched_val);
                    }
                    if (bind_ok) {
                        ParentMatch cont = { remainder, n_exprs - k, pats + 1, n_pats - 1,
                                             condition, pat_head, total_pats, parent };
                        if (match_internal(matched_val, inner_p, env, &cont)) {
                            expr_free(matched_val); { if (subset_owned) free(subset_owned); if (remainder_owned) free(remainder_owned); if (comb) free(comb); return true; }
                        }
                    }
                    env_rollback(env, saved_env_inner);
                    expr_free(matched_val);
                }

            }

            env_rollback(env, saved_env);
            if (subset_owned) free(subset_owned);
            if (remainder_owned) free(remainder_owned);

        } while (is_orderless && !lit_done && next_combination(comb, n_exprs, (int)k));
        
        if (comb) free(comb);
        }
    }

    if (is_optional && !is_shortest) {
        size_t saved_env_count = env->count;
        Expr* def_val = NULL;
        if (opt_container->data.function.arg_count == 2) {
            def_val = expr_copy(opt_container->data.function.args[1]);
        } else {
            def_val = get_default_value(pat_head, total_pats - n_pats + 1, total_pats);
        }

        if (def_val) {
            Expr* p_sym = bind_sym;
            Expr* inner_p = opt_pat;
            if (!p_sym) is_pattern(opt_pat, &p_sym, &inner_p);
            bool consistent = true;
            if (p_sym && p_sym->type == EXPR_SYMBOL) {
                Expr* existing = env_get(env, p_sym->data.symbol.name);
                if (existing) {
                    /* Don't overwrite an existing binding for this name; the
                     * Optional fallback can only commit to the default value if
                     * it is consistent with what is already bound. Otherwise the
                     * pattern would produce inconsistent bindings such as
                     * `a_. Sin[x_]^2 + a_. Cos[x_]^2` matching `Sin[x]^2 -
                     * Cos[x]^2` with one a=1 and another a=-1. */
                    if (!expr_eq(existing, def_val)) consistent = false;
                } else {
                    env_set(env, p_sym->data.symbol.name, def_val);
                }
            }
            expr_free(def_val);


            if (consistent && match_args_internal(exprs, n_exprs, pats + 1, n_pats - 1, env, condition, pat_head, total_pats, parent)) {
                return true;
            }
            env_rollback(env, saved_env_count);
        }
    }

    return false;
}

/* Collect the distinct symbol names occurring in `e` into out[0..*n) (bounded by
 * cap). Used to find which pattern-binding values carry a symbol that a scoping
 * construct in the RHS could capture (see rb_rec / scoping_capture_avoid). */
static void rb_collect_symbol_names(Expr* e, const char** out, size_t* n, size_t cap) {
    if (!e || *n >= cap) return;
    if (e->type == EXPR_SYMBOL) {
        /* Symbol names are interned (see intern_symbol in symtab.c), so the
         * canonical pointer IS the identity -- no strcmp fallback needed. */
        for (size_t i = 0; i < *n; i++)
            if (out[i] == e->data.symbol.name) return;
        out[(*n)++] = e->data.symbol.name;
        return;
    }
    if (e->type == EXPR_FUNCTION) {
        rb_collect_symbol_names(e->data.function.head, out, n, cap);
        for (size_t i = 0; i < e->data.function.arg_count; i++)
            rb_collect_symbol_names(e->data.function.args[i], out, n, cap);
    }
}

/* The capture-avoidance "danger set" of one replace_bindings call: the names a
 * scoping construct in the RHS could capture if a binding value carrying them
 * were substituted into its body (MATHILDA_DIVERGENCES A11).
 *
 * It is built LAZILY and at most ONCE per top-level call.  Building it walks
 * every binding value, and the overwhelming majority of rule applications have
 * no scoping construct in the RHS at all, so the walk is pure cost there; the
 * set is materialised only when rb_rec actually meets a Module/With/Block/
 * Function/Table-style binder.  (Before this, it was rebuilt on every NODE of
 * the RHS -- rb_rec recursed into arguments through the public entry point --
 * which is what made the small Charlwood rows double in cost when the A11 fix
 * landed in v0.212.) */
#define RB_DANGER_CAP 256
typedef struct {
    MatchEnv* env;
    const char* buf[RB_DANGER_CAP];
    size_t n;
    bool built;
} RbDanger;

static void rb_danger_build(RbDanger* d) {
    if (d->built) return;
    d->built = true;
    MatchEnv* env = d->env;
    if (!env) return;
    for (size_t i = 0; i < env->count && d->n < RB_DANGER_CAP; i++) {
        /* Free symbols of the value -- a local of that name would CAPTURE it. */
        rb_collect_symbol_names(env->values[i], d->buf, &d->n, RB_DANGER_CAP);
        /* The bound variable's own name -- a local of that name SHADOWS it, so
         * the binding must not be substituted into the (renamed) body:
         * sh[e_] := Module[{e = 1}, e + 1]; sh[99] must give 2, not 100. */
        if (d->n < RB_DANGER_CAP && env->symbols[i]) {
            bool seen = false;
            for (size_t k = 0; k < d->n; k++)
                if (d->buf[k] == env->symbols[i]) { seen = true; break; }
            if (!seen) d->buf[d->n++] = env->symbols[i];
        }
    }
}

static Expr* rb_rec(Expr* expr, MatchEnv* env, RbDanger* danger);

Expr* replace_bindings(Expr* expr, MatchEnv* env) {
    if (!expr) return NULL;
    /* No bindings: nothing to substitute, splice or rename, so the result is
     * the expression itself.  A literal rule (res[14] = <135135-term sum>)
     * otherwise walked its whole right-hand side on every application, only
     * to share every node back (issue #112: 20 ms per read of res[12]). */
    if (!env || env->count == 0) return expr_copy(expr);
    RbDanger danger;
    danger.env = env; danger.n = 0; danger.built = false;
    return rb_rec(expr, env, &danger);
}

static Expr* rb_rec(Expr* expr, MatchEnv* env, RbDanger* danger) {
    if (!expr) return NULL;

    if (expr->type == EXPR_SYMBOL) {
        Expr* bound = env_get(env, expr->data.symbol.name);
        if (bound) {
            return expr_copy(bound);
        }
    }

    if (expr->type == EXPR_FUNCTION) {
        /* Capture-avoiding substitution: if this is a scoping construct whose
         * bound locals collide with a name in a value we are about to inject,
         * alpha-rename those locals first, then substitute into the safe copy. */
        if (expr_is_binding_scope(expr)) {
            rb_danger_build(danger);                 /* only now is the set needed */
            if (danger->n > 0) {
                Expr* safe = scoping_capture_avoid(expr, danger->buf, danger->n);
                if (safe) {
                    Expr* result = rb_rec(safe, env, danger);
                    expr_free(safe);
                    return result;
                }
            }
        }

        Expr* new_head = rb_rec(expr->data.function.head, env, danger);
        
        bool skip_flattening = false;
        if (new_head->type == EXPR_SYMBOL && (
            new_head->data.symbol.name == SYM_Set ||
            new_head->data.symbol.name == SYM_SetDelayed ||
            new_head->data.symbol.name == SYM_Rule ||
            new_head->data.symbol.name == SYM_RuleDelayed)) {
            skip_flattening = true;
        }

        /* A Sequence is spliced only where a pattern variable's value is
         * substituted (f[x__] := HoldComplete[x] gives HoldComplete[1, 2]), as
         * in the reference. A literal Sequence written in the replacement is
         * left to the evaluator, which removes it except where the head keeps
         * it -- Unevaluated[Sequence[]] and HoldComplete[Sequence[], n] stay.
         * Splicing it here made Unevaluated[], so a Table whose body returns
         * Unevaluated[Sequence[]] gave {Unevaluated[], ...} instead of {}. */
        /* Small arities (nearly every node) use stack buffers: three heap
         * allocations per node were most of a rule application's cost
         * (f[x_, y___] := Module[...] applied in a loop, issue #112). */
        enum { RB_STACK_N = 8 };
        Expr* temp_stack[RB_STACK_N];
        bool spliced_stack[RB_STACK_N];
        size_t an = expr->data.function.arg_count;
        Expr** temp_args = an <= RB_STACK_N ? temp_stack : malloc(sizeof(Expr*) * an);
        bool* spliced = an <= RB_STACK_N ? spliced_stack : malloc(sizeof(bool) * (an ? an : 1));
        size_t new_count = 0;
        for (size_t i = 0; i < expr->data.function.arg_count; i++) {
            /* rb_rec, not replace_bindings: the danger set belongs to the whole
             * substitution and must not be rebuilt once per argument node. */
            Expr* a = expr->data.function.args[i];
            temp_args[i] = rb_rec(a, env, danger);
            spliced[i] = !skip_flattening && a->type == EXPR_SYMBOL &&
                env_get(env, a->data.symbol.name) != NULL &&
                temp_args[i]->type == EXPR_FUNCTION &&
                temp_args[i]->data.function.head->type == EXPR_SYMBOL &&
                temp_args[i]->data.function.head->data.symbol.name == SYM_Sequence;
            new_count += spliced[i] ? temp_args[i]->data.function.arg_count : 1;
        }

        /* Structural sharing: a subtree no binding reached comes back as the
         * very nodes it started from (leaves are returned by expr_copy, a
         * refcount bump), so share the original node instead of rebuilding
         * it.  A rebuilt copy is not merely an allocation per node: it loses
         * the evaluation stamp and cached hash, so the evaluator re-sorts and
         * re-evaluates a large stored right-hand side on every application --
         * res[n] := <135135-term sum> was rebuilt and re-canonicalised on each
         * read (issue #112).  Sharing is the OwnValue path's long-standing
         * behaviour (symtab apply_own_values returns expr_copy(replacement)). */
        if (new_head == expr->data.function.head &&
            new_count == expr->data.function.arg_count) {
            bool unchanged = true;
            for (size_t i = 0; i < expr->data.function.arg_count; i++) {
                if (spliced[i] || temp_args[i] != expr->data.function.args[i]) {
                    unchanged = false;
                    break;
                }
            }
            if (unchanged) {
                for (size_t i = 0; i < expr->data.function.arg_count; i++)
                    expr_free(temp_args[i]);
                expr_free(new_head);
                if (temp_args != temp_stack) free(temp_args);
                if (spliced != spliced_stack) free(spliced);
                return expr_copy(expr);
            }
        }

        Expr* final_stack[RB_STACK_N];
        Expr** final_args = new_count <= RB_STACK_N ? final_stack
                          : malloc(sizeof(Expr*) * (new_count ? new_count : 1));
        size_t idx = 0;
        for (size_t i = 0; i < expr->data.function.arg_count; i++) {
            if (spliced[i]) {
                for (size_t j = 0; j < temp_args[i]->data.function.arg_count; j++) {
                    final_args[idx++] = expr_copy(temp_args[i]->data.function.args[j]);
                }
                expr_free(temp_args[i]);
            } else {
                final_args[idx++] = temp_args[i];
            }
        }
        if (temp_args != temp_stack) free(temp_args);
        if (spliced != spliced_stack) free(spliced);

        Expr* result = expr_new_function(new_head, final_args, new_count);
        if (final_args != final_stack) free(final_args);
        return result;
    }
    
    return expr_copy(expr);
}

/* ------------------------------------------------------------------ *
 *  Flat-head "leftover" rule application
 *
 *  Mathematica applies a rewrite rule whose LHS head is Flat (Plus, Times,
 *  Dot, ...) to a SUBSET of a longer operand sequence: it rewrites the matched
 *  part and keeps the rest.  E.g. with
 *      mod /: mod[a_,p_] + mod[b_,p_] := mod[Mod[a+b,p],p]
 *  the sum  mod[1,5] + mod[2,5] + mod[3,5]  rewrites a matching pair and leaves
 *  the third term (which then re-fires to a fixed point).  The generic matcher
 *  has no such "leftover" path: its argument base case requires every operand
 *  to be consumed, and concrete-head pattern elements (mod[a_,p_]) cannot
 *  group-absorb the surplus the way a bare blank can.
 *
 *  We add this ONLY here, at the rule-application layer, never inside match()
 *  itself -- so MatchQ stays a whole-expression test (MatchQ[a+b+c, f[_]+f[_]]
 *  is still False).  The strategy is "normal match first, augment on failure":
 *
 *    1. Try the ordinary match.  This preserves today's Flat group-absorption,
 *       which is the WL-correct result whenever a trailing blank CAN swallow the
 *       surplus (pp /: Plus[pp,a_]:=... on pp+x+y binds a_ = x+y).  Trying this
 *       first is load-bearing for correctness, not just speed: augmenting that
 *       case would wrongly split a_ = x, leftover = y.  This step stays INLINE in
 *       each caller (a bare match()+replace_bindings, exactly as before), so the
 *       common path pays nothing for this feature.
 *    2. Only if that fails, the caller calls match_flat_leftover below: if the
 *       enclosing head is Flat with MORE operands than the pattern's flat-call
 *       arity, it retries against an augmented pattern
 *           head[p1, ..., pk, $FlatRest$___]
 *       with the replacement wrapped as  head[replacement, $FlatRest$].  The
 *       BlankNullSequence binds the surplus operands, which rb_rec then splices
 *       back into the (Flat) head.
 *
 *  Orderless heads (Plus, Times): the subset enumeration already reaches every
 *  position, so the one augmented match against the whole subject suffices.
 *
 *  Ordered Flat heads (Dot, **, And, Or, a user head with only Flat): the
 *  matcher consumes arguments left to right, so the augmented pattern alone
 *  finds only a PREFIX run.  The matched run may sit anywhere (b.c in a.b.c.d,
 *  y**z in x**y**z), so the augmented match is retried at each start offset
 *  i = 0, 1, ... against the suffix head[a_i, ..., a_n]; the first offset that
 *  matches wins -- the LEFTMOST run, as in WL -- and the operands before it are
 *  re-attached: head[a_1, ..., a_(i-1), <rewritten suffix>].  Offset 0 is
 *  exactly the prefix-run match.  This is equivalent to the pattern
 *  head[pre___, p1, ..., pk, rest___] with pre___ enumerated shortest-first,
 *  but cheaper, because offsets and spans that cannot succeed are skipped
 *  before any work is done (see the span-analysis helpers near the top):
 *    - the first ANCHOR (an element that consumes exactly one operand, e.g.
 *      op[i_] or Not[x_]) bounds the viable offsets: exactly when every
 *      element before it is single-operand, else by the last operand the
 *      anchor accepts; a failing rule whose anchor accepts nothing costs O(n)
 *      operand tests;
 *    - inside the matcher, a variable-span element followed by an anchor
 *      only tries spans whose next operand passes the anchor, and an
 *      already-bound name (ff[x_, x_]) rejects a group of the wrong size or
 *      contents without building it.
 *  What remains is the work the semantics requires: a test on a
 *  group-absorbing first element (ff[a_?Negative, b_]) is a distinct call for
 *  every (offset, span) group, so a failing rule of that shape still makes
 *  ~n^2/2 test calls.
 *
 *  One rewrite per application: like ReplaceAll at a single position, only the
 *  first (leftmost) run is rewritten -- ff[a,b,a,b] /. ff[a,b] -> w gives
 *  ff[w, a, b]; ReplaceRepeated (or a DownValue/UpValue re-firing to a fixed
 *  point) rewrites the rest.  This mirrors the Plus/Times leftover behaviour.
 * ------------------------------------------------------------------ */

/* Unwrap HoldPattern/Condition down to the inner flat-call node; return it when
 * its head is `flat_head`, else NULL. */
static Expr* flat_inner_call(Expr* pattern, const char* flat_head) {
    while (pattern && pattern->type == EXPR_FUNCTION &&
           pattern->data.function.head &&
           pattern->data.function.head->type == EXPR_SYMBOL) {
        const char* h = pattern->data.function.head->data.symbol.name;
        if (h == SYM_HoldPattern && pattern->data.function.arg_count == 1)
            pattern = pattern->data.function.args[0];
        else if (h == SYM_Condition && pattern->data.function.arg_count == 2)
            pattern = pattern->data.function.args[0];
        else
            return (h == flat_head) ? pattern : NULL;
    }
    return NULL;
}

/* Pattern[$FlatRest$, BlankNullSequence[]] -- the trailing leftover catcher.
 * `$FlatRest$` is a reserved, dollar-delimited internal name (cf. the
 * `$OptionsPattern$` convention above); it is transient -- only ever bound in a
 * MatchEnv and spliced out by rb_rec -- so it never reaches Names[] or output. */
static Expr* flat_rest_pattern(void) {
    Expr* bnull = expr_new_function(expr_new_symbol(SYM_BlankNullSequence), NULL, 0);
    Expr** pa = malloc(sizeof(Expr*) * 2);
    pa[0] = expr_new_symbol("$FlatRest$");
    pa[1] = bnull;
    Expr* p = expr_new_function(expr_new_symbol(SYM_Pattern), pa, 2);
    free(pa);
    return p;
}

/* Rebuild `pattern` with the leftover catcher appended to the inner flat-call's
 * argument list, preserving any HoldPattern / Condition wrappers (so a `/;`
 * guard, whose bound variables are untouched, still fires).  Returns a fully
 * owned copy, or NULL if `pattern` is not a (wrapped) `flat_head` call. */
static Expr* flat_augment_pattern(Expr* pattern, const char* flat_head) {
    if (!pattern || pattern->type != EXPR_FUNCTION ||
        !pattern->data.function.head ||
        pattern->data.function.head->type != EXPR_SYMBOL) return NULL;
    const char* h = pattern->data.function.head->data.symbol.name;
    if (h == SYM_HoldPattern && pattern->data.function.arg_count == 1) {
        Expr* inner = flat_augment_pattern(pattern->data.function.args[0], flat_head);
        if (!inner) return NULL;
        Expr** a = malloc(sizeof(Expr*));
        a[0] = inner;
        Expr* r = expr_new_function(expr_new_symbol(SYM_HoldPattern), a, 1);
        free(a);
        return r;
    }
    if (h == SYM_Condition && pattern->data.function.arg_count == 2) {
        Expr* inner = flat_augment_pattern(pattern->data.function.args[0], flat_head);
        if (!inner) return NULL;
        Expr** a = malloc(sizeof(Expr*) * 2);
        a[0] = inner;
        a[1] = expr_copy(pattern->data.function.args[1]);
        Expr* r = expr_new_function(expr_new_symbol(SYM_Condition), a, 2);
        free(a);
        return r;
    }
    if (h == flat_head) {
        size_t k = pattern->data.function.arg_count;
        Expr** a = malloc(sizeof(Expr*) * (k + 1));
        for (size_t i = 0; i < k; i++) a[i] = expr_copy(pattern->data.function.args[i]);
        a[k] = flat_rest_pattern();
        Expr* r = expr_new_function(expr_copy(pattern->data.function.head), a, k + 1);
        free(a);
        return r;
    }
    return NULL;
}

/* Occurrences of the symbol `name` anywhere in `e` (heads included). */
static size_t flat_symbol_count(const Expr* e, const char* name) {
    if (!e) return 0;
    if (e->type == EXPR_SYMBOL) return strcmp(e->data.symbol.name, name) == 0;
    if (e->type != EXPR_FUNCTION) return 0;
    size_t c = flat_symbol_count(e->data.function.head, name);
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        c += flat_symbol_count(e->data.function.args[i], name);
    return c;
}

/* True when `name` occurs in the guard of some Condition inside `e`. */
static bool flat_guard_mentions(const Expr* e, const char* name) {
    if (!e || e->type != EXPR_FUNCTION) return false;
    if (e->data.function.head->type == EXPR_SYMBOL
        && e->data.function.head->data.symbol.name == SYM_Condition
        && e->data.function.arg_count == 2
        && flat_symbol_count(e->data.function.args[1], name) > 0) return true;
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (flat_guard_mentions(e->data.function.args[i], name)) return true;
    return false;
}

/* True when the flat call `inner` (of `pattern`) has an element that is an
 * untyped, untested sequence blank -- `___`, `__`, or `r___` / `r__` whose name
 * occurs nowhere else in the pattern and in no guard of the replacement.  For
 * an Orderless head such an element makes the leftover retry redundant: any
 * augmented match {non-r elements -> their operands, r -> S, $FlatRest$ -> L}
 * is also a plain match with r -> S + L, because nothing reads r's binding to
 * accept or reject it.  The plain match has already failed, so the augmented
 * one must too -- and with r no longer last it would enumerate every subset of
 * the operands for r (2^n tries), e.g. the TrigReduce cancellation rule
 * c1_. Cosh[a_] + c2_. Cosh[b_] + r___ /; ... over a 16-term sum (20 s). */
static bool flat_has_free_rest(Expr* pattern, Expr* inner, Expr* replacement) {
    for (size_t i = 0; i < inner->data.function.arg_count; i++) {
        Expr* a = inner->data.function.args[i];
        const char* name = NULL;
        if (a->type == EXPR_FUNCTION && a->data.function.head->type == EXPR_SYMBOL
            && a->data.function.head->data.symbol.name == SYM_Pattern
            && a->data.function.arg_count == 2
            && a->data.function.args[0]->type == EXPR_SYMBOL) {
            name = a->data.function.args[0]->data.symbol.name;
            a = a->data.function.args[1];
        }
        if (a->type != EXPR_FUNCTION || a->data.function.head->type != EXPR_SYMBOL
            || a->data.function.arg_count != 0) continue;
        const char* bh = a->data.function.head->data.symbol.name;
        if (bh != SYM_BlankSequence && bh != SYM_BlankNullSequence) continue;
        if (!name) return true;
        if (flat_symbol_count(pattern, name) == 1
            && !flat_guard_mentions(replacement, name)) return true;
    }
    return false;
}

/* head[replacement, $FlatRest$] -- the matched RHS plus the spliced leftovers. */
static Expr* flat_augment_replacement(Expr* replacement, const char* flat_head) {
    Expr** a = malloc(sizeof(Expr*) * 2);
    a[0] = expr_copy(replacement);
    a[1] = expr_new_symbol("$FlatRest$");
    Expr* r = expr_new_function(expr_new_symbol(flat_head), a, 2);
    free(a);
    return r;
}

/* The augmented match at one start offset: `sub` is the subject (the whole
 * expression, or a suffix of it).  Returns the rewritten `sub`, or NULL. */
static Expr* flat_try_augmented(Expr* sub, Expr* aug_pat, Expr* replacement,
                                const char* flat_head, MatchEnv* env) {
    env_rollback(env, 0);       /* drop the previous attempt's partial bindings */
    if (!match(sub, aug_pat, env)) return NULL;
    Expr* aug_rep = flat_augment_replacement(replacement, flat_head);
    Expr* r = replace_bindings(aug_rep, env);      /* env left populated for caller */
    expr_free(aug_rep);
    return r;
}

Expr* match_flat_leftover(Expr* expr, Expr* pattern, Expr* replacement, MatchEnv* env) {
    /* Called by a rule-dispatch site AFTER its own inline match() has failed
     * (that failed match may have left partial bindings in `env`, which must be
     * a fresh per-rule env). Gate cheapest-first: the CHEAP, highly-selective
     * structural checks run before the (symbol-table) attribute lookup, so an
     * ordinary failed match on a non-Flat head costs only a few pointer derefs.
     * The dispatch prefilters relax the arity key ONLY for Flat heads, so for a
     * non-Flat head the pattern arity equals the input arity and the (b) test
     * below bails before get_attributes(). */
    if (!expr || expr->type != EXPR_FUNCTION || !env) return NULL;
    Expr* headnode = expr->data.function.head;
    if (!headnode || headnode->type != EXPR_SYMBOL) return NULL;
    const char* flat_head = headnode->data.symbol.name;

    /* (a) pattern must be a (HoldPattern/Condition-wrapped) call on this head */
    Expr* inner = flat_inner_call(pattern, flat_head);
    if (!inner) return NULL;
    /* (b) and the subject must carry strictly more operands than it specifies */
    size_t pat_arity = inner->data.function.arg_count;
    size_t n = expr->data.function.arg_count;
    if (pat_arity < 1 || n <= pat_arity) return NULL;
    /* (c) only now, the attribute lookup: the head must actually be Flat */
    uint32_t attrs = get_attributes(flat_head);
    if (!(attrs & ATTR_FLAT)) return NULL;
    /* (d) an Orderless head with a free trailing catcher cannot gain a match */
    if ((attrs & ATTR_ORDERLESS) && flat_has_free_rest(pattern, inner, replacement))
        return NULL;

    Expr* aug_pat = flat_augment_pattern(pattern, flat_head);
    if (!aug_pat) return NULL;

    /* Offset 0: the prefix run (the only offset an Orderless head needs). */
    Expr* result = flat_try_augmented(expr, aug_pat, replacement, flat_head, env);

    if (!result && !(attrs & ATTR_ORDERLESS)) {
        /* Ordered head: slide the start offset rightwards, leftmost first.
         *
         * Offset screen.  Let p_t be the first anchor of the pattern (an
         * element that consumes exactly one operand; for t > 0 it must also be
         * screenable, since p_0..p_(t-1) are bound before it).  At offset i
         * the anchor lands on operand i + t when every element before it is
         * single-span ("exact"), else at some operand >= i + minpre.  An
         * offset is viable only if that operand (exact) or SOME operand at or
         * after i + minpre (inexact: <= the last accepting operand, jmax)
         * passes the anchor standalone.  Independently, the suffix must hold
         * at least `minall` operands.  Every skipped offset would fail the full
         * match, so the leftmost viable match -- the result -- is unchanged;
         * a failing rule whose anchor accepts no operand now costs O(n)
         * operand tests instead of n full matches. */
        Expr** pargs = inner->data.function.args;
        size_t t = pat_arity, minpre = 0, minall = 0;
        bool exact = true, has_ops = false;
        for (size_t e = 0; e < pat_arity; e++) {
            Expr* bnd = NULL;
            if (ops_unwrap(pargs[e], &bnd)) has_ops = true;
            minall += flat_elem_min_span(pargs[e], flat_head);
        }
        if (!has_ops) {
            for (size_t e = 0; e < pat_arity && t == pat_arity; e++) {
                bool single = flat_elem_single_span(pargs[e], flat_head);
                if (single && (e == 0 || flat_elem_screenable(pargs[e], flat_head))) {
                    t = e;
                    break;
                }
                exact = exact && single;
                minpre += flat_elem_min_span(pargs[e], flat_head);
            }
        }
        Expr** args = expr->data.function.args;
        const char* thead = (t < pat_arity) ? flat_anchor_head(pargs[t]) : NULL;
        /* inexact anchor: the last operand it accepts bounds the viable offsets */
        size_t jmax = 0;
        bool have_jmax = false;
        if (t < pat_arity && !exact) {
            for (size_t j = n; j-- > minpre + 1; ) {
                if (flat_anchor_accepts(args[j], pargs[t], thead)) { jmax = j; have_jmax = true; break; }
            }
            if (!have_jmax) { expr_free(aug_pat); return NULL; }  /* no offset >= 1 viable */
        }
        for (size_t i = 1; i < n && !result; i++) {
            if (n - i < minall) break;            /* suffix too short, and shrinking */
            if (t < pat_arity) {
                if (exact) {
                    if (i + t >= n || !flat_anchor_accepts(args[i + t], pargs[t], thead)) continue;
                } else if (i + minpre > jmax) {
                    break;                        /* anchor cannot land at or after i+minpre */
                }
            }
            size_t m = n - i;
            Expr** sa = malloc(sizeof(Expr*) * m);
            for (size_t j = 0; j < m; j++) sa[j] = expr_copy(args[i + j]);
            Expr* sub = expr_new_function(expr_copy(headnode), sa, m);
            free(sa);
            Expr* rew = flat_try_augmented(sub, aug_pat, replacement, flat_head, env);
            expr_free(sub);
            if (!rew) continue;
            /* head[a_1, ..., a_(i-1), <rewritten suffix>]; the Flat head
             * re-flattens the nested call when the result is evaluated. */
            Expr** ra = malloc(sizeof(Expr*) * (i + 1));
            for (size_t j = 0; j < i; j++) ra[j] = expr_copy(args[j]);
            ra[i] = rew;
            result = expr_new_function(expr_copy(headnode), ra, i + 1);
            free(ra);
        }
    }
    expr_free(aug_pat);
    return result;
}

Expr* builtin_matchq(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 2) return NULL;
    
    Expr* expr = evaluate(res->data.function.args[0]);
    Expr* pattern = res->data.function.args[1];
    
    MatchEnv* env = env_new();
    bool result = match(expr, pattern, env);
    env_free(env);
    expr_free(expr);
    
    return expr_new_symbol(result ? "True" : "False");
}
