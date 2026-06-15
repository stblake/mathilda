#include "simp.h"
#include "simp_internal.h"
#include "arithmetic.h"
#include "attr.h"
#include "common.h"
#include "eval.h"
#include "expand.h"
#include "facpoly.h"
#include "numeric.h"
#include "parse.h"
#include "print.h"
#include "symtab.h"
#include "expr.h"
#include "rationalize.h"
#include "sym_names.h"
#include "sym_intern.h"
#include "trigrat.h"
#include "qa.h"
#include "qafactor.h"
#include "simp_log.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <gmp.h>
#ifdef USE_MPFR
#include <mpfr.h>
#endif


/* ----------------------------------------------------------------------- */
/* TanAddition: rewrite Tan[c]/Cot[c] when c = a+b and Tan[a],Tan[b] occur */
/* ----------------------------------------------------------------------- */

/* The standard angle-addition identities:
 *     Tan[a + b] = (Tan[a] + Tan[b]) / (1 - Tan[a] Tan[b])
 *     Cot[a + b] = (1 - Tan[a] Tan[b]) / (Tan[a] + Tan[b])
 * are always sound (away from the isolated singularities where the
 * denominator vanishes -- inside Simplify those are ignored).
 *
 * Mathilda's TrigExpand fires only when the argument is *literally* a
 * Plus, so it can rewrite Tan[x+y] but NOT Tan[5] -- even when the
 * surrounding expression also contains Tan[2] and Tan[3] making
 * Tan[5] = Tan[2 + 3] usable.  This transform performs that recognition:
 *
 *   1. Walk the input collecting every distinct expression that appears
 *      as the argument of any Tan, Cot, Sin, Cos, Sec, or Csc head.
 *   2. For each ordered pair (a, b) in that set with a != b, evaluate
 *      a + b.  If the sum is also in the set, the triple (a, b, c=a+b)
 *      witnesses an addition-formula opportunity.
 *   3. For each such triple, build a RuleDelayed list that rewrites the
 *      heads which actually appear at c (Tan[c], Cot[c], etc.) using the
 *      angle-addition formula in terms of Tan[a], Tan[b] (or Sin/Cos[a]
 *      and Sin/Cos[b], for Sin[c] / Cos[c] occurrences).
 *   4. Apply ReplaceAll with the constructed rules, then evaluate +
 *      Together + Cancel so the polynomial cancellation fires.
 *   5. Score-gated: keep only strict wins.  Inert when no pair-sum match
 *      exists in the input (the common case).
 *
 * Coverage of the case-6 shape:
 *     Tan[2] Tan[3] B A - (-Tan[2]/Tan[5] - Tan[3]/Tan[5] + 1) B A  ->  0
 * via the (2, 3, 5) triple's Cot[5] -> (1 - Tan[2] Tan[3]) / (Tan[2] +
 * Tan[3]) substitution.  Symbolic shapes (Tan[x], Tan[y], Tan[x+y])
 * also work; the integer-arg case is just the most surprising one.
 */

typedef struct {
    Expr** items;       /* Borrowed pointers into the input tree. */
    size_t count;
    size_t cap;
} TrigArgSet;

static void tas_init(TrigArgSet* s) {
    s->items = NULL; s->count = 0; s->cap = 0;
}
static void tas_free(TrigArgSet* s) {
    free(s->items);
    s->items = NULL; s->count = 0; s->cap = 0;
}
static bool tas_contains(const TrigArgSet* s, const Expr* e) {
    for (size_t i = 0; i < s->count; i++) {
        if (expr_eq(s->items[i], e)) return true;
    }
    return false;
}
static void tas_add_borrowed(TrigArgSet* s, Expr* e) {
    if (tas_contains(s, e)) return;
    if (s->count == s->cap) {
        s->cap = s->cap ? s->cap * 2 : 8;
        s->items = (Expr**)realloc(s->items, sizeof(Expr*) * s->cap);
    }
    s->items[s->count++] = e;
}

static bool is_addition_trig_head(const char* h) {
    return h == SYM_Tan || h == SYM_Cot || h == SYM_Sin
        || h == SYM_Cos || h == SYM_Sec || h == SYM_Csc;
}

static void collect_addition_trig_args(const Expr* e, TrigArgSet* set) {
    if (!e || e->type != EXPR_FUNCTION) return;
    Expr* head = e->data.function.head;
    if (head && head->type == EXPR_SYMBOL
        && is_addition_trig_head(head->data.symbol)
        && e->data.function.arg_count == 1) {
        tas_add_borrowed(set, e->data.function.args[0]);
    }
    if (head) collect_addition_trig_args(head, set);
    for (size_t i = 0; i < e->data.function.arg_count; i++) {
        collect_addition_trig_args(e->data.function.args[i], set);
    }
}

/* Returns true if `e` contains Power[head[c], _] or head[c] for the
 * given trig head symbol.  A cheap substring presence check without a
 * full pattern match. */
static bool tree_contains_trig_at(const Expr* e, const char* head_sym,
                                   const Expr* arg) {
    if (!e || e->type != EXPR_FUNCTION) return false;
    Expr* head = e->data.function.head;
    if (head && head->type == EXPR_SYMBOL
        && head->data.symbol == head_sym
        && e->data.function.arg_count == 1
        && expr_eq(e->data.function.args[0], arg)) {
        return true;
    }
    if (head && tree_contains_trig_at(head, head_sym, arg)) return true;
    for (size_t i = 0; i < e->data.function.arg_count; i++) {
        if (tree_contains_trig_at(e->data.function.args[i], head_sym, arg)) {
            return true;
        }
    }
    return false;
}

static Expr* mk_unary_call(const char* head, Expr* arg) {
    Expr* args[1] = { arg };
    return expr_new_function(expr_new_symbol(head), args, 1);
}

static Expr* mk_div(Expr* num, Expr* den) {
    Expr* inv_args[2] = { den, expr_new_integer(-1) };
    Expr* inv = expr_new_function(
        expr_new_symbol(SYM_Power), inv_args, 2);
    Expr* args[2] = { num, inv };
    return expr_new_function(expr_new_symbol(SYM_Times), args, 2);
}

/* Append rule list entries (RuleDelayed nodes) for the heads that
 * actually appear at `c` in `input`. */
static void append_addition_rules_for_triple(Expr* input,
                                              const Expr* a, const Expr* b,
                                              const Expr* c,
                                              Expr*** rules, size_t* rcount,
                                              size_t* rcap) {
    #define RAPPEND(node) do {                                \
        if (*rcount == *rcap) {                               \
            *rcap = *rcap ? *rcap * 2 : 8;                    \
            *rules = (Expr**)realloc(*rules,                  \
                                      sizeof(Expr*) * (*rcap)); \
        }                                                     \
        (*rules)[(*rcount)++] = (node);                       \
    } while (0)

    /* All six rules (Tan/Cot/Sin/Cos/Sec/Csc at c) are expressed in terms
     * of the same two Sin/Cos primitives:
     *     sin_rhs = Sin[a] Cos[b] + Cos[a] Sin[b]   (= Sin[a + b])
     *     cos_rhs = Cos[a] Cos[b] - Sin[a] Sin[b]   (= Cos[a + b])
     * Using the same primitives for Tan[c] = sin_rhs / cos_rhs and
     * Sec[c] = 1 / cos_rhs (rather than the Tan-based identity Tan[c] =
     * (Tan[a]+Tan[b])/(1 - Tan[a] Tan[b])) gives the post-substitution
     * expression a uniform Sin/Cos denominator structure, so a single
     * Together / Cancel pass collapses mixed Sec[c] + Tan[c] inputs that
     * the Tan-based formula leaves with mismatched (1 - Tan[a] Tan[b])
     * vs (Cos[a] Cos[b] - Sin[a] Sin[b]) denominators. */
    bool need_tan = tree_contains_trig_at(input, SYM_Tan, c);
    bool need_cot = tree_contains_trig_at(input, SYM_Cot, c);
    bool need_sin = tree_contains_trig_at(input, SYM_Sin, c);
    bool need_csc = tree_contains_trig_at(input, SYM_Csc, c);
    bool need_cos = tree_contains_trig_at(input, SYM_Cos, c);
    bool need_sec = tree_contains_trig_at(input, SYM_Sec, c);

    bool need_sin_rhs = need_sin || need_csc || need_tan || need_cot;
    bool need_cos_rhs = need_cos || need_sec || need_tan || need_cot;

    Expr* sin_rhs = NULL;
    if (need_sin_rhs) {
        Expr* sa_cb_args[2] = {
            mk_unary_call("Sin", expr_copy((Expr*)a)),
            mk_unary_call("Cos", expr_copy((Expr*)b))
        };
        Expr* sa_cb = expr_new_function(
            expr_new_symbol(SYM_Times), sa_cb_args, 2);
        Expr* ca_sb_args[2] = {
            mk_unary_call("Cos", expr_copy((Expr*)a)),
            mk_unary_call("Sin", expr_copy((Expr*)b))
        };
        Expr* ca_sb = expr_new_function(
            expr_new_symbol(SYM_Times), ca_sb_args, 2);
        Expr* sin_args[2] = { sa_cb, ca_sb };
        sin_rhs = expr_new_function(
            expr_new_symbol(SYM_Plus), sin_args, 2);
    }

    Expr* cos_rhs = NULL;
    if (need_cos_rhs) {
        Expr* ca_cb_args[2] = {
            mk_unary_call("Cos", expr_copy((Expr*)a)),
            mk_unary_call("Cos", expr_copy((Expr*)b))
        };
        Expr* ca_cb = expr_new_function(
            expr_new_symbol(SYM_Times), ca_cb_args, 2);
        Expr* neg_sa_sb_args[3] = {
            expr_new_integer(-1),
            mk_unary_call("Sin", expr_copy((Expr*)a)),
            mk_unary_call("Sin", expr_copy((Expr*)b))
        };
        Expr* neg_sa_sb = expr_new_function(
            expr_new_symbol(SYM_Times), neg_sa_sb_args, 3);
        Expr* cos_args[2] = { ca_cb, neg_sa_sb };
        cos_rhs = expr_new_function(
            expr_new_symbol(SYM_Plus), cos_args, 2);
    }

    if (need_sin) {
        Expr* lhs = mk_unary_call("Sin", expr_copy((Expr*)c));
        Expr* rule_args[2] = { lhs, expr_copy(sin_rhs) };
        Expr* rule = expr_new_function(
            expr_new_symbol(SYM_RuleDelayed), rule_args, 2);
        RAPPEND(rule);
    }
    if (need_csc) {
        Expr* lhs = mk_unary_call("Csc", expr_copy((Expr*)c));
        Expr* inv_args[2] = { expr_copy(sin_rhs), expr_new_integer(-1) };
        Expr* inv = expr_new_function(
            expr_new_symbol(SYM_Power), inv_args, 2);
        Expr* rule_args[2] = { lhs, inv };
        Expr* rule = expr_new_function(
            expr_new_symbol(SYM_RuleDelayed), rule_args, 2);
        RAPPEND(rule);
    }
    if (need_cos) {
        Expr* lhs = mk_unary_call("Cos", expr_copy((Expr*)c));
        Expr* rule_args[2] = { lhs, expr_copy(cos_rhs) };
        Expr* rule = expr_new_function(
            expr_new_symbol(SYM_RuleDelayed), rule_args, 2);
        RAPPEND(rule);
    }
    if (need_sec) {
        Expr* lhs = mk_unary_call("Sec", expr_copy((Expr*)c));
        Expr* inv_args[2] = { expr_copy(cos_rhs), expr_new_integer(-1) };
        Expr* inv = expr_new_function(
            expr_new_symbol(SYM_Power), inv_args, 2);
        Expr* rule_args[2] = { lhs, inv };
        Expr* rule = expr_new_function(
            expr_new_symbol(SYM_RuleDelayed), rule_args, 2);
        RAPPEND(rule);
    }
    if (need_tan) {
        Expr* lhs = mk_unary_call("Tan", expr_copy((Expr*)c));
        Expr* rhs = mk_div(expr_copy(sin_rhs), expr_copy(cos_rhs));
        Expr* rule_args[2] = { lhs, rhs };
        Expr* rule = expr_new_function(
            expr_new_symbol(SYM_RuleDelayed), rule_args, 2);
        RAPPEND(rule);
    }
    if (need_cot) {
        Expr* lhs = mk_unary_call("Cot", expr_copy((Expr*)c));
        Expr* rhs = mk_div(expr_copy(cos_rhs), expr_copy(sin_rhs));
        Expr* rule_args[2] = { lhs, rhs };
        Expr* rule = expr_new_function(
            expr_new_symbol(SYM_RuleDelayed), rule_args, 2);
        RAPPEND(rule);
    }

    if (sin_rhs) expr_free(sin_rhs);
    if (cos_rhs) expr_free(cos_rhs);
    #undef RAPPEND
}

static Expr* transform_tan_addition_impl(const Expr* e) {
    bool dbg = simp_debug_enabled();
    clock_t t0 = dbg ? clock() : 0;

    /* Cheap gate: if no trig head appears in the input, we cannot collect
     * any args and the rule list will be empty. */
    if (!has_pythag_head(e)) {
        Expr* out = expr_copy((Expr*)e);
        if (dbg) simp_debug_log("TanAddition", e, out,
                                simp_debug_elapsed_ms(t0));
        return out;
    }

    TrigArgSet args; tas_init(&args);
    collect_addition_trig_args(e, &args);
    if (args.count < 3) {
        /* Need at least three distinct args to have a (a, b, c=a+b) triple. */
        tas_free(&args);
        Expr* out = expr_copy((Expr*)e);
        if (dbg) simp_debug_log("TanAddition", e, out,
                                simp_debug_elapsed_ms(t0));
        return out;
    }

    Expr** rules = NULL;
    size_t rcount = 0, rcap = 0;
    Expr* input_copy = (Expr*)e;  /* tree_contains_trig_at takes const */
    for (size_t i = 0; i < args.count; i++) {
        for (size_t j = 0; j < args.count; j++) {
            if (i == j) continue;
            Expr* sum_args[2] = {
                expr_copy(args.items[i]), expr_copy(args.items[j])
            };
            Expr* sum_call = expr_new_function(
                expr_new_symbol(SYM_Plus), sum_args, 2);
            Expr* sum = eval_and_free(sum_call);
            if (!sum) continue;
            if (tas_contains(&args, sum)) {
                /* Use the canonical sum as `c` for the rule LHS.  The
                 * input might contain Tan[Plus[a, b]] (if a + b doesn't
                 * collapse to a literal) or Tan[5] (if a, b are
                 * integers), so we need to use the actually-occurring
                 * representative -- pull that out of the args set. */
                Expr* canonical_c = NULL;
                for (size_t k = 0; k < args.count; k++) {
                    if (expr_eq(args.items[k], sum)) {
                        canonical_c = args.items[k];
                        break;
                    }
                }
                if (canonical_c) {
                    append_addition_rules_for_triple(
                        input_copy,
                        args.items[i], args.items[j], canonical_c,
                        &rules, &rcount, &rcap);
                }
            }
            expr_free(sum);
        }
    }
    tas_free(&args);

    if (rcount == 0) {
        free(rules);
        Expr* out = expr_copy((Expr*)e);
        if (dbg) simp_debug_log("TanAddition", e, out,
                                simp_debug_elapsed_ms(t0));
        return out;
    }

    /* Assemble rule list { rule1, rule2, ... } and apply ReplaceAll. */
    Expr* rule_list = expr_new_function(
        expr_new_symbol(SYM_List), rules, rcount);
    free(rules);
    Expr* ra_args[2] = { expr_copy((Expr*)e), rule_list };
    Expr* ra_call = expr_new_function(
        expr_new_symbol(SYM_ReplaceAll), ra_args, 2);
    Expr* substituted = eval_and_free(ra_call);
    if (!substituted) {
        Expr* out = expr_copy((Expr*)e);
        if (dbg) simp_debug_log("TanAddition", e, out,
                                simp_debug_elapsed_ms(t0));
        return out;
    }

    /* Together + Cancel turns the rational expression into a single
     * fraction, then reduces by polynomial GCD.  When the substitution
     * was witnessing a true identity (the case-6 shape), this drops to
     * 0 cleanly. */
    Expr* tg_args[1] = { substituted };
    Expr* tg_call = expr_new_function(
        expr_new_symbol(SYM_Together), tg_args, 1);
    Expr* tg = eval_and_free(tg_call);
    Expr* result = tg ? tg : expr_copy((Expr*)e);

    Expr* cn_args[1] = { expr_copy(result) };
    Expr* cn_call = expr_new_function(
        expr_new_symbol(SYM_Cancel), cn_args, 1);
    Expr* cn = eval_and_free(cn_call);
    if (cn) { expr_free(result); result = cn; }

    if (dbg) simp_debug_log("TanAddition", e, result,
                            simp_debug_elapsed_ms(t0));
    return result;
}

Expr* transform_tan_addition(const Expr* e) {
    return simp_memo_wrap(e, "$TanAddition", transform_tan_addition_impl);
}

