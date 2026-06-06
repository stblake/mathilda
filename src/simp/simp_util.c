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
/* Helpers                                                                 */
/* ----------------------------------------------------------------------- */

/* Build f[arg], evaluate, and return the result. Takes ownership of `arg`. */
Expr* call_unary_owned(const char* head_name, Expr* arg) {
    Expr* a[1] = { arg };
    Expr* call = expr_new_function(expr_new_symbol(head_name), a, 1);
    Expr* r = evaluate(call);
    expr_free(call);
    return r;
}

Expr* call_unary_copy(const char* head_name, const Expr* arg) {
    return call_unary_owned(head_name, expr_copy((Expr*)arg));
}

/* ----------------------------------------------------------------------- */
/* $SimplifyDebug -- per-transform tracing                                 */
/* ----------------------------------------------------------------------- */

/*
 * When $SimplifyDebug is set to True, every transform invocation inside
 * simp_search emits one line on stderr in the format
 *   /<TransformName>/: <input> -> <output> [<elapsed> ms]
 * This is used to diagnose pathological inputs (Simplify hangs, runaway
 * candidate explosion, expensive single transforms). The check is read
 * directly off the OwnValue list -- evaluating $SimplifyDebug would
 * itself fire the OwnValue rule on every call. */
bool simp_debug_enabled(void) {
    Rule* r = symtab_get_own_values("$SimplifyDebug");
    if (!r || !r->replacement) return false;
    Expr* v = r->replacement;
    return v->type == EXPR_SYMBOL && v->data.symbol == SYM_True;
}

double simp_debug_elapsed_ms(clock_t t0) {
    return (double)(clock() - t0) * 1000.0 / (double)CLOCKS_PER_SEC;
}

void simp_debug_log(const char* xform, const Expr* in,
                           const Expr* out, double ms) {
    char* sin  = expr_to_string((Expr*)in);
    char* sout = out ? expr_to_string((Expr*)out) : NULL;
    fprintf(stderr, "/%s/: %s -> %s [%.2f ms]\n",
            xform,
            sin  ? sin  : "?",
            sout ? sout : "(no change)",
            ms);
    free(sin);
    free(sout);
    fflush(stderr);
}

/* Wrap call_unary_copy with tracing when $SimplifyDebug is True.
 *
 * Note: an experimental generic FactorMemo lookup at this layer was
 * tried (Phase 11 attempt) but reverted -- the per-transform memos
 * already in place (Factor, TrigFactor, TrigExpand, TrigRoundtrip,
 * PythagReduce, PythagSquareComplete, HalfAngle) cover the high-
 * volume duplicates, and the additional malloc/hash overhead at
 * every call exceeded the marginal gain on cheap transforms like
 * Together / Cancel / Apart that are individually fast and rarely
 * repeated. */
Expr* traced_call_unary(const char* xform, const Expr* in) {
    bool dbg = simp_debug_enabled();
    clock_t t0 = dbg ? clock() : 0;
    Expr* r = call_unary_copy(xform, in);
    if (dbg) simp_debug_log(xform, in, r, simp_debug_elapsed_ms(t0));
    return r;
}

bool is_rule_with_lhs(const Expr* e, const char* lhs_symbol) {
    if (!e || e->type != EXPR_FUNCTION) return false;
    if (e->data.function.arg_count != 2) return false;
    if (!e->data.function.head || e->data.function.head->type != EXPR_SYMBOL) return false;
    const char* h = e->data.function.head->data.symbol;
    if (h != SYM_Rule && h != SYM_RuleDelayed) return false;
    Expr* k = e->data.function.args[0];
    return k && k->type == EXPR_SYMBOL && strcmp(k->data.symbol, lhs_symbol) == 0;
}

bool head_threads_over(const char* h) {
    return strcmp(h, "List") == 0 ||
           strcmp(h, "Equal") == 0 ||
           strcmp(h, "Unequal") == 0 ||
           strcmp(h, "Less") == 0 ||
           strcmp(h, "LessEqual") == 0 ||
           strcmp(h, "Greater") == 0 ||
           strcmp(h, "GreaterEqual") == 0 ||
           strcmp(h, "And") == 0 ||
           strcmp(h, "Or") == 0 ||
           strcmp(h, "Not") == 0 ||
           strcmp(h, "Xor") == 0 ||
           strcmp(h, "Implies") == 0;
}

/* ----------------------------------------------------------------------- */
/* Candidate set                                                           */
/* ----------------------------------------------------------------------- */

/* SIMP_CAND_CAP, SIMP_ROUNDS, and the CandSet struct all live in
 * simp_internal.h since simp_search.c also instantiates a CandSet. */

void cs_init(CandSet* cs) {
    cs->items = NULL;
    cs->count = 0;
    cs->capacity = 0;
}

void cs_free(CandSet* cs) {
    for (size_t i = 0; i < cs->count; i++) expr_free(cs->items[i]);
    free(cs->items);
    cs->items = NULL;
    cs->count = 0;
    cs->capacity = 0;
}

bool cs_contains(const CandSet* cs, const Expr* e) {
    for (size_t i = 0; i < cs->count; i++) {
        if (expr_eq(cs->items[i], e)) return true;
    }
    return false;
}

/* Take ownership of `e`; free if duplicate or set is full. */
void cs_add_or_free(CandSet* cs, Expr* e) {
    if (!e) return;
    if (cs->count >= SIMP_CAND_CAP || cs_contains(cs, e)) {
        expr_free(e);
        return;
    }
    if (cs->count >= cs->capacity) {
        size_t new_cap = cs->capacity ? cs->capacity * 2 : 4;
        Expr** np = (Expr**)realloc(cs->items, new_cap * sizeof(Expr*));
        if (!np) { expr_free(e); return; }
        cs->items = np;
        cs->capacity = new_cap;
    }
    cs->items[cs->count++] = e;
}

/* ----------------------------------------------------------------------- */
/* Scoring                                                                 */
/* ----------------------------------------------------------------------- */

#define SIMP_SCORE_INF ((size_t)-1)

/* True iff a `Power[X, Rational[p, q]]` subtree with q >= 2 (a non-trivial
 * root) appears anywhere in e. Used by the nested-radical detector below and
 * by builtin_simplify's algebraic-top Together fast path -- it is the precise
 * "genuine radical" gate, excluding symbolic exponents (a^x) and real-power
 * exponents (x^1.5) that has_non_integer_power would also accept (and that
 * would otherwise be fed to Together, which can hang on symbolic powers). */
bool simp_has_rational_root(const Expr* e) {
    if (!e || e->type != EXPR_FUNCTION) return false;
    const Expr* head = e->data.function.head;
    size_t argc = e->data.function.arg_count;
    if (head && head->type == EXPR_SYMBOL && head->data.symbol == SYM_Power
        && argc == 2) {
        const Expr* exp = e->data.function.args[1];
        if (exp->type == EXPR_FUNCTION && exp->data.function.head
            && exp->data.function.head->type == EXPR_SYMBOL
            && exp->data.function.head->data.symbol == SYM_Rational
            && exp->data.function.arg_count == 2) {
            const Expr* qq = exp->data.function.args[1];
            if (qq->type == EXPR_INTEGER && qq->data.integer >= 2) return true;
        }
    }
    for (size_t i = 0; i < argc; i++) {
        if (simp_has_rational_root(e->data.function.args[i])) return true;
    }
    return false;
}

/* Penalty for *truly* nested radicals: a `Power[Compound, Rational[_, q]]`
 * (q >= 2) whose compound base itself contains another root. This lets
 * Simplify prefer denested forms like (1+Sqrt[5])/2 over the original
 * (2+Sqrt[5])^(1/3) -- their plain LeafCounts are tied or off by ~2,
 * but the denested form has no nested-radical structure, which is the
 * canonical preferred shape.
 *
 * Surcharge value (3) chosen as the smallest constant that lets the
 * cube-root denester's output win the user's (2+Sqrt[5])^(1/3) case
 * (LeafCount diff is 2) without dominating other transforms' choices.
 * Flat radicals like Sqrt[5] or Sqrt[x+y] are NOT penalised -- only
 * the truly-nested shape is. */
static size_t nested_radical_penalty(const Expr* e) {
    if (!e || e->type != EXPR_FUNCTION) return 0;
    size_t total = 0;
    const Expr* head = e->data.function.head;
    size_t argc = e->data.function.arg_count;
    if (head && head->type == EXPR_SYMBOL && head->data.symbol == SYM_Power
        && argc == 2) {
        const Expr* base = e->data.function.args[0];
        const Expr* exp  = e->data.function.args[1];
        if (base->type == EXPR_FUNCTION && base->data.function.head
            && base->data.function.head->type == EXPR_SYMBOL) {
            const char* bsym = base->data.function.head->data.symbol;
            if (bsym != SYM_Rational && bsym != SYM_Complex) {
                if (exp->type == EXPR_FUNCTION && exp->data.function.head
                    && exp->data.function.head->type == EXPR_SYMBOL
                    && exp->data.function.head->data.symbol == SYM_Rational
                    && exp->data.function.arg_count == 2) {
                    const Expr* qq = exp->data.function.args[1];
                    if (qq->type == EXPR_INTEGER && qq->data.integer >= 2
                        && simp_has_rational_root(base)) {
                        total += 3;
                    }
                }
            }
        }
    }
    for (size_t i = 0; i < argc; i++) {
        total += nested_radical_penalty(e->data.function.args[i]);
    }
    return total;
}

size_t score_with_func(const Expr* e, const Expr* complexity_func) {
    if (!complexity_func) {
        return simp_default_complexity(e) + nested_radical_penalty(e);
    }
    Expr* a[1] = { expr_copy((Expr*)e) };
    Expr* call = expr_new_function(expr_copy((Expr*)complexity_func), a, 1);
    Expr* result = evaluate(call);
    expr_free(call);
    size_t s;
    if (result->type == EXPR_INTEGER) {
        s = (result->data.integer < 0) ? 0 : (size_t)result->data.integer;
    } else if (result->type == EXPR_BIGINT) {
        s = SIMP_SCORE_INF;
    } else {
        s = simp_default_complexity(e);
    }
    expr_free(result);
    return s;
}

