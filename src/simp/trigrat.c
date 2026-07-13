/* `strdup` is POSIX, not C99.  glibc hides it under `-std=c99` unless
 * one of the standard feature-test macros is defined.  Darwin exposes
 * it implicitly, which masks the cross-platform issue; the rest of
 * Mathilda matches this by setting `_GNU_SOURCE` before any include
 * (see src/expr.c, src/match.c, src/modular.c).  Must precede every
 * other directive. */
#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif

/*
 * trigrat.c -- Fast algebraic normal form for rational functions of
 * trigonometric and hyperbolic functions.
 *
 * See trigrat.h for the contract and design overview.
 *
 * Pipeline (single static-storage entry point: simp_trig_rational):
 *
 *   1. Preprocess Tan/Cot/Sec/Csc/Tanh/Coth/Sech/Csch to Sin/Cos and
 *      Sinh/Cosh.
 *   2. Collect distinct (Sin|Cos, arg) pairs and (Sinh|Cosh, arg) pairs
 *      via a top-down walk that does NOT descend into a trig/hyp node
 *      once recorded (so Sin[Cos[x]] is opaque if there is no separate
 *      occurrence of Cos[x] outside it). Budget on count.
 *   3. Substitute each (head, arg) pair with fresh ground-field symbols
 *      s_i, c_i (or sigma_j, chi_j for hyperbolic).
 *   4. Together to combine into a single P/Q.
 *   5. Reduce P and Q modulo the trig ideal:
 *        s_i^2  -> 1 - c_i^2          (each trig pair)
 *        sigma_j^2 -> chi_j^2 - 1     (each hyp pair)
 *      Iterated with Expand until each s_i and sigma_j has degree <= 1.
 *   6. For each odd generator still in Q: rationalise by multiplying
 *      numerator and denominator by the conjugate (substitute generator
 *      with its negation in Q). Then re-reduce. After this loop, Q has
 *      no s_i and no sigma_j.
 *   7. Cancel the resulting rational function via the existing Cancel
 *      builtin -- this does the multivariate polynomial GCD in
 *      K[c_i, chi_j, other_vars].
 *   8. Substitute fresh symbols back to Sin/Cos/Sinh/Cosh and re-evaluate.
 *   9. Optional leaf-count-gated fold-back to Tan/Cot/Sec/Csc/Tanh/Coth/
 *      Sech/Csch.
 *  10. Final leaf-count gate: return only if STRICTLY smaller than the
 *      input. Otherwise return NULL.
 */

#include "trigrat.h"

#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "eval.h"
#include "expr.h"
#include "simp.h"
#include "sym_names.h"

/* ----------------------------------------------------------------------- */
/* Tunables                                                                */
/* ----------------------------------------------------------------------- */

#define TRIGRAT_BUDGET     4   /* max number of distinct (trig+hyp) args */
#define TRIGRAT_OPAQUE_MAX 24  /* max number of distinct opaque atoms    */
#define TRIGRAT_MAX_REDUCE 64  /* hard iteration cap for ideal reduction */
#define TRIGRAT_RADICAL_MAX 3  /* max number of distinct Sqrt[g] kernels */
/*
 * Below this leaf-count we always defer to simp_bottomup. simp_search's
 * heuristic transforms (TrigReduce, TrigFactor, HalfAngle, TanAddition,
 * etc.) produce more compact forms on small inputs by recognising
 * double-angle / half-angle / angle-addition identities, while the
 * algebraic normal form here is correctness-equivalent but typically
 * has the SAME leaf count -- so the strict-win gate downstream would
 * reject the win-by-one cases that the test corpus exercises. The
 * fast path's reason for existence is to rescue the cases where
 * simp_bottomup is impractical (hundreds of leaves with many opaque
 * Log/Exp atoms produced by integrators); for those, the floor is
 * easily exceeded.
 */
#define TRIGRAT_LEAF_FLOOR 60

/* ----------------------------------------------------------------------- */
/* Tiny helpers                                                            */
/* ----------------------------------------------------------------------- */

/* Build f[arg], evaluate, return the result. Takes ownership of arg. */
static Expr* tr_call_unary(const char* head_name, Expr* arg) {
    Expr* a[1] = { arg };
    Expr* call = expr_new_function(expr_new_symbol(head_name), a, 1);
    Expr* r = evaluate(call);
    expr_free(call);
    return r;
}

static Expr* tr_call_unary_copy(const char* head_name, const Expr* arg) {
    return tr_call_unary(head_name, expr_copy((Expr*)arg));
}

/* Check whether any descendant carries the given literal symbol name. */
static bool tr_has_symbol_name(const Expr* e, const char* name) {
    if (!e) return false;
    if (e->type == EXPR_SYMBOL) {
        return e->data.symbol.name && strcmp(e->data.symbol.name, name) == 0;
    }
    if (e->type != EXPR_FUNCTION) return false;
    if (tr_has_symbol_name(e->data.function.head, name)) return true;
    for (size_t i = 0; i < e->data.function.arg_count; i++) {
        if (tr_has_symbol_name(e->data.function.args[i], name)) return true;
    }
    return false;
}

/* Return a literal integer Expr. */
static Expr* tr_int(int64_t v) {
    return expr_new_integer(v);
}

/* Return Plus[a, b]; takes ownership of both. Evaluates. */
static Expr* tr_plus(Expr* a, Expr* b) {
    Expr* args[2] = { a, b };
    Expr* call = expr_new_function(expr_new_symbol(SYM_Plus), args, 2);
    Expr* r = evaluate(call);
    expr_free(call);
    return r;
}

/* Return Times[a, b]; takes ownership of both. Evaluates. */
static Expr* tr_times(Expr* a, Expr* b) {
    Expr* args[2] = { a, b };
    Expr* call = expr_new_function(expr_new_symbol(SYM_Times), args, 2);
    Expr* r = evaluate(call);
    expr_free(call);
    return r;
}

/* Return Power[a, b]; takes ownership of both. Evaluates. */
static Expr* tr_pow(Expr* a, Expr* b) {
    Expr* args[2] = { a, b };
    Expr* call = expr_new_function(expr_new_symbol(SYM_Power), args, 2);
    Expr* r = evaluate(call);
    expr_free(call);
    return r;
}

/* Return -e; takes ownership of e. */
static Expr* tr_neg(Expr* e) {
    return tr_times(tr_int(-1), e);
}

/* Raw (non-evaluating) Times/Plus/Power constructors used during the
 * top-down preprocessing rewrite (Tan -> Sin/Cos, etc.). Using the
 * evaluating versions would let the auto-canonicalisation rule
 * Sin[x]/Cos[x] -> Tan[x] fire and reproduce the very form we are
 * trying to eliminate, causing infinite recursion. The resulting raw
 * trees are immediately passed back through tr_walk_subst which
 * substitutes fresh symbols for Sin/Cos before any evaluator pass. */
static Expr* tr_raw_times(Expr* a, Expr* b) {
    Expr** args = (Expr**)calloc(2, sizeof(Expr*));
    args[0] = a; args[1] = b;
    /* expr_new_function copies the args array (it does not adopt the
     * buffer), so we own and must free it. */
    Expr* fn = expr_new_function(expr_new_symbol(SYM_Times), args, 2);
    free(args);
    return fn;
}

static Expr* tr_raw_pow(Expr* a, Expr* b) {
    Expr** args = (Expr**)calloc(2, sizeof(Expr*));
    args[0] = a; args[1] = b;
    Expr* fn = expr_new_function(expr_new_symbol(SYM_Power), args, 2);
    free(args);
    return fn;
}

/* ----------------------------------------------------------------------- */
/* Per-call binding table                                                  */
/* ----------------------------------------------------------------------- */

typedef struct {
    Expr*  arg;        /* the trig/hyp argument expression (owned)      */
    char*  s_name;     /* fresh symbol name for the "odd" generator     */
    char*  c_name;     /* fresh symbol name for the "even" generator    */
    bool   is_hyp;     /* false: trig (s^2+c^2=1); true: hyp (c^2-s^2=1)*/
} TRBind;

/* Opaque (transcendental) atom: a subtree that is not Plus/Times/Power/
 * Sin/Cos/Sinh/Cosh and cannot be broken down by rational-function
 * operations. Examples: Log[1+s+c], Exp[s], Sqrt[s+1], Erf[c], unknown
 * f[a,b,c]. Substituted with a fresh symbol so the algebraic core sees
 * a pure rational function in fresh symbols. */
typedef struct {
    Expr* subtree;   /* original subtree (owned)                        */
    char* name;      /* fresh symbol name                                */
} TROpaque;

/* Quadratic radical generator: a kernel Sqrt[g] (i.e. Power[g, p/2]).
 * Treated as an algebraic element l with the defining relation
 *   l^2 = radicand_subst   (the substituted value of g, rational in the
 *                           trig generators s_i, c_i and ground symbols).
 * This is structurally identical to a trig "odd" generator s_i, whose
 * relation is s_i^2 = 1 - c_i^2: in both cases the square of the generator
 * reduces to a known lower expression, so the existing ideal-reduction and
 * conjugate-rationalisation machinery applies verbatim once it learns the
 * generator's name and reduction RHS. All half-integer powers of the same g
 * (Sqrt[g], g^(3/2), 1/Sqrt[g], ...) collapse to integer powers of l. */
typedef struct {
    char* name;            /* fresh symbol name for l                    */
    Expr* radicand_subst;  /* l^2 = this (rational in s_i, c_i, ...) own */
    Expr* orig;            /* original Sqrt[g] subtree, for back-sub own */
} TRRadical;

typedef struct {
    TRBind*   items;
    size_t    count;
    size_t    capacity;
    TROpaque* opaques;
    size_t    opaque_count;
    size_t    opaque_capacity;
    TRRadical* radicals;
    size_t     radical_count;
    size_t     radical_capacity;
} TRBindings;

static void tr_bindings_init(TRBindings* b) {
    b->items = NULL; b->count = 0; b->capacity = 0;
    b->opaques = NULL; b->opaque_count = 0; b->opaque_capacity = 0;
    b->radicals = NULL; b->radical_count = 0; b->radical_capacity = 0;
}

static void tr_bindings_free(TRBindings* b) {
    if (!b) return;
    for (size_t i = 0; i < b->count; i++) {
        if (b->items[i].arg) expr_free(b->items[i].arg);
        free(b->items[i].s_name);
        free(b->items[i].c_name);
    }
    free(b->items);
    for (size_t i = 0; i < b->opaque_count; i++) {
        if (b->opaques[i].subtree) expr_free(b->opaques[i].subtree);
        free(b->opaques[i].name);
    }
    free(b->opaques);
    for (size_t i = 0; i < b->radical_count; i++) {
        if (b->radicals[i].radicand_subst) expr_free(b->radicals[i].radicand_subst);
        if (b->radicals[i].orig) expr_free(b->radicals[i].orig);
        free(b->radicals[i].name);
    }
    free(b->radicals);
    b->items = NULL; b->count = 0; b->capacity = 0;
    b->opaques = NULL; b->opaque_count = 0; b->opaque_capacity = 0;
    b->radicals = NULL; b->radical_count = 0; b->radical_capacity = 0;
}

/* Find a radical generator by its original Sqrt[g] subtree. */
static TRRadical* tr_radical_find(const TRBindings* b, const Expr* orig) {
    for (size_t i = 0; i < b->radical_count; i++) {
        if (expr_eq(b->radicals[i].orig, (Expr*)orig)) return &b->radicals[i];
    }
    return NULL;
}

/* True if `e` mentions any radical generator's fresh symbol name. */
static bool tr_uses_any_radical(const Expr* e, const TRBindings* b) {
    for (size_t i = 0; i < b->radical_count; i++) {
        if (tr_has_symbol_name(e, b->radicals[i].name)) return true;
    }
    return false;
}

/* Add a radical generator. Takes ownership of `radicand` (the l^2 value);
 * copies `orig`. Returns NULL on allocation failure. */
static TRRadical* tr_radical_add(TRBindings* b, const Expr* avoid,
                                 const Expr* orig, Expr* radicand,
                                 int* counter) {
    if (b->radical_count == b->radical_capacity) {
        size_t nc = b->radical_capacity == 0 ? 4 : b->radical_capacity * 2;
        TRRadical* nb = (TRRadical*)realloc(b->radicals, nc * sizeof(TRRadical));
        if (!nb) return NULL;
        b->radicals = nb;
        b->radical_capacity = nc;
    }
    char buf[64];
    int tries = 0;
    while (tries < 1000) {
        snprintf(buf, sizeof(buf), "$pc_rad_%d$", *counter);
        if (!tr_has_symbol_name(avoid, buf)) break;
        (*counter)++;
        tries++;
    }
    TRRadical* it = &b->radicals[b->radical_count++];
    it->name = strdup(buf);
    it->radicand_subst = radicand;
    it->orig = expr_copy((Expr*)orig);
    (*counter)++;
    return it;
}

static TRBind* tr_bindings_find(const TRBindings* b, const Expr* arg, bool is_hyp) {
    for (size_t i = 0; i < b->count; i++) {
        if (b->items[i].is_hyp == is_hyp && expr_eq(b->items[i].arg, arg)) {
            return &b->items[i];
        }
    }
    return NULL;
}

static TROpaque* tr_opaque_find(const TRBindings* b, const Expr* sub) {
    for (size_t i = 0; i < b->opaque_count; i++) {
        if (expr_eq(b->opaques[i].subtree, sub)) {
            return &b->opaques[i];
        }
    }
    return NULL;
}

static TROpaque* tr_opaque_add(TRBindings* b, const Expr* avoid,
                               const Expr* sub, int* counter) {
    if (b->opaque_count == b->opaque_capacity) {
        size_t nc = b->opaque_capacity == 0 ? 4 : b->opaque_capacity * 2;
        TROpaque* nb = (TROpaque*)realloc(b->opaques, nc * sizeof(TROpaque));
        if (!nb) return NULL;
        b->opaques = nb;
        b->opaque_capacity = nc;
    }
    char buf[64];
    int tries = 0;
    while (tries < 1000) {
        snprintf(buf, sizeof(buf), "$pc_opq_%d$", *counter);
        if (!tr_has_symbol_name(avoid, buf)) break;
        (*counter)++;
        tries++;
    }
    TROpaque* it = &b->opaques[b->opaque_count++];
    it->subtree = expr_copy((Expr*)sub);
    it->name = strdup(buf);
    (*counter)++;
    return it;
}

/* Allocate fresh s/c symbol names that do not clash with any symbol in
 * `avoid`. Names are stable across the call: trig pairs are
 *   $pc_trig_s_<N>$, $pc_trig_c_<N>$
 * hyp pairs are
 *   $pc_hyp_s_<N>$,  $pc_hyp_c_<N>$
 * with N a fresh counter per binding. The starting counter offset is
 * bumped if a generated name is already present in `avoid`. */
static TRBind* tr_bindings_add(TRBindings* b,
                               const Expr* avoid,
                               const Expr* arg,
                               bool is_hyp,
                               int* counter) {
    /* Grow */
    if (b->count == b->capacity) {
        size_t nc = b->capacity == 0 ? 4 : b->capacity * 2;
        TRBind* nb = (TRBind*)realloc(b->items, nc * sizeof(TRBind));
        if (!nb) return NULL;
        b->items = nb;
        b->capacity = nc;
    }
    char s_buf[64], c_buf[64];
    const char* fam_s = is_hyp ? "$pc_hyp_s_" : "$pc_trig_s_";
    const char* fam_c = is_hyp ? "$pc_hyp_c_" : "$pc_trig_c_";
    int tries = 0;
    while (tries < 1000) {
        snprintf(s_buf, sizeof(s_buf), "%s%d$", fam_s, *counter);
        snprintf(c_buf, sizeof(c_buf), "%s%d$", fam_c, *counter);
        if (!tr_has_symbol_name(avoid, s_buf) &&
            !tr_has_symbol_name(avoid, c_buf)) {
            break;
        }
        (*counter)++;
        tries++;
    }
    TRBind* it = &b->items[b->count++];
    it->arg = expr_copy((Expr*)arg);
    it->s_name = strdup(s_buf);
    it->c_name = strdup(c_buf);
    it->is_hyp = is_hyp;
    (*counter)++;
    return it;
}

/* ----------------------------------------------------------------------- */
/* Single-pass top-down walk: collect bindings AND produce a substituted   */
/* expression. Walks transparently through Plus/Times/Power[_, integer]    */
/* and Rational/Complex. At Tan/Cot/Sec/Csc/Tanh/Coth/Sech/Csch -- rewrites */
/* to a Sin/Cos (Sinh/Cosh) form and recurses on the rewrite. At           */
/* Sin/Cos/Sinh/Cosh with one arg -- substitutes with a fresh trig/hyp     */
/* generator symbol. Every other function call is recorded as an opaque    */
/* atom (its entire subtree stored verbatim) and substituted with a fresh  */
/* ground-field symbol. The original (un-substituted) subtree is what gets */
/* restored at back-substitution.                                          */
/* ----------------------------------------------------------------------- */

static Expr* tr_walk_subst(const Expr* e, TRBindings* b, int* counter,
                           const Expr* avoid, bool* over_budget) {
    if (!e || *over_budget) return e ? expr_copy((Expr*)e) : NULL;
    if (e->type != EXPR_FUNCTION) return expr_copy((Expr*)e);

    const Expr* h = e->data.function.head;
    size_t n = e->data.function.arg_count;
    const char* hs = (h && h->type == EXPR_SYMBOL) ? h->data.symbol.name : NULL;

    /* Transparent additive/multiplicative heads -- recurse into children. */
    if (hs == SYM_Plus || hs == SYM_Times) {
        Expr** args_new = (Expr**)calloc(n ? n : 1, sizeof(Expr*));
        for (size_t i = 0; i < n; i++) {
            args_new[i] = tr_walk_subst(e->data.function.args[i], b, counter,
                                        avoid, over_budget);
        }
        /* expr_new_function copies the args array; we own and must free it. */
        Expr* fn = expr_new_function(expr_new_symbol(hs), args_new, n);
        free(args_new);
        return fn;
    }

    /* Rational[n,d] / Complex[a,b] -- treat the wrapper as transparent
     * (number-like literal); just keep a deep copy. The integer args are
     * already canonical leaves. */
    if (hs == SYM_Rational || hs == SYM_Complex) {
        return expr_copy((Expr*)e);
    }

    /* Power[base, integer_exp] -- transparent in base. */
    if (hs == SYM_Power && n == 2 &&
        e->data.function.args[1] &&
        e->data.function.args[1]->type == EXPR_INTEGER) {
        Expr* base_new = tr_walk_subst(e->data.function.args[0], b, counter,
                                       avoid, over_budget);
        Expr* exp_new  = expr_copy(e->data.function.args[1]);
        Expr** args_new = (Expr**)calloc(2, sizeof(Expr*));
        args_new[0] = base_new;
        args_new[1] = exp_new;
        /* expr_new_function copies the args array; we own and must free it. */
        Expr* fn = expr_new_function(expr_new_symbol(hs), args_new, 2);
        free(args_new);
        return fn;
    }

    /* Power[base, Rational[p, 2]] -- a quadratic radical.  Introduce (or
     * reuse) a radical generator l = Sqrt[base] with l^2 = base-substituted,
     * and rewrite the whole power as the integer power l^p (p odd).  All
     * half-integer powers of the same base thereby collapse to powers of one
     * generator, and the algebraic relation l^2 = base ties the radical to
     * any plain (integer-power) occurrences of base elsewhere. */
    if (hs == SYM_Power && n == 2 && e->data.function.args[1] &&
        e->data.function.args[1]->type == EXPR_FUNCTION) {
        const Expr* exp_e = e->data.function.args[1];
        const Expr* eh = exp_e->data.function.head;
        if (eh && eh->type == EXPR_SYMBOL && eh->data.symbol.name == SYM_Rational &&
            exp_e->data.function.arg_count == 2 &&
            exp_e->data.function.args[0] &&
            exp_e->data.function.args[0]->type == EXPR_INTEGER &&
            exp_e->data.function.args[1] &&
            exp_e->data.function.args[1]->type == EXPR_INTEGER &&
            exp_e->data.function.args[1]->data.integer == 2) {
            int64_t p = exp_e->data.function.args[0]->data.integer;
            const Expr* base = e->data.function.args[0];

            /* Canonical key/back-sub form Sqrt[base] = Power[base, 1/2]. */
            Expr* half = expr_new_function(expr_new_symbol(SYM_Rational),
                            (Expr*[]){ tr_int(1), tr_int(2) }, 2);
            Expr* sqrt_key = expr_new_function(expr_new_symbol(SYM_Power),
                            (Expr*[]){ expr_copy((Expr*)base), half }, 2);

            TRRadical* rb = tr_radical_find(b, sqrt_key);
            if (!rb) {
                if (b->radical_count >= TRIGRAT_RADICAL_MAX) {
                    *over_budget = true;
                    expr_free(sqrt_key);
                    return expr_copy((Expr*)e);
                }
                /* Substitute the radicand in the SAME binding set so its
                 * trig generators are shared with the rest of the input. */
                Expr* g_sub = tr_walk_subst(base, b, counter, avoid, over_budget);
                /* Nested radicals (Sqrt of an expression that itself contains
                 * a radical) are not supported by the single-conjugate
                 * rationalisation: bail cleanly. */
                if (*over_budget || !g_sub || tr_uses_any_radical(g_sub, b)) {
                    *over_budget = true;
                    if (g_sub) expr_free(g_sub);
                    expr_free(sqrt_key);
                    return expr_copy((Expr*)e);
                }
                rb = tr_radical_add(b, avoid, sqrt_key, g_sub, counter);
            }
            expr_free(sqrt_key);
            if (!rb) { *over_budget = true; return expr_copy((Expr*)e); }
            if (p == 1) return expr_new_symbol(rb->name);
            return tr_raw_pow(expr_new_symbol(rb->name), tr_int(p));
        }
    }

    /* Tan/Cot/Sec/Csc and hyperbolic analogues -- rewrite to Sin/Cos
     * form, then recurse. Uses RAW (non-evaluating) constructors to
     * avoid the evaluator's Sin[x]/Cos[x] -> Tan[x] canonicalisation,
     * which would otherwise reproduce the input and infinite-loop. The
     * substitution at the recursive call replaces Sin/Cos with fresh
     * symbols before any subsequent evaluator pass. */
    if (hs && n == 1) {
        const Expr* x = e->data.function.args[0];
        Expr* rewritten = NULL;
        if (hs == SYM_Tan) {
            Expr* sn = expr_new_function(expr_new_symbol(SYM_Sin),
                                         (Expr*[]){expr_copy((Expr*)x)}, 1);
            Expr* cn = expr_new_function(expr_new_symbol(SYM_Cos),
                                         (Expr*[]){expr_copy((Expr*)x)}, 1);
            Expr* inv = tr_raw_pow(cn, tr_int(-1));
            rewritten = tr_raw_times(sn, inv);
        } else if (hs == SYM_Cot) {
            Expr* cn = expr_new_function(expr_new_symbol(SYM_Cos),
                                         (Expr*[]){expr_copy((Expr*)x)}, 1);
            Expr* sn = expr_new_function(expr_new_symbol(SYM_Sin),
                                         (Expr*[]){expr_copy((Expr*)x)}, 1);
            Expr* inv = tr_raw_pow(sn, tr_int(-1));
            rewritten = tr_raw_times(cn, inv);
        } else if (hs == SYM_Sec) {
            Expr* cn = expr_new_function(expr_new_symbol(SYM_Cos),
                                         (Expr*[]){expr_copy((Expr*)x)}, 1);
            rewritten = tr_raw_pow(cn, tr_int(-1));
        } else if (hs == SYM_Csc) {
            Expr* sn = expr_new_function(expr_new_symbol(SYM_Sin),
                                         (Expr*[]){expr_copy((Expr*)x)}, 1);
            rewritten = tr_raw_pow(sn, tr_int(-1));
        } else if (hs == SYM_Tanh) {
            Expr* sn = expr_new_function(expr_new_symbol(SYM_Sinh),
                                         (Expr*[]){expr_copy((Expr*)x)}, 1);
            Expr* cn = expr_new_function(expr_new_symbol(SYM_Cosh),
                                         (Expr*[]){expr_copy((Expr*)x)}, 1);
            Expr* inv = tr_raw_pow(cn, tr_int(-1));
            rewritten = tr_raw_times(sn, inv);
        } else if (hs == SYM_Coth) {
            Expr* cn = expr_new_function(expr_new_symbol(SYM_Cosh),
                                         (Expr*[]){expr_copy((Expr*)x)}, 1);
            Expr* sn = expr_new_function(expr_new_symbol(SYM_Sinh),
                                         (Expr*[]){expr_copy((Expr*)x)}, 1);
            Expr* inv = tr_raw_pow(sn, tr_int(-1));
            rewritten = tr_raw_times(cn, inv);
        } else if (hs == SYM_Sech) {
            Expr* cn = expr_new_function(expr_new_symbol(SYM_Cosh),
                                         (Expr*[]){expr_copy((Expr*)x)}, 1);
            rewritten = tr_raw_pow(cn, tr_int(-1));
        } else if (hs == SYM_Csch) {
            Expr* sn = expr_new_function(expr_new_symbol(SYM_Sinh),
                                         (Expr*[]){expr_copy((Expr*)x)}, 1);
            rewritten = tr_raw_pow(sn, tr_int(-1));
        }
        if (rewritten) {
            Expr* r = tr_walk_subst(rewritten, b, counter, avoid, over_budget);
            expr_free(rewritten);
            return r;
        }
    }

    /* Sin/Cos/Sinh/Cosh single-arg -- trig/hyp pair. */
    if ((hs == SYM_Sin || hs == SYM_Cos ||
         hs == SYM_Sinh || hs == SYM_Cosh) && n == 1) {
        bool is_hyp = (hs == SYM_Sinh || hs == SYM_Cosh);
        const Expr* arg = e->data.function.args[0];
        TRBind* it = tr_bindings_find(b, arg, is_hyp);
        if (!it) {
            if (b->count >= TRIGRAT_BUDGET) {
                *over_budget = true;
                return expr_copy((Expr*)e);
            }
            it = tr_bindings_add(b, avoid, arg, is_hyp, counter);
        }
        bool want_s = (hs == SYM_Sin || hs == SYM_Sinh);
        return expr_new_symbol(want_s ? it->s_name : it->c_name);
    }

    /* Everything else: opaque atom. Store the original subtree, return a
     * fresh symbol. */
    TROpaque* op = tr_opaque_find(b, e);
    if (!op) {
        if (b->opaque_count >= TRIGRAT_OPAQUE_MAX) {
            *over_budget = true;
            return expr_copy((Expr*)e);
        }
        op = tr_opaque_add(b, avoid, e, counter);
    }
    return expr_new_symbol(op->name);
}

/* ----------------------------------------------------------------------- */
/* Backward substitution: replace fresh symbols with Sin[a_i] etc.         */
/* ----------------------------------------------------------------------- */

static Expr* tr_subst_back(const Expr* e, const TRBindings* b) {
    if (!e) return NULL;
    if (e->type == EXPR_SYMBOL && e->data.symbol.name) {
        for (size_t i = 0; i < b->count; i++) {
            if (strcmp(e->data.symbol.name, b->items[i].s_name) == 0) {
                const char* head = b->items[i].is_hyp ? SYM_Sinh : SYM_Sin;
                return expr_new_function(expr_new_symbol(head),
                    (Expr*[]){ expr_copy(b->items[i].arg) }, 1);
            }
            if (strcmp(e->data.symbol.name, b->items[i].c_name) == 0) {
                const char* head = b->items[i].is_hyp ? SYM_Cosh : SYM_Cos;
                return expr_new_function(expr_new_symbol(head),
                    (Expr*[]){ expr_copy(b->items[i].arg) }, 1);
            }
        }
        for (size_t i = 0; i < b->opaque_count; i++) {
            if (strcmp(e->data.symbol.name, b->opaques[i].name) == 0) {
                return expr_copy(b->opaques[i].subtree);
            }
        }
        for (size_t i = 0; i < b->radical_count; i++) {
            if (strcmp(e->data.symbol.name, b->radicals[i].name) == 0) {
                return expr_copy(b->radicals[i].orig);
            }
        }
        return expr_copy((Expr*)e);
    }
    if (e->type != EXPR_FUNCTION) return expr_copy((Expr*)e);
    Expr* head_new = tr_subst_back(e->data.function.head, b);
    size_t n = e->data.function.arg_count;
    Expr** args_new = (Expr**)calloc(n ? n : 1, sizeof(Expr*));
    for (size_t i = 0; i < n; i++) {
        args_new[i] = tr_subst_back(e->data.function.args[i], b);
    }
    /* expr_new_function copies the args array; we own and must free it. */
    Expr* fn = expr_new_function(head_new, args_new, n);
    free(args_new);
    return fn;
}

/* ----------------------------------------------------------------------- */
/* Reduce a polynomial modulo the trig ideals.                             */
/*                                                                         */
/* For every Power[symbol_named_s_i, n] with n>=2:                         */
/*   trig pair: rewrite as Power[s_i, n-2] * (1 - Power[c_i, 2])           */
/*   hyp  pair: rewrite as Power[sigma_j, n-2] * (Power[chi_j, 2] - 1)     */
/* (We also handle the rare case where s_i appears directly squared via    */
/*  Times[s_i, s_i] without a Power node: a final Expand call normalises   */
/*  these into Power before we walk.)                                      */
/*                                                                         */
/* Returns a freshly allocated, Expand'd, reduced expression.              */
/* ----------------------------------------------------------------------- */

/* Look up a binding whose s_name matches the given symbol name. */
static const TRBind* tr_find_by_s_name(const TRBindings* b, const char* sym) {
    if (!sym) return NULL;
    for (size_t i = 0; i < b->count; i++) {
        if (strcmp(b->items[i].s_name, sym) == 0) return &b->items[i];
    }
    return NULL;
}

/* Look up a radical generator by its fresh symbol name. */
static const TRRadical* tr_find_radical_by_name(const TRBindings* b,
                                                const char* sym) {
    if (!sym) return NULL;
    for (size_t i = 0; i < b->radical_count; i++) {
        if (strcmp(b->radicals[i].name, sym) == 0) return &b->radicals[i];
    }
    return NULL;
}

/* Build the replacement for a radical generator l^n (n >= 2):
 *   l^(n-2) * radicand   (since l^2 = radicand). */
static Expr* tr_build_radical_rhs(const TRRadical* rb, int64_t n) {
    Expr* base_pow = (n == 2)
        ? tr_int(1)
        : tr_pow(expr_new_symbol(rb->name), tr_int(n - 2));
    return tr_times(base_pow, expr_copy(rb->radicand_subst));
}

/* Build the replacement expression for a generator s_i^n with n >= 2. */
static Expr* tr_build_reduce_rhs(const TRBind* it, int64_t n) {
    /* s_i^(n-2) * (1 - c_i^2)        for trig
     * sigma^(n-2) * (chi^2 - 1)      for hyp                            */
    Expr* base_pow;
    if (n == 2) {
        base_pow = tr_int(1);
    } else {
        base_pow = tr_pow(expr_new_symbol(it->s_name), tr_int(n - 2));
    }
    Expr* c_sq = tr_pow(expr_new_symbol(it->c_name), tr_int(2));
    Expr* factor;
    if (it->is_hyp) {
        /* chi^2 - 1 */
        factor = tr_plus(c_sq, tr_int(-1));
    } else {
        /* 1 - c^2 */
        factor = tr_plus(tr_int(1), tr_neg(c_sq));
    }
    return tr_times(base_pow, factor);
}

/* Walk the tree and rewrite Power[s_i, n] with n>=2 once. Returns a
 * fresh tree; sets *changed if any rewrite happened. Does NOT recurse
 * into a rewritten node (the caller will Expand+iterate). */
static Expr* tr_rewrite_pass(const Expr* e, const TRBindings* b, bool* changed) {
    if (!e) return NULL;
    if (e->type == EXPR_FUNCTION) {
        const Expr* h = e->data.function.head;
        /* Power[s_i, n] with n >= 2  -> rewrite */
        if (h && h->type == EXPR_SYMBOL && h->data.symbol.name == SYM_Power &&
            e->data.function.arg_count == 2) {
            const Expr* base = e->data.function.args[0];
            const Expr* exp_e = e->data.function.args[1];
            if (base && base->type == EXPR_SYMBOL &&
                exp_e && exp_e->type == EXPR_INTEGER &&
                exp_e->data.integer >= 2) {
                const TRBind* it = tr_find_by_s_name(b, base->data.symbol.name);
                if (it) {
                    *changed = true;
                    return tr_build_reduce_rhs(it, exp_e->data.integer);
                }
                const TRRadical* rb =
                    tr_find_radical_by_name(b, base->data.symbol.name);
                if (rb) {
                    *changed = true;
                    return tr_build_radical_rhs(rb, exp_e->data.integer);
                }
            }
        }
        /* Recurse into head and children */
        Expr* head_new = tr_rewrite_pass(h, b, changed);
        size_t n = e->data.function.arg_count;
        Expr** args_new = (Expr**)calloc(n ? n : 1, sizeof(Expr*));
        for (size_t i = 0; i < n; i++) {
            args_new[i] = tr_rewrite_pass(e->data.function.args[i], b, changed);
        }
        /* expr_new_function copies the args array; we own and must free it. */
        Expr* fn = expr_new_function(head_new, args_new, n);
        free(args_new);
        return fn;
    }
    return expr_copy((Expr*)e);
}

static Expr* tr_reduce_mod_ideal(const Expr* poly, const TRBindings* b) {
    /* First, ensure the polynomial is in expanded form so that
     * s_i * s_i becomes Power[s_i, 2] for our pattern to match. */
    Expr* cur = tr_call_unary_copy("Expand", poly);
    if (!cur) cur = expr_copy((Expr*)poly);
    int iters = 0;
    while (iters++ < TRIGRAT_MAX_REDUCE) {
        bool changed = false;
        Expr* nx = tr_rewrite_pass(cur, b, &changed);
        if (!changed) {
            expr_free(nx);
            break;
        }
        expr_free(cur);
        cur = tr_call_unary("Expand", nx);
        if (!cur) {
            /* Defensive: should never happen, but bail safely. */
            return expr_copy((Expr*)poly);
        }
    }
    return cur;
}

/* ----------------------------------------------------------------------- */
/* Detect whether a tree contains a particular fresh symbol name.          */
/* ----------------------------------------------------------------------- */

static bool tr_uses_sname(const Expr* e, const char* name) {
    return tr_has_symbol_name(e, name);
}

/* ----------------------------------------------------------------------- */
/* Numerator / Denominator split.                                          */
/* Returns fresh trees. Numerator/Denominator are core builtins.           */
/* ----------------------------------------------------------------------- */

static void tr_split_frac(const Expr* tg, Expr** num_out, Expr** den_out) {
    *num_out = tr_call_unary_copy("Numerator", tg);
    *den_out = tr_call_unary_copy("Denominator", tg);
    if (!*num_out) *num_out = expr_copy((Expr*)tg);
    if (!*den_out) *den_out = tr_int(1);
}

/* ----------------------------------------------------------------------- */
/* Compute Q with s -> 0 (replace a specific symbol by 0 everywhere).      */
/* ----------------------------------------------------------------------- */

static Expr* tr_subst_sym_zero(const Expr* e, const char* name) {
    if (!e) return NULL;
    if (e->type == EXPR_SYMBOL && e->data.symbol.name &&
        strcmp(e->data.symbol.name, name) == 0) {
        return tr_int(0);
    }
    if (e->type != EXPR_FUNCTION) return expr_copy((Expr*)e);
    Expr* head_new = tr_subst_sym_zero(e->data.function.head, name);
    size_t n = e->data.function.arg_count;
    Expr** args_new = (Expr**)calloc(n ? n : 1, sizeof(Expr*));
    for (size_t i = 0; i < n; i++) {
        args_new[i] = tr_subst_sym_zero(e->data.function.args[i], name);
    }
    Expr* call = expr_new_function(head_new, args_new, n);
    free(args_new);  /* expr_new_function copies the array; free our buffer. */
    /* Evaluate to fold the new constants. */
    Expr* r = evaluate(call);
    expr_free(call);
    return r;
}

/* ----------------------------------------------------------------------- */
/* Inverse odd-generator powers.                                           */
/*                                                                         */
/* A radical Sqrt[g] whose radicand g is *rational* in the generators with */
/* an odd generator s_i in its denominator -- Cot = Cos/Sin, Csc = 1/Sin   */
/* -- carries the relation l^2 = g = .../s_i.  Reducing l^2 then injects   */
/* Power[s_i, -1] into the working num/den.  The conjugate rationalisation  */
/* below clears s_i via the substitution den|_{s_i -> 0}, which would       */
/* evaluate Power[0, -1] -- a literal 1/0 -- on such a term.  We detect and */
/* clear these inverse powers first.  (For Tan = Sin/Cos the inverse        */
/* generator is the *even* c_i, which is never substitute-zeroed, so that   */
/* path never triggers this and is left exactly as before.)                */
/* ----------------------------------------------------------------------- */

/* True if `e` contains Power[s_i, k] with k a negative integer for some odd
 * trig generator s_i. */
static bool tr_has_neg_sgen_power(const Expr* e, const TRBindings* b) {
    if (!e || e->type != EXPR_FUNCTION) return false;
    const Expr* h = e->data.function.head;
    if (h && h->type == EXPR_SYMBOL && h->data.symbol.name == SYM_Power &&
        e->data.function.arg_count == 2) {
        const Expr* base = e->data.function.args[0];
        const Expr* ex   = e->data.function.args[1];
        if (base && base->type == EXPR_SYMBOL &&
            ex && ex->type == EXPR_INTEGER && ex->data.integer < 0 &&
            tr_find_by_s_name(b, base->data.symbol.name)) {
            return true;
        }
    }
    if (tr_has_neg_sgen_power(e->data.function.head, b)) return true;
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (tr_has_neg_sgen_power(e->data.function.args[i], b)) return true;
    return false;
}

/* If num/den carry inverse odd-generator powers, recombine over a common
 * denominator and re-split so both are polynomials in every s_i, then
 * re-reduce mod the ideal.  This is only reached after the radical
 * generators are cleared, so the re-reduction cannot reintroduce a fresh
 * inverse power (no l^2 remains to expand into .../s_i).  Borrows nothing;
 * rewrites *num and *den in place. */
static void tr_clear_neg_sgen(Expr** num, Expr** den, const TRBindings* b) {
    if (!tr_has_neg_sgen_power(*num, b) && !tr_has_neg_sgen_power(*den, b))
        return;
    Expr* frac = tr_times(expr_copy(*num),
                          tr_pow(expr_copy(*den), tr_int(-1)));
    Expr* tg = tr_call_unary("Together", frac);   /* consumes frac */
    if (!tg) return;
    Expr* n2; Expr* d2;
    tr_split_frac(tg, &n2, &d2);
    expr_free(tg);
    Expr* nr = tr_reduce_mod_ideal(n2, b);
    Expr* dr = tr_reduce_mod_ideal(d2, b);
    expr_free(n2); expr_free(d2);
    expr_free(*num); expr_free(*den);
    *num = nr; *den = dr;
}

/* ----------------------------------------------------------------------- */
/* Rationalise the denominator: clear every s_i / sigma_j from den.        */
/* Returns new (num, den) by output pointers, both freshly allocated.      */
/* Inputs are borrowed.                                                    */
/* ----------------------------------------------------------------------- */

static void tr_rationalise_denom(const Expr* num_in, const Expr* den_in,
                                 const TRBindings* b,
                                 Expr** num_out, Expr** den_out) {
    Expr* num = expr_copy((Expr*)num_in);
    Expr* den = expr_copy((Expr*)den_in);

    /* Clear radical generators FIRST. Clearing l reintroduces its radicand
     * (rational in s_i, c_i), so the trig generators must be cleared after.
     * The reverse order would re-dirty the denominator with the radical. */
    for (size_t i = 0; i < b->radical_count; i++) {
        const char* lname = b->radicals[i].name;
        if (!tr_uses_sname(den, lname)) continue;
        /* den = den_0 + l * den_1 (l degree <= 1 after reduction). Conjugate
         * den_0 - l * den_1 = 2*den_0 - den; product den_0^2 - l^2 den_1^2
         * is l-free after l^2 -> radicand reduction. */
        Expr* den_0 = tr_subst_sym_zero(den, lname);
        Expr* conj = tr_plus(tr_times(tr_int(2), expr_copy(den_0)),
                             tr_neg(expr_copy(den)));
        expr_free(den_0);

        Expr* new_num = tr_call_unary("Expand", tr_times(num, expr_copy(conj)));
        Expr* new_den = tr_call_unary("Expand", tr_times(den, conj));

        Expr* num_red = tr_reduce_mod_ideal(new_num, b);
        Expr* den_red = tr_reduce_mod_ideal(new_den, b);
        expr_free(new_num);
        expr_free(new_den);
        num = num_red;
        den = den_red;
    }

    /* Clearing the radicals above can leave inverse odd-generator powers
     * (s_i^(-1)) when a radicand was rational with s_i in its denominator
     * (Cot, Csc).  Clear them now so the conjugate substitution den|_{s->0}
     * below never evaluates a literal 1/0. */
    tr_clear_neg_sgen(&num, &den, b);

    for (size_t i = 0; i < b->count; i++) {
        const char* sname = b->items[i].s_name;
        if (!tr_uses_sname(den, sname)) continue;
        /* den = den_0 + s * den_1, with s having degree <= 1 (after
         * prior reduction). The conjugate is den_0 - s * den_1
         * = 2*den_0 - den. */
        Expr* den_0 = tr_subst_sym_zero(den, sname);
        Expr* two_d0 = tr_times(tr_int(2), expr_copy(den_0));
        Expr* den_neg = tr_neg(expr_copy(den));
        Expr* conj = tr_plus(two_d0, den_neg);
        expr_free(den_0);

        /* Multiply numerator and denominator by conj. */
        Expr* new_num_raw = tr_times(num, expr_copy(conj));
        Expr* new_den_raw = tr_times(den, conj);

        Expr* new_num = tr_call_unary("Expand", new_num_raw);
        Expr* new_den = tr_call_unary("Expand", new_den_raw);

        /* Reduce mod ideal to push s^2 out of den. */
        Expr* num_red = tr_reduce_mod_ideal(new_num, b);
        Expr* den_red = tr_reduce_mod_ideal(new_den, b);
        expr_free(new_num);
        expr_free(new_den);
        num = num_red;
        den = den_red;
    }

    *num_out = num;
    *den_out = den;
}

/* ----------------------------------------------------------------------- */
/* Fold-back to Tan/Sec/Csc/Cot/Tanh/Sech/Csch/Coth.                       */
/*                                                                         */
/* The evaluator already auto-canonicalises Sin[x]/Cos[x] -> Tan[x] and    */
/* 1/Cos[x] -> Sec[x] (and the hyperbolic analogues) on re-evaluation, so  */
/* a separate explicit fold-back pass is unnecessary -- the `evaluate`    */
/* call at the end of simp_trig_rational picks up these forms. Any        */
/* additional rewrites belong inside the evaluator's canonicalisation,   */
/* not as an inline post-pass here.                                       */
/* ----------------------------------------------------------------------- */

/* ----------------------------------------------------------------------- */
/* Pure (s, c, sigma, chi)-rational simplification helper.                 */
/*                                                                         */
/* Together -> ideal reduction -> denominator rationalisation -> Cancel.   */
/* Borrows the input, returns a freshly-allocated expression.              */
/*                                                                         */
/* Used by the opaque-monomial grouping path so each per-group simplify    */
/* runs on a small bivariate (or trivariate-with-hyp) expression -- the   */
/* multivariate GCD inside Cancel is then O(d^2) instead of O(d^(2+k))     */
/* with k opaque-atom generators.                                         */
/* ----------------------------------------------------------------------- */

static Expr* tr_simp_pure_sc(const Expr* e, const TRBindings* b) {
    /* Step a: Together (combine over common denominator). */
    Expr* tg = tr_call_unary_copy("Together", e);
    if (!tg) return expr_copy((Expr*)e);

    /* Step b: split. */
    Expr* num0; Expr* den0;
    tr_split_frac(tg, &num0, &den0);
    expr_free(tg);

    /* Step c: reduce both modulo the trig/hyp ideal. */
    Expr* num_r = tr_reduce_mod_ideal(num0, b);
    Expr* den_r = tr_reduce_mod_ideal(den0, b);
    expr_free(num0);
    expr_free(den0);

    /* Step d: rationalise the denominator. */
    Expr* num_rt; Expr* den_rt;
    tr_rationalise_denom(num_r, den_r, b, &num_rt, &den_rt);
    expr_free(num_r);
    expr_free(den_r);

    /* Step e: build num/den, Cancel. */
    Expr* inv_den = tr_pow(expr_copy(den_rt), tr_int(-1));
    Expr* frac = tr_times(num_rt, inv_den);
    expr_free(den_rt);

    Expr* cancelled = tr_call_unary("Cancel", frac);
    if (!cancelled) cancelled = frac;
    return cancelled;
}

/* ----------------------------------------------------------------------- */
/* Opaque-monomial grouping.                                               */
/*                                                                         */
/* After substitution and Expand, the substituted expression is a Plus    */
/* of Times terms. Each Times factor is either:                            */
/*   (a) a positive-integer power of an opaque atom L_k (or L_k itself);  */
/*   (b) something built from s_i, c_i, sigma_j, chi_j, other ground vars.*/
/* For each Plus addend, partition the factors and form                    */
/*   addend = opaque_monomial * trig_rational                              */
/* Group addends with structurally-equal opaque_monomial, sum the         */
/* trig_rational coefficients, and apply tr_simp_pure_sc to each sum.     */
/*                                                                         */
/* This sidesteps the big multivariate GCD that would otherwise run        */
/* inside Cancel on the full P/Q form. Each per-group Cancel is bivariate */
/* (or trivariate with hyperbolic), which is fast.                        */
/*                                                                         */
/* Returns NULL when any addend has an opaque atom appearing in a way     */
/* that doesn't factor cleanly: opaque atom in a denominator               */
/* (Power[L_k, negative_int]) or opaque atom inside a Plus or other       */
/* non-multiplicative structure that Expand didn't distribute.  The      */
/* caller then falls back to the existing single-Cancel pipeline.        */
/* ----------------------------------------------------------------------- */

static bool tr_is_opaque_name(const char* sym, const TRBindings* b) {
    if (!sym) return false;
    for (size_t i = 0; i < b->opaque_count; i++) {
        if (strcmp(b->opaques[i].name, sym) == 0) return true;
    }
    return false;
}

static bool tr_contains_opaque(const Expr* e, const TRBindings* b) {
    if (!e) return false;
    if (e->type == EXPR_SYMBOL) return tr_is_opaque_name(e->data.symbol.name, b);
    if (e->type != EXPR_FUNCTION) return false;
    if (tr_contains_opaque(e->data.function.head, b)) return true;
    for (size_t i = 0; i < e->data.function.arg_count; i++) {
        if (tr_contains_opaque(e->data.function.args[i], b)) return true;
    }
    return false;
}

/* Classify a single multiplicative factor (a Times child or a standalone
 * term):
 *   0  -- factor has no opaque atom
 *   1  -- factor IS a positive-integer power of an opaque atom (clean)
 *  -1  -- factor IS a negative-integer power of an opaque atom (bail)
 *   2  -- factor contains an opaque atom in a non-separable way (bail) */
static int tr_classify_factor(const Expr* f, const TRBindings* b) {
    if (!f) return 0;
    if (f->type == EXPR_SYMBOL) {
        return tr_is_opaque_name(f->data.symbol.name, b) ? 1 : 0;
    }
    if (f->type != EXPR_FUNCTION) return 0;
    const Expr* h = f->data.function.head;
    if (h && h->type == EXPR_SYMBOL && h->data.symbol.name == SYM_Power &&
        f->data.function.arg_count == 2) {
        const Expr* base = f->data.function.args[0];
        const Expr* exp_e = f->data.function.args[1];
        bool base_opq_sym = (base && base->type == EXPR_SYMBOL &&
                              tr_is_opaque_name(base->data.symbol.name, b));
        if (base_opq_sym && exp_e && exp_e->type == EXPR_INTEGER) {
            return (exp_e->data.integer >= 1) ? 1 : -1;
        }
        /* Power but base contains an opaque (e.g. Power[Plus[L,1], 2])
         * -- non-separable. */
        if (tr_contains_opaque(f, b)) return 2;
        return 0;
    }
    /* Any other function call: if it has an opaque inside, it's a
     * non-separable mix (e.g. Plus[L, 1] as a multiplicative factor). */
    return tr_contains_opaque(f, b) ? 2 : 0;
}

/* Recursively flatten a term into its multiplicative factors. Times is
 * supposed to be Flat (so any Times-of-Times should already be one big
 * Times by evaluator canonicalisation), but Expand can produce
 * unflattened intermediate results, so we defensively walk through
 * nested Times here. */
static void tr_collect_factors(const Expr* e, const Expr*** out, size_t* n, size_t* cap) {
    if (e->type == EXPR_FUNCTION && e->data.function.head &&
        e->data.function.head->type == EXPR_SYMBOL &&
        e->data.function.head->data.symbol.name == SYM_Times) {
        for (size_t i = 0; i < e->data.function.arg_count; i++) {
            tr_collect_factors(e->data.function.args[i], out, n, cap);
        }
        return;
    }
    if (*n == *cap) {
        size_t nc = *cap == 0 ? 8 : *cap * 2;
        *out = (const Expr**)realloc((void*)*out, nc * sizeof(const Expr*));
        *cap = nc;
    }
    (*out)[(*n)++] = e;
}

/* Take a Plus addend; partition its multiplicative factors into an
 * opaque-monomial and a trig-rational rest. Returns true on success. */
static bool tr_split_addend(const Expr* term, const TRBindings* b,
                            Expr** opq_out, Expr** rest_out) {
    const Expr** factors = NULL;
    size_t n_factors = 0;
    size_t cap = 0;
    tr_collect_factors(term, &factors, &n_factors, &cap);
    /* tr_collect_factors guarantees n_factors >= 1 for non-NULL term. */

    Expr** opq_factors  = (Expr**)calloc(n_factors ? n_factors : 1, sizeof(Expr*));
    Expr** rest_factors = (Expr**)calloc(n_factors ? n_factors : 1, sizeof(Expr*));
    size_t opq_n = 0, rest_n = 0;
    bool ok = true;

    for (size_t i = 0; i < n_factors; i++) {
        int k = tr_classify_factor(factors[i], b);
        if (k == -1 || k == 2) {
            ok = false; break;
        }
        if (k == 1) {
            opq_factors[opq_n++] = expr_copy((Expr*)factors[i]);
        } else {
            rest_factors[rest_n++] = expr_copy((Expr*)factors[i]);
        }
    }
    free((void*)factors);

    if (!ok) {
        for (size_t j = 0; j < opq_n;  j++) expr_free(opq_factors[j]);
        for (size_t j = 0; j < rest_n; j++) expr_free(rest_factors[j]);
        free(opq_factors);
        free(rest_factors);
        return false;
    }

    Expr* opq;
    if (opq_n == 0) {
        opq = tr_int(1);
        free(opq_factors);
    } else if (opq_n == 1) {
        opq = opq_factors[0];
        free(opq_factors);
    } else {
        Expr* call = expr_new_function(expr_new_symbol(SYM_Times),
                                       opq_factors, opq_n);
        free(opq_factors);
        opq = evaluate(call);
        expr_free(call);
    }

    Expr* rest;
    if (rest_n == 0) {
        rest = tr_int(1);
        free(rest_factors);
    } else if (rest_n == 1) {
        rest = rest_factors[0];
        free(rest_factors);
    } else {
        Expr* call = expr_new_function(expr_new_symbol(SYM_Times),
                                       rest_factors, rest_n);
        free(rest_factors);
        rest = evaluate(call);
        expr_free(call);
    }

    *opq_out = opq;
    *rest_out = rest;
    return true;
}

typedef struct {
    Expr*  key;        /* opaque-monomial (owned) */
    Expr** terms;      /* trig-rational pieces to be summed (owned) */
    size_t n;
    size_t cap;
} TROpqGroup;

typedef struct {
    TROpqGroup* items;
    size_t      count;
    size_t      cap;
} TROpqGroups;

static void tr_groups_init(TROpqGroups* g) {
    g->items = NULL; g->count = 0; g->cap = 0;
}

static void tr_groups_free(TROpqGroups* g) {
    for (size_t i = 0; i < g->count; i++) {
        if (g->items[i].key) expr_free(g->items[i].key);
        for (size_t j = 0; j < g->items[i].n; j++) expr_free(g->items[i].terms[j]);
        free(g->items[i].terms);
    }
    free(g->items);
    g->items = NULL; g->count = 0; g->cap = 0;
}

static TROpqGroup* tr_groups_find(TROpqGroups* g, const Expr* key) {
    for (size_t i = 0; i < g->count; i++) {
        if (expr_eq(g->items[i].key, key)) return &g->items[i];
    }
    return NULL;
}

static TROpqGroup* tr_groups_add(TROpqGroups* g, Expr* key /* takes ownership */) {
    if (g->count == g->cap) {
        size_t nc = g->cap == 0 ? 4 : g->cap * 2;
        g->items = (TROpqGroup*)realloc(g->items, nc * sizeof(TROpqGroup));
        g->cap = nc;
    }
    TROpqGroup* it = &g->items[g->count++];
    it->key = key;
    it->terms = NULL; it->n = 0; it->cap = 0;
    return it;
}

static void tr_group_add_term(TROpqGroup* gr, Expr* term /* takes ownership */) {
    if (gr->n == gr->cap) {
        size_t nc = gr->cap == 0 ? 4 : gr->cap * 2;
        gr->terms = (Expr**)realloc(gr->terms, nc * sizeof(Expr*));
        gr->cap = nc;
    }
    gr->terms[gr->n++] = term;
}

static Expr* tr_simp_l_grouped(const Expr* substed, const TRBindings* b) {
    /* Grouping is only useful when there are at least two distinct
     * opaque atoms; with zero opaques every addend lives in the same
     * "key = 1" bucket so grouping is a no-op (and the Expand cost is
     * a strict regression vs. the single-fraction fallback). With one
     * opaque, grouping is still profitable iff the opaque appears
     * multiple times -- but the caller has no cheap way to count
     * occurrences, so bail unless we have multiple distinct opaques. */
    if (b->opaque_count == 0) return NULL;

    /* Step a: try without Expand first. If the substituted form is
     * already a flat Plus of multiplicative addends (the common case
     * after evaluation), no distribution is needed and we save the
     * Expand cost entirely. If any addend contains a non-separable
     * opaque-bearing factor (FACTOR_OPQ_MIXED), retry with Expand. */
    const Expr* source = substed;
    Expr* expanded = NULL;

    bool top_is_plus = (source->type == EXPR_FUNCTION &&
                        source->data.function.head &&
                        source->data.function.head->type == EXPR_SYMBOL &&
                        source->data.function.head->data.symbol.name == SYM_Plus);

    size_t n_addends = top_is_plus ? source->data.function.arg_count : 1;

    /* Step b: split each addend into (opaque_monomial, trig_rest), group.
     * If we hit a non-separable factor, retry with Expand once; if it
     * still fails after expansion, bail to the caller. */
    TROpqGroups groups;
    tr_groups_init(&groups);
    bool ok = true;
    bool retried = false;
retry:
    for (size_t i = 0; i < n_addends; i++) {
        const Expr* term = top_is_plus ? source->data.function.args[i]
                                       : source;
        Expr* opq; Expr* rest;
        if (!tr_split_addend(term, b, &opq, &rest)) {
            ok = false;
            break;
        }
        TROpqGroup* gr = tr_groups_find(&groups, opq);
        if (gr) {
            expr_free(opq);
        } else {
            gr = tr_groups_add(&groups, opq);
        }
        tr_group_add_term(gr, rest);
    }

    if (!ok && !retried) {
        /* Distribute products-of-sums and retry once. */
        tr_groups_free(&groups);
        tr_groups_init(&groups);
        expanded = tr_call_unary_copy("Expand", substed);
        if (!expanded) return NULL;
        source = expanded;
        top_is_plus = (source->type == EXPR_FUNCTION &&
                       source->data.function.head &&
                       source->data.function.head->type == EXPR_SYMBOL &&
                       source->data.function.head->data.symbol.name == SYM_Plus);
        n_addends = top_is_plus ? source->data.function.arg_count : 1;
        ok = true;
        retried = true;
        goto retry;
    }

    if (!ok) {
        tr_groups_free(&groups);
        if (expanded) expr_free(expanded);
        return NULL;
    }

    /* Step c: for each group, sum the trig-rest pieces and apply
     * tr_simp_pure_sc on the sum (small, fast). */
    Expr** result_terms = (Expr**)calloc(groups.count ? groups.count : 1,
                                          sizeof(Expr*));
    size_t n_result = 0;
    for (size_t i = 0; i < groups.count; i++) {
        TROpqGroup* gr = &groups.items[i];
        Expr* sum;
        if (gr->n == 1) {
            sum = expr_copy(gr->terms[0]);
        } else {
            Expr** args = (Expr**)calloc(gr->n, sizeof(Expr*));
            for (size_t j = 0; j < gr->n; j++) args[j] = expr_copy(gr->terms[j]);
            Expr* call = expr_new_function(expr_new_symbol(SYM_Plus), args, gr->n);
            free(args);
            sum = evaluate(call);
            expr_free(call);
        }
        Expr* simp = tr_simp_pure_sc(sum, b);
        expr_free(sum);

        /* Multiply by the opaque-monomial key (which may be 1 for the
         * "constant" bucket). evaluate folds Times[1, x] -> x and
         * Times[k, 0] -> 0 so trivial groups vanish naturally. */
        Expr* args2[2] = { expr_copy(gr->key), simp };
        Expr* call2 = expr_new_function(expr_new_symbol(SYM_Times), args2, 2);
        Expr* term = evaluate(call2);
        expr_free(call2);
        result_terms[n_result++] = term;
    }
    tr_groups_free(&groups);
    if (expanded) expr_free(expanded);

    /* Step d: sum all per-group results. */
    Expr* result;
    if (n_result == 0) {
        result = tr_int(0);
        free(result_terms);
    } else if (n_result == 1) {
        result = result_terms[0];
        free(result_terms);
    } else {
        Expr* call = expr_new_function(expr_new_symbol(SYM_Plus),
                                       result_terms, n_result);
        free(result_terms);
        result = evaluate(call);
        expr_free(call);
    }
    return result;
}

/* ----------------------------------------------------------------------- */
/* The entry point.                                                        */
/* ----------------------------------------------------------------------- */

Expr* simp_trig_rational(const Expr* input,
                         const AssumeCtx* ctx,
                         const Expr* complexity_func) {
    (void)ctx; /* Unused for now; the algorithm is correctness-preserving
                  without assumptions. */
    if (!input) return NULL;

    size_t in_score = complexity_func
        ? simp_default_complexity(input)
        : simp_default_complexity(input);
    /* complexity_func support deferred: this fast path uses the default
     * leaf count for its strict-win gate. If the caller provided a
     * custom complexity, we'd still compute the result here using the
     * default and let simp_search's loop score with the custom one. */
    (void)complexity_func;

    if (in_score < TRIGRAT_LEAF_FLOOR) return NULL;

    /* Steps 1-3 combined: top-down walk that rewrites Tan/Cot/etc. to
     * Sin/Cos form, collects trig/hyp pairs and opaque atoms, and
     * substitutes everything into fresh ground-field symbols. */
    TRBindings b;
    tr_bindings_init(&b);
    int counter = 0;
    bool over_budget = false;
    Expr* substed = tr_walk_subst(input, &b, &counter, input, &over_budget);
    if (over_budget || (b.count == 0 && b.radical_count == 0) || !substed) {
        tr_bindings_free(&b);
        if (substed) expr_free(substed);
        return NULL;
    }

    /* Step 3.5: If any trig argument is non-atomic (e.g. Sin[x + Pi/6]
     * or Sin[2 x]), retry after TrigExpand'ing the input. TrigExpand
     * rewrites Sin[x + Pi/6] -> (Sqrt[3]/2) Sin[x] + (1/2) Cos[x] and
     * Sin[2 x] -> 2 Sin[x] Cos[x], canonicalising every trig argument
     * into a sum of atomic-argument trig calls so the algebraic
     * normal-form algorithm can recognise the angle-addition /
     * multiple-angle relationships. TrigExpand is cheap when arguments
     * are already atomic (a no-op tree walk), so checking first avoids
     * paying its cost on inputs like the user's derivative-of-Risch-
     * Norman cases where every Sin/Cos has argument `x`. */
    bool any_compound = false;
    for (size_t i = 0; i < b.count; i++) {
        if (b.items[i].arg && b.items[i].arg->type != EXPR_SYMBOL) {
            any_compound = true; break;
        }
    }
    if (any_compound) {
        Expr* trigexp = tr_call_unary_copy("TrigExpand", input);
        if (!trigexp || expr_eq(trigexp, input)) {
            /* TrigExpand didn't change the form -- our compound args
             * weren't of a kind TrigExpand can unfold (e.g. Sin[f[x]]).
             * Defer to simp_search rather than risk producing a partial
             * normal form that locks out the canonical answer. */
            if (trigexp) expr_free(trigexp);
            tr_bindings_free(&b);
            expr_free(substed);
            return NULL;
        }
        /* Redo the walk on the expanded input. Reset bindings. */
        expr_free(substed);
        tr_bindings_free(&b);
        tr_bindings_init(&b);
        counter = 0;
        over_budget = false;
        substed = tr_walk_subst(trigexp, &b, &counter, trigexp, &over_budget);
        expr_free(trigexp);
        if (over_budget || b.count == 0 || !substed) {
            tr_bindings_free(&b);
            if (substed) expr_free(substed);
            return NULL;
        }
        /* If TrigExpand still left compound args, give up cleanly. */
        for (size_t i = 0; i < b.count; i++) {
            if (b.items[i].arg && b.items[i].arg->type != EXPR_SYMBOL) {
                tr_bindings_free(&b);
                expr_free(substed);
                return NULL;
            }
        }
    }

    /* Steps 4-7: simplify the substituted form to a fraction in
     * (Sin/Cos) generators. Two paths:
     *
     *   - Opaque-monomial grouping (fast): when each addend separates
     *     cleanly into (opaque_monomial * trig_rational_rest), we
     *     simplify each group's trig_rational sum in (s, c) only --
     *     Cancel runs on a bivariate polynomial, not on the full
     *     multivariate (s, c, L_1, ..., L_k) ring. This is the main
     *     reason Simplify can return in ~0.1 s on derivative-of-Risch-
     *     Norman inputs whose direct Together blows up because of the
     *     Log/Exp atoms.
     *
     *   - Fallback single-fraction pipeline: when the grouping cannot
     *     factor cleanly (opaque atom in a denominator, opaque atom
     *     fused with trig inside a Power[Plus[...], _], etc.), fall
     *     back to one big Together / ideal-reduce / rationalise /
     *     Cancel. Slower in the worst case but always correct. */
    Expr* consolidated = tr_simp_l_grouped(substed, &b);
    if (!consolidated) {
        /* Fallback path. */
        Expr* tg = tr_call_unary_copy("Together", substed);
        if (!tg) tg = expr_copy((Expr*)input);

        Expr* num0; Expr* den0;
        tr_split_frac(tg, &num0, &den0);
        expr_free(tg);

        Expr* num_red = tr_reduce_mod_ideal(num0, &b);
        Expr* den_red = tr_reduce_mod_ideal(den0, &b);
        expr_free(num0);
        expr_free(den0);

        Expr* num_r; Expr* den_r;
        tr_rationalise_denom(num_red, den_red, &b, &num_r, &den_r);
        expr_free(num_red);
        expr_free(den_red);

        Expr* inv_den = tr_pow(expr_copy(den_r), tr_int(-1));
        Expr* frac = tr_times(num_r, inv_den);
        expr_free(den_r);

        Expr* cancelled = tr_call_unary("Cancel", frac);
        if (!cancelled) cancelled = frac;

        consolidated = tr_call_unary_copy("Together", cancelled);
        expr_free(cancelled);
        if (!consolidated) consolidated = expr_copy((Expr*)input);
    }
    expr_free(substed);

    /* Step 8: substitute back to Sin/Cos/Sinh/Cosh. */
    Expr* back = tr_subst_back(consolidated, &b);
    expr_free(consolidated);
    tr_bindings_free(&b);
    if (!back) return NULL;

    /* Re-evaluate to fold any new constants (Sin[0], etc.). */
    Expr* evald = evaluate(back);
    expr_free(back);
    if (!evald) return NULL;

    /* Step 9 (fold-back) is handled implicitly by `evaluate` above:
     * Sin[x]/Cos[x] -> Tan[x] and 1/Cos[x] -> Sec[x] (plus hyperbolic
     * analogues) are auto-canonicalisations in the evaluator. */

    /* Step 10: strict leaf-count gate. */
    size_t out_score = simp_default_complexity(evald);
    if (out_score < in_score) {
        return evald;
    }
    expr_free(evald);
    return NULL;
}
