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
/* Assumption-driven seed rewriters                                        */
/* ----------------------------------------------------------------------- */

/* For each direct EXPR_SYMBOL fact-target, generate context-specific
 * rewrite rules and apply them via ReplaceRepeated. The rules are
 * unconditional in pattern form: their conditional nature is captured by
 * the choice of the rule's free symbol -- e.g., we only emit
 *   Power[Power[<x>, 2], Rational[1, 2]] :> <x>
 * when <x> is the literal symbol that the assumption set says is
 * positive. So the rules are valid by construction whenever applied.
 *
 * The generated rule list is built as a string and parsed; this is
 * cheaper to maintain than constructing the AST by hand and matches the
 * style used in trigsimp.c.
 */

static bool sym_already_listed(char** list, size_t n, const char* s) {
    for (size_t i = 0; i < n; i++) if (strcmp(list[i], s) == 0) return true;
    return false;
}

/* Walk the assumption fact list and collect every EXPR_SYMBOL that the
 * context proves positive, real, integer, or even. The caller passes
 * pre-sized arrays plus the maximum count. */
static void collect_known_symbols(const AssumeCtx* ctx,
                                  char** positives, size_t* npos,
                                  char** reals,     size_t* nreal,
                                  char** integers,  size_t* nint,
                                  char** negatives, size_t* nneg,
                                  char** evens,     size_t* neven,
                                  size_t cap) {
    *npos = *nreal = *nint = *nneg = *neven = 0;
    if (!ctx) return;
    /* Mod[m, 2] == 0 hides `m` inside a function argument, so a top-level
     * "scan operands of facts" alone misses it. We additionally walk
     * Mod[s, _] for any symbol s that appears under an even-type fact. */
    for (size_t i = 0; i < ctx->count; i++) {
        const Expr* f = ctx->facts[i];
        if (f->type != EXPR_FUNCTION) continue;
        for (size_t j = 0; j < f->data.function.arg_count; j++) {
            Expr* a = f->data.function.args[j];
            if (a->type == EXPR_SYMBOL) {
                const char* nm = a->data.symbol.name;
                if (assume_known_positive(ctx, a) && *npos < cap && !sym_already_listed(positives, *npos, nm)) {
                    positives[(*npos)++] = (char*)nm;
                }
                if (assume_known_negative(ctx, a) && *nneg < cap && !sym_already_listed(negatives, *nneg, nm)) {
                    negatives[(*nneg)++] = (char*)nm;
                }
                if (assume_known_real(ctx, a) && *nreal < cap && !sym_already_listed(reals, *nreal, nm)) {
                    reals[(*nreal)++] = (char*)nm;
                }
                if (assume_known_integer(ctx, a) && *nint < cap && !sym_already_listed(integers, *nint, nm)) {
                    integers[(*nint)++] = (char*)nm;
                }
                if (assume_known_even(ctx, a) && *neven < cap && !sym_already_listed(evens, *neven, nm)) {
                    evens[(*neven)++] = (char*)nm;
                }
            } else if (a->type == EXPR_FUNCTION &&
                       a->data.function.head &&
                       a->data.function.head->type == EXPR_SYMBOL &&
                       a->data.function.head->data.symbol.name == SYM_Mod &&
                       a->data.function.arg_count == 2 &&
                       a->data.function.args[0]->type == EXPR_SYMBOL) {
                Expr* sym = a->data.function.args[0];
                const char* nm = sym->data.symbol.name;
                if (assume_known_even(ctx, sym) && *neven < cap && !sym_already_listed(evens, *neven, nm)) {
                    evens[(*neven)++] = (char*)nm;
                }
            }
        }
    }
}

/* ----------------------------------------------------------------------- */
/* Structural assumption rewrites                                          */
/*                                                                         */
/* Some assumption-driven simplifications cannot be expressed as a rule    */
/* keyed on a pre-enumerated symbol: Floor[e]->e needs to know e is an     */
/* integer for an ARBITRARY subexpression e; Ceiling[x]->3 needs interval  */
/* reasoning; Mod[a,4]->1 needs modular reasoning from an Element fact;    */
/* Re[a+b I]->a needs every free symbol to be real. These are handled by a */
/* bottom-up structural walk that queries the assume_known_* API at each   */
/* node. The walk runs inside apply_assumption_rules, so Simplify (which   */
/* calls apply_assumption_rules) benefits from it too.                     */
/* ----------------------------------------------------------------------- */

static bool sr_head(const Expr* e, const char* name) {
    return e && e->type == EXPR_FUNCTION && e->data.function.head &&
           e->data.function.head->type == EXPR_SYMBOL &&
           strcmp(e->data.function.head->data.symbol.name, name) == 0;
}

/* Numeric value of an Integer/Real/Bigint/Rational literal. */
static bool sr_num(const Expr* e, double* v) {
    if (!e) return false;
    if (e->type == EXPR_INTEGER) { *v = (double)e->data.integer; return true; }
    if (e->type == EXPR_REAL)    { *v = e->data.real; return true; }
    if (e->type == EXPR_BIGINT)  { *v = mpz_get_d(e->data.bigint); return true; }
    if (sr_head(e, "Rational") && e->data.function.arg_count == 2) {
        const Expr* p = e->data.function.args[0];
        const Expr* q = e->data.function.args[1];
        if (p->type == EXPR_INTEGER && q->type == EXPR_INTEGER && q->data.integer != 0) {
            *v = (double)p->data.integer / (double)q->data.integer; return true;
        }
    }
    return false;
}

/* Does ctx prove the strict order A < B by a direct inequality fact? */
static bool sr_proves_slt(const AssumeCtx* ctx, const Expr* A, const Expr* B) {
    if (!ctx) return false;
    for (size_t i = 0; i < ctx->count; i++) {
        const Expr* f = ctx->facts[i];
        if (f->type != EXPR_FUNCTION || !f->data.function.head ||
            f->data.function.head->type != EXPR_SYMBOL ||
            f->data.function.arg_count != 2) continue;
        const char* h = f->data.function.head->data.symbol.name;
        const Expr* x = f->data.function.args[0];
        const Expr* y = f->data.function.args[1];
        if (h == SYM_Less && expr_eq((Expr*)x, (Expr*)A) && expr_eq((Expr*)y, (Expr*)B))
            return true;
        if (h == SYM_Greater && expr_eq((Expr*)x, (Expr*)B) && expr_eq((Expr*)y, (Expr*)A))
            return true;
    }
    return false;
}

/* Tightest numeric bounds on the bare symbol `sym` from direct inequality
 * facts. Returns whether a lower / upper bound was found; *lo_strict indicates
 * a strict `>`. */
static void sr_bounds(const AssumeCtx* ctx, const Expr* sym,
                      bool* has_lo, double* lo, bool* lo_strict,
                      bool* has_hi, double* hi, bool* hi_strict) {
    *has_lo = *has_hi = false;
    if (!ctx) return;
    for (size_t i = 0; i < ctx->count; i++) {
        const Expr* f = ctx->facts[i];
        if (f->type != EXPR_FUNCTION || !f->data.function.head ||
            f->data.function.head->type != EXPR_SYMBOL ||
            f->data.function.arg_count != 2) continue;
        const char* h = f->data.function.head->data.symbol.name;
        const Expr* A = f->data.function.args[0];
        const Expr* B = f->data.function.args[1];
        double k;
        bool sym_left, strict; int kind; /* kind: -1 lower, +1 upper */
        if (expr_eq((Expr*)A, (Expr*)sym) && sr_num(B, &k)) { sym_left = true; }
        else if (expr_eq((Expr*)B, (Expr*)sym) && sr_num(A, &k)) { sym_left = false; }
        else continue;
        /* Normalise to "sym <op> k". */
        if (h == SYM_Less)        { kind = sym_left ? +1 : -1; strict = true;  }
        else if (h == SYM_LessEqual)   { kind = sym_left ? +1 : -1; strict = false; }
        else if (h == SYM_Greater)     { kind = sym_left ? -1 : +1; strict = true;  }
        else if (h == SYM_GreaterEqual){ kind = sym_left ? -1 : +1; strict = false; }
        else continue;
        if (kind < 0) { /* lower bound sym > k or sym >= k */
            if (!*has_lo || k > *lo || (k == *lo && strict && !*lo_strict)) {
                *has_lo = true; *lo = k; *lo_strict = strict;
            }
        } else {        /* upper bound sym < k or sym <= k */
            if (!*has_hi || k < *hi || (k == *hi && strict && !*hi_strict)) {
                *has_hi = true; *hi = k; *hi_strict = strict;
            }
        }
    }
}

/* Collect distinct bare non-constant symbols appearing in e. */
static void sr_collect_syms(const Expr* e, const Expr*** out, size_t* n, size_t* cap) {
    if (!e) return;
    if (e->type == EXPR_SYMBOL) {
        if (is_real_constant_symbol(e->data.symbol.name)) return;
        for (size_t i = 0; i < *n; i++)
            if ((*out)[i]->data.symbol.name == e->data.symbol.name) return;
        if (*n == *cap) { *cap = *cap ? *cap * 2 : 8; *out = realloc(*out, *cap * sizeof(Expr*)); }
        (*out)[(*n)++] = e;
        return;
    }
    if (e->type == EXPR_FUNCTION) {
        /* Descend into arguments only: a function's head (Plus, Times, Re, ...)
         * is a structural operator, not a variable whose reality we require. */
        for (size_t i = 0; i < e->data.function.arg_count; i++)
            sr_collect_syms(e->data.function.args[i], out, n, cap);
    }
}

/* Bottom-up walk applying the structural assumption rewrites. Always returns a
 * freshly owned tree; sets *changed when any node was rewritten. */
static Expr* assume_structural_rewrite(const Expr* e, const AssumeCtx* ctx, int* changed) {
    if (!e) return NULL;
    if (e->type != EXPR_FUNCTION) return expr_copy((Expr*)e);

    size_t n = e->data.function.arg_count;
    Expr** args = (Expr**)calloc(n ? n : 1, sizeof(Expr*));
    for (size_t i = 0; i < n; i++)
        args[i] = assume_structural_rewrite(e->data.function.args[i], ctx, changed);
    Expr* node = expr_new_function(expr_copy(e->data.function.head), args, n);
    free(args);

    if (node->data.function.head->type != EXPR_SYMBOL) return node;
    const char* h = node->data.function.head->data.symbol.name;
    size_t nn = node->data.function.arg_count;

    if (nn == 1) {
        Expr* a0 = node->data.function.args[0];

        /* Integer-valued roundings of a provably-integer argument. */
        if ((strcmp(h, "Floor") == 0 || strcmp(h, "Ceiling") == 0 ||
             strcmp(h, "Round") == 0 || strcmp(h, "IntegerPart") == 0) &&
            assume_known_integer(ctx, a0)) {
            Expr* out = expr_copy(a0); expr_free(node); *changed = 1; return out;
        }
        if (strcmp(h, "FractionalPart") == 0 && assume_known_integer(ctx, a0)) {
            expr_free(node); *changed = 1; return expr_new_integer(0);
        }

        /* Ceiling[x] / Floor[x] pinned to a single integer by numeric bounds. */
        if ((strcmp(h, "Ceiling") == 0 || strcmp(h, "Floor") == 0) &&
            a0->type == EXPR_SYMBOL) {
            bool hl, hu, ls = false, us = false; double lo = 0, hi = 0;
            sr_bounds(ctx, a0, &hl, &lo, &ls, &hu, &hi, &us);
            if (strcmp(h, "Ceiling") == 0 && hu && hl) {
                long c = (long)ceil(hi);
                if ((double)c >= hi &&                       /* s <= c */
                    (lo > (double)(c - 1) || (lo == (double)(c - 1) && ls))) { /* s > c-1 */
                    expr_free(node); *changed = 1; return expr_new_integer(c);
                }
            }
            if (strcmp(h, "Floor") == 0 && hl && hu) {
                long fl = (long)floor(lo);
                if ((double)fl <= lo &&                      /* s >= f */
                    (hi < (double)(fl + 1) || (hi == (double)(fl + 1) && us))) { /* s < f+1 */
                    expr_free(node); *changed = 1; return expr_new_integer(fl);
                }
            }
        }

        /* FractionalPart[a] from a Mod[a,1] fact plus the sign of a. */
        if (strcmp(h, "FractionalPart") == 0 && a0->type == EXPR_SYMBOL && ctx) {
            for (size_t i = 0; i < ctx->count; i++) {
                const Expr* f = ctx->facts[i];
                if (!fact_is_function(f, "Equal", 2)) continue;
                const Expr* L = f->data.function.args[0];
                const Expr* R = f->data.function.args[1];
                const Expr* modv = NULL; const Expr* rv = NULL;
                if (sr_head(L, "Mod") && L->data.function.arg_count == 2) { modv = L; rv = R; }
                else if (sr_head(R, "Mod") && R->data.function.arg_count == 2) { modv = R; rv = L; }
                if (!modv) continue;
                if (!expr_eq((Expr*)modv->data.function.args[0], a0)) continue;
                double one; const Expr* m1 = modv->data.function.args[1];
                if (!(sr_num(m1, &one) && one == 1.0)) continue;
                double r;
                if (!sr_num(rv, &r) || !(r > 0.0 && r < 1.0)) continue;
                Expr* out = NULL;
                if (assume_known_negative(ctx, a0))
                    out = eval_and_free(expr_new_function(expr_new_symbol(SYM_Plus),
                              (Expr*[]){ expr_copy((Expr*)rv), expr_new_integer(-1) }, 2));
                else if (assume_known_nonneg(ctx, a0))
                    out = expr_copy((Expr*)rv);
                if (out) { expr_free(node); *changed = 1; return out; }
            }
        }

        /* Re/Im/Conjugate/Arg/Abs of an expression whose every free symbol is
         * real: ComplexExpand assumes reality, so it yields the refined form
         * (e.g. Abs[a + b I] -> Sqrt[a^2 + b^2]). When ComplexExpand cannot
         * reduce (unknown sign, as in Abs[x]) it returns the input unchanged,
         * which the expr_eq guard below discards. */
        if ((strcmp(h, "Re") == 0 || strcmp(h, "Im") == 0 ||
             strcmp(h, "Conjugate") == 0 || strcmp(h, "Arg") == 0 ||
             strcmp(h, "Abs") == 0)) {
            const Expr** syms = NULL; size_t ns = 0, cap = 0;
            sr_collect_syms(a0, &syms, &ns, &cap);
            bool all_real = (ns > 0);
            for (size_t i = 0; i < ns; i++)
                if (!prov_re(ctx, syms[i])) { all_real = false; break; }
            free(syms);
            if (all_real) {
                Expr* ce = eval_and_free(expr_new_function(expr_new_symbol("ComplexExpand"),
                               (Expr*[]){ expr_copy(node) }, 1));
                if (ce && !expr_eq(ce, node)) { expr_free(node); *changed = 1; return ce; }
                if (ce) expr_free(ce);
            }
        }

        /* ArcTan[Tan[e]] -> e when Re[e] (or e itself) lies in (-Pi/2, Pi/2). */
        if (strcmp(h, "ArcTan") == 0 && sr_head(a0, "Tan") &&
            a0->data.function.arg_count == 1) {
            const Expr* inner = a0->data.function.args[0];
            Expr* pihalf = eval_and_free(parse_expression("Pi/2"));
            Expr* neghalf = eval_and_free(parse_expression("-Pi/2"));
            Expr* u_re = eval_and_free(expr_new_function(expr_new_symbol(SYM_Re),
                             (Expr*[]){ expr_copy((Expr*)inner) }, 1));
            bool ok = (sr_proves_slt(ctx, neghalf, u_re) && sr_proves_slt(ctx, u_re, pihalf)) ||
                      (sr_proves_slt(ctx, neghalf, inner) && sr_proves_slt(ctx, inner, pihalf));
            expr_free(pihalf); expr_free(neghalf); expr_free(u_re);
            if (ok) { Expr* out = expr_copy((Expr*)inner); expr_free(node); *changed = 1; return out; }
        }
    }

    /* Mod[a, m] -> r from an Element[(a + c)/m, Integers] fact. */
    if (strcmp(h, "Mod") == 0 && nn == 2 && ctx) {
        Expr* a = node->data.function.args[0];
        Expr* m = node->data.function.args[1];
        if (m->type == EXPR_INTEGER && m->data.integer > 0) {
            long mm = m->data.integer;
            for (size_t i = 0; i < ctx->count; i++) {
                const Expr* f = ctx->facts[i];
                if (!fact_is_function(f, "Element", 2)) continue;
                const Expr* E = f->data.function.args[0];
                const Expr* dom = f->data.function.args[1];
                if (dom->type != EXPR_SYMBOL || dom->data.symbol.name != SYM_Integers) continue;
                /* d = m*E - a; valid iff d is a constant integer (a cancels). */
                Expr* d = eval_and_free(expr_new_function(expr_new_symbol(SYM_Plus),
                    (Expr*[]){ expr_new_function(expr_new_symbol(SYM_Times),
                                   (Expr*[]){ expr_new_integer(mm), expr_copy((Expr*)E) }, 2),
                               expr_new_function(expr_new_symbol(SYM_Times),
                                   (Expr*[]){ expr_new_integer(-1), expr_copy(a) }, 2) }, 2));
                if (d->type == EXPR_INTEGER) {
                    long dv = d->data.integer;
                    long r = ((-dv) % mm + mm) % mm;
                    expr_free(d); expr_free(node); *changed = 1; return expr_new_integer(r);
                }
                expr_free(d);
            }
        }
    }

    return node;
}

/* Produce a rewritten expression by applying assumption-derived rules via
 * ReplaceRepeated. Returns a newly owned expression, or NULL if no rules
 * were generated. The input is not consumed. */
/* True if `name` (an interned symbol pointer) occurs anywhere in `e`. */
static bool ar_expr_has_symbol(const Expr* e, const char* name) {
    if (!e) return false;
    if (e->type == EXPR_SYMBOL) return e->data.symbol.name == name;
    if (e->type == EXPR_FUNCTION) {
        if (ar_expr_has_symbol(e->data.function.head, name)) return true;
        for (size_t i = 0; i < e->data.function.arg_count; i++)
            if (ar_expr_has_symbol(e->data.function.args[i], name)) return true;
    }
    return false;
}

/* Compact a symbol bucket down to those symbols that actually appear in the
 * target expression. A per-symbol rewrite rule for a symbol absent from `input`
 * can never fire, so dropping it loses nothing -- and it keeps the O(n^2)
 * pairwise positive rules and the 8192-byte rule buffer from overflowing when
 * many symbols are assumed but few are used (e.g. Abs[a1] under 20 positivity
 * facts, which otherwise dropped ALL rules on overflow). */
static void ar_filter_to_input(char** arr, size_t* n, const Expr* input) {
    size_t w = 0;
    for (size_t i = 0; i < *n; i++)
        if (ar_expr_has_symbol(input, arr[i])) arr[w++] = arr[i];
    *n = w;
}

Expr* apply_assumption_rules(const Expr* input, const AssumeCtx* ctx) {
    if (!ctx) return NULL;

    /* Conservative caps for the per-symbol rule synthesis. */
    enum { MAX_SYM = 16 };
    char* positives[MAX_SYM]; size_t npos;
    char* reals    [MAX_SYM]; size_t nreal;
    char* integers [MAX_SYM]; size_t nint;
    char* negatives[MAX_SYM]; size_t nneg;
    char* evens    [MAX_SYM]; size_t neven;
    collect_known_symbols(ctx, positives, &npos, reals, &nreal,
                          integers, &nint, negatives, &nneg,
                          evens, &neven, MAX_SYM);

    /* Three further buckets not gathered by collect_known_symbols:
     *   nonnegs -- x >= 0  (for (x^m)^r -> x^(m r) and Sqrt[x^2] -> x at x = 0)
     *   nonpos  -- x <= 0  (for Sqrt[x^2] -> -x and Abs[x] -> -x at x = 0)
     *   abslt1  -- -1 < b < 1  (for (a^b)^c -> a^(b c)). */
    char* nonnegs[MAX_SYM]; size_t nnn = 0;
    char* nonpos [MAX_SYM]; size_t nnp = 0;
    char* abslt1 [MAX_SYM]; size_t nab = 0;
    if (ctx) {
        Expr* neg1 = expr_new_integer(-1);
        Expr* pos1 = expr_new_integer(1);
        for (size_t i = 0; i < ctx->count; i++) {
            const Expr* f = ctx->facts[i];
            if (f->type != EXPR_FUNCTION) continue;
            for (size_t j = 0; j < f->data.function.arg_count; j++) {
                Expr* a = f->data.function.args[j];
                if (a->type != EXPR_SYMBOL) continue;
                const char* nm = a->data.symbol.name;
                if (assume_known_nonneg(ctx, a) && nnn < MAX_SYM &&
                    !sym_already_listed(nonnegs, nnn, nm)) nonnegs[nnn++] = (char*)nm;
                if (assume_known_nonpos(ctx, a) && nnp < MAX_SYM &&
                    !sym_already_listed(nonpos, nnp, nm)) nonpos[nnp++] = (char*)nm;
                if (assume_known_gt(ctx, a, neg1) && assume_known_lt(ctx, a, pos1) &&
                    nab < MAX_SYM && !sym_already_listed(abslt1, nab, nm)) abslt1[nab++] = (char*)nm;
            }
        }
        expr_free(neg1);
        expr_free(pos1);
    }

    /* Keep only symbols the target actually mentions: rules for the rest can
     * never match, and pruning them bounds the rule-string size (see the helper). */
    ar_filter_to_input(positives, &npos,  input);
    ar_filter_to_input(reals,     &nreal, input);
    ar_filter_to_input(integers,  &nint,  input);
    ar_filter_to_input(negatives, &nneg,  input);
    ar_filter_to_input(evens,     &neven, input);
    ar_filter_to_input(nonnegs,   &nnn,   input);
    ar_filter_to_input(nonpos,    &nnp,   input);
    ar_filter_to_input(abslt1,    &nab,   input);

    /* Build a single rule list "{r1, r2, ...}" as a string, then parse. */
    char buf[8192];
    size_t off = 0;
    int wrote_any = 0;

    #define EMIT(...) do { \
        int w = snprintf(buf + off, sizeof(buf) - off, __VA_ARGS__); \
        if (w < 0 || (size_t)w >= sizeof(buf) - off) goto overflow; \
        off += (size_t)w; \
    } while (0)
    #define SEP() do { if (wrote_any) EMIT(", "); wrote_any = 1; } while (0)

    EMIT("{");

    for (size_t i = 0; i < npos; i++) {
        const char* x = positives[i];
        /* Sqrt[x^2] forms */
        SEP(); EMIT("Power[Power[%s, 2], Rational[1, 2]] :> %s", x, x);
        SEP(); EMIT("Power[Power[%s, -1], Rational[1, 2]] :> Power[%s, Rational[-1, 2]]", x, x);
        SEP(); EMIT("Power[Power[%s, -2], Rational[1, 2]] :> Power[%s, -1]", x, x);
        /* General (x^2)^r -> x^(2r) for x > 0 (valid for any rational r, since
         * x^2 > 0): covers 1/Sqrt[x^2] = (x^2)^(-1/2) -> 1/x and the
         * (x^2)^(3/2) -> x^3 that arise from Laurent/fractional radical powers. */
        SEP(); EMIT("Power[Power[%s, 2], r_] :> Power[%s, 2 r]", x, x);
        /* Sqrt[x^2 * rest] -> x * Sqrt[rest] for x > 0; lets multi-factor
         * radicals like Sqrt[x^2 y^2] reduce one symbol at a time. */
        SEP(); EMIT("Power[Times[Power[%s, 2], rest___], Rational[1, 2]] :> %s Power[Times[rest], Rational[1, 2]]", x, x);
        /* Abs[x] -> x  for x > 0 */
        SEP(); EMIT("Abs[%s] :> %s", x, x);
        /* Sign[x] -> 1  for x > 0 */
        SEP(); EMIT("Sign[%s] :> 1", x);
        /* Arg[x] -> 0  for x > 0 (positive real) */
        SEP(); EMIT("Arg[%s] :> 0", x);
        /* Conjugate[x] -> x  for x > 0 (x is real) */
        SEP(); EMIT("Conjugate[%s] :> %s", x, x);
        /* Log[x^p] -> p Log[x]  for x > 0 (any real p; v1 accepts symbolic p too) */
        SEP(); EMIT("Log[Power[%s, p_]] :> p Log[%s]", x, x);
        /* Log[x rest] -> Log[x] + Log[rest]  for x > 0: a positive real factor
         * never shifts the branch (arg[x rest] = arg[rest]), so it splits off
         * exactly. Peels one positive factor at a time; Log[x] (bare) does not
         * re-match this Times pattern, so it terminates. */
        SEP(); EMIT("Log[Times[%s, rest___]] :> Log[%s] + Log[Times[rest]]", x, x);
        /* Inverse-trig sum identity: ArcTan[x] + ArcTan[1/x] -> Pi/2  for x > 0.
         * Mathilda's matcher does NOT perform orderless-Plus subset matching
         * out of the box (unlike Mathematica), so the rule must explicitly
         * absorb the trailing terms via `+ rest___` and re-emit them. The
         * BlankNullSequence pattern matches 0 or more remaining terms, so
         * the bare two-term sum reduces too. */
        SEP(); EMIT("ArcTan[%s] + ArcTan[Power[%s, -1]] + rest___ :> Pi/2 + rest", x, x);
        /* ArcCosh double-angle reduction: ArcCosh[2x^2 - 1] -> 2 ArcCosh[x]
         * for x > 0. The identity holds for any complex x with the
         * principal branch (verified for x in [0,1] via the i*ArcCos bridge
         * and for x >= 1 directly), so the x > 0 condition is the weakest
         * sufficient assumption. The Plus arg list is canonical
         * (Plus[-1, Times[2, x^2]]). */
        SEP(); EMIT("ArcCosh[Plus[-1, Times[2, Power[%s, 2]]]] :> 2 ArcCosh[%s]", x, x);
        /* General (x^m)^r -> x^(m r) for x > 0 and any m, r (covers
         * (x^3)^(1/3) -> x and the like). */
        SEP(); EMIT("Power[Power[%s, m_], r_] :> Power[%s, m r]", x, x);
        /* x^p y^p -> (x y)^p for two positive bases x, y (any common
         * exponent p). Emitted once per unordered pair. */
        for (size_t j = i + 1; j < npos; j++) {
            const char* y = positives[j];
            SEP(); EMIT("Times[Power[%s, p_], Power[%s, p_], rest___] :> "
                        "Power[Times[%s, %s], p] rest", x, y, x, y);
        }
    }

    for (size_t i = 0; i < nneg; i++) {
        const char* x = negatives[i];
        /* Abs[x] -> -x  for x < 0 */
        SEP(); EMIT("Abs[%s] :> -%s", x, x);
        /* Sign[x] -> -1  for x < 0 */
        SEP(); EMIT("Sign[%s] :> -1", x);
        /* Arg[x] -> Pi  for x < 0 (negative real, principal branch) */
        SEP(); EMIT("Arg[%s] :> Pi", x);
        /* Conjugate[x] -> x  for x < 0 (x is real) */
        SEP(); EMIT("Conjugate[%s] :> %s", x, x);
        /* Sqrt[x^2] -> -x  for x < 0  (Power[Power[x,2], 1/2]) */
        SEP(); EMIT("Power[Power[%s, 2], Rational[1, 2]] :> -%s", x, x);
        /* General (x^2)^r -> (-x)^(2r) for x < 0 (x^2 = (-x)^2, -x > 0):
         * covers 1/Sqrt[x^2] = (x^2)^(-1/2) -> 1/(-x) and (x^2)^(3/2) -> -x^3. */
        SEP(); EMIT("Power[Power[%s, 2], r_] :> Power[Times[-1, %s], 2 r]", x, x);
        /* Sqrt[x^2 * rest] -> -x * Sqrt[rest] for x < 0. */
        SEP(); EMIT("Power[Times[Power[%s, 2], rest___], Rational[1, 2]] :> -%s Power[Times[rest], Rational[1, 2]]", x, x);
        /* Mirror of the x > 0 inverse-trig sum: ArcTan[x] + ArcTan[1/x]
         * -> -Pi/2 for x < 0. Same `+ rest___` trick as above to handle
         * the embedded-in-larger-sum case. */
        SEP(); EMIT("ArcTan[%s] + ArcTan[Power[%s, -1]] + rest___ :> -Pi/2 + rest", x, x);
        /* Log[x] -> I Pi + Log[-x]  for x < 0 (principal branch). */
        SEP(); EMIT("Log[%s] :> I Pi + Log[-%s]", x, x);
    }

    /* Non-positive (x <= 0, but not proven strictly negative): Abs[x] -> -x,
     * Sqrt[x^2] -> -x, (x^2)^r -> (-x)^(2r) all hold at x <= 0 including x = 0
     * (where -x = 0 = |x|). Sign/Arg are NOT constant on x <= 0 (0 at x = 0),
     * so they are deliberately omitted. Skip symbols already strictly negative. */
    for (size_t i = 0; i < nnp; i++) {
        const char* x = nonpos[i];
        if (sym_already_listed(negatives, nneg, x)) continue;
        if (sym_already_listed(positives, npos, x)) continue;   /* 0 only: handled by nonneg too */
        SEP(); EMIT("Abs[%s] :> -%s", x, x);
        SEP(); EMIT("Conjugate[%s] :> %s", x, x);
        SEP(); EMIT("Power[Power[%s, 2], Rational[1, 2]] :> -%s", x, x);
        SEP(); EMIT("Power[Power[%s, 2], r_] :> Power[Times[-1, %s], 2 r]", x, x);
        SEP(); EMIT("Power[Times[Power[%s, 2], rest___], Rational[1, 2]] :> -%s Power[Times[rest], Rational[1, 2]]", x, x);
    }

    /* For real-but-unknown-sign, Sqrt[x^2] -> Abs[x]. Skip symbols already
     * proven positive, negative, or nonpositive (their stronger rule above wins). */
    for (size_t i = 0; i < nreal; i++) {
        const char* x = reals[i];
        if (sym_already_listed(positives, npos, x)) continue;
        if (sym_already_listed(negatives, nneg, x)) continue;
        if (sym_already_listed(nonpos, nnp, x)) continue;
        SEP(); EMIT("Power[Power[%s, 2], Rational[1, 2]] :> Abs[%s]", x, x);
        /* Sqrt[x^2 * rest] -> Abs[x] * Sqrt[rest] for real x. */
        SEP(); EMIT("Power[Times[Power[%s, 2], rest___], Rational[1, 2]] :> Abs[%s] Power[Times[rest], Rational[1, 2]]", x, x);
        /* Conjugate[x] -> x  for real x. */
        SEP(); EMIT("Conjugate[%s] :> %s", x, x);
        /* Log[x^2] -> 2 Log[Abs[x]] for real x (|x| keeps the argument real
         * and positive across the sign of x). */
        SEP(); EMIT("Log[Power[%s, 2]] :> 2 Log[Abs[%s]]", x, x);
        /* Log[E^x] -> x for real x (E^x is a positive real). */
        SEP(); EMIT("Log[Power[E, %s]] :> %s", x, x);
    }

    /* Sin[n Pi] -> 0, Cos[n Pi] -> (-1)^n, Tan[n Pi] -> 0 for integer n.
     * Plus: (-1)^(even_int * n) -> 1 and ((-1)^n)^even_int -> 1, so the
     * Cos rule can collapse all the way (the standalone Cos result
     * Cos[k Pi]^4 -> Power[-1, 4 k], for instance). */
    for (size_t i = 0; i < nint; i++) {
        const char* n = integers[i];
        SEP(); EMIT("Sin[%s Pi] :> 0", n);
        SEP(); EMIT("Sin[Pi %s] :> 0", n);
        SEP(); EMIT("Cos[%s Pi] :> Power[-1, %s]", n, n);
        SEP(); EMIT("Cos[Pi %s] :> Power[-1, %s]", n, n);
        SEP(); EMIT("Tan[%s Pi] :> 0", n);
        SEP(); EMIT("Tan[Pi %s] :> 0", n);
        SEP(); EMIT("Power[-1, Times[m_Integer /; EvenQ[m], %s]] :> 1", n);
        SEP(); EMIT("Power[Power[-1, %s], m_Integer /; EvenQ[m]] :> 1", n);
        /* Shift identities: Cos[x + n Pi] -> (-1)^n Cos[x], likewise Sin, and
         * Tan[x + n Pi] -> Tan[x] (period Pi), for integer n. The rest___
         * absorbs the remaining Plus terms (Plus[] -> 0 gives the bare
         * Cos[n Pi] case too). */
        SEP(); EMIT("Cos[Plus[Times[%s, Pi], rest___]] :> Power[-1, %s] Cos[Plus[rest]]", n, n);
        SEP(); EMIT("Sin[Plus[Times[%s, Pi], rest___]] :> Power[-1, %s] Sin[Plus[rest]]", n, n);
        SEP(); EMIT("Tan[Plus[Times[%s, Pi], rest___]] :> Tan[Plus[rest]]", n);
    }

    /* Even-exponent identities: (-1)^m = 1 when m is even, and the
     * lifted forms ((-1)^k)^m and (-1)^(k m) when additionally k is an
     * integer (so k m is also even). The pair-wise integer/even rules
     * cover the common Cos[k Pi]^m -> 1 path -- the existing Cos rule
     * above rewrites Cos[k Pi] to Power[-1, k], so the Power surface form
     * we land on is Power[Power[-1, k], m]. */
    for (size_t i = 0; i < neven; i++) {
        const char* m = evens[i];
        SEP(); EMIT("Power[-1, %s] :> 1", m);
        /* Literal-integer multiplier (handles concrete numbers without
         * needing them in the integer set). */
        SEP(); EMIT("Power[-1, k_Integer %s] :> 1", m);
        SEP(); EMIT("Power[Power[-1, k_Integer], %s] :> 1", m);
        for (size_t j = 0; j < nint; j++) {
            const char* k = integers[j];
            SEP(); EMIT("Power[-1, %s %s] :> 1", k, m);
            SEP(); EMIT("Power[Power[-1, %s], %s] :> 1", k, m);
        }
    }

    /* Non-negative (but not proven strictly positive) symbols: (x^m)^r ->
     * x^(m r), Sqrt[x^2] -> x, Abs[x] -> x all hold at x >= 0 including x = 0.
     * Strictly-positive symbols already received the stronger rule set above. */
    for (size_t i = 0; i < nnn; i++) {
        const char* x = nonnegs[i];
        if (sym_already_listed(positives, npos, x)) continue;
        if (sym_already_listed(negatives, nneg, x)) continue;
        SEP(); EMIT("Abs[%s] :> %s", x, x);
        SEP(); EMIT("Power[Power[%s, 2], Rational[1, 2]] :> %s", x, x);
        SEP(); EMIT("Power[Power[%s, m_], r_] :> Power[%s, m r]", x, x);
        SEP(); EMIT("Power[Times[Power[%s, 2], rest___], Rational[1, 2]] :> %s Power[Times[rest], Rational[1, 2]]", x, x);
    }

    /* Symbols with -1 < b < 1: (a^b)^c -> a^(b c) on the principal branch,
     * for any base a. */
    for (size_t i = 0; i < nab; i++) {
        const char* b = abslt1[i];
        SEP(); EMIT("Power[Power[base_, %s], c_] :> Power[base, %s c]", b, b);
    }

    /* Equal[u, v] facts -> two-way substitution rules. We use immediate
     * Rule (->) so the pattern uses exact structural matching. */
    for (size_t i = 0; i < ctx->count; i++) {
        const Expr* f = ctx->facts[i];
        if (!fact_is_function(f, "Equal", 2)) continue;
        /* We can't easily re-emit arbitrary Expr* into our string buffer;
         * instead, build these rules as Expr* and merge them in below. */
        (void)f; /* handled in the Expr* merge step below */
    }

    EMIT("}");

    if (!wrote_any) {
        /* No string-built rules. We may still have Equal substitutions. */
    }

    Expr* string_rules = wrote_any ? parse_expression(buf) : NULL;

    /* Now build Equal-substitution rules.
     *
     * Two complementary rules per equation:
     *   1. The direct rule heavier(lhs,rhs) -> lighter (catches cases
     *      where the equation's LHS appears verbatim as a subterm).
     *   2. ONE monomial-isolation rule when diff = lhs - rhs is a Plus
     *      with >= 3 terms: pick the heaviest non-numeric term t and
     *      emit t -> -(other terms). Polynomial relations like
     *      a^2 + b^2 == 1 then rewrite occurrences of a^2 even when
     *      the full "a^2 + b^2" sum is not present in the input.
     *
     * Emitting only one monomial rule (instead of one per term) avoids
     * the bidirectional cycle a^2 -> 1-b^2 ; b^2 -> 1-a^2 that
     * ReplaceRepeated would chase up to its 65536 iteration cap. */
    Expr** eq_diffs = (Expr**)calloc(ctx->count, sizeof(Expr*));
    size_t eq_count = 0;
    for (size_t i = 0; i < ctx->count; i++) {
        const Expr* f = ctx->facts[i];
        if (!fact_is_function(f, "Equal", 2)) continue;
        Expr* lhs = f->data.function.args[0];
        Expr* rhs = f->data.function.args[1];
        Expr* sub_args[2] = { expr_copy(lhs),
                              expr_new_function(expr_new_symbol(SYM_Times),
                                  (Expr*[]){ expr_new_integer(-1), expr_copy(rhs) }, 2) };
        Expr* sum = expr_new_function(expr_new_symbol(SYM_Plus), sub_args, 2);
        Expr* diff = eval_and_free(sum);
        eq_diffs[i] = diff;
        /* Always emit the direct heavier->lighter rule. */
        eq_count++;
        /* Plus extra monomial-isolation rule for polynomial relations. */
        if (diff->type == EXPR_FUNCTION &&
            diff->data.function.head &&
            diff->data.function.head->type == EXPR_SYMBOL &&
            diff->data.function.head->data.symbol.name == SYM_Plus &&
            diff->data.function.arg_count >= 3) {
            for (size_t j = 0; j < diff->data.function.arg_count; j++) {
                Expr* term = diff->data.function.args[j];
                if (term->type == EXPR_INTEGER || term->type == EXPR_BIGINT ||
                    term->type == EXPR_REAL) continue;
                eq_count++;
                break; /* one monomial rule per equation */
            }
        }
    }

    /* Apply the synthesized string + Equal-substitution rules (if any) to
     * obtain base_result; the structural pass below runs regardless. */
    Expr* base_result = NULL;
    if (string_rules || eq_count > 0) {
    size_t string_len = 0;
    if (string_rules && string_rules->type == EXPR_FUNCTION) {
        string_len = string_rules->data.function.arg_count;
    }
    size_t total = string_len + eq_count;
    Expr** all = (Expr**)calloc(total, sizeof(Expr*));
    size_t fill = 0;
    if (string_rules && string_rules->type == EXPR_FUNCTION) {
        for (size_t i = 0; i < string_len; i++) {
            all[fill++] = expr_copy(string_rules->data.function.args[i]);
        }
    }
    for (size_t i = 0; i < ctx->count; i++) {
        const Expr* f = ctx->facts[i];
        if (!fact_is_function(f, "Equal", 2)) continue;
        Expr* lhs = f->data.function.args[0];
        Expr* rhs = f->data.function.args[1];
        Expr* diff = eq_diffs[i];

        /* Direct heavier->lighter rule. */
        Expr *src, *dst;
        if (simp_default_complexity(lhs) >= simp_default_complexity(rhs)) {
            src = lhs; dst = rhs;
        } else {
            src = rhs; dst = lhs;
        }
        Expr* direct[2] = { expr_copy(src), expr_copy(dst) };
        all[fill++] = expr_new_function(expr_new_symbol(SYM_Rule), direct, 2);

        /* Polynomial-relation monomial-isolation rule (one per fact). */
        if (diff->type == EXPR_FUNCTION &&
            diff->data.function.head &&
            diff->data.function.head->type == EXPR_SYMBOL &&
            diff->data.function.head->data.symbol.name == SYM_Plus &&
            diff->data.function.arg_count >= 3) {
            size_t n = diff->data.function.arg_count;
            /* Pick the first non-numeric term, breaking ties by canonical
             * (Plus-Orderless) order which is already applied by the
             * evaluator. */
            size_t pick = (size_t)-1;
            size_t pick_score = 0;
            for (size_t j = 0; j < n; j++) {
                Expr* term = diff->data.function.args[j];
                if (term->type == EXPR_INTEGER || term->type == EXPR_BIGINT ||
                    term->type == EXPR_REAL) continue;
                size_t s = simp_default_complexity(term);
                if (pick == (size_t)-1 || s > pick_score) {
                    pick = j;
                    pick_score = s;
                }
            }
            if (pick != (size_t)-1) {
                Expr* term = diff->data.function.args[pick];
                Expr** other_args = (Expr**)calloc(n - 1, sizeof(Expr*));
                size_t oi = 0;
                for (size_t k = 0; k < n; k++) {
                    if (k == pick) continue;
                    other_args[oi++] = expr_new_function(expr_new_symbol(SYM_Times),
                        (Expr*[]){ expr_new_integer(-1),
                                   expr_copy(diff->data.function.args[k]) }, 2);
                }
                Expr* iso_rhs;
                if (n - 1 == 1) {
                    iso_rhs = other_args[0];
                    free(other_args);
                } else {
                    iso_rhs = expr_new_function(expr_new_symbol(SYM_Plus), other_args, n - 1);
                    free(other_args);
                }
                Expr* iso[2] = { expr_copy(term), iso_rhs };
                all[fill++] = expr_new_function(expr_new_symbol(SYM_Rule), iso, 2);
            }
        }
    }
    for (size_t i = 0; i < ctx->count; i++) if (eq_diffs[i]) expr_free(eq_diffs[i]);
    free(eq_diffs);
    if (string_rules) expr_free(string_rules);

    Expr* rules_list = expr_new_function(expr_new_symbol(SYM_List), all, fill);
    free(all);

    Expr* call_args[2] = { expr_copy((Expr*)input), rules_list };
    Expr* call = expr_new_function(expr_new_symbol(SYM_ReplaceRepeated), call_args, 2);
    Expr* out = evaluate(call);
    expr_free(call);
    base_result = out;
    } else {
        free(eq_diffs);   /* eq_count == 0: no entries were populated */
    }

    /* Structural assumption rewrites (Floor/Ceiling/Round/IntegerPart/
     * FractionalPart/Mod, Re/Im/Conjugate/Arg, ArcTan[Tan]) run whether or not
     * any string/Equal rule was generated, so a lone structural fact still
     * fires and Simplify inherits the same rewrites. */
    {
        const Expr* walk_in = base_result ? base_result : input;
        int schanged = 0;
        Expr* walked = assume_structural_rewrite(walk_in, ctx, &schanged);
        if (schanged && walked) {
            Expr* ev = eval_and_free(walked);
            if (base_result) expr_free(base_result);
            return ev;
        }
        if (walked) expr_free(walked);
    }
    return base_result;

overflow:
    /* Buffer was too small; bail out, no rules applied. */
    return NULL;

    #undef EMIT
    #undef SEP
}

