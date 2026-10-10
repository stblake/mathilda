/* integrate_linarg.c
 *
 * Linear-argument substitution stage.  See integrate_linarg.h for the rationale.
 *
 * Reduces  Integrate[f(a*x+b), x]  to  (1/a) * (Integrate[f(u), u] /. u->a*x+b)
 * when every occurrence of x in f sits inside a trig/hyperbolic kernel sharing
 * one linear argument w = a*x + b (a a non-zero number, a != 1 or b != 0), and
 * at least one such kernel appears in a denominator.  The substitution is an
 * exact identity, so the result is correct by construction whenever the bare
 * sub-integral closes; otherwise the stage declines and the cascade continues to
 * Weierstrass.
 *
 * Memory: a builtin stage takes no ownership of its arguments -- f and x are
 * borrowed.  Every Expr built here is freed explicitly; eval_take consumes its
 * argument.  We never expr_free(f) or expr_free(x).
 */

#include "integrate_linarg.h"

#include "expr.h"
#include "eval.h"
#include "symtab.h"
#include "common.h"
#include "internal.h"
#include "sym_names.h"

#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

/* ---------------------------------------------------------------------- */
/* Small builders / evaluation helpers (mirror integrate_jeffrey.c)       */
/* ---------------------------------------------------------------------- */

static Expr* mk_int(int64_t v) { return expr_new_integer(v); }

static Expr* mk_fn2(const char* name, Expr* a, Expr* b) {
    Expr* args[2] = { a, b };
    return expr_new_function(expr_new_symbol(name), args, 2);
}
static Expr* mk_plus2(Expr* a, Expr* b)  { return mk_fn2(SYM_Plus,  a, b); }
static Expr* mk_times2(Expr* a, Expr* b) { return mk_fn2(SYM_Times, a, b); }
static Expr* mk_pow(Expr* a, Expr* b)    { return mk_fn2(SYM_Power, a, b); }

/* Evaluate `call` to a fixed point, freeing `call`. */
static Expr* eval_take(Expr* call) {
    Expr* r = evaluate(call);
    expr_free(call);
    return r;
}

/* True if `f` contains no subexpression structurally equal to `x`. */
static bool expr_free_of(const Expr* f, const Expr* x) {
    if (expr_eq((Expr*)f, (Expr*)x)) return false;
    if (f->type != EXPR_FUNCTION) return true;
    if (!expr_free_of(f->data.function.head, x)) return false;
    for (size_t i = 0; i < f->data.function.arg_count; i++)
        if (!expr_free_of(f->data.function.args[i], x)) return false;
    return true;
}

/* ---------------------------------------------------------------------- */
/* Kernel classification                                                  */
/* ---------------------------------------------------------------------- */

static bool is_trig_head(const char* h) {
    return h == SYM_Sin || h == SYM_Cos || h == SYM_Tan ||
           h == SYM_Cot || h == SYM_Sec || h == SYM_Csc;
}
static bool is_hyp_head(const char* h) {
    return h == SYM_Sinh || h == SYM_Cosh || h == SYM_Tanh ||
           h == SYM_Coth || h == SYM_Sech || h == SYM_Csch;
}
static bool is_kernel_head(const char* h) {
    return is_trig_head(h) || is_hyp_head(h);
}
/* A kernel that is itself a quotient (Tan/Cot/Sec/Csc + hyperbolics): its
 * presence means a trig kernel sits in a denominator. */
static bool is_recip_kernel_head(const char* h) {
    return h == SYM_Tan || h == SYM_Cot || h == SYM_Sec || h == SYM_Csc ||
           h == SYM_Tanh || h == SYM_Coth || h == SYM_Sech || h == SYM_Csch;
}

/* True if `e` is a unary trig/hyp kernel H[arg] with arg depending on x. */
static bool is_xkernel(const Expr* e, const Expr* x) {
    if (!e || e->type != EXPR_FUNCTION || e->data.function.arg_count != 1)
        return false;
    const Expr* h = e->data.function.head;
    if (h->type != EXPR_SYMBOL || !is_kernel_head(h->data.symbol.name))
        return false;
    return !expr_free_of(e->data.function.args[0], x);
}

/* Does `e` carry a trig/hyp kernel of x in a denominator?  Either a reciprocal
 * kernel (Sec/Csc/Tan/Cot[...x...]) anywhere, or a negative integer power of any
 * kernel of x.  This restricts the stage to the genuine rational-trig domain and
 * leaves polynomial trig (Sin[2x]^3) to the cleaner table/linearity stages. */
static bool has_denominator_trig(const Expr* e, const Expr* x) {
    if (!e || e->type != EXPR_FUNCTION) return false;
    const Expr* h = e->data.function.head;
    if (h->type == EXPR_SYMBOL) {
        const char* n = h->data.symbol.name;
        if (e->data.function.arg_count == 1 && is_recip_kernel_head(n) &&
            !expr_free_of(e->data.function.args[0], x))
            return true;
        if (n == SYM_Power && e->data.function.arg_count == 2) {
            const Expr* base = e->data.function.args[0];
            const Expr* ex   = e->data.function.args[1];
            if (is_xkernel(base, x) && ex->type == EXPR_INTEGER &&
                ex->data.integer < 0)
                return true;
        }
    }
    if (has_denominator_trig(h, x)) return true;
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (has_denominator_trig(e->data.function.args[i], x)) return true;
    return false;
}

/* Walk `e`, requiring every occurrence of x to sit inside a trig/hyp kernel, and
 * all such kernels to share ONE argument `w`.  *w points at the (borrowed) common
 * argument once found; *found becomes true; *bad becomes true on a bare x outside
 * a kernel or a second, different kernel argument. */
static void collect_linarg(const Expr* e, const Expr* x,
                           const Expr** w, bool* found, bool* bad) {
    if (*bad) return;
    if (expr_free_of(e, x)) return;                 /* no x here: nothing to do  */
    if (is_xkernel(e, x)) {                          /* shields its argument      */
        const Expr* arg = e->data.function.args[0];
        if (!*found) { *w = arg; *found = true; }
        else if (!expr_eq((Expr*)*w, (Expr*)arg)) *bad = true;
        return;                                      /* do NOT recurse into arg   */
    }
    if (e->type == EXPR_SYMBOL) {                    /* a bare x outside a kernel */
        if (e->data.symbol.name == x->data.symbol.name) *bad = true;
        return;
    }
    if (e->type != EXPR_FUNCTION) return;
    collect_linarg(e->data.function.head, x, w, found, bad);
    for (size_t i = 0; i < e->data.function.arg_count && !*bad; i++)
        collect_linarg(e->data.function.args[i], x, w, found, bad);
}

/* Is `e` a non-zero numeric constant (Integer / Real / BigInt / Rational)? */
static bool is_number(const Expr* e) {
    if (!e) return false;
    if (e->type == EXPR_INTEGER || e->type == EXPR_REAL ||
        e->type == EXPR_BIGINT) return true;
    return (e->type == EXPR_FUNCTION && e->data.function.head->type == EXPR_SYMBOL
            && e->data.function.head->data.symbol.name == SYM_Rational);
}
static bool num_equals_int(const Expr* e, int64_t v) {
    if (e->type == EXPR_INTEGER) return e->data.integer == v;
    if (e->type == EXPR_REAL)    return e->data.real == (double)v;
    return false;
}

/* True if `r` still contains an unevaluated Integrate head anywhere. */
static bool contains_integrate(const Expr* r) {
    if (!r || r->type != EXPR_FUNCTION) return false;
    const Expr* h = r->data.function.head;
    if (h->type == EXPR_SYMBOL && h->data.symbol.name == SYM_Integrate)
        return true;
    if (contains_integrate(h)) return true;
    for (size_t i = 0; i < r->data.function.arg_count; i++)
        if (contains_integrate(r->data.function.args[i])) return true;
    return false;
}

/* A fresh symbol guaranteed absent from the current symbol table. */
static Expr* fresh_symbol(void) {
    static uint64_t ctr = 0;
    char buf[64];
    for (;;) {
        snprintf(buf, sizeof(buf), "Integrate$lin$%llu", (unsigned long long)ctr++);
        if (!symtab_lookup(buf)) break;
    }
    return expr_new_symbol(buf);
}

/* ReplaceAll[expr, Rule[lhs, rhs]]; consumes expr, lhs, rhs; returns owned. */
static Expr* subst_take(Expr* expr, Expr* lhs, Expr* rhs) {
    Expr* rule = mk_fn2(SYM_Rule, lhs, rhs);
    return eval_take(internal_replace_all((Expr*[]){ expr, rule }, 2));
}

/* ---------------------------------------------------------------------- */

#define LINARG_MAX_DEPTH 8
static int linarg_depth = 0;

Expr* integrate_linarg_try(Expr* f, Expr* x) {
    if (!f || !x || x->type != EXPR_SYMBOL) return NULL;
    if (linarg_depth >= LINARG_MAX_DEPTH) return NULL;

    /* 1. Gate on the genuine rational-trig-with-denominator domain. */
    if (!has_denominator_trig(f, x)) return NULL;

    /* 2. Every x inside one shared kernel argument w. */
    const Expr* w = NULL; bool found = false, bad = false;
    collect_linarg(f, x, &w, &found, &bad);
    if (!found || bad || !w) return NULL;

    /* 3. w = a*x + b, a a non-zero number; non-trivial (a != 1 or b != 0). */
    Expr* a = eval_take(mk_fn2(SYM_D, expr_copy((Expr*)w), expr_copy(x)));
    if (!a || !is_number(a) || num_equals_int(a, 0)) { expr_free(a); return NULL; }
    Expr* b = subst_take(expr_copy((Expr*)w), expr_copy(x), mk_int(0));
    if (!b) { expr_free(a); return NULL; }
    bool trivial = num_equals_int(a, 1) &&
                   (is_number(b) ? num_equals_int(b, 0)
                                 : false);
    if (trivial) { expr_free(a); expr_free(b); return NULL; }

    /* 4. Substitute x -> (u - b)/a, making the integrand a function of bare u. */
    Expr* u = fresh_symbol();
    Expr* repl = mk_times2(mk_plus2(expr_copy(u),
                                    mk_times2(mk_int(-1), expr_copy(b))),
                           mk_pow(expr_copy(a), mk_int(-1)));
    Expr* h = subst_take(expr_copy(f), expr_copy(x), repl);
    if (!h || !expr_free_of(h, x)) {
        expr_free(a); expr_free(b); expr_free(u); expr_free(h); return NULL;
    }

    /* 5. Close the bare-argument sub-integral, recursing into the full cascade. */
    linarg_depth++;
    Expr* Hu = eval_take(mk_fn2(SYM_Integrate, h, expr_copy(u)));  /* consumes h */
    linarg_depth--;
    if (!Hu || contains_integrate(Hu)) {
        expr_free(a); expr_free(b); expr_free(u); expr_free(Hu); return NULL;
    }

    /* 6. Back-substitute u -> a*x + b and divide by a. */
    Expr* axb = mk_plus2(mk_times2(expr_copy(a), expr_copy(x)), expr_copy(b));
    Expr* Hx  = subst_take(Hu, expr_copy(u), axb);                 /* consumes Hu */
    Expr* result = eval_take(mk_times2(mk_pow(expr_copy(a), mk_int(-1)), Hx));

    expr_free(a); expr_free(b); expr_free(u);
    return result;
}
