/* product_util.c -- shared helpers for the Product sub-algorithms. */

#include "product_internal.h"
#include "sym_names.h"
#include "eval.h"
#include "expr.h"
#include "poly.h"
#include <string.h>
#include <stdlib.h>

Expr* prod_eval(const char* head, Expr** args, size_t n) {
    Expr* call = expr_new_function(expr_new_symbol(head), args, n);
    Expr* r = evaluate(call);
    expr_free(call);
    return r;
}

Expr* prod_int(int64_t v) { return expr_new_integer(v); }

Expr* prod_subst(Expr* e, Expr* var, Expr* val) {
    Expr* rule = expr_new_function(expr_new_symbol(SYM_Rule),
                    (Expr*[]){ expr_copy(var), expr_copy(val) }, 2);
    Expr* args[2] = { expr_copy(e), rule };
    return prod_eval("ReplaceAll", args, 2);
}

Expr* prod_factor(Expr* e) {
    Expr* arg[1] = { expr_copy(e) };
    Expr* r = prod_eval("Factor", arg, 1);
    /* If Factor came back unevaluated (head Factor), keep the input form. */
    if (r && r->type == EXPR_FUNCTION
          && r->data.function.head->type == EXPR_SYMBOL
          && strcmp(r->data.function.head->data.symbol, "Factor") == 0) {
        expr_free(r);
        return expr_copy(e);
    }
    return r;
}

bool prod_free_of(Expr* e, Expr* var) {
    Expr* args[2] = { expr_copy(e), expr_copy(var) };
    Expr* r = prod_eval("FreeQ", args, 2);
    bool yes = (r && r->type == EXPR_SYMBOL && strcmp(r->data.symbol, "True") == 0);
    if (r) expr_free(r);
    return yes;
}

/* a / b  =  Times[a, Power[b, -1]], evaluated. */
Expr* prod_div(Expr* a, Expr* b) {
    Expr* inv = expr_new_function(expr_new_symbol(SYM_Power),
                    (Expr*[]){ expr_copy(b), expr_new_integer(-1) }, 2);
    Expr* args[2] = { expr_copy(a), inv };
    return prod_eval("Times", args, 2);
}

bool product_stage_args(Expr* res, Expr** f, Expr** var,
                        Expr** imin, Expr** imax, bool* definite) {
    if (res->type != EXPR_FUNCTION) return false;
    size_t argc = res->data.function.arg_count;
    Expr** a = res->data.function.args;
    if (argc != 2 && argc != 4) return false;
    if (a[1]->type != EXPR_SYMBOL) return false;   /* product variable */
    *f = a[0];
    *var = a[1];
    if (argc == 4) {
        *imin = a[2];
        *imax = a[3];
        *definite = true;
    } else {
        *imin = NULL;
        *imax = NULL;
        *definite = false;
    }
    return true;
}

/* ------------------------------------------------------------------ */
/*  Symbolic-power guard                                               */
/* ------------------------------------------------------------------ */

/* Structural test: does the interned symbol `v` appear anywhere in e? */
static bool contains_sym(const Expr* e, const char* v) {
    if (!e) return false;
    if (e->type == EXPR_SYMBOL) return e->data.symbol == v;
    if (e->type == EXPR_FUNCTION) {
        if (contains_sym(e->data.function.head, v)) return true;
        for (size_t i = 0; i < e->data.function.arg_count; i++)
            if (contains_sym(e->data.function.args[i], v)) return true;
    }
    return false;
}

bool prod_has_symbolic_power(Expr* e, Expr* var) {
    if (!e || e->type != EXPR_FUNCTION) return false;
    const char* v = var->data.symbol;
    Expr* head = e->data.function.head;
    if (head->type == EXPR_SYMBOL && head->data.symbol == SYM_Power
            && e->data.function.arg_count == 2) {
        if (contains_sym(e->data.function.args[1], v)) return true;
    }
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (prod_has_symbolic_power(e->data.function.args[i], var)) return true;
    return false;
}

/* ------------------------------------------------------------------ */
/*  Linear-factor extraction over Q                                    */
/* ------------------------------------------------------------------ */

/* Evaluate Expand[e]. */
static Expr* prod_expand(Expr* e) {
    Expr* arg[1] = { expr_copy(e) };
    return prod_eval("Expand", arg, 1);
}

/* lead <- evaluate(Times[lead, factor]); consumes neither argument's copy
 * semantics beyond the standard prod_eval (lead is freed and replaced). */
static void lead_mul(Expr** lead, Expr* factor) {
    Expr* args[2] = { *lead, expr_copy(factor) };  /* adopt *lead directly */
    Expr* r = prod_eval("Times", args, 2);
    *lead = r;
}

/* Append a (root, mult) pair, growing the arrays. */
static void push_root(Expr*** roots, int** mults, size_t* n, size_t* cap,
                      Expr* root, int mult) {
    if (*n == *cap) {
        *cap = *cap ? *cap * 2 : 8;
        *roots = realloc(*roots, sizeof(Expr*) * (*cap));
        *mults = realloc(*mults, sizeof(int) * (*cap));
    }
    (*roots)[*n] = root;
    (*mults)[*n] = mult;
    (*n)++;
}

/* Process one factor `base` raised to integer multiplicity `m`, updating the
 * leading constant, root list, and all_linear flag. */
static void process_factor(Expr* base, Expr* var, int m,
                           Expr** lead, Expr*** roots, int** mults,
                           size_t* n, size_t* cap, bool* all_linear) {
    Expr* eb = prod_expand(base);
    int dg = get_degree_poly(eb, var);
    if (dg == 0) {
        /* Constant in var (may be a symbolic constant): fold into lead^m. */
        Expr* p = expr_new_function(expr_new_symbol(SYM_Power),
                      (Expr*[]){ expr_copy(eb), expr_new_integer(m) }, 2);
        Expr* pe = evaluate(p);
        expr_free(p);
        lead_mul(lead, pe);
        expr_free(pe);
    } else if (dg == 1) {
        Expr* c1 = get_coeff(eb, var, 1);
        Expr* c0 = get_coeff(eb, var, 0);
        /* root r = -c0 / c1 */
        Expr* negc0 = expr_new_function(expr_new_symbol(SYM_Times),
                          (Expr*[]){ expr_new_integer(-1), expr_copy(c0) }, 2);
        Expr* r = prod_div(negc0, c1);
        expr_free(negc0);
        /* leading-coefficient contribution c1^m */
        Expr* p = expr_new_function(expr_new_symbol(SYM_Power),
                      (Expr*[]){ expr_copy(c1), expr_new_integer(m) }, 2);
        Expr* pe = evaluate(p);
        expr_free(p);
        lead_mul(lead, pe);
        expr_free(pe);
        expr_free(c1); expr_free(c0);
        push_root(roots, mults, n, cap, r, m);
    } else {
        /* Irreducible quadratic+ : not representable as a linear chain. */
        *all_linear = false;
    }
    expr_free(eb);
}

/* Process one factor `base` raised to integer multiplicity `m`, allowing
 * irreducible QUADRATICS (complex-conjugate roots via the quadratic formula) in
 * addition to linear factors.  Degree >= 3 irreducible factors set *ok=false. */
static void process_factor_c(Expr* base, Expr* var, int m,
                             Expr** lead, Expr*** roots, int** mults,
                             size_t* n, size_t* cap, bool* ok) {
    Expr* eb = prod_expand(base);
    int dg = get_degree_poly(eb, var);
    if (dg == 0) {
        Expr* p = expr_new_function(expr_new_symbol(SYM_Power),
                      (Expr*[]){ expr_copy(eb), expr_new_integer(m) }, 2);
        Expr* pe = evaluate(p); expr_free(p);
        lead_mul(lead, pe); expr_free(pe);
    } else if (dg == 1) {
        Expr* c1 = get_coeff(eb, var, 1);
        Expr* c0 = get_coeff(eb, var, 0);
        Expr* negc0 = expr_new_function(expr_new_symbol(SYM_Times),
                          (Expr*[]){ expr_new_integer(-1), expr_copy(c0) }, 2);
        Expr* r = prod_div(negc0, c1); expr_free(negc0);
        Expr* p = expr_new_function(expr_new_symbol(SYM_Power),
                      (Expr*[]){ expr_copy(c1), expr_new_integer(m) }, 2);
        Expr* pe = evaluate(p); expr_free(p);
        lead_mul(lead, pe); expr_free(pe);
        expr_free(c1); expr_free(c0);
        push_root(roots, mults, n, cap, r, m);
    } else if (dg == 2) {
        /* c2 k^2 + c1 k + c0 = c2 (k - r+)(k - r-),
         * r± = (-c1 ± Sqrt[c1^2 - 4 c2 c0]) / (2 c2). */
        Expr* c2 = get_coeff(eb, var, 2);
        Expr* c1 = get_coeff(eb, var, 1);
        Expr* c0 = get_coeff(eb, var, 0);
        /* disc = c1^2 - 4 c2 c0 */
        Expr* c1sq = expr_new_function(expr_new_symbol(SYM_Power),
                         (Expr*[]){ expr_copy(c1), expr_new_integer(2) }, 2);
        Expr* four_c2c0 = expr_new_function(expr_new_symbol(SYM_Times),
                              (Expr*[]){ expr_new_integer(-4), expr_copy(c2),
                                         expr_copy(c0) }, 3);
        Expr* disc = expr_new_function(expr_new_symbol(SYM_Plus),
                         (Expr*[]){ c1sq, four_c2c0 }, 2);         /* adopts both */
        Expr* sq = expr_new_function(expr_new_symbol("Sqrt"),
                       (Expr*[]){ disc }, 1);                      /* adopts disc */
        Expr* sqe = evaluate(sq); expr_free(sq);
        Expr* negc1 = expr_new_function(expr_new_symbol(SYM_Times),
                          (Expr*[]){ expr_new_integer(-1), expr_copy(c1) }, 2);
        Expr* twoc2 = expr_new_function(expr_new_symbol(SYM_Times),
                          (Expr*[]){ expr_new_integer(2), expr_copy(c2) }, 2);
        Expr* twoc2e = evaluate(twoc2); expr_free(twoc2);
        /* r+ = (negc1 + sqe)/twoc2, r- = (negc1 - sqe)/twoc2 */
        Expr* np = expr_new_function(expr_new_symbol(SYM_Plus),
                       (Expr*[]){ expr_copy(negc1), expr_copy(sqe) }, 2);
        Expr* rp = prod_div(np, twoc2e); expr_free(np);
        Expr* nm = expr_new_function(expr_new_symbol(SYM_Plus),
                       (Expr*[]){ expr_copy(negc1),
                                  expr_new_function(expr_new_symbol(SYM_Times),
                                      (Expr*[]){ expr_new_integer(-1), expr_copy(sqe) }, 2) }, 2);
        Expr* rm = prod_div(nm, twoc2e); expr_free(nm);
        expr_free(negc1); expr_free(sqe); expr_free(twoc2e);
        /* lead contribution c2^m */
        Expr* p = expr_new_function(expr_new_symbol(SYM_Power),
                      (Expr*[]){ expr_copy(c2), expr_new_integer(m) }, 2);
        Expr* pe = evaluate(p); expr_free(p);
        lead_mul(lead, pe); expr_free(pe);
        expr_free(c2); expr_free(c1); expr_free(c0);
        push_root(roots, mults, n, cap, rp, m);
        push_root(roots, mults, n, cap, rm, m);
    } else {
        *ok = false;   /* irreducible cubic+ : not handled */
    }
    expr_free(eb);
}

bool prod_rational_roots(Expr* e, Expr* var,
                         Expr** lead_out,
                         Expr*** roots_out, int** mults_out, size_t* n_out,
                         bool* ok_out) {
    Expr* vars[1] = { var };
    if (!is_polynomial(e, vars, 1)) return false;

    Expr* lead = expr_new_integer(1);
    Expr** roots = NULL;
    int* mults = NULL;
    size_t n = 0, cap = 0;
    bool ok = true;

    Expr* F = prod_factor(e);
    bool is_times = (F->type == EXPR_FUNCTION
                     && F->data.function.head->type == EXPR_SYMBOL
                     && F->data.function.head->data.symbol == SYM_Times);
    size_t fc = is_times ? F->data.function.arg_count : 1;
    for (size_t i = 0; i < fc; i++) {
        Expr* fac = is_times ? F->data.function.args[i] : F;
        Expr* base = fac;
        int m = 1;
        if (fac->type == EXPR_FUNCTION
                && fac->data.function.head->type == EXPR_SYMBOL
                && fac->data.function.head->data.symbol == SYM_Power
                && fac->data.function.arg_count == 2
                && fac->data.function.args[1]->type == EXPR_INTEGER) {
            base = fac->data.function.args[0];
            m = (int)fac->data.function.args[1]->data.integer;
        }
        if (m <= 0) { ok = false; continue; }
        if (base->type == EXPR_INTEGER || base->type == EXPR_BIGINT
                || base->type == EXPR_REAL
                || (base->type == EXPR_FUNCTION
                    && base->data.function.head->type == EXPR_SYMBOL
                    && base->data.function.head->data.symbol == SYM_Rational)) {
            Expr* p = expr_new_function(expr_new_symbol(SYM_Power),
                          (Expr*[]){ expr_copy(base), expr_new_integer(m) }, 2);
            Expr* pe = evaluate(p); expr_free(p);
            lead_mul(&lead, pe); expr_free(pe);
            continue;
        }
        process_factor_c(base, var, m, &lead, &roots, &mults, &n, &cap, &ok);
    }
    expr_free(F);

    *lead_out = lead;
    *roots_out = roots;
    *mults_out = mults;
    *n_out = n;
    *ok_out = ok;
    return true;
}

bool prod_linear_factors(Expr* e, Expr* var,
                         Expr** lead_out,
                         Expr*** roots_out, int** mults_out, size_t* n_out,
                         bool* all_linear_out) {
    Expr* vars[1] = { var };
    if (!is_polynomial(e, vars, 1)) return false;

    Expr* lead = expr_new_integer(1);
    Expr** roots = NULL;
    int* mults = NULL;
    size_t n = 0, cap = 0;
    bool all_linear = true;

    Expr* F = prod_factor(e);   /* product of irreducible factors over Q */

    /* Iterate the top-level factors of F (Times[...] or a single factor). */
    bool is_times = (F->type == EXPR_FUNCTION
                     && F->data.function.head->type == EXPR_SYMBOL
                     && F->data.function.head->data.symbol == SYM_Times);
    size_t fc = is_times ? F->data.function.arg_count : 1;
    for (size_t i = 0; i < fc; i++) {
        Expr* fac = is_times ? F->data.function.args[i] : F;
        Expr* base = fac;
        int m = 1;
        if (fac->type == EXPR_FUNCTION
                && fac->data.function.head->type == EXPR_SYMBOL
                && fac->data.function.head->data.symbol == SYM_Power
                && fac->data.function.arg_count == 2
                && fac->data.function.args[1]->type == EXPR_INTEGER) {
            base = fac->data.function.args[0];
            m = (int)fac->data.function.args[1]->data.integer;
        }
        if (m <= 0) { all_linear = false; continue; }  /* shouldn't happen for a poly */
        /* A bare number factor (no var) folds into lead. */
        if (base->type == EXPR_INTEGER || base->type == EXPR_BIGINT
                || base->type == EXPR_REAL
                || (base->type == EXPR_FUNCTION
                    && base->data.function.head->type == EXPR_SYMBOL
                    && base->data.function.head->data.symbol == SYM_Rational)) {
            Expr* p = expr_new_function(expr_new_symbol(SYM_Power),
                          (Expr*[]){ expr_copy(base), expr_new_integer(m) }, 2);
            Expr* pe = evaluate(p);
            expr_free(p);
            lead_mul(&lead, pe);
            expr_free(pe);
            continue;
        }
        process_factor(base, var, m, &lead, &roots, &mults, &n, &cap, &all_linear);
    }
    expr_free(F);

    *lead_out = lead;
    *roots_out = roots;
    *mults_out = mults;
    *n_out = n;
    *all_linear_out = all_linear;
    return true;
}
