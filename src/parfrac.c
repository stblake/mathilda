#include "common.h"      /* alg_thread_call */
#include "expr.h"
#include "pack.h"
#include "symtab.h"
#include "eval.h"
#include "parse.h"
#include "print.h"
#include "poly.h"
#include "attr.h"
#include "arithmetic.h"
#include "expand.h"
#include "core.h"
#include "rationalize.h"
#include "sym_names.h"
#include "options.h"
#include "qafactor.h"
#include "flint_bridge.h"
#include "ratfun.h"
#include "sym_intern.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <stdbool.h>

/* Inner implementation of Apart.  Receives the option-stripped arg count
 * `apart_argc` ∈ {1, 2} and a borrowed `apart_alpha` (possibly NULL, possibly
 * coming from extension_autodetect).  Defined below `builtin_apart`. */
static Expr* apart_impl(Expr* res, size_t apart_argc, const Expr* apart_alpha);

/* True if `e` contains a symbolic radical Power[base, Rational[p,q]] (q >= 2)
 * whose base is not an integer — i.e. a Sqrt[k]-style radical over a parameter.
 * These are the coefficients that make the classical partial-fraction RowReduce
 * (symbolic Gaussian elimination) blow up; the residue-formula fast path is
 * gated on this so plain-rational Apart keeps its existing output. */
static bool apart_has_symbolic_radical(const Expr* e) {
    if (!e || e->type != EXPR_FUNCTION) return false;
    const Expr* h = e->data.function.head;
    if (h && h->type == EXPR_SYMBOL && h->data.symbol.name == SYM_Power
        && e->data.function.arg_count == 2) {
        const Expr* base = e->data.function.args[0];
        const Expr* exp  = e->data.function.args[1];
        if (exp->type == EXPR_FUNCTION && exp->data.function.head
            && exp->data.function.head->type == EXPR_SYMBOL
            && exp->data.function.head->data.symbol.name == SYM_Rational
            && base->type != EXPR_INTEGER && base->type != EXPR_BIGINT) {
            return true;
        }
    }
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (apart_has_symbolic_radical(e->data.function.args[i])) return true;
    return false;
}

/* ------------------------------------------------------------------ */
/*  Radical coefficients                                               */
/* ------------------------------------------------------------------ */

/* A denominator already split over radicals keeps its factors, as in the
 * reference: Apart[c/((z - (-1)^(1/3)) (z + (-1)^(1/3)) ...), z] is a sum of
 * linear terms over Q((-1)^(1/3)), and Apart[1/((z - Sqrt[2])(z + Sqrt[2])), z]
 * splits in two.  The rational engine expands the denominator over Q and
 * refactors it, merging the factors back (1/(-2 + z^2)).  Each radical atom
 * -- a number to a non-integer rational power -- is replaced by a fresh
 * symbol, the expression is apart-ed as if the radicals were parameters, and
 * they are put back: a partial-fraction identity in independent symbols
 * stays one under any substitution that keeps its denominators nonzero. */
typedef struct { Expr** from; Expr** to; size_t n, cap; } ApRad;

static bool ap_is_radical(const Expr* e) {
    if (e->type != EXPR_FUNCTION || e->data.function.head->type != EXPR_SYMBOL
        || e->data.function.head->data.symbol.name != SYM_Power || e->data.function.arg_count != 2)
        return false;
    const Expr* b = e->data.function.args[0];
    const Expr* x = e->data.function.args[1];
    bool num_base = b->type == EXPR_INTEGER || b->type == EXPR_BIGINT
        || (b->type == EXPR_FUNCTION && b->data.function.head->type == EXPR_SYMBOL
            && b->data.function.head->data.symbol.name == SYM_Rational);
    bool frac_exp = x->type == EXPR_FUNCTION && x->data.function.head->type == EXPR_SYMBOL
        && x->data.function.head->data.symbol.name == SYM_Rational;
    return num_base && frac_exp;
}

static Expr* ap_abstract(const Expr* e, ApRad* t) {
    if (ap_is_radical(e)) {
        for (size_t i = 0; i < t->n; i++)
            if (expr_eq(t->from[i], (Expr*)e)) return expr_copy(t->to[i]);
        if (t->n == t->cap) {
            t->cap = t->cap ? 2 * t->cap : 8;
            t->from = realloc(t->from, sizeof(Expr*) * t->cap);
            t->to = realloc(t->to, sizeof(Expr*) * t->cap);
        }
        char name[48];
        snprintf(name, sizeof name, "Apart`Private`r%zu", t->n + 1);
        t->from[t->n] = expr_copy((Expr*)e);
        t->to[t->n] = expr_new_symbol(intern_symbol(name));
        return expr_copy(t->to[t->n++]);
    }
    if (e->type != EXPR_FUNCTION) return expr_copy((Expr*)e);
    size_t n = e->data.function.arg_count;
    Expr** a = malloc(sizeof(Expr*) * (n ? n : 1));
    for (size_t i = 0; i < n; i++) a[i] = ap_abstract(e->data.function.args[i], t);
    Expr* out = expr_new_function(expr_copy(e->data.function.head), a, n);
    free(a);
    return out;
}

static Expr* ap_restore(const Expr* e, const ApRad* t) {
    if (e->type == EXPR_SYMBOL) {
        for (size_t i = 0; i < t->n; i++)
            if (t->to[i]->data.symbol.name == e->data.symbol.name) return expr_copy(t->from[i]);
        return expr_copy((Expr*)e);
    }
    if (e->type != EXPR_FUNCTION) return expr_copy((Expr*)e);
    size_t n = e->data.function.arg_count;
    Expr** a = malloc(sizeof(Expr*) * (n ? n : 1));
    for (size_t i = 0; i < n; i++) a[i] = ap_restore(e->data.function.args[i], t);
    Expr* out = expr_new_function(ap_restore(e->data.function.head, t), a, n);
    free(a);
    return out;
}

/* Apart[expr, x] with radicals in a denominator of two or more x-factors:
 * the result, or NULL to take the ordinary path. */
static Expr* apart_with_radical_coefficients(Expr* res) {
    if (res->data.function.arg_count != 2) return NULL;
    Expr* x = res->data.function.args[1];
    if (x->type != EXPR_SYMBOL) return NULL;
    Expr* den = eval_and_free(expr_new_function(expr_new_symbol(SYM_Denominator),
                                                (Expr*[]){ expr_copy(res->data.function.args[0]) }, 1));
    size_t nx = 0;
    bool any_rad = false;
    if (den->type == EXPR_FUNCTION && den->data.function.head->type == EXPR_SYMBOL
        && den->data.function.head->data.symbol.name == SYM_Times) {
        for (size_t i = 0; i < den->data.function.arg_count; i++)
            if (contains_any_symbol_from(den->data.function.args[i], x)) nx++;
    }
    ApRad t = { NULL, NULL, 0, 0 };
    Expr* abs = NULL;
    if (nx >= 2) {
        abs = ap_abstract(res->data.function.args[0], &t);
        any_rad = t.n > 0;
    }
    expr_free(den);
    Expr* out = NULL;
    if (any_rad) {
        Expr* call = expr_new_function(expr_new_symbol(SYM_Apart), (Expr*[]){ abs, expr_copy(x) }, 2);
        abs = NULL;
        Expr* r = eval_and_free(call);
        if (!(r->type == EXPR_FUNCTION && r->data.function.head->type == EXPR_SYMBOL
              && r->data.function.head->data.symbol.name == SYM_Apart))
            out = eval_and_free(ap_restore(r, &t));
        expr_free(r);
    }
    if (abs) expr_free(abs);
    for (size_t i = 0; i < t.n; i++) { expr_free(t.from[i]); expr_free(t.to[i]); }
    free(t.from); free(t.to);
    return out;
}

/* Apart[expr, v] for a compound variable v (x[3], as FeynCalc's Feynman
 * parameters): v is replaced by a fresh symbol, decomposed, and put back, as
 * the reference treats it.  NULL when v is a symbol or the call declines. */
static Expr* apart_compound_variable(Expr* res) {
    if (res->data.function.arg_count != 2) return NULL;
    Expr* v = res->data.function.args[1];
    if (v->type != EXPR_FUNCTION) return NULL;
    if (v->data.function.head->type == EXPR_SYMBOL) {
        const char* h = v->data.function.head->data.symbol.name;
        if (h == SYM_Rule || h == SYM_RuleDelayed || h == SYM_Plus || h == SYM_Times
            || h == SYM_Power || h == SYM_List) return NULL;
    }
    ApRad t = { NULL, NULL, 0, 0 };
    t.cap = 1;
    t.from = malloc(sizeof(Expr*)); t.to = malloc(sizeof(Expr*));
    t.from[0] = expr_copy(v);
    t.to[0] = expr_new_symbol(intern_symbol("Apart`Private`v"));
    t.n = 1;
    /* structural substitution of v */
    Expr* body = res->data.function.args[0];
    Expr* sub = eval_and_free(expr_new_function(expr_new_symbol(SYM_ReplaceAll),
                    (Expr*[]){ expr_copy(body),
                               expr_new_function(expr_new_symbol(SYM_Rule),
                                   (Expr*[]){ expr_copy(v), expr_copy(t.to[0]) }, 2) }, 2));
    Expr* r = eval_and_free(expr_new_function(expr_new_symbol(SYM_Apart),
                                              (Expr*[]){ sub, expr_copy(t.to[0]) }, 2));
    Expr* out = NULL;
    if (!(r->type == EXPR_FUNCTION && r->data.function.head->type == EXPR_SYMBOL
          && r->data.function.head->data.symbol.name == SYM_Apart))
        out = eval_and_free(ap_restore(r, &t));
    expr_free(r);
    expr_free(t.from[0]); expr_free(t.to[0]); free(t.from); free(t.to);
    return out;
}

Expr* builtin_apart(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count < 1) return NULL;
    /* A Rule, relation or logical connective threads, as in the reference
     * (Apart[x -> 1/(a (a + 1)), a] is x -> 1/a - 1/(1 + a)). */
    {
        Expr* t = alg_thread_call(res);
        if (t) return t;
    }

    {
        Expr* ar = apart_compound_variable(res);
        if (ar) return ar;
    }
    /* Apart[expr, x] with expr free of x is expr, as in the reference
     * (Apart[(a + b)/(c d), y] is not split into a/(c d) + b/(c d)). */
    if (res->data.function.arg_count == 2) {
        Expr* x = res->data.function.args[1];
        if (x->type == EXPR_SYMBOL && !contains_any_symbol_from(res->data.function.args[0], x)) {
            Expr* e0 = res->data.function.args[0];
            res->data.function.args[0] = NULL;
            return e0;
        }
    }
    {
        Expr* ar = apart_with_radical_coefficients(res);
        if (ar) return ar;
    }

    /* Modulus -> p: partial fractions over GF(p) (one variable). */
    {
        mpz_t p; mpz_init(p);
        Expr* stripped = NULL;
        int mr = ratfun_modulus_option(res, SYM_Apart, p, &stripped);
        Expr* out = NULL;
        if (mr == 1) {
            const Expr* e0 = res->data.function.args[0];
            const Expr* v = NULL;
            if (res->data.function.arg_count >= 3 || (res->data.function.arg_count == 2 &&
                res->data.function.args[1]->type == EXPR_SYMBOL))
                v = res->data.function.args[1];
            if (v && v->type == EXPR_SYMBOL) out = ratfun_apart_mod(e0, v, p);
            else if (!v) {
                size_t vc = 0, vcap = 16;
                Expr** vs = malloc(sizeof(Expr*) * vcap);
                collect_variables((Expr*)e0, &vs, &vc, &vcap);
                if (vc == 1) out = ratfun_apart_mod(e0, vs[0], p);
                for (size_t i = 0; i < vc; i++) expr_free(vs[i]);
                free(vs);
            }
        }
        mpz_clear(p);
        if (mr != 0) return out;
        if (stripped) {
            out = builtin_apart(stripped);
            expr_free(stripped);
            return out;
        }
    }

    /* Strip a trailing Extension -> α option, if any.  When present,
     * the option is propagated to the Factor[D, ...] call inside apart_impl
     * so the denominator is split over Q(α) before partial-fraction
     * decomposition runs.  Extension -> Automatic triggers
     * extension_autodetect on the input. */
    size_t apart_argc = res->data.function.arg_count;
    bool auto_flag = false;
    const Expr* apart_alpha = extract_extension_option_full(res, &apart_argc, &auto_flag);
    if (apart_argc != 1 && apart_argc != 2) return NULL;

    /* Auto-detect: when Extension -> Automatic was given but no explicit α,
     * scan the input for algebraic generators.  A single generator (n == 1)
     * is passed as the bare α; a multi-generator tower (n >= 2) is passed as
     * the List {α_1, ..., α_n}, which the downstream Together / Factor calls
     * honour (Factor via qa_factor_with_extension_tower, Together via
     * qa_cancel_with_tower). This lets Apart split a denominator like
     * (x - Sqrt[2])(x + Sqrt[3]) over the whole tower Q(Sqrt[2], Sqrt[3]). */
    Expr* apart_alpha_auto = NULL;
    QATower* apart_auto_tower = NULL;
    if (!apart_alpha && auto_flag) {
        apart_auto_tower = extension_autodetect(res->data.function.args[0]);
        if (apart_auto_tower && apart_auto_tower->n == 1) {
            apart_alpha_auto = expr_copy(apart_auto_tower->alpha_renders[0]);
            apart_alpha = apart_alpha_auto;
        } else if (apart_auto_tower && apart_auto_tower->n >= 2) {
            Expr** gens = malloc(sizeof(Expr*) * apart_auto_tower->n);
            for (int i = 0; i < apart_auto_tower->n; i++)
                gens[i] = expr_copy(apart_auto_tower->alpha_renders[i]);
            apart_alpha_auto = expr_new_function(
                expr_new_symbol(SYM_List), gens, (size_t)apart_auto_tower->n);
            free(gens);
            apart_alpha = apart_alpha_auto;
        }
    }

    Expr* result = apart_impl(res, apart_argc, apart_alpha);

    if (apart_alpha_auto) expr_free(apart_alpha_auto);
    if (apart_auto_tower) qa_tower_free(apart_auto_tower);
    return result;
}

static Expr* apart_impl(Expr* res, size_t apart_argc, const Expr* apart_alpha) {

    /* Apart relies on PolynomialQuotient/Remainder and partial-fraction
     * decomposition over rationals — inexact coefficients break both. */
    if (internal_args_contain_inexact_coefficient(res)) {
        return internal_rationalize_then_numericalize(res, builtin_apart);
    }

    Expr* expr = res->data.function.args[0];

    /* Build a fresh Apart[<inner>, <var?>, Extension -> α?] call.  Used
     * to recursively dispatch under threading and the radical-gen
     * substitution path.  `inner_var` may be NULL (1-arg form). */
    /* Forward declaration not needed; we inline construct via lambdas. */

    // Check if it's an equation, inequality, list, or logical function to thread over
    if (expr->type == EXPR_FUNCTION) {
        const char* head = expr->data.function.head->type == EXPR_SYMBOL ? expr->data.function.head->data.symbol.name : "";
        if (head == SYM_List || head == SYM_Equal || head == SYM_Unequal ||
            head == SYM_Less || head == SYM_LessEqual || head == SYM_Greater ||
            head == SYM_GreaterEqual || head == SYM_And || head == SYM_Or ||
            head == SYM_Not) {
            Expr** args = malloc(sizeof(Expr*) * expr->data.function.arg_count);
            for (size_t i = 0; i < expr->data.function.arg_count; i++) {
                /* Build Apart[args[i], var?, Extension -> α?]. */
                size_t inner_argc = apart_argc + (apart_alpha ? 1 : 0);
                Expr** ap_args = malloc(sizeof(Expr*) * inner_argc);
                ap_args[0] = expr_copy(expr->data.function.args[i]);
                if (apart_argc == 2) ap_args[1] = expr_copy(res->data.function.args[1]);
                if (apart_alpha) {
                    ap_args[apart_argc] = expr_new_function(
                        expr_new_symbol(SYM_Rule),
                        (Expr*[]){
                            expr_new_symbol(SYM_Extension),
                            expr_copy((Expr*)apart_alpha)
                        }, 2);
                }
                args[i] = eval_and_free(expr_new_function(
                    expr_new_symbol(SYM_Apart), ap_args, inner_argc));
                free(ap_args);
            }
            Expr* ret = eval_and_free(expr_new_function(expr_copy(expr->data.function.head), args, expr->data.function.arg_count));
            free(args);
            return ret;
        }
    }

    /* Factored rational-function engine (src/poly/ratfun.c): one conversion,
     * partial fractions natively over Q or Q(params).  Declines (NULL) on
     * algebraic coefficients, radicals and everything else outside
     * Q(vars, kernels); those keep the path below.  MATHILDA_RATFUN=off
     * bypasses it. */
    if (!apart_alpha) {
        static int en = -1;
        if (en < 0) { const char* ev = getenv("MATHILDA_RATFUN"); en = !(ev && strcmp(ev, "off") == 0); }
        Expr* v = NULL;
        if (en && apart_argc == 2) v = expr_copy(res->data.function.args[1]);
        else if (en) {
            /* the variable the classical path picks: the last in canonical order */
            size_t vc = 0, vcap = 16;
            Expr** vs = malloc(sizeof(Expr*) * vcap);
            collect_variables(expr, &vs, &vc, &vcap);
            if (vc > 0) {
                qsort(vs, vc, sizeof(Expr*), compare_expr_ptrs);
                v = expr_copy(vs[vc - 1]);
            }
            for (size_t i = 0; i < vc; i++) expr_free(vs[i]);
            free(vs);
        }
        /* An auto-picked kernel variable (y[3], Sin[z], Log[x]): the engine
         * takes symbol generators only, so map every kernel to a fresh
         * symbol, decompose in the one standing for v, and map back.  The
         * fresh names keep the kernels' canonical order, so v's symbol is
         * still the last variable. */
        if (v && v->type != EXPR_SYMBOL && apart_argc == 1) {
            PolyKernelMap km;
            Expr* ks = poly_kernels_to_symbols(expr, &km);
            Expr* rf = NULL;
            if (ks) {
                for (size_t i = 0; i < km.n; i++) {
                    if (expr_eq(km.kern[i], v)) { rf = ratfun_apart(ks, km.sym[i]); break; }
                }
                if (rf) {
                    Expr* back = eval_and_free(poly_symbols_to_kernels(rf, &km));
                    expr_free(rf);
                    rf = back;
                }
                expr_free(ks);
                poly_kernel_map_free(&km);
            }
            if (rf) { expr_free(v); return rf; }
        }
        if (v) {
            Expr* rf = ratfun_apart(expr, v);
            /* one variable over Q(i) / one number field */
            if (!rf && v->type == EXPR_SYMBOL) rf = flint_field_apart(expr, v);
            expr_free(v);
            if (rf) return rf;
        }
    }

    /* Algebraic-generator pass: if the input has a sub-expression u with
     * fractional rational exponents, substitute u -> g^m so the rational
     * function becomes polynomial in g, run Apart in g recursively, then
     * back-substitute.  This converts Apart[1/(-1+r^(3/7))] into partial
     * fractions in r^(1/7).
     *
     * The 1-arg form runs this unconditionally (the fresh generator g
     * becomes the auto-detected variable).  The explicit 2-arg form
     * Apart[expr, v] runs it too, but only when the radical is a genuine
     * radical of v itself (base == v, atom == 1): then g = v^(1/m) is the
     * right generator and Apart[f(Sqrt[x]), x] decomposes in g exactly like
     * the 1-arg form (issue: the 2-arg form used to return the Together'd
     * fraction because Sqrt[x] is not polynomial in x).  A radical of some
     * *other* symbol — Apart[f(Sqrt[y]), x] — is left as an opaque constant
     * in v, so the substitution is skipped there. */
    {
        Expr* base = NULL;
        Expr* atom = NULL;
        int64_t m = 1;
        if (poly_find_radical_gen(expr, &base, &atom, &m)) {
            bool run_subst = (apart_argc == 1);
            if (apart_argc == 2) {
                const Expr* the_var = res->data.function.args[1];
                run_subst = (atom->type == EXPR_INTEGER && atom->data.integer == 1)
                            && expr_eq(base, (Expr*)the_var);
            }
            if (run_subst) {
                char* gen = poly_make_fresh_gen(expr);
                Expr* substituted = poly_subst_radical_to_gen(expr, base, atom, m, gen);
                /* Recursive Apart call, propagating the extension if any.
                 * The 2-arg form pins the decomposition variable to the fresh
                 * generator g, so it is not swayed by whatever other symbol
                 * collect_variables would otherwise auto-pick. */
                size_t inner_argc = (apart_argc == 2 ? 2 : 1) + (apart_alpha ? 1 : 0);
                Expr** ap_args = malloc(sizeof(Expr*) * inner_argc);
                ap_args[0] = substituted;
                if (apart_argc == 2) ap_args[1] = expr_new_symbol(gen);
                if (apart_alpha) {
                    ap_args[inner_argc - 1] = expr_new_function(
                        expr_new_symbol(SYM_Rule),
                        (Expr*[]){
                            expr_new_symbol(SYM_Extension),
                            expr_copy((Expr*)apart_alpha)
                        }, 2);
                }
                Expr* call = expr_new_function(expr_new_symbol(SYM_Apart),
                                  ap_args, inner_argc);
                free(ap_args);
                Expr* result_in_g = evaluate(call);
                expr_free(call);
                Expr* final = poly_subst_radical_from_gen(result_in_g, base, atom, m, gen);
                expr_free(result_in_g);
                expr_free(base);
                expr_free(atom);
                free(gen);
                return final;
            }
            expr_free(base);
            expr_free(atom);
        }
    }

    /* Pre-combine into a single fraction (with Extension if supplied,
     * so any algebraic-number cancellations fire before partial-fraction
     * decomposition starts). */
    Expr* together;
    if (apart_alpha) {
        together = eval_and_free(expr_new_function(
            expr_new_symbol(SYM_Together),
            (Expr*[]){
                expr_copy(expr),
                expr_new_function(
                    expr_new_symbol(SYM_Rule),
                    (Expr*[]){
                        expr_new_symbol(SYM_Extension),
                        expr_copy((Expr*)apart_alpha)
                    }, 2)
            }, 2));
    } else {
        together = eval_and_free(expr_new_function(
            expr_new_symbol(SYM_Together), (Expr*[]){expr_copy(expr)}, 1));
    }

    Expr* var = NULL;
    if (apart_argc == 2) {
        var = expr_copy(res->data.function.args[1]);
    } else {
        size_t v_count = 0, v_cap = 16;
        Expr** vars = malloc(sizeof(Expr*) * v_cap);
        collect_variables(together, &vars, &v_count, &v_cap);
        if (v_count > 0) {
            qsort(vars, v_count, sizeof(Expr*), compare_expr_ptrs);
            var = expr_copy(vars[v_count - 1]);
        }
        for (size_t i = 0; i < v_count; i++) expr_free(vars[i]);
        free(vars);
    }
    
    if (!var) {
        return eval_and_free(expr_new_function(expr_new_symbol(SYM_Expand), (Expr*[]){together}, 1));
    }
    
    Expr* N = eval_and_free(expr_new_function(expr_new_symbol(SYM_Numerator), (Expr*[]){expr_copy(together)}, 1));
    Expr* D = eval_and_free(expr_new_function(expr_new_symbol(SYM_Denominator), (Expr*[]){expr_copy(together)}, 1));

    /* Apart's matrix-construction step assumes both N and D are polynomials
     * in `var`. When Together produces a numerator or denominator with
     * non-integer powers of var (e.g. Together[1/(y^(2/3)-1/y^(1/3))] =
     * y^(1/3)/(y-1)), get_coeff(y^(1/3), y, 0) returns 0 for every row,
     * yielding a zero solution and a wrong answer of 0. Detect this case
     * up front and return the Together'd form unchanged -- partial-fraction
     * decomposition is undefined when the numerator is not polynomial in
     * the chosen variable. */
    Expr* npq_args[2] = { expr_copy(N), expr_copy(var) };
    Expr* npq = eval_and_free(expr_new_function(expr_new_symbol(SYM_PolynomialQ), npq_args, 2));
    Expr* dpq_args[2] = { expr_copy(D), expr_copy(var) };
    Expr* dpq = eval_and_free(expr_new_function(expr_new_symbol(SYM_PolynomialQ), dpq_args, 2));
    bool n_poly = (npq->type == EXPR_SYMBOL && npq->data.symbol.name == SYM_True);
    bool d_poly = (dpq->type == EXPR_SYMBOL && dpq->data.symbol.name == SYM_True);
    expr_free(npq); expr_free(dpq);
    if (!n_poly || !d_poly) {
        expr_free(N); expr_free(D); expr_free(var);
        return together;
    }

    expr_free(together);

    if (get_degree_poly(D, var) == 0) {
        Expr* expanded = eval_and_free(expr_new_function(expr_new_symbol(SYM_Expand), (Expr*[]){eval_and_free(expr_new_function(expr_new_symbol(SYM_Times), (Expr*[]){N, eval_and_free(expr_new_function(expr_new_symbol(SYM_Power), (Expr*[]){D, expr_new_integer(-1)}, 2))}, 2))}, 1));
        expr_free(var);
        return expanded;
    }
    
    Expr* Q = eval_and_free(expr_new_function(expr_new_symbol(SYM_PolynomialQuotient), (Expr*[]){expr_copy(N), expr_copy(D), expr_copy(var)}, 3));
    Expr* R = eval_and_free(expr_new_function(expr_new_symbol(SYM_PolynomialRemainder), (Expr*[]){expr_copy(N), expr_copy(D), expr_copy(var)}, 3));

    /* Factor the denominator over Q (default) or Q(α) when an Extension
     * option was supplied to Apart.  Splitting D over the extension lets
     * partial-fraction decomposition produce e.g. linear factors in
     * Sqrt[2] for `Apart[1/(x^2 - 2), x, Extension -> Sqrt[2]]`. */
    Expr* f_den;
    if (apart_alpha) {
        f_den = eval_and_free(expr_new_function(
            expr_new_symbol(SYM_Factor),
            (Expr*[]){
                expr_copy(D),
                expr_new_function(
                    expr_new_symbol(SYM_Rule),
                    (Expr*[]){
                        expr_new_symbol(SYM_Extension),
                        expr_copy((Expr*)apart_alpha)
                    }, 2)
            }, 2));
    } else {
        f_den = eval_and_free(expr_new_function(
            expr_new_symbol(SYM_Factor), (Expr*[]){expr_copy(D)}, 1));
    }
    expr_free(N); expr_free(D);

    size_t num_args = 1;
    Expr** args = &f_den;
    if (f_den->type == EXPR_FUNCTION && f_den->data.function.head->data.symbol.name == SYM_Times) {
        num_args = f_den->data.function.arg_count;
        args = f_den->data.function.args;
    }
    
    Expr* C = expr_new_integer(1);
    size_t m = 0;
    Expr** bases = malloc(sizeof(Expr*) * num_args);
    int64_t* ks = malloc(sizeof(int64_t) * num_args);
    
    for (size_t i = 0; i < num_args; i++) {
        Expr* factor = args[i];
        Expr* base = factor;
        int64_t k = 1;
        if (factor->type == EXPR_FUNCTION && factor->data.function.head->data.symbol.name == SYM_Power) {
            base = factor->data.function.args[0];
            if (factor->data.function.args[1]->type == EXPR_INTEGER) {
                k = factor->data.function.args[1]->data.integer;
            }
        }
        
        if (get_degree_poly(base, var) == 0) {
            Expr* term = (factor == base && k != 1) ? eval_and_free(expr_new_function(expr_new_symbol(SYM_Power), (Expr*[]){expr_copy(base), expr_new_integer(k)}, 2)) : expr_copy(factor);
            C = eval_and_free(expr_new_function(expr_new_symbol(SYM_Times), (Expr*[]){C, term}, 2));
        } else {
            bases[m] = expr_copy(base);
            ks[m] = k;
            m++;
        }
    }
    
    if (m == 0) {
        Expr* C_inv = eval_and_free(expr_new_function(expr_new_symbol(SYM_Power), (Expr*[]){C, expr_new_integer(-1)}, 2));
        Expr* R_term = eval_and_free(expr_new_function(expr_new_symbol(SYM_Expand), (Expr*[]){eval_and_free(expr_new_function(expr_new_symbol(SYM_Times), (Expr*[]){R, C_inv}, 2))}, 1));
        Expr* res = eval_and_free(expr_new_function(expr_new_symbol(SYM_Plus), (Expr*[]){Q, R_term}, 2));
        expr_free(f_den); free(bases); free(ks); expr_free(var);
        return res;
    }
    
    /* Fast path — Bézout modular-inverse partial fraction for a squarefree
     * denominator over a radical field.  The classical solve below builds an
     * S×(S+1) matrix and RowReduces it; over Sqrt[k]-laden coefficients that
     * symbolic Gaussian elimination blows up (the parametric rational-
     * integrator's Apart blocker), and even after the Sqrt[k]->S collapse a
     * multivariate-rational RowReduce is exponential.  Instead, for distinct
     * irreducible factors p_i (multiplicity 1, ANY degree — so it also covers
     * the irreducible quadratics that yield ArcTan), the partial fraction is
     *     R/D = sum_i A_i / p_i,   A_i = R * (D/p_i)^{-1} mod p_i,
     * where (D/p_i)^{-1} mod p_i comes from PolynomialExtendedGCD(D/p_i, p_i)
     * = {g, {u, v}} with u*(D/p_i) + v*p_i = g (a unit, since the factors are
     * coprime): A_i = PolynomialRemainder[R*u, p_i, x] / g.  Every step is a
     * low-degree operation over the field, each a fast Sqrt[k] Cancel via the
     * FLINT parametric path.  Gated on a symbolic radical so plain-rational
     * Apart keeps its existing (RowReduce) output. */
    {
        bool squarefree = (m > 0);
        bool bases_have_radical = false;
        for (size_t i = 0; i < m; i++) {
            if (ks[i] != 1) { squarefree = false; break; }
            if (apart_has_symbolic_radical(bases[i])) bases_have_radical = true;
        }
        if (squarefree && (bases_have_radical || apart_has_symbolic_radical(R))) {
            Expr* sum = Q;   /* take ownership of the polynomial part */
            for (size_t i = 0; i < m; i++) {
                /* q_i = D / p_i = C * prod_{j != i} p_j. */
                Expr* qi = expr_copy(C);
                for (size_t j = 0; j < m; j++) {
                    if (j == i) continue;
                    qi = eval_and_free(expr_new_function(expr_new_symbol(SYM_Times),
                        (Expr*[]){qi, expr_copy(bases[j])}, 2));
                }
                /* {g, {u, v}} = PolynomialExtendedGCD[q_i, p_i, x]. */
                Expr* egcd = eval_and_free(expr_new_function(
                    expr_new_symbol(SYM_PolynomialExtendedGCD),
                    (Expr*[]){ qi, expr_copy(bases[i]), expr_copy(var) }, 3));
                Expr* g_unit = NULL; Expr* u = NULL;
                if (egcd->type == EXPR_FUNCTION && egcd->data.function.arg_count == 2
                    && egcd->data.function.args[1]->type == EXPR_FUNCTION
                    && egcd->data.function.args[1]->data.function.arg_count >= 1) {
                    g_unit = expr_copy(egcd->data.function.args[0]);
                    u = expr_copy(egcd->data.function.args[1]->data.function.args[0]);
                }
                expr_free(egcd);
                if (!g_unit || !u) {   /* unexpected shape: abandon fast path */
                    if (g_unit) expr_free(g_unit);
                    if (u) expr_free(u);
                    expr_free(sum); sum = NULL; break;
                }
                /* A_i = PolynomialRemainder[R*u, p_i, x] / g. */
                Expr* Ru = eval_and_free(expr_new_function(expr_new_symbol(SYM_Times),
                    (Expr*[]){expr_copy(R), u}, 2));
                Expr* rem = eval_and_free(expr_new_function(
                    expr_new_symbol(SYM_PolynomialRemainder),
                    (Expr*[]){ Ru, expr_copy(bases[i]), expr_copy(var) }, 3));
                Expr* g_inv = eval_and_free(expr_new_function(expr_new_symbol(SYM_Power),
                    (Expr*[]){g_unit, expr_new_integer(-1)}, 2));
                Expr* Ai = eval_and_free(expr_new_function(expr_new_symbol(SYM_Cancel),
                    (Expr*[]){ eval_and_free(expr_new_function(expr_new_symbol(SYM_Times),
                        (Expr*[]){rem, g_inv}, 2)) }, 1));
                /* term = A_i / p_i. */
                Expr* pi_inv = eval_and_free(expr_new_function(expr_new_symbol(SYM_Power),
                    (Expr*[]){expr_copy(bases[i]), expr_new_integer(-1)}, 2));
                Expr* term = eval_and_free(expr_new_function(expr_new_symbol(SYM_Times),
                    (Expr*[]){Ai, pi_inv}, 2));
                sum = eval_and_free(expr_new_function(expr_new_symbol(SYM_Plus),
                    (Expr*[]){sum, term}, 2));
            }
            if (sum) {
                expr_free(R); expr_free(C); expr_free(f_den);
                for (size_t i = 0; i < m; i++) expr_free(bases[i]);
                free(bases); free(ks); expr_free(var);
                return sum;
            }
            /* fast path bailed (Q consumed into a freed sum): rebuild Q=0 so the
             * classical path below still has a valid polynomial part. */
            Q = expr_new_integer(0);
        }
    }

    /* FLINT-native partial fraction for a plain rational function over Q:
     * distinct-factor CRT split + p-adic expansion, entirely in fmpq_poly.
     * Replaces the O(S^2+) symbolic RowReduce below for the pure-Q case; the
     * multivariate / radical case (any base or coefficient carrying another
     * symbol) declines here and keeps the classical elimination path. */
    if (var->type == EXPR_SYMBOL) {
        Expr* frac = flint_apart_over_q(R, (const Expr* const*)bases, ks,
                                        (int)m, C, var->data.symbol.name);
        if (frac) {
            Expr* res_pf = eval_and_free(expr_new_function(
                expr_new_symbol(SYM_Plus), (Expr*[]){ Q, frac }, 2));
            expr_free(R); expr_free(C); expr_free(f_den);
            for (size_t i = 0; i < m; i++) expr_free(bases[i]);
            free(bases); free(ks); expr_free(var);
            return res_pf;
        }
    }

    Expr* D_main = expr_new_integer(1);
    for (size_t i = 0; i < m; i++) {
        Expr* term = ks[i] == 1 ? expr_copy(bases[i]) : eval_and_free(expr_new_function(expr_new_symbol(SYM_Power), (Expr*[]){expr_copy(bases[i]), expr_new_integer(ks[i])}, 2));
        D_main = eval_and_free(expr_new_function(expr_new_symbol(SYM_Times), (Expr*[]){D_main, term}, 2));
    }
    
    Expr* exp_D_main = expr_expand(D_main);
    int S = get_degree_poly(exp_D_main, var);
    expr_free(exp_D_main);
    expr_free(D_main);
    
    Expr*** M = malloc(sizeof(Expr**) * S);
    for (int row = 0; row < S; row++) {
        M[row] = malloc(sizeof(Expr*) * (S + 1));
        for (int col = 0; col <= S; col++) M[row][col] = expr_new_integer(0);
    }
    
    int col = 0;
    for (size_t i = 0; i < m; i++) {
        Expr* exp_base = expr_expand(bases[i]);
        int d_B = get_degree_poly(exp_base, var);
        expr_free(exp_base);
        
        for (int j = 1; j <= ks[i]; j++) {
            for (int r = 0; r < d_B; r++) {
                Expr* M_ij = expr_new_integer(1);
                for (size_t p = 0; p < m; p++) {
                    int exp = (p == i) ? ks[p] - j : ks[p];
                    if (exp > 0) {
                        Expr* term = (exp == 1) ? expr_copy(bases[p]) : eval_and_free(expr_new_function(expr_new_symbol(SYM_Power), (Expr*[]){expr_copy(bases[p]), expr_new_integer(exp)}, 2));
                        M_ij = eval_and_free(expr_new_function(expr_new_symbol(SYM_Times), (Expr*[]){M_ij, term}, 2));
                    }
                }
                
                Expr* var_pow = (r == 0) ? expr_new_integer(1) : ((r == 1) ? expr_copy(var) : eval_and_free(expr_new_function(expr_new_symbol(SYM_Power), (Expr*[]){expr_copy(var), expr_new_integer(r)}, 2)));
                Expr* term_mult = eval_and_free(expr_new_function(expr_new_symbol(SYM_Times), (Expr*[]){M_ij, var_pow}, 2));
                Expr* expanded_term = expr_expand(term_mult);
                expr_free(term_mult);
                
                for (int row = 0; row < S; row++) {
                    expr_free(M[row][col]);
                    M[row][col] = get_coeff(expanded_term, var, row);
                }
                expr_free(expanded_term);
                col++;
            }
        }
    }
    
    Expr* C_inv = eval_and_free(expr_new_function(expr_new_symbol(SYM_Power), (Expr*[]){C, expr_new_integer(-1)}, 2));
    Expr* tmp_R = eval_and_free(expr_new_function(expr_new_symbol(SYM_Times), (Expr*[]){R, C_inv}, 2));
    Expr* R_rhs = expr_expand(tmp_R);
    expr_free(tmp_R);
    for (int row = 0; row < S; row++) {
        expr_free(M[row][S]);
        M[row][S] = get_coeff(R_rhs, var, row);
    }
    expr_free(R_rhs);
    
    Expr** row_args = malloc(sizeof(Expr*) * S);
    for (int row = 0; row < S; row++) {
        row_args[row] = expr_new_function(expr_new_symbol(SYM_List), M[row], S + 1);
        free(M[row]);
    }
    free(M);
    
    Expr* matrix_expr = expr_new_function(expr_new_symbol(SYM_List), row_args, S);
    free(row_args);
    
    /* pack_eval_plain, not evaluate: RowReduce now answers with a BUFFER for a
     * uniform machine matrix, and the loop below walks the result's args. See
     * pack.h -- the transparency gate protects a builtin's arguments, never the
     * value it gets back from its own internal evaluate(). */
    Expr* rr_call = expr_new_function(expr_new_symbol(SYM_RowReduce), (Expr*[]){matrix_expr}, 1);
    Expr* reduced = pack_eval_plain(rr_call);
    expr_free(rr_call);
    
    Expr* pfrac_sum = expr_new_integer(0);
    if (reduced->type == EXPR_FUNCTION && reduced->data.function.head->data.symbol.name == SYM_List && reduced->data.function.arg_count == (size_t)S) {
        col = 0;
        for (size_t i = 0; i < m; i++) {
            Expr* exp_base = expr_expand(bases[i]);
            int d_B = get_degree_poly(exp_base, var);
            expr_free(exp_base);
            
            for (int j = 1; j <= ks[i]; j++) {
                Expr* A_ij = expr_new_integer(0);
                for (int r = 0; r < d_B; r++) {
                    Expr* coeff = expr_copy(reduced->data.function.args[col]->data.function.args[S]);
                    
                    Expr* fact_coeff = eval_and_free(expr_new_function(expr_new_symbol(SYM_Factor), (Expr*[]){coeff}, 1));
                    
                    Expr* var_pow = (r == 0) ? expr_new_integer(1) : ((r == 1) ? expr_copy(var) : eval_and_free(expr_new_function(expr_new_symbol(SYM_Power), (Expr*[]){expr_copy(var), expr_new_integer(r)}, 2)));
                    Expr* term = eval_and_free(expr_new_function(expr_new_symbol(SYM_Times), (Expr*[]){fact_coeff, var_pow}, 2));
                    
                    A_ij = eval_and_free(expr_new_function(expr_new_symbol(SYM_Plus), (Expr*[]){A_ij, term}, 2));
                    col++;
                }
                
                Expr* denom = (j == 1) ? expr_copy(bases[i]) : eval_and_free(expr_new_function(expr_new_symbol(SYM_Power), (Expr*[]){expr_copy(bases[i]), expr_new_integer(j)}, 2));
                Expr* fraction = eval_and_free(expr_new_function(expr_new_symbol(SYM_Times), (Expr*[]){A_ij, eval_and_free(expr_new_function(expr_new_symbol(SYM_Power), (Expr*[]){denom, expr_new_integer(-1)}, 2))}, 2));
                
                pfrac_sum = eval_and_free(expr_new_function(expr_new_symbol(SYM_Plus), (Expr*[]){pfrac_sum, fraction}, 2));
            }
        }
    } else {
        pfrac_sum = expr_new_integer(0); 
    }
    expr_free(reduced);
    
    Expr* final_res = eval_and_free(expr_new_function(expr_new_symbol(SYM_Plus), (Expr*[]){Q, pfrac_sum}, 2));
    
    for (size_t i = 0; i < m; i++) expr_free(bases[i]);
    free(bases); free(ks);
    expr_free(f_den); expr_free(var);
    
    return final_res;
}

void parfrac_init(void) {
    symtab_add_builtin("Apart", builtin_apart);
    symtab_get_def("Apart")->attributes |= ATTR_PROTECTED | ATTR_LISTABLE;
}
