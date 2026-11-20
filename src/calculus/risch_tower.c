/* risch_tower.c — differential transcendental-tower construction for the
 * recursive Risch integrator (Bronstein/Roach lineage).
 *
 * This is the operational Risch STRUCTURE-THEOREM layer, split out of
 * integrate_risch_transcendental.c: it collects the log / exp / tangent kernels
 * of an integrand, decides their commensurability and (triangular)
 * independence, and assembles the ordered differential tower `RtTower`
 * (rt_tower_build_min) together with the tower's structural operations —
 * lifecycle (rt_tower_free), monomial derivation (rt_tower_deriv / rt_dt_i /
 * rt_build_deriv_rules), structural kernel substitution (rt_subst_kernels), and
 * the pre-normalisations that expose the minimal generator set
 * (rt_powers_to_exp, rt_expand_logs, rt_expand_exp_sums).
 *
 * The abstract structure-theorem DECISION builtins (Risch`RationalSpan /
 * LogReducible / ExpReducible, Bronstein §9.3) live separately in
 * risch_structure.c; this file is the concrete tower the recursive integrator
 * and its Risch-DE layer run over.
 *
 * Layering: the field integrator (rt_field_integrate, in the transcendental
 * file) and the Risch DE (integrate_risch_rde.c) sit ON TOP of this tower layer.
 * The one upward call here — rt_tower_build_min's structure-soundness probe via
 * rt_field_integrate — uses the shared entry point already declared in
 * integrate_risch_rde.h, so there is no dependency cycle.  A handful of small
 * eval/predicate helpers remain DEFINED in integrate_risch_transcendental.c and
 * are shared back in via risch_tower.h.
 *
 * Memory contract follows the standard Mathilda ownership rule (see the
 * per-function comments): every builder returns freshly-owned trees.
 */

#include "risch_tower.h"
#include "risch_structure.h"
#include "integrate_risch_rde.h"

#include "expr.h"
#include "eval.h"
#include "parse.h"
#include "print.h"
#include "symtab.h"
#include "attr.h"
#include "sym_intern.h"
#include "sym_names.h"
#include "arithmetic.h"

#include <stdbool.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>

/* Structural (non-evaluating) rewrite exposing a transcendental general power as
 * an explicit base-e exponential:  b^e  ->  Power[E, e Log[b]],  for every
 * x-dependent power whose exponent `e` is NOT a rational-number constant (so the
 * genuine algebraic / polynomial powers x^2, x^(1/2), x^(-1) are left intact).
 *
 * The evaluator canonicalizes E^(c Log[b]) straight back to b^c (power.c
 * simplify_exp_log), so an exponential whose exponent carries a Log — precisely
 * the image of a circular/inverse trig of a logarithm under TrigToExp, e.g.
 * Cos[x Log x] -> (x^(I x) + x^(-I x))/2, whose kernels are E^(±I x Log x) —
 * is STORED as the general power b^e and hides from every exponential recognizer
 * (all keyed to Exp / Power[E, .]).  This rebuilds the raw Power[E, .] spelling
 * WITHOUT passing through the evaluator, so the tower machinery
 * (rt_collect_exp_exponents / rt_subst_kernels, which match Power[E, .]
 * structurally) sees the hyperexponential monomial.  Here Log[b] is a GENUINE
 * logarithmic sub-kernel (D[Log x] = 1/x), collected and differentiated as such —
 * unlike the constant-base debasing (rt_debase_exponentials), whose opaque L_p
 * stands for a CONSTANT Log[p].  The result is a raw tree that MUST be
 * kernel-substituted (rt_subst_kernels) before any evaluation; back-substituting
 * the tower kernels then lets the evaluator re-collapse E^(e Log b) -> b^e, so the
 * antiderivative is rendered in the original power form.  Applying it to an
 * already base-e form is a no-op (Power[E, .] is skipped). */
Expr* rt_powers_to_exp(Expr* e, Expr* x) {
    if (!e) return NULL;
    if (e->type != EXPR_FUNCTION) return expr_copy(e);
    if (rt_head_is(e, "Power") && e->data.function.arg_count == 2) {
        Expr* b  = e->data.function.args[0];
        Expr* ex = e->data.function.args[1];
        bool b_is_E = (b->type == EXPR_SYMBOL && b->data.symbol.name == intern_symbol("E"));
        /* Transcendental exponential kernel: non-e base, non-rational exponent,
         * and the whole power depends on x (so a constant power like 2^Sqrt[2] is
         * left alone — nothing to integrate there). */
        if (!b_is_E && !rt_is_rat_const(ex) && !rt_free_of_x(e, x)) {
            Expr* nb  = rt_powers_to_exp(b, x);
            /* Split the exponent's additive terms: rational-constant terms stay as
             * an ALGEBRAIC power b^(alg) (a base-field factor), the rest form the
             * transcendental kernel E^(trans Log b).  This is essential because the
             * evaluator MERGES a polynomial coefficient into the exponent — x^2 x^x
             * -> x^(2+x), x^2 x^(I x) -> x^(2+I x) — which would otherwise read as a
             * kernel independent of x^(3+x)/x^(3+I x) (exponent ratio non-constant),
             * defeating the commensurability reduction.  Peeling the constant part
             * back off restores the shared primitive kernel x^x / x^(I x). */
            Expr** terms; size_t nt; Expr* single[1];
            if (rt_head_is(ex, "Plus")) {
                terms = ex->data.function.args; nt = ex->data.function.arg_count;
            } else { single[0] = ex; terms = single; nt = 1; }
            Expr** algt = malloc(nt * sizeof(Expr*)); size_t na = 0;
            Expr** trt  = malloc(nt * sizeof(Expr*)); size_t ntr = 0;
            for (size_t i = 0; i < nt; i++) {
                if (rt_is_rat_const(terms[i])) algt[na++] = expr_copy(terms[i]);
                else trt[ntr++] = rt_powers_to_exp(terms[i], x);
            }
            /* trt is never empty: ex is non-rational overall, so it has a term that
             * is not a rational constant. */
            Expr* trans = (ntr == 1) ? trt[0]
                : expr_new_function(expr_new_symbol("Plus"), trt, ntr);
            free(trt);
            Expr* ker = expr_new_function(expr_new_symbol("Power"),
                (Expr*[]){ expr_new_symbol("E"),
                    expr_new_function(expr_new_symbol("Times"),
                        (Expr*[]){ trans, expr_new_function(expr_new_symbol("Log"),
                            (Expr*[]){ expr_copy(nb) }, 1) }, 2) }, 2);
            Expr* result;
            if (na == 0) {
                expr_free(nb);
                result = ker;
            } else {
                Expr* algsum = (na == 1) ? algt[0]
                    : expr_new_function(expr_new_symbol("Plus"), algt, na);
                Expr* algpow = expr_new_function(expr_new_symbol("Power"),
                    (Expr*[]){ nb, algsum }, 2);
                result = expr_new_function(expr_new_symbol("Times"),
                    (Expr*[]){ algpow, ker }, 2);
            }
            free(algt);
            return result;
        }
    }
    Expr* nh = rt_powers_to_exp(e->data.function.head, x);
    size_t k = e->data.function.arg_count;
    Expr** na = malloc((k ? k : 1) * sizeof(Expr*));
    for (size_t i = 0; i < k; i++)
        na[i] = rt_powers_to_exp(e->data.function.args[i], x);
    Expr* r = expr_new_function(nh, na, k);
    free(na);
    return r;
}
/* Collect (as owned copies, deduplicated) the exponents of every E^w / Exp[w]
 * kernel of `e` whose exponent depends on x. */
static bool rt_gauss_const_exponent(Expr* w, int64_t* a, int64_t* b);
void rt_collect_exp_exponents(Expr* e, Expr* x,
                                     Expr*** arr, size_t* n, size_t* cap) {
    if (!e || e->type != EXPR_FUNCTION) return;
    const char* h = (e->data.function.head->type == EXPR_SYMBOL)
        ? e->data.function.head->data.symbol.name : NULL;
    Expr* w = NULL;
    if (h == intern_symbol("Exp") && e->data.function.arg_count == 1
        && !rt_free_of_x(e->data.function.args[0], x))
        w = e->data.function.args[0];
    else if (h == intern_symbol("Power") && e->data.function.arg_count == 2
        && e->data.function.args[0]->type == EXPR_SYMBOL
        && e->data.function.args[0]->data.symbol.name == intern_symbol("E")
        && !rt_free_of_x(e->data.function.args[1], x))
        w = e->data.function.args[1];
    /* a constant exponent spelled with x (I u - I (1 + u)) is no monomial:
     * rt_subst_kernels maps E^(a + b i) to E^a Z^b */
    if (w && rt_head_is(w, "Plus")) {
        int64_t ga, gb;
        if (rt_gauss_const_exponent(w, &ga, &gb)) w = NULL;
    }
    if (w) {
        bool dup = false;
        for (size_t i = 0; i < *n; i++) if (expr_eq((*arr)[i], w)) { dup = true; break; }
        if (!dup) {
            if (*n == *cap) {
                *cap = *cap ? *cap * 2 : 4;
                *arr = realloc(*arr, *cap * sizeof(Expr*));
            }
            (*arr)[(*n)++] = expr_copy(w);
        }
    }
    rt_collect_exp_exponents(e->data.function.head, x, arr, n, cap);
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        rt_collect_exp_exponents(e->data.function.args[i], x, arr, n, cap);
}
/* Collect (owned, deduplicated) every Log[u] kernel of `e` whose argument
 * depends on x. */
void rt_collect_logs(Expr* e, Expr* x, Expr*** arr, size_t* n, size_t* cap) {
    if (!e || e->type != EXPR_FUNCTION) return;
    if (e->data.function.head->type == EXPR_SYMBOL
        && e->data.function.head->data.symbol.name == intern_symbol("Log")
        && e->data.function.arg_count == 1
        && !rt_free_of_x(e->data.function.args[0], x)) {
        bool dup = false;
        for (size_t i = 0; i < *n; i++) if (expr_eq((*arr)[i], e)) { dup = true; break; }
        if (!dup) {
            if (*n == *cap) { *cap = *cap ? *cap * 2 : 4;
                              *arr = realloc(*arr, *cap * sizeof(Expr*)); }
            (*arr)[(*n)++] = expr_copy(e);
        }
    }
    rt_collect_logs(e->data.function.head, x, arr, n, cap);
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        rt_collect_logs(e->data.function.args[i], x, arr, n, cap);
}

/* Non-elementary primitive special functions Knowles admits as tower monomials:
 * each is a single-argument primitive theta = SF[u] with theta' = Dcoef a field
 * element of the lower tower (Erf: (2/Sqrt[Pi])E^(-u^2)u'; Erfi: (2/Sqrt[Pi])E^(u^2)u';
 * Erfc: -(2/Sqrt[Pi])E^(-u^2)u'; ExpIntegralEi: E^u u'/u; LogIntegral: u'/Log[u]).
 * All have D-rules in deriv.c.  (SinIntegral/CosIntegral/Fresnel are reducible but
 * their derivatives introduce trig kernels — deferred; see KNOWLES_DESIGN.md §2.1.) */
bool rt_is_primitive_head(const char* h) {
    return h == intern_symbol("Erf")  || h == intern_symbol("Erfi")
        || h == intern_symbol("Erfc") || h == intern_symbol("ExpIntegralEi")
        || h == intern_symbol("LogIntegral");
}

/* Collect (owned, deduplicated) every primitive-SF kernel of `e` whose single
 * argument depends on x. */
void rt_collect_primitives(Expr* e, Expr* x, Expr*** arr, size_t* n, size_t* cap) {
    if (!e || e->type != EXPR_FUNCTION) return;
    if (e->data.function.head->type == EXPR_SYMBOL
        && rt_is_primitive_head(e->data.function.head->data.symbol.name)
        && e->data.function.arg_count == 1
        && !rt_free_of_x(e->data.function.args[0], x)) {
        bool dup = false;
        for (size_t i = 0; i < *n; i++) if (expr_eq((*arr)[i], e)) { dup = true; break; }
        if (!dup) {
            if (*n == *cap) { *cap = *cap ? *cap * 2 : 4;
                              *arr = realloc(*arr, *cap * sizeof(Expr*)); }
            (*arr)[(*n)++] = expr_copy(e);
        }
    }
    rt_collect_primitives(e->data.function.head, x, arr, n, cap);
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        rt_collect_primitives(e->data.function.args[i], x, arr, n, cap);
}

/* Structural containment: does `big` contain `small` as a subexpression? */
bool rt_contains(Expr* big, Expr* small) {
    if (!big) return false;
    if (expr_eq(big, small)) return true;
    if (big->type != EXPR_FUNCTION) return false;
    if (rt_contains(big->data.function.head, small)) return true;
    for (size_t i = 0; i < big->data.function.arg_count; i++)
        if (rt_contains(big->data.function.args[i], small)) return true;
    return false;
}

/* Does an exponential E^w with w a nonzero rational multiple of `a` occur in
 * e?  (An exponential level is stored as its class primitive E^a, its members
 * E^(m a) as aliases t^m, so a kernel that needs the level may spell it E^(2 a).) */
static bool rt_has_exp_multiple(const Expr* e, Expr* a) {
    if (!e || e->type != EXPR_FUNCTION) return false;
    if (e->data.function.head->type == EXPR_SYMBOL
        && e->data.function.head->data.symbol.name == intern_symbol("Power")
        && e->data.function.arg_count == 2
        && e->data.function.args[0]->type == EXPR_SYMBOL
        && e->data.function.args[0]->data.symbol.name == intern_symbol("E")) {
        Expr* q = rt_eval1("Cancel", expr_new_function(expr_new_symbol("Times"),
            (Expr*[]){ expr_copy(e->data.function.args[1]), expr_new_function(expr_new_symbol("Power"),
                (Expr*[]){ expr_copy(a), expr_new_integer(-1) }, 2) }, 2));
        bool mult = q && (q->type == EXPR_INTEGER || rt_head_is(q, "Rational")) && !rt_is_zero(q);
        if (q) expr_free(q);
        if (mult) return true;
    }
    if (rt_has_exp_multiple(e->data.function.head, a)) return true;
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (rt_has_exp_multiple(e->data.function.args[i], a)) return true;
    return false;
}
/* Does kernel K need level j below it: it contains that level's kernel, or
 * (an exponential level) an alias E^(m a) of it? */
static bool rt_kernel_needs(Expr* K, RtTower* T, size_t j) {
    if (rt_contains(K, T->kernel[j])) return true;
    return T->kind[j] == RT_EXP && T->arg[j] && rt_has_exp_multiple(K, T->arg[j]);
}

/* Build the monomial  prod_j lv[j]^e[j]  (a Times of Powers; owned). */
Expr* rt_build_monomial(Expr** lv, const long* e, size_t nlv) {
    Expr** fs = malloc((nlv ? nlv : 1) * sizeof(Expr*));
    size_t nf = 0;
    for (size_t j = 0; j < nlv; j++) {
        if (e[j] == 0) continue;
        fs[nf++] = expr_new_function(expr_new_symbol("Power"),
            (Expr*[]){ expr_copy(lv[j]), expr_new_integer(e[j]) }, 2);
    }
    Expr* m;
    if (nf == 0) m = expr_new_integer(1);
    else if (nf == 1) m = fs[0];
    else m = expr_new_function(expr_new_symbol("Times"), fs, nf);
    free(fs);
    return m;
}

/* Decode the flat monomial index `idx` into the exponent vector `e` over the
 * per-variable degree bounds `bd` (odometer, lowest var fastest). */
void rt_decode_mono(long idx, const long* bd, size_t nlv, long* e) {
    for (size_t j = 0; j < nlv; j++) {
        long base = bd[j] + 1;
        e[j] = idx % base;
        idx /= base;
    }
}

/* ------------------------------------------------------------------ */
/* Logarithm-argument normalization for the tower builder.             */
/* ------------------------------------------------------------------ */
/* Rewrite  Log[a b ...] -> Log[a] + Log[b] + ...  and  Log[b^p] -> p Log[b]
 * recursively, so a nested-log integrand presents the MINIMAL set of
 * multiplicatively-independent Log generators to rt_collect_logs.  Without
 * this, a composite kernel like Log[x/Log[x]] is treated as an INDEPENDENT
 * transcendental on top of Log[x] and Log[Log[x]], inflating the tower depth
 * (and hence the undetermined-coefficient ansatz — the dominant cost) by a
 * spurious, functionally-redundant generator.
 *
 * The rewrite is exact for the DERIVATIVE (d/dx Log[u] = u'/u regardless of
 * branch), which is all the tower's correct-by-construction certificate and
 * the final diff-back gate (rt_verify_antideriv, against the ORIGINAL f)
 * depend on; the branch-cut constants it drops are absorbed by the constant
 * of integration.  Sums inside a Log are left intact (Log[a+b] does not
 * split). */

/* Expanded form of Log[a] (the argument `a` is borrowed). */
Expr* rt_log_of(Expr* a) {
    if (a->type == EXPR_FUNCTION && a->data.function.head->type == EXPR_SYMBOL) {
        const char* h = a->data.function.head->data.symbol.name;
        if (h == intern_symbol("Times")) {
            size_t m = a->data.function.arg_count;
            Expr** ts = malloc((m ? m : 1) * sizeof(Expr*));
            for (size_t i = 0; i < m; i++) ts[i] = rt_log_of(a->data.function.args[i]);
            Expr* s = expr_new_function(expr_new_symbol("Plus"), ts, m);
            free(ts);
            return s;
        }
        if (h == intern_symbol("Power") && a->data.function.arg_count == 2) {
            Expr* lb = rt_log_of(a->data.function.args[0]);
            return expr_new_function(expr_new_symbol("Times"),
                (Expr*[]){ rt_expand_logs(a->data.function.args[1]), lb }, 2);
        }
    }
    /* Atomic (or Plus) argument: keep the Log, but still expand logs nested
     * inside the argument (e.g. the Log[x] inside Log[Log[x]]). */
    return expr_new_function(expr_new_symbol("Log"),
        (Expr*[]){ rt_expand_logs(a) }, 1);
}

Expr* rt_expand_logs(Expr* e) {
    if (!e || e->type != EXPR_FUNCTION) return expr_copy(e);
    if (e->data.function.head->type == EXPR_SYMBOL
        && e->data.function.head->data.symbol.name == intern_symbol("Log")
        && e->data.function.arg_count == 1) {
        return rt_log_of(e->data.function.args[0]);
    }
    size_t n = e->data.function.arg_count;
    Expr** args = malloc((n ? n : 1) * sizeof(Expr*));
    for (size_t i = 0; i < n; i++) args[i] = rt_expand_logs(e->data.function.args[i]);
    Expr* r = expr_new_function(expr_copy(e->data.function.head), args, n);
    free(args);
    return r;
}
/* Leaf count (orders the composite-class elimination, largest first). */
static size_t rt_leaves(const Expr* e) {
    if (!e || e->type != EXPR_FUNCTION) return 1;
    size_t n = rt_leaves(e->data.function.head);
    for (size_t i = 0; i < e->data.function.arg_count; i++) n += rt_leaves(e->data.function.args[i]);
    return n;
}

void rt_tower_free(RtTower* T) {
    for (size_t i = 0; i < T->n; i++) {
        if (T->kernel && T->kernel[i]) expr_free(T->kernel[i]);
        if (T->arg && T->arg[i]) expr_free(T->arg[i]);
        if (T->t && T->t[i]) expr_free(T->t[i]);
        if (T->Dcoef && T->Dcoef[i]) expr_free(T->Dcoef[i]);
    }
    if (T->subrules) expr_free(T->subrules);
    for (size_t i = 0; i < T->nm; i++) {
        if (T->marg && T->marg[i]) expr_free(T->marg[i]);
        if (T->mval && T->mval[i]) expr_free(T->mval[i]);
    }
    free(T->mval); T->mval = NULL;
    for (size_t i = 0; i < T->nlm; i++) { expr_free(T->lmk[i]); expr_free(T->lmv[i]); }
    free(T->lmk); free(T->lmv); T->lmk = NULL; T->lmv = NULL; T->nlm = 0;
    free(T->kind); free(T->kernel); free(T->arg); free(T->t); free(T->Dcoef);
    free(T->tsg); free(T->marg); free(T->mprim); free(T->mmult);
    T->kind = NULL; T->kernel = NULL; T->arg = NULL; T->t = NULL; T->Dcoef = NULL;
    T->tsg = NULL; T->marg = NULL; T->mprim = NULL; T->mmult = NULL;
    T->n = 0; T->nm = 0; T->subrules = NULL;
}
/* Collect distinct x-dependent tangent-family monomials: each is (argument u, sigma)
 * with sigma = +1 for Tan/Cot (special polynomial t^2+1) and -1 for Tanh/Coth
 * (t^2-1).  Cot[u]/Coth[u] are reciprocals of the Tan[u]/Tanh[u] tower variable, so
 * they yield the same monomial (deduped by (u, sigma)).  The tower variable is
 * always the Tan/Tanh spelling; rt_tower_build adds subst rules for all four. */
void rt_collect_tangents(Expr* e, Expr* x, Expr*** args, long** sigs,
                                size_t* n, size_t* cap) {
    if (!e || e->type != EXPR_FUNCTION) return;
    const char* h = (e->data.function.head->type == EXPR_SYMBOL)
        ? e->data.function.head->data.symbol.name : NULL;
    if (h && e->data.function.arg_count == 1
        && !rt_free_of_x(e->data.function.args[0], x)) {
        long sg = 0;
        if (h == intern_symbol("Tan") || h == intern_symbol("Cot")) sg = 1;
        else if (h == intern_symbol("Tanh") || h == intern_symbol("Coth")) sg = -1;
        if (sg) {
            Expr* u = e->data.function.args[0];
            bool dup = false;
            for (size_t i = 0; i < *n; i++)
                if ((*sigs)[i] == sg && expr_eq((*args)[i], u)) { dup = true; break; }
            if (!dup) {
                if (*n == *cap) { *cap = *cap ? *cap * 2 : 4;
                    *args = realloc(*args, *cap * sizeof(Expr*));
                    *sigs = realloc(*sigs, *cap * sizeof(long)); }
                (*args)[*n] = expr_copy(u); (*sigs)[*n] = sg; (*n)++;
            }
        }
    }
    rt_collect_tangents(e->data.function.head, x, args, sigs, n, cap);
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        rt_collect_tangents(e->data.function.args[i], x, args, sigs, n, cap);
}
/* True when k occurs in e as a literal subtree.  Not FreeQ: a member value
 * is built unevaluated, and evaluating it would merge (E^w)^2 into E^(2 w). */
static bool rt_occurs_literally(const Expr* e, const Expr* k) {
    if (!e) return false;
    if (expr_eq((Expr*)e, (Expr*)k)) return true;
    if (e->type != EXPR_FUNCTION) return false;
    if (rt_occurs_literally(e->data.function.head, k)) return true;
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (rt_occurs_literally(e->data.function.args[i], k)) return true;
    return false;
}
/* Build the ordered differential tower of f over C(x).  Collect every
 * x-dependent Log and E^ kernel, order them innermost-first (deepest at index 0)
 * by structural containment — tie-breaking independent kernels EXP-before-LOG so
 * the primitive (logarithmic) recursion sits on top and the exponential Risch DEs
 * bottom out in C(x) — assign tower variables, and compute each derivation
 * coefficient Dcoef_i (log: u_i'/u_i ; exp: w_i').  The structure-theorem
 * soundness check requires every Dcoef_i to lie in K_{i-1} = C(x, t_1..t_{i-1})
 * (triangular: free of t_i..t_n and of any residual foreign kernel).  Returns
 * true with T populated (n >= min_n; no upper depth cap — the arrays are
 * heap-allocated to the actual kernel count); false otherwise (caller still calls
 * rt_tower_free). */
static bool rt_tower_build_core(Expr* f, Expr* x, RtTower* T, size_t min_n,
                                Expr** xk, Expr** xr, size_t nx) {
    T->kind = NULL; T->kernel = NULL; T->arg = NULL; T->t = NULL; T->Dcoef = NULL;
    T->tsg = NULL; T->marg = NULL; T->mprim = NULL; T->mmult = NULL; T->mval = NULL;
    T->subrules = NULL; T->n = 0; T->nm = 0; T->dependent = false;
    T->nlm = 0; T->lmk = NULL; T->lmv = NULL;

    /* Knowles Liouvillian primitives (Erf/Erfi/Erfc/ExpIntegralEi/LogIntegral).
     * Collect them FIRST: when present, a merged prim-bearing exponential must be
     * split (E^(-x^2 - Erf[x]^2) -> E^(-x^2) E^(-Erf[x]^2)) so each factor is a
     * proper tower monomial — otherwise the merged kernel's Dcoef would reference a
     * lower exponential that is not itself a tower generator and the structure check
     * would reject it.  We therefore collect the log/exp/tan kernels from the
     * prim-split form fx.  Elementary integrands carry no primitive kernel: fx == f,
     * so this is byte-identical to before (rt_expand_exp_sums is not even called). */
    Expr** prims = NULL; size_t npr = 0, prc = 0;
    rt_collect_primitives(f, x, &prims, &npr, &prc);
    Expr* fx = (npr > 0) ? rt_expand_exp_sums(f) : f;

    Expr** logs = NULL; size_t nl = 0, lc = 0; rt_collect_logs(fx, x, &logs, &nl, &lc);
    Expr** exps = NULL; size_t ne = 0, ec = 0;
    rt_collect_exp_exponents(fx, x, &exps, &ne, &ec);
    Expr** tans = NULL; long* tsigs = NULL; size_t nt = 0, ntc = 0;
    rt_collect_tangents(fx, x, &tans, &tsigs, &nt, &ntc);   /* tangent-family monomials */
    if (npr > 0) expr_free(fx);

    /* CLOSE the tower under each primitive's derivative — theta' introduces exp/log
     * kernels (Erf -> E^(-u^2), Ei -> E^u, li -> Log[u]) that must themselves be
     * lower tower monomials.  Split the derivative's prim-bearing exponentials the
     * same way and seed the exps/logs/prims arrays BEFORE the commensurability pass
     * (npr may grow via nested SF like Erf[Erf[x]] — the loop condition re-reads it). */
    for (size_t i = 0; i < npr; i++) {
        Expr* dP = rt_eval2("D", expr_copy(prims[i]), expr_copy(x));
        if (dP) {
            Expr* dPx = rt_expand_exp_sums(dP);
            rt_collect_logs(dPx, x, &logs, &nl, &lc);
            rt_collect_exp_exponents(dPx, x, &exps, &ne, &ec);
            rt_collect_primitives(dPx, x, &prims, &npr, &prc);
            expr_free(dPx); expr_free(dP);
        }
    }

    /* the dependent logarithms (members, not levels: see rt_tower_build_min) */
    if (nx) {
        size_t k = 0;
        for (size_t i = 0; i < nl; i++) {
            bool ex = false;
            for (size_t j = 0; j < nx && !ex; j++) if (expr_eq(logs[i], xk[j])) ex = true;
            if (ex) expr_free(logs[i]); else logs[k++] = logs[i];
        }
        nl = k;
        /* ... and the dependent exponentials (rt_dependent_exp_rule): xk = E^w */
        k = 0;
        for (size_t i = 0; i < ne; i++) {
            bool ex = false;
            for (size_t j = 0; j < nx && !ex; j++)
                if (rt_head_is(xk[j], "Power") && xk[j]->data.function.arg_count == 2
                    && expr_eq(exps[i], xk[j]->data.function.args[1])) ex = true;
            if (ex) expr_free(exps[i]); else exps[k++] = exps[i];
        }
        ne = k;
    }
    Expr** mprim_pexp = malloc((ne ? ne : 1) * sizeof(Expr*)); /* per-member class primitive (borrowed) */

    /* --- Multiplicatively commensurate reduction of the exponential kernels. ---
     * Collected exponents w_i, w_j define algebraically DEPENDENT kernels
     * E^w_i, E^w_j when w_i/w_j is a nonzero rational: then E^w = (E^prim)^k for
     * a class primitive exponent `prim` and integer k (e.g. E^(2 E^x) =
     * (E^(E^x))^2).  Partition the exponents into such commensurability classes,
     * SYNTHESIZE one primitive exponent per class (prim = member/lcm(ratio
     * denominators), possibly NOT itself a member — e.g. E^x for {2 E^x, 3 E^x},
     * or E^(x/6) for {E^(x/2), E^(x/3)}), keep ONLY the synthesized primitives as
     * tower extensions, and record each dependent member as the integer power
     * E^w -> t[prim]^k of its primitive's tower variable (T->m*).  Without this a
     * dependent kernel would spuriously add an extension, breaking independence.
     * A class whose members are not all rational multiples of one another has no
     * common primitive and declines the whole tower (never wrong). */
    long* clsrep = malloc((ne ? ne : 1) * sizeof(long));
    long* multof = malloc((ne ? ne : 1) * sizeof(long));
    for (size_t i = 0; i < ne; i++) { clsrep[i] = -1; multof[i] = 1; }
    bool okc = true;
    for (size_t i = 0; i < ne && okc; i++) {        /* group by commensurability */
        if (clsrep[i] != -1) continue;
        clsrep[i] = (long)i;
        for (size_t j = i + 1; j < ne; j++) {
            if (clsrep[j] != -1) continue;
            Expr* r = rt_eval1("Cancel", expr_new_function(expr_new_symbol("Times"),
                (Expr*[]){ expr_copy(exps[j]), expr_new_function(expr_new_symbol("Power"),
                    (Expr*[]){ expr_copy(exps[i]), expr_new_integer(-1) }, 2) }, 2));
            if (rt_is_rat_const(r)) clsrep[j] = (long)i;
            if (r) expr_free(r);
        }
    }
    /* Synthesize one primitive exponent per class (clsprim[rep], owned) and the
     * per-member integer multiplier multof[i] (w_i = multof[i] * clsprim[rep]). */
    Expr** clsprim = calloc((ne ? ne : 1), sizeof(Expr*));
    for (size_t rep = 0; rep < ne && okc; rep++) {
        if (clsrep[rep] != (long)rep) continue;
        Expr** cls = malloc((ne ? ne : 1) * sizeof(Expr*));
        size_t* idxs = malloc((ne ? ne : 1) * sizeof(size_t));
        size_t cm = 0;
        for (size_t m = 0; m < ne; m++)
            if (clsrep[m] == (long)rep) { cls[cm] = exps[m]; idxs[cm] = m; cm++; }
        long* kk = malloc((cm ? cm : 1) * sizeof(long));
        Expr* p = rt_class_primitive(cls, cm, kk);   /* owned primitive, or NULL */
        if (!p) okc = false;
        else {
            clsprim[rep] = p;
            for (size_t c = 0; c < cm; c++) multof[idxs[c]] = kk[c];
        }
        free(cls); free(idxs); free(kk);
    }
    /* Composite classes (review F12).  The commensurability classes catch
     * rational multiples of ONE exponent; a primitive that is an integer
     * combination of several others (2u + u^2 beside 2u and u^2 -- the
     * evaluator re-merges E^(2u) E^(u^2) into E^(2u + u^2)) made the tower
     * dependent, and every certificate over it was demoted.  Such a class gets
     * no level: its members become the monomials Prod t_s^(m k_s) of the
     * surviving primitives.  Largest exponent first; a primitive another
     * combination uses is never eliminated. */
    bool* elim = calloc((ne ? ne : 1), sizeof(bool));
    bool* used = calloc((ne ? ne : 1), sizeof(bool));
    long** ccoef = calloc((ne ? ne : 1), sizeof(long*));   /* per eliminated rep: k over reps */
    if (okc) {
        size_t* order = malloc((ne ? ne : 1) * sizeof(size_t));
        size_t no = 0;
        for (size_t rep = 0; rep < ne; rep++) if (clsrep[rep] == (long)rep) order[no++] = rep;
        for (size_t a = 1; a < no; a++)                  /* by descending size */
            for (size_t b = a; b > 0 && rt_leaves(clsprim[order[b]]) > rt_leaves(clsprim[order[b - 1]]); b--) {
                size_t tmp = order[b]; order[b] = order[b - 1]; order[b - 1] = tmp;
            }
        Expr* xl = expr_new_function(expr_new_symbol("List"), (Expr*[]){ expr_copy(x) }, 1);
        for (size_t a = 0; a < no; a++) {
            size_t r = order[a];
            if (used[r]) continue;
            Expr** gens = malloc((no ? no : 1) * sizeof(Expr*));
            size_t* greps = malloc((no ? no : 1) * sizeof(size_t));
            size_t ng = 0;
            for (size_t b = 0; b < no; b++)
                if (order[b] != r && !elim[order[b]]) { gens[ng] = clsprim[order[b]]; greps[ng] = order[b]; ng++; }
            if (ng >= 2) {
                bool failed = false;
                Expr* co = risch_rational_span(clsprim[r], gens, ng, xl, &failed);
                bool integral = co && co->type == EXPR_FUNCTION && co->data.function.arg_count == ng;
                for (size_t g = 0; g < ng && integral; g++)
                    if (co->data.function.args[g]->type != EXPR_INTEGER) integral = false;
                size_t nz = 0;
                for (size_t g = 0; g < ng && integral; g++)
                    if (co->data.function.args[g]->data.integer != 0) nz++;
                if (integral && nz >= 2) {      /* nz == 1 is a commensurate multiple (handled) */
                    elim[r] = true;
                    ccoef[r] = calloc(ne, sizeof(long));
                    for (size_t g = 0; g < ng; g++) {
                        long k = (long)co->data.function.args[g]->data.integer;
                        ccoef[r][greps[g]] = k;
                        if (k) used[greps[g]] = true;
                    }
                }
                if (co) expr_free(co);
            }
            free(gens); free(greps);
        }
        expr_free(xl);
        free(order);
    }
    size_t np = 0;
    for (size_t rep = 0; rep < ne; rep++) if (clsrep[rep] == (long)rep && !elim[rep]) np++;
    size_t n = nl + np + nt + npr;
    if (!okc || n < min_n) {                         /* no upper depth cap */
        if (getenv("MATHILDA_RDE_CERT_DEBUG"))
            fprintf(stderr, "[rt-tower] not built: %s (n=%zu)\n", okc ? "fewer kernels than min_n" : "kernel collection failed", n);
        for (size_t rep = 0; rep < ne; rep++) if (clsprim[rep]) expr_free(clsprim[rep]);
        free(clsprim); free(clsrep); free(multof); free(mprim_pexp);
        for (size_t rep = 0; rep < ne; rep++) free(ccoef[rep]);
        free(ccoef); free(elim); free(used);
        for (size_t i = 0; i < nl; i++) expr_free(logs[i]);
        for (size_t i = 0; i < ne; i++) expr_free(exps[i]);
        for (size_t i = 0; i < nt; i++) expr_free(tans[i]);
        for (size_t i = 0; i < npr; i++) expr_free(prims[i]);
        free(logs); free(exps); free(tans); free(tsigs); free(prims);
        return false;
    }
    /* Allocate the tower arrays now the depth n and member bound ne are known. */
    T->kind   = calloc(n ? n : 1, sizeof(RtKind));
    T->kernel = calloc(n ? n : 1, sizeof(Expr*));
    T->arg    = calloc(n ? n : 1, sizeof(Expr*));
    T->t      = calloc(n ? n : 1, sizeof(Expr*));
    T->Dcoef  = calloc(n ? n : 1, sizeof(Expr*));
    T->tsg    = calloc(n ? n : 1, sizeof(long));
    T->marg   = calloc(ne ? ne : 1, sizeof(Expr*));
    T->mprim  = calloc(ne ? ne : 1, sizeof(long));
    T->mmult  = calloc(ne ? ne : 1, sizeof(long));
    T->mval   = calloc(ne ? ne : 1, sizeof(Expr*));
    long* mrep = calloc(ne ? ne : 1, sizeof(long));   /* per member: its class rep */

    size_t idx = 0;
    for (size_t i = 0; i < nl; i++) {
        T->kind[idx] = RT_LOG;
        T->kernel[idx] = logs[i];                                  /* adopt Log[u] */
        T->arg[idx] = expr_copy(logs[i]->data.function.args[0]);
        idx++;
    }
    free(logs);
    /* One EXP tower kernel per synthesized class primitive. */
    for (size_t rep = 0; rep < ne; rep++) {
        if (clsrep[rep] != (long)rep || elim[rep]) continue;
        T->kind[idx] = RT_EXP;
        T->arg[idx] = expr_copy(clsprim[rep]);                     /* primitive exponent */
        T->kernel[idx] = expr_new_function(expr_new_symbol("Power"),
            (Expr*[]){ expr_new_symbol("E"), expr_copy(clsprim[rep]) }, 2);
        idx++;
    }
    /* One RT_TAN tower kernel per collected tangent monomial (Tan/Tanh spelling). */
    for (size_t i = 0; i < nt; i++) {
        T->kind[idx] = RT_TAN;
        T->tsg[idx] = tsigs[i];
        T->arg[idx] = expr_copy(tans[i]);
        T->kernel[idx] = expr_new_function(
            expr_new_symbol(tsigs[i] > 0 ? "Tan" : "Tanh"),
            (Expr*[]){ expr_copy(tans[i]) }, 1);
        idx++;
    }
    for (size_t i = 0; i < nt; i++) expr_free(tans[i]);
    free(tans); free(tsigs);
    /* One RT_PRIM tower kernel per collected primitive SF.  theta' = D[kernel]
     * is computed in the Dcoef pass below (expressed in the lower tower vars). */
    for (size_t i = 0; i < npr; i++) {
        T->kind[idx] = RT_PRIM;
        T->kernel[idx] = prims[i];                                 /* adopt SF[u] */
        T->arg[idx] = expr_copy(prims[i]->data.function.args[0]);  /* u */
        idx++;
    }
    free(prims);
    /* Every exp member that is not itself the class primitive is an alias
     * E^w -> t[prim]^k; store its primitive exponent for the post-reorder remap. */
    for (size_t i = 0; i < ne; i++) {
        if (!elim[clsrep[i]] && expr_eq(exps[i], clsprim[clsrep[i]]))
            continue;                                         /* realized by the primitive kernel */
        T->marg[T->nm] = expr_copy(exps[i]);                  /* member exponent w */
        T->mmult[T->nm] = multof[i];                          /* w = k * clsprim[rep]  */
        mprim_pexp[T->nm] = clsprim[clsrep[i]];               /* borrowed; remapped below */
        mrep[T->nm] = clsrep[i];
        T->nm++;
    }
    T->n = n;

    /* Order innermost-first (deepest at index 0); tie-break EXP/TAN before LOG.
     * A LOG monomial's derivation coefficient may be a rational function of an
     * EXP or TAN monomial (e.g. D[Log[Cos x]] = -Tan[x], D[Log[1+E^x]] =
     * E^x/(1+E^x)), so those non-log monomials must sit BELOW the log for the
     * triangular structure-theorem check to pass.  Conversely an exponential /
     * tangent kernel's own coefficient lies in C(x) (free of logs), so it can
     * always be pushed deeper.  When the LOG is genuinely nested inside the
     * EXP/TAN argument the rt_contains test above keeps the correct order. */
    /* Containment first, as a topological order: every kernel above every
     * kernel it contains.  The bubble below compares ADJACENT pairs only, and
     * containment is not a total order, so a kernel separated from one it
     * contains by an unrelated level never moved past it: E^(u Tan[1 + u])
     * stayed below Tan[1 + u] across an E^u, its Dcoef then named a higher
     * level, and the tower was refused (a silent decline of the decision,
     * found by the f8 hunt, 2026-11-19).  From a topological order, the
     * bubble's swaps of adjacent INDEPENDENT pairs keep it topological. */
    for (bool moved = true; moved; ) {
        moved = false;
        for (size_t i = 0; i < n && !moved; i++)
            for (size_t j = i + 1; j < n && !moved; j++) {
                if (!rt_kernel_needs(T->kernel[i], T, j)) continue;
                /* rotate level j down to position i */
                RtKind kk = T->kind[j]; long sg = T->tsg[j];
                Expr* kj = T->kernel[j]; Expr* aj = T->arg[j];
                for (size_t m = j; m > i; m--) {
                    T->kind[m] = T->kind[m - 1]; T->tsg[m] = T->tsg[m - 1];
                    T->kernel[m] = T->kernel[m - 1]; T->arg[m] = T->arg[m - 1];
                }
                T->kind[i] = kk; T->tsg[i] = sg; T->kernel[i] = kj; T->arg[i] = aj;
                moved = true;
            }
    }
    for (size_t pass = 0; pass < n; pass++)
        for (size_t i = 0; i + 1 < n; i++) {
            bool swap = false;
            if (rt_contains(T->kernel[i], T->kernel[i + 1])) swap = true;
            else if (!rt_contains(T->kernel[i + 1], T->kernel[i])
                     && (T->kind[i] == RT_LOG || T->kind[i] == RT_PRIM)
                     && (T->kind[i + 1] == RT_EXP || T->kind[i + 1] == RT_TAN)) swap = true;
            /* Knowles towers only (npr>0): two monomials independent by structural
             * containment can still be dependency-ordered when one's argument nests a
             * transcendental kernel and the other's is base-field.  E.g. E^(-Erf[x]^2)
             * and E^(-x^2) do not contain each other, yet the former's derivation
             * coefficient needs Erf[x] (hence E^(-x^2)) below it.  Sink the base one
             * deeper.  Gated on npr>0 so elementary towers are byte-identical. */
            else if (npr > 0
                     && !rt_contains(T->kernel[i + 1], T->kernel[i])
                     && !rt_contains(T->kernel[i], T->kernel[i + 1])
                     && rt_has_explog_kernel(T->arg[i])
                     && !rt_has_explog_kernel(T->arg[i + 1])) swap = true;
            if (swap) {
                RtKind kk = T->kind[i]; T->kind[i] = T->kind[i + 1]; T->kind[i + 1] = kk;
                long sg = T->tsg[i]; T->tsg[i] = T->tsg[i + 1]; T->tsg[i + 1] = sg;
                Expr* a = T->kernel[i]; T->kernel[i] = T->kernel[i + 1]; T->kernel[i + 1] = a;
                a = T->arg[i]; T->arg[i] = T->arg[i + 1]; T->arg[i + 1] = a;
            }
        }

    /* Tower variables t_i and the combined substitution rule list.  Each member
     * kernel E^(marg) contributes a rule E^(marg) -> t[prim]^mmult (both the
     * Exp[] and Power[E,] spellings); mprim is remapped from an exps[] index to
     * the primitive's post-ordering tower index by matching its exponent. */
    Expr** rules = malloc((6 * n + 2 * T->nm + nx + 1) * sizeof(Expr*)); size_t nr = 0;
    for (size_t i = 0; i < n; i++) {
        char nm[32]; snprintf(nm, sizeof(nm), "rmR%zu", i);
        T->t[i] = expr_new_symbol(nm);
        if (T->kind[i] == RT_LOG || T->kind[i] == RT_PRIM) {
            rules[nr++] = expr_new_function(expr_new_symbol("Rule"),
                (Expr*[]){ expr_copy(T->kernel[i]), expr_copy(T->t[i]) }, 2);
        } else if (T->kind[i] == RT_TAN) {
            /* Rationalise the whole trig family of the tangent argument to the tower
             * variable: Tan/Tanh -> t, Cot/Coth -> 1/t, and Sin=t/rad, Cos=1/rad,
             * Sec=rad, Csc=rad/t with rad = Sqrt[1 + sigma t^2] (sigma = T->tsg).
             * The evaluator-canonical form Csc*Sec of (1+Tan^2)/Tan then substitutes. */
            long sg = T->tsg[i]; Expr* u = T->arg[i]; Expr* ti = T->t[i];
            Expr* rad = expr_new_function(expr_new_symbol("Power"),
                (Expr*[]){ expr_new_function(expr_new_symbol("Plus"),
                    (Expr*[]){ expr_new_integer(1), expr_new_function(expr_new_symbol("Times"),
                        (Expr*[]){ expr_new_integer(sg), expr_new_function(expr_new_symbol("Power"),
                            (Expr*[]){ expr_copy(ti), expr_new_integer(2) }, 2) }, 2) }, 2),
                  expr_new_function(expr_new_symbol("Rational"),
                    (Expr*[]){ expr_new_integer(1), expr_new_integer(2) }, 2) }, 2);
            Expr* invt = expr_new_function(expr_new_symbol("Power"),
                (Expr*[]){ expr_copy(ti), expr_new_integer(-1) }, 2);
            #define RT_TRULE(HD, RHS) rules[nr++] = expr_new_function(expr_new_symbol("Rule"), \
                (Expr*[]){ expr_new_function(expr_new_symbol(HD), (Expr*[]){ expr_copy(u) }, 1), (RHS) }, 2)
            rules[nr++] = expr_new_function(expr_new_symbol("Rule"),   /* Tan/Tanh kernel -> t */
                (Expr*[]){ expr_copy(T->kernel[i]), expr_copy(ti) }, 2);
            RT_TRULE(sg > 0 ? "Cot" : "Coth", expr_copy(invt));
            RT_TRULE(sg > 0 ? "Sin" : "Sinh", expr_new_function(expr_new_symbol("Times"),
                (Expr*[]){ expr_copy(ti), expr_new_function(expr_new_symbol("Power"),
                    (Expr*[]){ expr_copy(rad), expr_new_integer(-1) }, 2) }, 2));
            RT_TRULE(sg > 0 ? "Cos" : "Cosh", expr_new_function(expr_new_symbol("Power"),
                (Expr*[]){ expr_copy(rad), expr_new_integer(-1) }, 2));
            RT_TRULE(sg > 0 ? "Sec" : "Sech", expr_copy(rad));
            RT_TRULE(sg > 0 ? "Csc" : "Csch", expr_new_function(expr_new_symbol("Times"),
                (Expr*[]){ expr_copy(rad), expr_copy(invt) }, 2));
            #undef RT_TRULE
            expr_free(rad); expr_free(invt);
        } else {
            rules[nr++] = expr_new_function(expr_new_symbol("Rule"),
                (Expr*[]){ expr_new_function(expr_new_symbol("Exp"),
                    (Expr*[]){ expr_copy(T->arg[i]) }, 1), expr_copy(T->t[i]) }, 2);
            rules[nr++] = expr_new_function(expr_new_symbol("Rule"),
                (Expr*[]){ expr_copy(T->kernel[i]), expr_copy(T->t[i]) }, 2);
        }
    }
    for (size_t m = 0; m < T->nm; m++) {
        Expr* tk;
        if (!elim[mrep[m]]) {
            Expr* pw = mprim_pexp[m];              /* class primitive exponent value */
            long ti = -1;
            for (size_t i = 0; i < n; i++)
                if (T->kind[i] == RT_EXP && expr_eq(T->arg[i], pw)) { ti = (long)i; break; }
            T->mprim[m] = ti;                      /* now a tower index (>= 0) */
            tk = expr_new_function(expr_new_symbol("Power"),
                (Expr*[]){ expr_copy(T->t[ti]), expr_new_integer(T->mmult[m]) }, 2);
        } else {
            /* composite: w = multof * Sum_s k_s prim_s  ->  Prod t_s^(multof k_s) */
            long* kc = ccoef[mrep[m]];
            Expr** fs = malloc((ne ? ne : 1) * sizeof(Expr*));
            size_t nf = 0;
            T->mprim[m] = -1;
            for (size_t srep = 0; srep < ne; srep++) {
                if (!kc[srep]) continue;
                long ti = -1;
                for (size_t i = 0; i < n; i++)
                    if (T->kind[i] == RT_EXP && expr_eq(T->arg[i], clsprim[srep])) { ti = (long)i; break; }
                if (T->mprim[m] < 0) T->mprim[m] = ti;
                fs[nf++] = expr_new_function(expr_new_symbol("Power"),
                    (Expr*[]){ expr_copy(T->t[ti]), expr_new_integer(T->mmult[m] * kc[srep]) }, 2);
            }
            tk = expr_new_function(expr_new_symbol("Times"), fs, nf);
            free(fs);
        }
        T->mval[m] = expr_copy(tk);
        rules[nr++] = expr_new_function(expr_new_symbol("Rule"),
            (Expr*[]){ expr_new_function(expr_new_symbol("Exp"),
                (Expr*[]){ expr_copy(T->marg[m]) }, 1), expr_copy(tk) }, 2);
        rules[nr++] = expr_new_function(expr_new_symbol("Rule"),
            (Expr*[]){ expr_new_function(expr_new_symbol("Power"),
                (Expr*[]){ expr_new_symbol("E"), expr_copy(T->marg[m]) }, 2), tk }, 2);
    }
    if (nx) {                                      /* dependent Log members */
        T->lmk = calloc(nx, sizeof(Expr*));
        T->lmv = calloc(nx, sizeof(Expr*));
        for (size_t j = 0; j < nx; j++) {
            T->lmk[j] = expr_copy(xk[j]);
            T->lmv[j] = rt_subst_kernels(xr[j], T);   /* sees members 0..j-1 */
            T->nlm = j + 1;
        }
        /* A member's value is written in the levels of the tower it was found
         * in, and a LATER round may eliminate one of those levels in turn
         * (gauss-elem-27: E^(-I + 2 I x) -> Z^-3 (E^(-I (1 + x)))^-2, then
         * E^(-I (1 + x)) -> Z^-3 / E^(-2 I + I x)): the earlier value still
         * names the now-eliminated kernel, a foreign kernel after
         * substitution.  Resolve it latest-first, each value through the
         * members after it (each round's value names only the levels of its
         * own round, so the chain has no cycle). */
        for (size_t j = nx; j-- > 1; ) {
            for (size_t m = j; m < nx; m++) {
                if (!rt_occurs_literally(T->lmv[j - 1], T->lmk[m])) continue;
                Expr* r = rt_subst_kernels(T->lmv[j - 1], T);
                expr_free(T->lmv[j - 1]);
                T->lmv[j - 1] = r;
                break;
            }
        }
        for (size_t j = 0; j < nx; j++)
            rules[nr++] = expr_new_function(expr_new_symbol("Rule"),
                (Expr*[]){ expr_copy(T->lmk[j]), expr_copy(T->lmv[j]) }, 2);
    }
    T->subrules = expr_new_function(expr_new_symbol("List"), rules, nr);
    free(rules);
    for (size_t rep = 0; rep < ne; rep++) if (clsprim[rep]) expr_free(clsprim[rep]);
    free(clsprim); free(clsrep); free(multof); free(mprim_pexp);
    for (size_t rep = 0; rep < ne; rep++) free(ccoef[rep]);
    free(ccoef); free(elim); free(used); free(mrep);
    for (size_t i = 0; i < ne; i++) expr_free(exps[i]);
    free(exps);

    /* Derivation coefficients + structure-theorem (triangularity) soundness. */
    bool ok = true;
    for (size_t i = 0; i < n && ok; i++) {
        /* Compute the RAW derivative (Log: u'/u, Exp: w') WITHOUT normalising,
         * substitute the kernels structurally, and ONLY THEN Cancel in the tower
         * variables.  Cancelling first expands a factored kernel power such as
         * (x + E^x)^2 into x^2 + 2 x E^x + E^(2 x); the E^(2 x) term is not the
         * literal kernel E^x, so the substitution misses it and the structure
         * check spuriously rejects a valid tower (e.g. Bronstein's own example
         * E^((x^2-1)/x + 1/(x + E^x))).  Substituting on the unexpanded form maps
         * (x + E^x)^2 -> (x + t_0)^2 cleanly. */
        Expr* d;
        if (T->kind[i] == RT_LOG) {
            Expr* du = rt_eval2("D", expr_copy(T->arg[i]), expr_copy(x));
            d = expr_new_function(expr_new_symbol("Times"),
                (Expr*[]){ du, expr_new_function(expr_new_symbol("Power"),
                    (Expr*[]){ expr_copy(T->arg[i]), expr_new_integer(-1) }, 2) }, 2);
        } else if (T->kind[i] == RT_TAN) {
            /* Dcoef = Dt/(t^2 + sigma) = sigma * u': D[Tan[u]]=(t^2+1)u' (sigma=+1),
             * D[Tanh[u]]=(1-t^2)u' = -(t^2-1)u' (sigma=-1). */
            Expr* du = rt_eval2("D", expr_copy(T->arg[i]), expr_copy(x));
            d = expr_new_function(expr_new_symbol("Times"),
                (Expr*[]){ expr_new_integer(T->tsg[i]), du }, 2);
        } else if (T->kind[i] == RT_PRIM) {
            /* theta' = D[SF[u]] straight from deriv.c — this introduces the lower
             * exp/log kernel (Erf -> E^(-u^2), Ei -> E^u, li -> Log[u]) which the
             * kernel substitution below folds back to the tower variable, so Dcoef
             * lands in K_{i-1}. */
            d = rt_eval2("D", expr_copy(T->kernel[i]), expr_copy(x));
        } else {
            d = rt_eval2("D", expr_copy(T->arg[i]), expr_copy(x));   /* RT_EXP: w' */
        }
        if (!d) { ok = false; break; }
        /* Rewrite the tangent-derivative trig squares so the tower substitution
         * catches them (D[Tan[a]] = Sec[a]^2 would otherwise strand a foreign
         * kernel): Sec[a]^2 -> 1+Tan[a]^2, Csc -> 1+Cot^2, Sech -> 1-Tanh^2,
         * Csch -> Coth^2-1. */
        /* The right-hand sides name the pattern variable, NOT the pattern: with
         * `Tan[rmTA_]` on the right the binding substituted into Pattern[rmTA, _],
         * so D[Tan[u]] gave 1 + Tan[u_]^2, a foreign kernel that no tower variable
         * matches, and every exponential of a tangent argument (E^Tan[u] Sec[u]^2)
         * carried it in its Dcoef (review F8, v0.6601).  Private context, so a
         * user symbol of the same name is not captured. */
        Expr* trig = parse_expression(
            "{Sec[Risch`Private`rmTA_]^2 :> 1 + Tan[Risch`Private`rmTA]^2, "
            "Csc[Risch`Private`rmTA_]^2 :> 1 + Cot[Risch`Private`rmTA]^2, "
            "Sech[Risch`Private`rmTA_]^2 :> 1 - Tanh[Risch`Private`rmTA]^2, "
            "Csch[Risch`Private`rmTA_]^2 :> Coth[Risch`Private`rmTA]^2 - 1}");
        if (trig) d = rt_eval_own(expr_new_function(expr_new_symbol("ReplaceRepeated"),
            (Expr*[]){ d, trig }, 2));                             /* adopts d, trig */
        if (!d) { ok = false; break; }
        /* Split any exponential product the D-evaluation MERGED back into a summed
         * exponent (D[E^(E^x)] * x = x E^(E^x) E^x, which the evaluator canonicalises
         * to x E^(x + E^x)); without this the merged E^(x+E^x) is not the literal
         * kernel E^x / E^(E^x) that T->subrules aliases, so it reads as a foreign
         * exponential and the structure check spuriously rejects a valid tower — the
         * depth-2 nested-exp-trig case Cos[x E^(E^x)] after TrigToExp. */
        { Expr* de = rt_expand_exp_sums(d); expr_free(d); d = de; }
        if (!d) { ok = false; break; }
        /* Substitute the kernels to tower variables STRUCTURALLY (rt_subst_kernels),
         * NOT via an evaluated ReplaceAll: the latter re-evaluates and re-merges the
         * exp product just split above (E^(E^x) E^x -> E^(x+E^x)), stranding a foreign
         * kernel.  rt_subst_kernels also carries the tangent-trig rationalisation, so
         * the evaluator-canonicalised (1+Tan^2)/Tan = Csc*Sec still maps cleanly. */
        Expr* dsub = rt_subst_kernels(d, T);
        expr_free(d);
        Expr* ds = dsub ? rt_eval1("Cancel", dsub) : NULL;
        if (!ds) { ok = false; break; }
        T->Dcoef[i] = ds;
        if (rt_find_exp_of_x(ds, x) != NULL || rt_find_log_of_x(ds, x) != NULL) ok = false;
        for (size_t j = i; j < n && ok; j++)
            if (!rt_free_of_x(ds, T->t[j])) ok = false;
        if (!ok && getenv("MATHILDA_RDE_CERT_DEBUG")) {
            char* a = expr_to_string(T->kernel[i]); char* b = expr_to_string(ds);
            fprintf(stderr, "[rt-tower] not built: level %zu %s has Dcoef %.300s outside the levels below\n", i, a, b);
            free(a); free(b);
        }
    }
    /* Algebraic independence of the Log/Exp monomials (Risch structure
     * theorem, Bronstein Thm 9.3.x): E^w_i or Log[u_i] is a new monomial over
     * the earlier ones iff Dcoef_i (w_i' or u_i'/u_i) is NOT in the Q-span of
     * their Dcoef.  The builder's own commensurability classes catch rational
     * multiples of one exponent only; a sum of exponents (E^(x + x^2) beside
     * E^x and E^(x^2)) slips through.  Checking each monomial against those
     * before it finds any dependence in the set. */
    if (ok && n >= 2) {
        size_t nv = n + 1;
        Expr** vv = malloc(nv * sizeof(Expr*));
        vv[0] = expr_copy(x);
        for (size_t j2 = 0; j2 < n; j2++) vv[j2 + 1] = expr_copy(T->t[j2]);
        Expr* vars = expr_new_function(expr_new_symbol("List"), vv, nv);
        free(vv);
        Expr** gens = malloc(n * sizeof(Expr*));
        size_t ng = 0;
        for (size_t i2 = 0; i2 < n && !T->dependent; i2++) {
            if (T->kind[i2] != RT_LOG && T->kind[i2] != RT_EXP) continue;
            if (ng > 0) {
                bool failed = false;
                Expr* co = risch_rational_span(T->Dcoef[i2], gens, ng, vars, &failed);
                if (co || failed) T->dependent = true;
                if (co) expr_free(co);
            }
            gens[ng++] = T->Dcoef[i2];                 /* borrowed */
        }
        free(gens);
        expr_free(vars);
    }
    return ok;
}

/* ---- Dependent logarithms (C completeness campaign, step 8) ---------------- */

/* Transcendental number atoms of e (Pi, E, EulerGamma, Log[<number>]), collected
 * distinct into a[] (borrowed); x-dependent E^w contributes its exponent only. */
static void rt_tconst_atoms(Expr* e, Expr* x, Expr** a, size_t* n, size_t cap) {
    if (!e) return;
    if (e->type == EXPR_SYMBOL) {
        const char* h = e->data.symbol.name;
        /* E^I (RT_EI_SYMBOL) is counted as E: the two are algebraically
         * independent (Lindemann-Weierstrass), one group of known constants */
        if (h == intern_symbol(RT_EI_SYMBOL)) {
            static Expr* esym = NULL;              /* process-lifetime, never freed */
            if (!esym) esym = expr_new_symbol("E");
            e = esym;
        }
        else if (h != intern_symbol("Pi") && h != intern_symbol("E") && h != intern_symbol("EulerGamma")
            && h != intern_symbol("Catalan")) return;
    } else if (e->type != EXPR_FUNCTION) return;
    else if (rt_head_is(e, "Power") && e->data.function.arg_count == 2
             && e->data.function.args[0]->type == EXPR_SYMBOL
             && e->data.function.args[0]->data.symbol.name == intern_symbol("E")) {
        rt_tconst_atoms(e->data.function.args[1], x, a, n, cap);
        if (!rt_free_of_x(e->data.function.args[1], x)) return;
        e = e->data.function.args[0];
    } else if (!(rt_head_is(e, "Log") && e->data.function.arg_count == 1
                 && rt_free_of_x(e->data.function.args[0], x))) {
        for (size_t i = 0; i < e->data.function.arg_count; i++)
            rt_tconst_atoms(e->data.function.args[i], x, a, n, cap);
        return;
    }
    for (size_t i = 0; i < *n; i++) if (expr_eq(a[i], e)) return;
    if (*n < cap) a[(*n)++] = e;
    else (*n)++;                                   /* over capacity: still counted */
}

/* For the first Log level of T that the structure theorem finds DEPENDENT,
 * Log[v] = Sum r_k K_k + c over the lower Log/Exp levels (K_k = Log[u_k], or
 * the exponent w_k of E^w_k), set *kern = Log[v] and *val = that value in
 * kernel form, with c exact: v^d / (Prod u_k^(d r_k) E^(d Sum s_k w_k)) = kappa
 * is a constant (d the lcm of the denominators), so
 *     Log[v] - Sum r_k K_k = (Log[kappa] + 2 Pi I m) / d
 * on each region where the principal branches do not jump; m is read off
 * three real sample points and must agree on all of them (else no rule:
 * the integrand is genuinely piecewise there).  The certificate is then
 * about the germ at those points, which is what a non-elementarity proof
 * needs.  False when no such rule is found. */
static bool rt_dependent_log_rule(RtTower* T, Expr* x, bool wide, Expr** kern, Expr** val) {
    *kern = NULL; *val = NULL;
    size_t n = T->n;
    if (n < 2) return false;
    Expr** vv = malloc((n + 1) * sizeof(Expr*));
    vv[0] = expr_copy(x);
    for (size_t j = 0; j < n; j++) vv[j + 1] = expr_copy(T->t[j]);
    Expr* vars = expr_new_function(expr_new_symbol("List"), vv, n + 1);
    free(vv);
    Expr** gens = malloc(n * sizeof(Expr*));
    size_t* gi = malloc(n * sizeof(size_t));
    size_t ng = 0;
    bool found = false;
    for (size_t i2 = 0; i2 < n && !found; i2++) {
        if (T->kind[i2] != RT_LOG && T->kind[i2] != RT_EXP) continue;
        if (ng > 0 && T->kind[i2] == RT_LOG) {
            bool failed = false;
            Expr* co = risch_rational_span(T->Dcoef[i2], gens, ng, vars, &failed);
            if (co && !failed && co->type == EXPR_FUNCTION && co->data.function.arg_count == ng) {
                bool okr = true;
                long d = 1;
                for (size_t k = 0; k < ng && okr; k++) {
                    Expr* c = co->data.function.args[k];
                    if (c->type == EXPR_INTEGER) continue;
                    if (rt_head_is(c, "Rational") && c->data.function.args[1]->type == EXPR_INTEGER) {
                        long q = (long)c->data.function.args[1]->data.integer, a = d, b = q;
                        while (b) { long t = a % b; a = b; b = t; }
                        if (q > 64 || d / a * q > 4096) okr = false; else d = d / a * q;
                    } else okr = false;
                }
                /* Sum r_k K_k (kernel form) and the ratio kappa */
                Expr** sm = okr ? malloc((ng + 1) * sizeof(Expr*)) : NULL;
                Expr** pf = okr ? malloc((ng + 2) * sizeof(Expr*)) : NULL;
                size_t ns = 0, np = 0;
                Expr* v = T->arg[i2];
                if (okr) pf[np++] = expr_new_function(expr_new_symbol("Power"),
                    (Expr*[]){ expr_copy(v), expr_new_integer(d) }, 2);
                for (size_t k = 0; okr && k < ng; k++) {
                    Expr* c = co->data.function.args[k];
                    if (rt_is_zero(c)) continue;
                    size_t j = gi[k];
                    Expr* K = T->kind[j] == RT_LOG ? expr_copy(T->kernel[j]) : expr_copy(T->arg[j]);
                    sm[ns++] = expr_new_function(expr_new_symbol("Times"), (Expr*[]){ expr_copy(c), K }, 2);
                    /* 1/u^(d r) or E^(-d s w) */
                    Expr* ex = rt_eval_own(expr_new_function(expr_new_symbol("Times"),
                        (Expr*[]){ expr_new_integer(-d), expr_copy(c) }, 2));
                    pf[np++] = T->kind[j] == RT_LOG
                        ? expr_new_function(expr_new_symbol("Power"), (Expr*[]){ expr_copy(T->arg[j]), ex }, 2)
                        : expr_new_function(expr_new_symbol("Power"), (Expr*[]){ expr_new_symbol("E"),
                              expr_new_function(expr_new_symbol("Times"), (Expr*[]){ ex, expr_copy(T->arg[j]) }, 2) }, 2);
                }
                Expr* kap = okr ? rt_eval1("Together", rt_eval_own(expr_new_function(expr_new_symbol("Times"), pf, np))) : NULL;
                if (pf) free(pf);
                Expr* S = okr ? (ns ? rt_eval_own(expr_new_function(expr_new_symbol("Plus"), sm, ns)) : expr_new_integer(0)) : NULL;
                if (sm) free(sm);
                if (kap && S && rt_free_of_x(kap, x) && !rt_is_zero(kap)) {
                    /* m from (d (Log[v] - S) - Log[kappa]) / (2 Pi I) at three points */
                    Expr* lk = rt_eval1("Log", expr_copy(kap));
                    Expr* del = expr_new_function(expr_new_symbol("Times"), (Expr*[]){
                        expr_new_function(expr_new_symbol("Plus"), (Expr*[]){
                            expr_new_function(expr_new_symbol("Times"), (Expr*[]){ expr_new_integer(d),
                                expr_new_function(expr_new_symbol("Plus"), (Expr*[]){ expr_copy(T->kernel[i2]),
                                    expr_new_function(expr_new_symbol("Times"), (Expr*[]){ expr_new_integer(-1), expr_copy(S) }, 2) }, 2) }, 2),
                            expr_new_function(expr_new_symbol("Times"), (Expr*[]){ expr_new_integer(-1), expr_copy(lk) }, 2) }, 2),
                        expr_new_function(expr_new_symbol("Power"), (Expr*[]){
                            expr_new_function(expr_new_symbol("Times"), (Expr*[]){ expr_new_integer(2), expr_new_symbol("Pi"),
                                expr_new_function(expr_new_symbol("Complex"), (Expr*[]){ expr_new_integer(0), expr_new_integer(1) }, 2) }, 3),
                            expr_new_integer(-1) }, 2) }, 2);
                    /* {re num, re den, im num, im den}: three positive reals;
                     * `wide` adds negative reals and both half-planes */
                    static const long pts[8][4] = { { 7, 19, 0, 1 }, { 13, 8, 0, 1 }, { 29, 11, 0, 1 },
                        { -7, 19, 0, 1 }, { -13, 8, 0, 1 }, { 1, 3, 1, 2 }, { -2, 3, -3, 5 }, { 5, 4, -7, 9 } };
                    bool agree = true; long m = 0;
                    for (int p = 0; p < (wide ? 8 : 3) && agree; p++) {
                        Expr* x0 = rt_eval_own(expr_new_function(expr_new_symbol("Plus"), (Expr*[]){
                            expr_new_function(expr_new_symbol("Rational"),
                                (Expr*[]){ expr_new_integer(pts[p][0]), expr_new_integer(pts[p][1]) }, 2),
                            expr_new_function(expr_new_symbol("Times"), (Expr*[]){
                                expr_new_function(expr_new_symbol("Rational"),
                                    (Expr*[]){ expr_new_integer(pts[p][2]), expr_new_integer(pts[p][3]) }, 2),
                                expr_new_function(expr_new_symbol("Complex"),
                                    (Expr*[]){ expr_new_integer(0), expr_new_integer(1) }, 2) }, 2) }, 2));
                        Expr* dn = rt_eval1("N", rt_eval2("ReplaceAll", expr_copy(del),
                            expr_new_function(expr_new_symbol("Rule"), (Expr*[]){ expr_copy(x), x0 }, 2)));
                        Expr* rd = dn ? rt_eval1("Round", rt_eval1("Re", expr_copy(dn))) : NULL;
                        Expr* off = (dn && rd) ? rt_eval1("Abs", rt_eval_own(expr_new_function(expr_new_symbol("Plus"),
                            (Expr*[]){ expr_copy(dn), expr_new_function(expr_new_symbol("Times"),
                                (Expr*[]){ expr_new_integer(-1), expr_copy(rd) }, 2) }, 2))) : NULL;
                        Expr* small = off ? rt_eval2("Less", off, expr_new_real(1e-8)) : NULL;
                        if (!(small && rt_is_true(small) && rd->type == EXPR_INTEGER)) agree = false;
                        else if (p == 0) m = (long)rd->data.integer;
                        else if ((long)rd->data.integer != m) agree = false;
                        if (small) expr_free(small);
                        if (rd) expr_free(rd);
                        if (dn) expr_free(dn);
                    }
                    expr_free(del);
                    if (agree) {
                        /* c = (Log[kappa] + 2 Pi I m) / d */
                        Expr* c = rt_eval1("Expand", rt_eval_own(expr_new_function(expr_new_symbol("Times"), (Expr*[]){
                            expr_new_function(expr_new_symbol("Rational"), (Expr*[]){ expr_new_integer(1), expr_new_integer(d) }, 2),
                            expr_new_function(expr_new_symbol("Plus"), (Expr*[]){ expr_copy(lk),
                                expr_new_function(expr_new_symbol("Times"), (Expr*[]){ expr_new_integer(2 * m), expr_new_symbol("Pi"),
                                    expr_new_function(expr_new_symbol("Complex"), (Expr*[]){ expr_new_integer(0), expr_new_integer(1) }, 2) }, 3) }, 2) }, 2)));
                        if (c) {
                            *kern = expr_copy(T->kernel[i2]);
                            *val = rt_eval_own(expr_new_function(expr_new_symbol("Plus"), (Expr*[]){ expr_copy(S), c }, 2));
                            found = *val != NULL;
                            if (!found) { expr_free(*kern); *kern = NULL; }
                        }
                    }
                    expr_free(lk);
                }
                if (kap) expr_free(kap);
                if (S) expr_free(S);
            }
            if (co) expr_free(co);
        }
        gi[ng] = i2;
        gens[ng++] = T->Dcoef[i2];                     /* borrowed */
    }
    free(gens); free(gi);
    expr_free(vars);
    return found;
}

/* Build the tower; a dependent logarithm (Log[1 + I x] beside Log[x - I]) is
 * then eliminated -- it becomes a member Log[v] = Sum r_j t_j + c with the
 * exact constant c (rt_dependent_log_rule), not a level -- and the tower is
 * rebuilt without it, so the structure theorem's independence, and every
 * certificate it carries, holds again (campaign step 8; before, T->dependent
 * demoted every verdict).  At most one transcendental number may enter the
 * constant field this way (Pi from Log[I], or one Log[2]), so the constants
 * stay algebraically independent; otherwise the dependent tower is kept. */
/* For the first Exp level of T that the structure theorem finds DEPENDENT,
 * w' = Sum k_j Dcoef_j over the lower Log/Exp levels, so
 *     w = Sum_{exp j} k_j w_j + Sum_{log j} k_j Log[u_j] + c
 * with c a constant.  When every k_j is an integer and c is an integer, E^w is
 * no new monomial: set *kern = E^w and *val = E^c Prod (E^w_j)^k_j Prod u_j^k_j,
 * built unevaluated (the evaluator would merge the product back into E^w).
 * c is read exactly in the tower variables, where the lower kernels are
 * independent symbols.  The case it closes: TrigToExp[E^Tanh[1 + u]] carries
 * E^(-r^2/(1 + r^2)) and E^(1/(1 + r^2)) over r = E^(-1 - u), whose exponents
 * differ by exactly 1, and E^(2 + 2 u) beside E^u.  A non-integer k_j makes E^w
 * algebraic over the lower field (out of scope).  A Gaussian-integer
 * c = a + b i is E^a Z^b with Z the symbol RT_EI_SYMBOL for E^I: E and E^I are
 * algebraically independent (Lindemann-Weierstrass), so both are parameters
 * and their powers are exact polynomial arithmetic.  (Left as atoms, E^I and
 * E^(2 I) are unrelated to Together.)  Any other c would bring in a constant
 * whose relations rational arithmetic does not see (E^(1/2), E^(1/3)).
 * False when no such rule is found. */
/* n/d of an Integer or Rational expression; false otherwise. */
static bool rt_rat_parts(const Expr* e, int64_t* n, int64_t* d) {
    if (e->type == EXPR_INTEGER) { *n = e->data.integer; *d = 1; return true; }
    if (e->type == EXPR_FUNCTION && e->data.function.head->type == EXPR_SYMBOL
        && e->data.function.head->data.symbol.name == intern_symbol("Rational")
        && e->data.function.arg_count == 2
        && e->data.function.args[0]->type == EXPR_INTEGER && e->data.function.args[1]->type == EXPR_INTEGER) {
        *n = e->data.function.args[0]->data.integer; *d = e->data.function.args[1]->data.integer;
        return *d != 0;
    }
    return false;
}

static bool rt_dependent_exp_rule(RtTower* T, Expr* x, Expr** kern, Expr** val) {
    *kern = NULL; *val = NULL;
    size_t n = T->n;
    if (n < 2) return false;
    Expr** vv = malloc((n + 1) * sizeof(Expr*));
    vv[0] = expr_copy(x);
    for (size_t j = 0; j < n; j++) vv[j + 1] = expr_copy(T->t[j]);
    Expr* vars = expr_new_function(expr_new_symbol("List"), vv, n + 1);
    free(vv);
    Expr** gens = malloc(n * sizeof(Expr*));
    size_t* gi = malloc(n * sizeof(size_t));
    size_t ng = 0;
    bool found = false;
    for (size_t i2 = 0; i2 < n && !found; i2++) {
        if (T->kind[i2] != RT_LOG && T->kind[i2] != RT_EXP) continue;
        if (ng > 0 && T->kind[i2] == RT_EXP) {
            bool failed = false;
            Expr* co = risch_rational_span(T->Dcoef[i2], gens, ng, vars, &failed);
            bool okk = co && !failed && co->type == EXPR_FUNCTION && co->data.function.arg_count == ng;
            /* the relation w_i2 = Sum k_j w_j + c, with rational k_j */
            int64_t* kn = calloc(ng + 1, sizeof(int64_t));
            int64_t* kd = calloc(ng + 1, sizeof(int64_t));
            bool allint = okk;
            for (size_t k = 0; okk && k < ng; k++) {
                if (!rt_rat_parts(co->data.function.args[k], &kn[k], &kd[k])) okk = false;
                else if (kd[k] != 1) allint = false;
            }
            /* the level eliminated, and the integer multipliers of the others:
             *   allint:  target i2,  w_i2 = Sum k_j w_j + c;
             *   inverted (C2.1): the lower EXPONENTIAL level j0 with k_j0 = +-1/m,
             *     every m k_i an integer -- w_j0 = (1/k_j0) w_i2 - Sum_{i != j0}
             *     (k_i/k_j0) w_i + c'.  The tower orders E^(1 + 2 x) below E^x,
             *     and E^x = (E^(1 + 2 x))^(1/2) E^(-1/2) has no integer form,
             *     while E^(1 + 2 x) = E (E^x)^2 does (benchmark 259: f12/gauss
             *     cases were undecided, "dependent monomials"). */
            long target = -1;                 /* tower index of the eliminated level */
            int64_t* mult = calloc(ng + 1, sizeof(int64_t));   /* per gen; [ng] = i2's */
            if (okk && allint) {
                target = (long)i2;
                for (size_t k = 0; k < ng; k++) mult[k] = kn[k];
                mult[ng] = 0;
            } else if (okk) {
                for (size_t j0 = 0; j0 < ng && target < 0; j0++) {
                    if (T->kind[gi[j0]] != RT_EXP || kn[j0] == 0) continue;
                    if (kn[j0] != 1 && kn[j0] != -1) continue;
                    int64_t m = kd[j0], s = kn[j0];          /* 1/k_j0 = s m */
                    bool ok = true;
                    for (size_t k = 0; k < ng && ok; k++) {
                        if (k == j0) continue;
                        /* -k_k/k_j0 = -k_k s m */
                        if ((kn[k] * m) % kd[k] != 0) ok = false;
                        else mult[k] = -(kn[k] * m / kd[k]) * s;
                    }
                    if (ok) { target = (long)gi[j0]; mult[j0] = 0; mult[ng] = s * m; }
                }
            }
            if (target >= 0) {
                /* c = w_target - Sum mult w, exactly, in the tower variables */
                Expr** tm = malloc((ng + 2) * sizeof(Expr*));
                size_t nt = 0;
                tm[nt++] = rt_subst_kernels(T->arg[target], T);
                for (size_t k = 0; k <= ng; k++) {
                    if (!mult[k]) continue;
                    size_t j = k < ng ? gi[k] : i2;
                    Expr* wj = T->kind[j] == RT_EXP ? rt_subst_kernels(T->arg[j], T) : expr_copy(T->t[j]);
                    tm[nt++] = expr_new_function(expr_new_symbol("Times"),
                        (Expr*[]){ expr_new_integer(-mult[k]), wj }, 2);
                }
                Expr* c = rt_eval1("Together", rt_eval_own(expr_new_function(expr_new_symbol("Plus"), tm, nt)));
                free(tm);
                /* c = ca + cb i, Gaussian integer */
                int64_t ca = 0, cb = 0;
                bool gint = false;
                if (c && c->type == EXPR_INTEGER) { ca = c->data.integer; gint = true; }
                else if (c && rt_head_is(c, "Complex") && c->data.function.arg_count == 2
                         && c->data.function.args[0]->type == EXPR_INTEGER
                         && c->data.function.args[1]->type == EXPR_INTEGER) {
                    ca = c->data.function.args[0]->data.integer;
                    cb = c->data.function.args[1]->data.integer;
                    gint = true;
                }
                if (gint) {
                    Expr** pf = malloc((ng + 3) * sizeof(Expr*));
                    size_t np = 0;
                    if (ca != 0)
                        pf[np++] = expr_new_function(expr_new_symbol("Power"),
                            (Expr*[]){ expr_new_symbol("E"), expr_new_integer(ca) }, 2);
                    /* E^(cb i) = Z^cb, Z the symbol for E^I (RT_EI_SYMBOL) */
                    if (cb != 0)
                        pf[np++] = expr_new_function(expr_new_symbol("Power"),
                            (Expr*[]){ expr_new_symbol(RT_EI_SYMBOL), expr_new_integer(cb) }, 2);
                    for (size_t k = 0; k <= ng; k++) {
                        if (!mult[k]) continue;
                        size_t j = k < ng ? gi[k] : i2;
                        Expr* b = T->kind[j] == RT_EXP ? expr_copy(T->kernel[j]) : expr_copy(T->arg[j]);
                        pf[np++] = expr_new_function(expr_new_symbol("Power"),
                            (Expr*[]){ b, expr_new_integer(mult[k]) }, 2);
                    }
                    /* the integrand's own spelling of E^w: when the class primitive
                     * was re-spelled (E^(2 u)/(1 + E^(2 u)) for E^u/(E^-u + E^u)),
                     * the collected exponent is a multiplicity-1 member of level
                     * target, and it is that spelling the rebuild must exclude and map */
                    Expr* kw = T->kernel[target];
                    for (size_t m = 0; m < T->nm; m++)
                        if (T->mprim[m] == target && T->mmult[m] == 1 && T->mval[m]
                            && rt_head_is(T->mval[m], "Power") && T->mval[m]->data.function.arg_count == 2
                            && expr_eq(T->mval[m]->data.function.args[0], T->t[target])) {   /* a simple t^1 member */
                            kw = NULL;
                            *kern = expr_new_function(expr_new_symbol("Power"),
                                (Expr*[]){ expr_new_symbol("E"), expr_copy(T->marg[m]) }, 2);
                            break;
                        }
                    if (kw) *kern = expr_copy(kw);
                    *val = np == 1 ? pf[0] : expr_new_function(expr_new_symbol("Times"), pf, np);
                    free(pf);
                    found = true;
                }
                if (c) expr_free(c);
            }
            free(kn); free(kd); free(mult);
            if (co) expr_free(co);
        }
        gens[ng] = T->Dcoef[i2]; gi[ng] = i2; ng++;      /* borrowed */
    }
    free(gens); free(gi);
    expr_free(vars);
    return found;
}

/* True when e carries a constant E^c with c not an integer (E^I, E^(1/2)):
 * a constant the symbol RT_EI_SYMBOL would be treated as independent of.
 * Constant means D[c] = 0, not only x-free: TrigToExp[E^(I u) Tan[1 + u]]
 * spells E^(-I) as E^(I u - I (1 + u)), and the tower made that a level of
 * derivative 0 -- an independent "constant" t = E^(-I) beside Z = E^I, whose
 * relation t Z = 1 went unseen: a false "no" (2026-11-20). */
static bool rt_has_nonint_econst(const Expr* e, Expr* x) {
    if (!e || e->type != EXPR_FUNCTION) return false;
    if (rt_head_is(e, "Power") && e->data.function.arg_count == 2
        && e->data.function.args[0]->type == EXPR_SYMBOL
        && e->data.function.args[0]->data.symbol.name == intern_symbol("E")) {
        Expr* w = e->data.function.args[1];
        bool cst = rt_free_of_x(w, x);
        if (!cst) {
            Expr* dw = rt_eval1("Together", rt_eval2("D", expr_copy(w), expr_copy(x)));
            cst = dw && rt_is_zero(dw);
            if (dw) expr_free(dw);
        }
        if (cst) {
            int64_t ga, gb;   /* a Gaussian integer maps to E^a Z^b (rt_subst_kernels) */
            if (!rt_gauss_const_exponent(w, &ga, &gb)) return true;
        }
    }
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (rt_has_nonint_econst(e->data.function.args[i], x)) return true;
    return false;
}

/* e with every subexpression equal to k replaced by a copy of v, without
 * evaluation (owned result). */
static Expr* rt_struct_subst(const Expr* e, const Expr* k, const Expr* v) {
    if (expr_eq((Expr*)e, (Expr*)k)) return expr_copy((Expr*)v);
    if (e->type != EXPR_FUNCTION) return expr_copy((Expr*)e);
    size_t n = e->data.function.arg_count;
    Expr** a = malloc((n ? n : 1) * sizeof(Expr*));
    for (size_t i = 0; i < n; i++) a[i] = rt_struct_subst(e->data.function.args[i], k, v);
    Expr* out = expr_new_function(rt_struct_subst(e->data.function.head, k, v), a, n);
    free(a);
    return out;
}

static bool rt_tower_build_elim(Expr* f, Expr* x, RtTower* T, size_t min_n, bool wide);
/* The integrator's tower: a dependent logarithm is eliminated only when its
 * constant holds on the wide sample set (Log[1 + I x] = Log[x - I] + I Pi/2
 * on the whole line).  A constant of one region only (Log[-x] = Log[x] + I Pi
 * for x > 0) is left dependent: an answer built on it fails the
 * differentiation check elsewhere (dependent-log hunt: True 78 -> 67, 58). */
bool rt_tower_build_min(Expr* f, Expr* x, RtTower* T, size_t min_n) {
    return rt_tower_build_elim(f, x, T, min_n, true);
}
/* The decision's tower (rt_decide_tower): also eliminates a region constant;
 * the certificate is then about the germ on the positive sample points. */
bool rt_tower_build_decide(Expr* f, Expr* x, RtTower* T, size_t min_n) {
    return rt_tower_build_elim(f, x, T, min_n, false);
}
static bool rt_tower_build_elim(Expr* f, Expr* x, RtTower* T, size_t min_n, bool wide) {
    bool ok = rt_tower_build_core(f, x, T, min_n, NULL, NULL, 0);
    if (!ok || !T->dependent) return ok;
    Expr* tc[8]; size_t ntc = 0;
    rt_tconst_atoms(f, x, tc, &ntc, 8);
    Expr** xk = NULL; Expr** xr = NULL; size_t nx = 0;
    for (int round = 0; round < 8; round++) {
        Expr* k = NULL; Expr* v = NULL;
        if (!rt_dependent_log_rule(T, x, wide, &k, &v)
            && !rt_dependent_exp_rule(T, x, &k, &v)) break;
        size_t nt2 = ntc;
        rt_tconst_atoms(v, x, tc, &nt2, 8);
        /* Z = E^I next to a literal E^I (or E^(1/2)) of the integrand would be
         * treated as independent of a constant it is a power of */
        if (!rt_free_of_head(v, RT_EI_SYMBOL) && rt_has_nonint_econst(f, x)) { expr_free(k); expr_free(v); break; }
        if (getenv("MATHILDA_RDE_CERT_DEBUG")) { char* ks = expr_to_string(k); char* vs = expr_to_string(v);   /* trace */
            fprintf(stderr, "[rt-elim] %.200s -> %.200s (atoms %zu)\n", ks, vs, nt2); free(ks); free(vs); }
        if (nt2 > 1) { expr_free(k); expr_free(v); break; }
        xk = realloc(xk, (nx + 1) * sizeof(Expr*));
        xr = realloc(xr, (nx + 1) * sizeof(Expr*));
        /* compose: an earlier rule's value may hold the kernel eliminated now
         * (E^(1 + 3 I x) -> E/(E^(-I x))^3, then E^(-I x) -> E/E^(1 + I x)):
         * substitute it there too, structurally -- the values are built
         * unevaluated, and evaluating would merge them back into one E^w */
        for (size_t i = 0; i < nx; i++) {
            Expr* r = rt_struct_subst(xr[i], k, v);
            expr_free(xr[i]); xr[i] = r;
        }
        xk[nx] = k; xr[nx] = v; nx++;
        RtTower T2;
        bool ok2 = rt_tower_build_core(f, x, &T2, min_n, xk, xr, nx);
        if (!ok2 && getenv("MATHILDA_RDE_CERT_DEBUG")) fprintf(stderr, "[rt-elim] rebuild failed\n");   /* trace */
        if (!ok2) { rt_tower_free(&T2); nx--; expr_free(xk[nx]); expr_free(xr[nx]); break; }
        ntc = nt2;
        rt_tower_free(T);
        *T = T2;
        if (!T->dependent) break;
    }
    for (size_t i = 0; i < nx; i++) { expr_free(xk[i]); expr_free(xr[i]); }
    free(xk); free(xr);
    return true;
}
Expr* rt_tower_back_rules(const RtTower* T) {
    Expr** rs = malloc((T->n + 1) * sizeof(Expr*));
    for (size_t i = 0; i < T->n; i++)
        rs[i] = expr_new_function(expr_new_symbol("Rule"),
            (Expr*[]){ expr_copy(T->t[i]), expr_copy(T->kernel[i]) }, 2);
    rs[T->n] = expr_new_function(expr_new_symbol("Rule"), (Expr*[]){ expr_new_symbol(RT_EI_SYMBOL),
        expr_new_function(expr_new_symbol("Power"), (Expr*[]){ expr_new_symbol("E"),
            expr_new_function(expr_new_symbol("Complex"), (Expr*[]){ expr_new_integer(0), expr_new_integer(1) }, 2) }, 2) }, 2);
    Expr* l = expr_new_function(expr_new_symbol("List"), rs, T->n + 1);
    free(rs);
    return l;
}
/* Dt_i: the derivation of the tower variable t_i — Dcoef (log), Dcoef*t (exp),
 * or Dcoef*(t^2 + sigma) (tan, sigma = +-1).  Owned. */
Expr* rt_dt_i(RtTower* T, size_t i) {
    if (T->kind[i] == RT_LOG || T->kind[i] == RT_PRIM) return expr_copy(T->Dcoef[i]);
    if (T->kind[i] == RT_TAN)
        return expr_new_function(expr_new_symbol("Times"),
            (Expr*[]){ expr_copy(T->Dcoef[i]),
                expr_new_function(expr_new_symbol("Plus"),
                    (Expr*[]){ expr_new_function(expr_new_symbol("Power"),
                        (Expr*[]){ expr_copy(T->t[i]), expr_new_integer(2) }, 2),
                      expr_new_integer(T->tsg[i]) }, 2) }, 2);
    return expr_new_function(expr_new_symbol("Times"),   /* RT_EXP: Dcoef * t */
        (Expr*[]){ expr_copy(T->Dcoef[i]), expr_copy(T->t[i]) }, 2);
}
bool rt_tower_mobius(RtTower* T, size_t m, Expr* alpha, RtTower* T2,
                     Expr** to_theta, Expr** from_theta) {
    Expr* t = T->t[m];
    Expr* den = expr_new_function(expr_new_symbol("Plus"), (Expr*[]){
        expr_new_function(expr_new_symbol("Power"), (Expr*[]){ expr_copy(alpha), expr_new_integer(2) }, 2),
        expr_new_function(expr_new_symbol("Times"), (Expr*[]){ expr_new_integer(-1),
            expr_new_function(expr_new_symbol("Power"), (Expr*[]){ expr_copy(t), expr_new_integer(2) }, 2) }, 2) }, 2);
    Expr* a = rt_eval1("Together", rt_eval_own(expr_new_function(expr_new_symbol("Times"), (Expr*[]){
        expr_new_integer(2), expr_copy(alpha), rt_dt_i(T, m),
        expr_new_function(expr_new_symbol("Power"), (Expr*[]){ den, expr_new_integer(-1) }, 2) }, 4)));
    if (!a || !rt_free_of_x(a, t)) { if (a) expr_free(a); return false; }
    size_t n = m + 1;
    *T2 = *T;
    T2->n = n;
    T2->kind = malloc(n * sizeof(RtKind));
    T2->t = malloc(n * sizeof(Expr*));
    T2->Dcoef = malloc(n * sizeof(Expr*));
    T2->tsg = malloc(n * sizeof(long));
    T2->kernel = malloc(n * sizeof(Expr*));
    T2->arg = malloc(n * sizeof(Expr*));
    for (size_t i = 0; i < n; i++) {
        T2->kind[i] = T->kind[i]; T2->t[i] = T->t[i]; T2->Dcoef[i] = T->Dcoef[i];
        T2->tsg[i] = T->tsg[i]; T2->kernel[i] = T->kernel[i]; T2->arg[i] = T->arg[i];
    }
    char nm[48];
    snprintf(nm, sizeof nm, "Risch`Private`rmMob%zu", m);
    Expr* th = expr_new_symbol(nm);
    /* theta's exponent 2 v/alpha: 2 v for Tanh[v], -2 I v for Tan[v] */
    Expr* argm = rt_eval_own(expr_new_function(expr_new_symbol("Times"), (Expr*[]){
        expr_new_integer(2), expr_new_function(expr_new_symbol("Power"),
            (Expr*[]){ expr_copy(alpha), expr_new_integer(-1) }, 2), expr_copy(T->arg[m]) }, 3));
    T2->kind[m] = RT_EXP; T2->tsg[m] = 0; T2->t[m] = th; T2->Dcoef[m] = a; T2->arg[m] = argm;
    T2->kernel[m] = expr_new_function(expr_new_symbol("Power"),
        (Expr*[]){ expr_new_symbol("E"), expr_copy(argm) }, 2);
    Expr* tin = expr_new_function(expr_new_symbol("Times"), (Expr*[]){ expr_copy(alpha),
        expr_new_function(expr_new_symbol("Plus"), (Expr*[]){ expr_copy(th), expr_new_integer(-1) }, 2),
        expr_new_function(expr_new_symbol("Power"), (Expr*[]){
            expr_new_function(expr_new_symbol("Plus"), (Expr*[]){ expr_copy(th), expr_new_integer(1) }, 2),
            expr_new_integer(-1) }, 2) }, 3);
    /* subrules: the tangent kernel now reads as a rational function of theta */
    Expr* sr = expr_copy(T->subrules);
    if (sr && sr->type == EXPR_FUNCTION) {
        for (size_t i = 0; i < sr->data.function.arg_count; i++) {
            Expr* r = sr->data.function.args[i];
            if (r && r->type == EXPR_FUNCTION && r->data.function.arg_count == 2
                && expr_eq(r->data.function.args[1], t)) {
                Expr* nr = expr_new_function(expr_new_symbol("Rule"),
                    (Expr*[]){ expr_copy(r->data.function.args[0]), expr_copy(tin) }, 2);
                expr_free(r);
                sr->data.function.args[i] = nr;
            }
        }
    }
    T2->subrules = sr;
    *to_theta = expr_new_function(expr_new_symbol("Rule"), (Expr*[]){ expr_copy(t), tin }, 2);
    *from_theta = expr_new_function(expr_new_symbol("Rule"), (Expr*[]){ expr_copy(th),
        expr_new_function(expr_new_symbol("Times"), (Expr*[]){
            expr_new_function(expr_new_symbol("Plus"), (Expr*[]){ expr_copy(alpha), expr_copy(t) }, 2),
            expr_new_function(expr_new_symbol("Power"), (Expr*[]){
                expr_new_function(expr_new_symbol("Plus"), (Expr*[]){ expr_copy(alpha),
                    expr_new_function(expr_new_symbol("Times"), (Expr*[]){ expr_new_integer(-1), expr_copy(t) }, 2) }, 2),
                expr_new_integer(-1) }, 2) }, 2) }, 2);
    return true;
}
void rt_tower_mobius_free(RtTower* T2) {
    size_t m = T2->n - 1;
    expr_free(T2->kernel[m]); expr_free(T2->arg[m]); expr_free(T2->Dcoef[m]); expr_free(T2->t[m]);
    if (T2->subrules) expr_free(T2->subrules);
    free(T2->kind); free(T2->t); free(T2->Dcoef); free(T2->tsg); free(T2->kernel); free(T2->arg);
    T2->n = 0; T2->subrules = NULL;
}
/* Tower derivation D_tower[e] = D[e,x] + sum_i Dt_i D[e,t_i], with Dt_i from
 * rt_dt_i (log / exp / tan).  Owned result. */
Expr* rt_tower_deriv(Expr* e, RtTower* T, Expr* x) {
    Expr** terms = malloc((T->n + 1) * sizeof(Expr*));
    terms[0] = rt_eval2("D", expr_copy(e), expr_copy(x));
    for (size_t i = 0; i < T->n; i++) {
        Expr* dei = rt_eval2("D", expr_copy(e), expr_copy(T->t[i]));
        Expr* dti = rt_dt_i(T, i);
        terms[i + 1] = expr_new_function(expr_new_symbol("Times"),
            (Expr*[]){ dti, dei }, 2);
    }
    Expr* r = rt_eval_own(expr_new_function(expr_new_symbol("Plus"), terms, T->n + 1));
    free(terms);
    return r;
}
/* Build the tower derivation as a List[Rule[var, Dvar], ...] for risch_field /
 * risch_hermite: x -> 1, and each tower variable t_i -> Dt_i (Dcoef_i for a log,
 * Dcoef_i * t_i for an exp).  Owned; matches rt_tower_deriv's D_tower exactly. */
Expr* rt_build_deriv_rules(RtTower* T, Expr* x) {
    size_t n = T->n;
    Expr** rules = malloc((n + 1) * sizeof(Expr*));
    rules[0] = expr_new_function(expr_new_symbol("Rule"),
        (Expr*[]){ expr_copy(x), expr_new_integer(1) }, 2);
    for (size_t i = 0; i < n; i++) {
        Expr* dti = rt_dt_i(T, i);
        rules[i + 1] = expr_new_function(expr_new_symbol("Rule"),
            (Expr*[]){ expr_copy(T->t[i]), dti }, 2);
    }
    Expr* r = expr_new_function(expr_new_symbol("List"), rules, n + 1);
    free(rules);
    return r;
}
/* True if `e` contains an exponential (Exp[...] / E^...), logarithm (Log[...]), or
 * a Knowles primitive SF kernel (Erf/Erfi/Erfc/ExpIntegralEi/LogIntegral) anywhere.
 * Used to tell a "base-field" exponent term (rational in x, no nested transcendental
 * — e.g. 1, 1/x) from a term carrying a nested kernel (x E^(1+1/x), or the Erf-bearing
 * -Erf[x]^2 term of exp(-x^2 - Erf[x]^2)) when deciding how to split a merged
 * exponential.  Including the primitive kernels lets rt_expand_exp_sums split
 * E^(-x^2 - Erf[x]^2) -> E^(-x^2) E^(-Erf[x]^2) so each factor is a proper tower
 * monomial (KNOWLES_DESIGN.md §2.1).  Elementary integrands carry no primitive
 * kernel, so this predicate is unchanged for them (byte-identical towers). */
bool rt_has_explog_kernel(Expr* e) {
    if (!e || e->type != EXPR_FUNCTION) return false;
    if (e->data.function.head->type == EXPR_SYMBOL) {
        const char* h = e->data.function.head->data.symbol.name;
        if (h == intern_symbol("Exp") || h == intern_symbol("Log")) return true;
        if (rt_is_primitive_head(h) && e->data.function.arg_count == 1) return true;
        if (h == intern_symbol("Power") && e->data.function.arg_count == 2
            && e->data.function.args[0]->type == EXPR_SYMBOL
            && e->data.function.args[0]->data.symbol.name == intern_symbol("E"))
            return true;
    }
    if (rt_has_explog_kernel(e->data.function.head)) return true;
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (rt_has_explog_kernel(e->data.function.args[i])) return true;
    return false;
}
Expr* rt_expand_exp_sums(Expr* e) {
    if (!e) return NULL;
    if (e->type != EXPR_FUNCTION) return expr_copy(e);
    const char* h = (e->data.function.head->type == EXPR_SYMBOL)
        ? e->data.function.head->data.symbol.name : NULL;
    Expr* expo = NULL;
    if (h == intern_symbol("Power") && e->data.function.arg_count == 2
        && e->data.function.args[0]->type == EXPR_SYMBOL
        && e->data.function.args[0]->data.symbol.name == intern_symbol("E"))
        expo = e->data.function.args[1];
    else if (h == intern_symbol("Exp") && e->data.function.arg_count == 1)
        expo = e->data.function.args[0];
    if (expo && rt_head_is(expo, "Plus") && expo->data.function.arg_count >= 2) {
        /* Split E^(a1 + a2 + ...) only to ISOLATE nested-kernel terms (a term
         * carrying an Exp/Log, e.g. x E^(1+1/x)), which would otherwise strand a
         * foreign kernel in the tower derivation.  All base-field terms (rational
         * in x — 1, 1/x, ...) are grouped into ONE factor E^(sum of them).  This
         * avoids peeling a bare constant E^c off a Plus like 1 + 1/x: E^1 · E^(1/x)
         * re-merges under evaluation to E^(1+1/x), so the tower kernel E^(1/x) would
         * never match its E^(1+1/x) occurrences (the In[14] nested-exp failure).
         * Keeping 1+1/x together gives the valid single kernel E^(1+1/x). */
        size_t m = expo->data.function.arg_count;
        Expr** base = malloc(m * sizeof(Expr*)); size_t nb = 0;
        Expr** facs = malloc((m + 1) * sizeof(Expr*)); size_t nf = 0;
        for (size_t i = 0; i < m; i++) {
            Expr* term = expo->data.function.args[i];
            if (rt_has_explog_kernel(term))
                facs[nf++] = expr_new_function(expr_new_symbol("Power"),
                    (Expr*[]){ expr_new_symbol("E"), rt_expand_exp_sums(term) }, 2);
            else
                base[nb++] = expr_copy(term);
        }
        if (nf == 0) {                 /* no nested term: splitting would only peel a
                                        * constant — keep the exponential whole. */
            for (size_t i = 0; i < nb; i++) expr_free(base[i]);
            free(base); free(facs);
        } else {
            if (nb > 0) {
                Expr* basesum = (nb == 1) ? base[0]
                    : expr_new_function(expr_new_symbol("Plus"), base, nb);
                facs[nf++] = expr_new_function(expr_new_symbol("Power"),
                    (Expr*[]){ expr_new_symbol("E"), basesum }, 2);
            }
            free(base);
            Expr* prod = expr_new_function(expr_new_symbol("Times"), facs, nf);
            free(facs);
            return prod;
        }
    }
    Expr* nh = rt_expand_exp_sums(e->data.function.head);
    size_t k = e->data.function.arg_count;
    Expr** na = malloc((k ? k : 1) * sizeof(Expr*));
    for (size_t i = 0; i < k; i++)
        na[i] = rt_expand_exp_sums(e->data.function.args[i]);
    Expr* r = expr_new_function(nh, na, k);
    free(na);
    return r;
}
/* Structural (non-evaluating) kernel substitution: replace each tower-kernel
 * subtree (E^arg_i / Log[arg_i], plus the Exp[arg_i] spelling) by its tower
 * variable t_i, top-down (a matched node is not descended into).  Unlike an
 * evaluated ReplaceAll this never lets the evaluator re-merge a split
 * exponential product (E^x E^(E^x) -> E^(x+E^x)) before the kernels are
 * aliased; because the kernels were collected from the very tree being
 * substituted, structural equality holds by construction.  Returns a
 * freshly-owned tree (caller frees). */
/* When w (borrowed) is a constant exponent -- x-free, or a sum whose
 * derivative is 0 (TrigToExp spells E^(-I) as E^(I u - I (1 + u))) -- equal to
 * a Gaussian integer a + b i, set *a, *b and return true.  Expand decides it:
 * a constant sum such as I u - I (1 + u) expands to its value. */
static bool rt_gauss_const_exponent(Expr* w, int64_t* a, int64_t* b) {
    Expr* we = rt_eval1("Expand", expr_copy(w));
    if (!we) return false;
    bool ok = false;
    if (we->type == EXPR_INTEGER) { *a = we->data.integer; *b = 0; ok = true; }
    else if (rt_head_is(we, "Complex") && we->data.function.arg_count == 2
             && we->data.function.args[0]->type == EXPR_INTEGER
             && we->data.function.args[1]->type == EXPR_INTEGER) {
        *a = we->data.function.args[0]->data.integer;
        *b = we->data.function.args[1]->data.integer;
        ok = true;
    }
    expr_free(we);
    return ok;
}

Expr* rt_subst_kernels(Expr* e, RtTower* T) {
    if (!e) return NULL;
    /* a constant E^(a + b i), however spelled, is E^a Z^b with Z the symbol for
     * E^I (RT_EI_SYMBOL): every spelling of the same constant becomes the same
     * polynomial in E and Z, so a relation between them is an identity the
     * arithmetic sees.  (Left alone, E^(I u - I (1 + u)) was a level of
     * derivative 0 -- an independent "constant" beside Z = E^I.) */
    if (e->type == EXPR_FUNCTION && rt_head_is(e, "Power") && e->data.function.arg_count == 2
        && e->data.function.args[0]->type == EXPR_SYMBOL
        && e->data.function.args[0]->data.symbol.name == intern_symbol("E")) {
        Expr* w = e->data.function.args[1];
        int64_t ga, gb;
        if ((w->type == EXPR_FUNCTION && (rt_head_is(w, "Complex") || rt_head_is(w, "Plus")))
            && rt_gauss_const_exponent(w, &ga, &gb) && gb != 0) {
            Expr* zb = expr_new_function(expr_new_symbol("Power"),
                (Expr*[]){ expr_new_symbol(RT_EI_SYMBOL), expr_new_integer(gb) }, 2);
            if (ga == 0) return zb;
            return expr_new_function(expr_new_symbol("Times"), (Expr*[]){
                expr_new_function(expr_new_symbol("Power"),
                    (Expr*[]){ expr_new_symbol("E"), expr_new_integer(ga) }, 2), zb }, 2);
        }
    }
    for (size_t m = 0; m < T->nlm; m++)               /* dependent Log member */
        if (expr_eq(e, T->lmk[m])) return expr_copy(T->lmv[m]);
    for (size_t i = 0; i < T->n; i++) {
        if (expr_eq(e, T->kernel[i])) return expr_copy(T->t[i]);
        if (T->kind[i] == RT_EXP
            && e->type == EXPR_FUNCTION
            && e->data.function.head->type == EXPR_SYMBOL
            && e->data.function.head->data.symbol.name == intern_symbol("Exp")
            && e->data.function.arg_count == 1
            && expr_eq(e->data.function.args[0], T->arg[i]))
            return expr_copy(T->t[i]);
        /* E^(q a) = t^q for any INTEGER q, a the level's exponent: the build
         * records only the member exponents it collected, so an E^(4 u) that
         * first appears in a derivative after expansion ((1 + E^(2 u))^2) was a
         * foreign kernel and the tower was refused (the Tanh hunt, through
         * TrigToExp, 2026-11-19) */
        if (T->kind[i] == RT_EXP && T->arg[i] && e->type == EXPR_FUNCTION
            && e->data.function.head->type == EXPR_SYMBOL
            && ((e->data.function.head->data.symbol.name == intern_symbol("Power")
                 && e->data.function.arg_count == 2
                 && e->data.function.args[0]->type == EXPR_SYMBOL
                 && e->data.function.args[0]->data.symbol.name == intern_symbol("E"))
                || (e->data.function.head->data.symbol.name == intern_symbol("Exp")
                    && e->data.function.arg_count == 1))) {
            Expr* w = e->data.function.arg_count == 2 ? e->data.function.args[1] : e->data.function.args[0];
            if (!rt_free_of_x(w, T->arg[i]) || T->arg[i]->type != EXPR_SYMBOL) {
                Expr* q = rt_eval1("Cancel", expr_new_function(expr_new_symbol("Times"),
                    (Expr*[]){ expr_copy(w), expr_new_function(expr_new_symbol("Power"),
                        (Expr*[]){ expr_copy(T->arg[i]), expr_new_integer(-1) }, 2) }, 2));
                bool integer = q && q->type == EXPR_INTEGER && q->data.integer != 0;
                if (integer) {
                    Expr* r = expr_new_function(expr_new_symbol("Power"),
                        (Expr*[]){ expr_copy(T->t[i]), q }, 2);
                    return r;
                }
                if (q) expr_free(q);
            }
        }
        /* RT_TAN: rationalise the OTHER circular/hyperbolic trig of the tangent
         * argument to the tower variable (the fresh symbol keeps the evaluator from
         * canonicalising the result back to Csc/Sec).  Tan/Tanh already matched the
         * kernel above; here Sin=t/rad, Cos=1/rad, Sec=rad, Csc=rad/t, Cot=1/t with
         * rad = Sqrt[1 + sigma t^2]. */
        if (T->kind[i] == RT_TAN && e->type == EXPR_FUNCTION
            && e->data.function.head->type == EXPR_SYMBOL
            && e->data.function.arg_count == 1
            && expr_eq(e->data.function.args[0], T->arg[i])) {
            const char* h = e->data.function.head->data.symbol.name;
            long sg = T->tsg[i]; Expr* ti = T->t[i];
            Expr* rad = expr_new_function(expr_new_symbol("Power"),
                (Expr*[]){ expr_new_function(expr_new_symbol("Plus"),
                    (Expr*[]){ expr_new_integer(1), expr_new_function(expr_new_symbol("Times"),
                        (Expr*[]){ expr_new_integer(sg), expr_new_function(expr_new_symbol("Power"),
                            (Expr*[]){ expr_copy(ti), expr_new_integer(2) }, 2) }, 2) }, 2),
                  expr_new_function(expr_new_symbol("Rational"),
                    (Expr*[]){ expr_new_integer(1), expr_new_integer(2) }, 2) }, 2);
            if (h == intern_symbol(sg > 0 ? "Sin" : "Sinh"))
                return expr_new_function(expr_new_symbol("Times"),
                    (Expr*[]){ expr_copy(ti), expr_new_function(expr_new_symbol("Power"),
                        (Expr*[]){ rad, expr_new_integer(-1) }, 2) }, 2);
            if (h == intern_symbol(sg > 0 ? "Cos" : "Cosh"))
                return expr_new_function(expr_new_symbol("Power"),
                    (Expr*[]){ rad, expr_new_integer(-1) }, 2);
            if (h == intern_symbol(sg > 0 ? "Sec" : "Sech")) return rad;
            if (h == intern_symbol(sg > 0 ? "Csc" : "Csch"))
                return expr_new_function(expr_new_symbol("Times"),
                    (Expr*[]){ rad, expr_new_function(expr_new_symbol("Power"),
                        (Expr*[]){ expr_copy(ti), expr_new_integer(-1) }, 2) }, 2);
            if (h == intern_symbol(sg > 0 ? "Cot" : "Coth")) {
                expr_free(rad);
                return expr_new_function(expr_new_symbol("Power"),
                    (Expr*[]){ expr_copy(ti), expr_new_integer(-1) }, 2);
            }
            expr_free(rad);
        }
    }
    /* Multiplicatively commensurate member kernel E^(marg) = t[mprim]^mmult. */
    for (size_t m = 0; m < T->nm; m++) {
        if (e->type == EXPR_FUNCTION && e->data.function.arg_count >= 1
            && e->data.function.head->type == EXPR_SYMBOL) {
            const char* h = e->data.function.head->data.symbol.name;
            Expr* w = NULL;
            if (h == intern_symbol("Exp") && e->data.function.arg_count == 1)
                w = e->data.function.args[0];
            else if (h == intern_symbol("Power") && e->data.function.arg_count == 2
                     && e->data.function.args[0]->type == EXPR_SYMBOL
                     && e->data.function.args[0]->data.symbol.name == intern_symbol("E"))
                w = e->data.function.args[1];
            if (w && expr_eq(w, T->marg[m]))
                return expr_copy(T->mval[m]);
        }
    }
    if (e->type != EXPR_FUNCTION) return expr_copy(e);
    Expr* nh = rt_subst_kernels(e->data.function.head, T);
    size_t k = e->data.function.arg_count;
    Expr** na = malloc((k ? k : 1) * sizeof(Expr*));
    for (size_t i = 0; i < k; i++)
        na[i] = rt_subst_kernels(e->data.function.args[i], T);
    Expr* r = expr_new_function(nh, na, k);
    free(na);
    return r;
}
