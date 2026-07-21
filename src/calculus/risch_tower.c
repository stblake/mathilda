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
void rt_tower_free(RtTower* T) {
    for (size_t i = 0; i < T->n; i++) {
        if (T->kernel && T->kernel[i]) expr_free(T->kernel[i]);
        if (T->arg && T->arg[i]) expr_free(T->arg[i]);
        if (T->t && T->t[i]) expr_free(T->t[i]);
        if (T->Dcoef && T->Dcoef[i]) expr_free(T->Dcoef[i]);
    }
    if (T->subrules) expr_free(T->subrules);
    for (size_t i = 0; i < T->nm; i++)
        if (T->marg && T->marg[i]) expr_free(T->marg[i]);
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
bool rt_tower_build_min(Expr* f, Expr* x, RtTower* T, size_t min_n) {
    T->kind = NULL; T->kernel = NULL; T->arg = NULL; T->t = NULL; T->Dcoef = NULL;
    T->tsg = NULL; T->marg = NULL; T->mprim = NULL; T->mmult = NULL;
    T->subrules = NULL; T->n = 0; T->nm = 0;

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
    size_t np = 0;
    for (size_t rep = 0; rep < ne; rep++) if (clsrep[rep] == (long)rep) np++;
    size_t n = nl + np + nt + npr;
    if (!okc || n < min_n) {                         /* no upper depth cap */
        for (size_t rep = 0; rep < ne; rep++) if (clsprim[rep]) expr_free(clsprim[rep]);
        free(clsprim); free(clsrep); free(multof); free(mprim_pexp);
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
        if (clsrep[rep] != (long)rep) continue;
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
        if (expr_eq(exps[i], clsprim[clsrep[i]])) continue;   /* realized by the primitive kernel */
        T->marg[T->nm] = expr_copy(exps[i]);                  /* member exponent w */
        T->mmult[T->nm] = multof[i];                          /* w = k * clsprim[rep]  */
        mprim_pexp[T->nm] = clsprim[clsrep[i]];               /* borrowed; remapped below */
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
    Expr** rules = malloc((6 * n + 2 * T->nm) * sizeof(Expr*)); size_t nr = 0;
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
        Expr* pw = mprim_pexp[m];                  /* class primitive exponent value */
        long ti = -1;
        for (size_t i = 0; i < n; i++)
            if (T->kind[i] == RT_EXP && expr_eq(T->arg[i], pw)) { ti = (long)i; break; }
        T->mprim[m] = ti;                          /* now a tower index (>= 0) */
        Expr* tk = expr_new_function(expr_new_symbol("Power"),
            (Expr*[]){ expr_copy(T->t[ti]), expr_new_integer(T->mmult[m]) }, 2);
        rules[nr++] = expr_new_function(expr_new_symbol("Rule"),
            (Expr*[]){ expr_new_function(expr_new_symbol("Exp"),
                (Expr*[]){ expr_copy(T->marg[m]) }, 1), expr_copy(tk) }, 2);
        rules[nr++] = expr_new_function(expr_new_symbol("Rule"),
            (Expr*[]){ expr_new_function(expr_new_symbol("Power"),
                (Expr*[]){ expr_new_symbol("E"), expr_copy(T->marg[m]) }, 2), tk }, 2);
    }
    T->subrules = expr_new_function(expr_new_symbol("List"), rules, nr);
    free(rules);
    for (size_t rep = 0; rep < ne; rep++) if (clsprim[rep]) expr_free(clsprim[rep]);
    free(clsprim); free(clsrep); free(multof); free(mprim_pexp);
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
        Expr* trig = parse_expression(
            "{Sec[rmTA_]^2 -> 1 + Tan[rmTA_]^2, Csc[rmTA_]^2 -> 1 + Cot[rmTA_]^2, "
            "Sech[rmTA_]^2 -> 1 - Tanh[rmTA_]^2, Csch[rmTA_]^2 -> Coth[rmTA_]^2 - 1}");
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
    }
    return ok;
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
Expr* rt_subst_kernels(Expr* e, RtTower* T) {
    if (!e) return NULL;
    for (size_t i = 0; i < T->n; i++) {
        if (expr_eq(e, T->kernel[i])) return expr_copy(T->t[i]);
        if (T->kind[i] == RT_EXP
            && e->type == EXPR_FUNCTION
            && e->data.function.head->type == EXPR_SYMBOL
            && e->data.function.head->data.symbol.name == intern_symbol("Exp")
            && e->data.function.arg_count == 1
            && expr_eq(e->data.function.args[0], T->arg[i]))
            return expr_copy(T->t[i]);
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
                return expr_new_function(expr_new_symbol("Power"),
                    (Expr*[]){ expr_copy(T->t[T->mprim[m]]),
                               expr_new_integer(T->mmult[m]) }, 2);
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
