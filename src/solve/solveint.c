/*
 * solveint.c
 *
 * Integer-domain (Diophantine) solving pre-pass for
 * Solve[eqns && constraints, vars, Integers].  See solveint.h for the
 * contract.
 *
 * Phase 1 engine:
 *   Stage A  separate equations from inequality / ordering / disequation
 *            constraints; convert each equation residual to a sparse
 *            integer MPoly (denominators cleared).
 *   Stage B  derive a finite integer box [lo_i, hi_i] per variable by a
 *            fixpoint of: explicit bounds, ordering propagation, and an
 *            interval-positivity rule (a sign-definite term is bounded by
 *            the rest of its (in)equality).  Decline if any variable stays
 *            unbounded (that is a later phase: linear-parametric / Pell /
 *            research-grade forms).
 *   Stage C  recursive elimination over the search variables with an exact
 *            univariate leaf (integer k-th root / quadratic discriminant /
 *            rational-root), every candidate re-verified against the
 *            original conjunction before it is emitted.
 *
 * Only necessary conditions are ever used to tighten a bound, so an
 * exhausted finite search that finds nothing returns the empty List {}
 * (a proof of no integer solutions); an input we cannot bound or evaluate
 * returns NULL (Solve stays unevaluated).
 */

#include "solveint.h"

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gmp.h>

#include "attr.h"
#include "eval.h"
#include "expr.h"
#include "internal.h"
#include "sym_names.h"
#include "symtab.h"
#include "checked_int.h"
#include "poly/mpoly.h"
#include "numbertheory/numbertheory_internal.h"
#include "linalg/hnf.h"
#include "solvethue.h"
#include "solveint_internal.h"


/* Dispatch the special forms.  Returns true if one handled the input
 * (candidates emitted into st). */
static bool si_try_special_forms(SICtx* c, SearchState* st) {
    if (si_solve_three_cubes_booker(c, st)) return true;
    if (si_solve_pell(c, st)) return true;
    if (si_solve_conic(c, st)) return true;
    if (si_solve_factorable_conic(c, st)) return true;
    if (si_solve_elliptic_bqf(c, st)) return true;
    if (si_solve_reciprocal(c, st)) return true;
    if (si_solve_linelim_bilinear(c, st)) return true;
    if (si_solve_biquadrate_frye(c, st)) return true;
    if (si_solve_separable_mitm(c, st)) return true;
    return false;
}


/* Fermat's Last Theorem short-circuit.  a*x^n + a*y^n - a*z^n == 0 (equal
 * coefficient magnitudes, exponent n >= 3, no constant term) with x, y, z all
 * strictly positive has NO solutions (Wiles 1995) -- return {} immediately
 * rather than search the box, and independently of any upper bound so the
 * unbounded "indefinite" form is decided too.  Returns {} on a match, else NULL
 * (fall through to the ordinary machinery). */
static Expr* si_solve_fermat(SICtx* c) {
    if (c->neq != 1 || c->n != 3) return NULL;
    const MPoly* eq = c->eq[0];
    int vexp[3] = {0, 0, 0}, vsgn[3] = {0, 0, 0}, seen[3] = {0, 0, 0};
    mpz_t mag; mpz_init(mag); bool have_mag = false, ok = true;
    for (size_t t = 0; t < eq->n_terms && ok; t++) {
        const int* ex = eq->exps + t * 3;
        int nz = -1, cnt = 0;
        for (int v = 0; v < 3; v++) if (ex[v] > 0) { nz = v; cnt++; }
        if (cnt == 0) { ok = false; break; }          /* constant term: k != 0 */
        if (cnt > 1 || seen[nz]) { ok = false; break; }
        int s = mpz_sgn(eq->coefs[t]);
        if (s == 0) { ok = false; break; }
        mpz_t a; mpz_init(a); mpz_abs(a, eq->coefs[t]);
        if (!have_mag) { mpz_set(mag, a); have_mag = true; }
        else if (mpz_cmp(mag, a) != 0) ok = false;
        mpz_clear(a);
        seen[nz] = 1; vexp[nz] = ex[nz]; vsgn[nz] = s;
    }
    mpz_clear(mag);
    if (!ok || !seen[0] || !seen[1] || !seen[2]) return NULL;
    int n = vexp[0];
    if (n < 3 || vexp[1] != n || vexp[2] != n) return NULL;   /* FLT is n >= 3 */
    int ssum = vsgn[0] + vsgn[1] + vsgn[2];
    if (ssum != 1 && ssum != -1) return NULL;         /* must be a^n + b^n == c^n */
    for (int i = 0; i < 3; i++)
        if (!(c->has_lo[i] && c->lo[i] >= 1)) return NULL;    /* strictly positive */
    return mk_list(NULL, 0);                          /* provably no solutions */
}

Expr* solveint_solve_integer(Expr* expr, Expr* vars, Expr* dom) {
    if (!expr || !vars || !dom) return NULL;
    if (!(dom->type == EXPR_SYMBOL && dom->data.symbol.name == SYM_Integers)) return NULL;

    /* Parse the variable list into a symbol array. */
    SICtx c;
    memset(&c, 0, sizeof(c));
    Expr* var_storage[SI_MAX_VARS];
    if (vars->type == EXPR_SYMBOL) {
        var_storage[0] = vars; c.n = 1;
    } else if (vars->type == EXPR_FUNCTION
        && vars->data.function.head->type == EXPR_SYMBOL
        && vars->data.function.head->data.symbol.name == SYM_List) {
        if (vars->data.function.arg_count < 1
            || vars->data.function.arg_count > SI_MAX_VARS) return NULL;
        c.n = (int)vars->data.function.arg_count;
        for (int i = 0; i < c.n; i++) {
            if (vars->data.function.args[i]->type != EXPR_SYMBOL) return NULL;
            var_storage[i] = vars->data.function.args[i];
        }
    } else return NULL;
    c.var = var_storage;
    c.original = expr;
    c.all_captured = true;      /* cleared by any constraint the store can't hold */

    /* This pre-pass engages when there is at least one inequality / ordering /
     * disequation constraint, OR the equation is multivariable (so the
     * parametric linear path can produce a solution family).  A bare
     * single-variable equation with no constraints is left to the ordinary
     * polynomial dispatch (which, with the Integers reality + integer filter,
     * already handles x^2 == 4). */
    Expr** conj; int ncj;
    flatten_conjuncts(expr, &conj, &ncj);
    si_warn_free_symbols(&c, conj, ncj);   /* Solve::svars for stray symbols */
    bool has_constraint = false, has_equation = false;
    for (int i = 0; i < ncj; i++) {
        Expr* e = conj[i];
        if (is_fun(e, SYM_Equal, 2)) has_equation = true;
        else if (is_fun(e, SYM_Less, 2) || is_fun(e, SYM_LessEqual, 2)
              || is_fun(e, SYM_Greater, 2) || is_fun(e, SYM_GreaterEqual, 2)
              || is_fun(e, SYM_Unequal, 2)
              || (e->type == EXPR_FUNCTION && e->data.function.head->type == EXPR_SYMBOL
                  && e->data.function.head->data.symbol.name == SYM_Inequality))
            has_constraint = true;
    }
    if (!has_equation || (!has_constraint && c.n < 2)) return NULL;

    /* Exponential Diophantine (variable exponents) is handled before the MPoly
     * stage, which cannot represent x^a. */
    {
        Expr* ex = si_solve_exponential(expr, c.var, c.n);
        if (ex) return ex;
    }

    /* Unbounded Ramanujan-Nagell-type  x^2 + D == 2^n  (also a variable-exponent
     * shape, so before the MPoly stage): the complete finite set via the
     * class-number-1 factorisation + BHV bound, or a decline outside the gate. */
    {
        Expr* rn = si_solve_ramanujan_nagell(expr, c.var, c.n);
        if (rn) return rn;
    }

    /* Non-polynomial bounded power-leaf (n! + 1 == m^2, ...): also before the
     * MPoly stage, which cannot represent Factorial / Binomial / ... */
    {
        Expr* pl = si_solve_bounded_powerleaf(expr, c.var, c.n);
        if (pl) return pl;
    }

    /* Stage A. */
    for (int i = 0; i < ncj; i++) {
        if (!classify_conjunct(&c, conj[i])) { ctx_free(&c); return NULL; }
    }
    if (c.neq == 0) { ctx_free(&c); return NULL; }

    /* Stage B. */
    derive_bounds(&c);
    derive_even_only_bounds(&c);   /* sign-symmetric even-power vars -> [-B, B] */

    /* Fermat's Last Theorem: x^n + y^n == z^n (n >= 3, all positive) has no
     * solutions -- decide it instantly, before any (possibly unbounded) search. */
    { Expr* flt = si_solve_fermat(&c); if (flt) { ctx_free(&c); return flt; } }

    /* Per-variable degree (max over equations) and whether it is solvable as
     * an exact leaf. */
    int maxdeg[SI_MAX_VARS];
    for (int i = 0; i < c.n; i++) {
        maxdeg[i] = 0;
        for (int q = 0; q < c.neq; q++) {
            int dg = mpoly_deg_var(c.eq[q], i);
            if (dg > maxdeg[i]) maxdeg[i] = dg;
        }
    }

    /* An unbounded variable is admissible ONLY as the leaf (it is solved
     * exactly, never enumerated).  Count them: two or more and the box is not
     * finite -> decline (later phases handle parametric / Pell). */
    int n_unbounded = 0, unbounded_var = -1;
    for (int i = 0; i < c.n; i++)
        if (!(c.has_lo[i] && c.has_hi[i])) { n_unbounded++; unbounded_var = i; }

    SearchState st; memset(&st, 0, sizeof(st));
    st.ctx = &c;

    /* Unconstrained single linear equation -> parametric family (symbolic,
     * so it bypasses the numeric candidate machinery). */
    {
        Expr* lin = si_solve_linear_parametric(&c);
        if (lin) { ctx_free(&c); return lin; }
    }

    /* Unbounded positive Pell -> parametric fundamental-unit family (symbolic). */
    {
        Expr* pell = si_solve_pell_parametric(&c);
        if (pell) { ctx_free(&c); return pell; }
    }

    /* Unbounded positive generalised Pell  x^2 - D y^2 == N  (|N| >= 2) -> one
     * parametric family per solution class (Nagell fundamentals + orbit). */
    {
        Expr* gp = si_solve_genpell_parametric(&c);
        if (gp) { ctx_free(&c); return gp; }
    }

    /* Homogeneous ternary quadratic  a x^2 + b y^2 + c z^2 == 0  (unbounded):
     * Legendre solvability (a proof) -> the complete 2-parameter integer family,
     * or the trivial-only {{x->0,y->0,z->0}} when no nontrivial solution exists. */
    {
        Expr* tq = si_solve_ternary_quadratic(&c);
        if (tq) { ctx_free(&c); return tq; }
    }

    /* Homogeneous linear system with positivity -> parametric ray (symbolic). */
    {
        Expr* ray = si_solve_linear_system_ray(&c);
        if (ray) { ctx_free(&c); return ray; }
    }

    /* General unconstrained linear system (m >= 2 equations, unbounded) -> the
     * complete integer solution via HNF: particular solution + kernel lattice
     * as C[k], or {} when provably unsolvable.  (Replaces the silent wrong `{}`
     * the Complexes-oriented linear-system dispatch produced for underdetermined
     * integer systems.) */
    {
        Expr* hs = si_solve_linear_system_hnf(&c);
        if (hs) { ctx_free(&c); return hs; }
    }

    /* Prouhet-Tarry-Escott equal-power-sum system with a forcing disequation
     * -> provably {} via Newton's identities (even when unbounded). */
    {
        Expr* pte = si_solve_power_sum_equal(&c);
        if (pte) { ctx_free(&c); return pte; }
    }

    /* Unbounded Mordell y^2 == x^3 + k (k = -1, -2): the complete integer-point
     * set via factorisation in Z[sqrt k]. */
    {
        Expr* mor = si_solve_mordell(&c);
        if (mor) { ctx_free(&c); return mor; }
    }

    /* Thue equation F(x,y) == m (irreducible homogeneous binary form,
     * deg >= 3): the complete finite solution set via Tzanakis-de Weger,
     * or a decline (the engine never returns an unproven set). */
    {
        Expr* th = si_solve_thue(&c);
        if (th) { ctx_free(&c); return th; }
    }

    /* Special forms first: divisor-factoring bilinear and unit-fraction
     * recursion are exact and O(#divisors) / O(bounded), so they beat the
     * enumerative fallback whenever they match -- including fully bounded
     * systems whose box would otherwise force a large leaf search (e.g. the
     * Pythagorean-perimeter case, bounded by its linear equation). */
    if (si_try_special_forms(&c, &st)) {
        if (st.overflow) { free(st.sols); ctx_free(&c); return NULL; }
        Expr* result = build_result(&st);
        free(st.sols); ctx_free(&c);
        return result;
    }
    /* Multi-leaf staged elimination: a system whose "determined" variables each
     * fall out of a single equation (e.g. the Euler brick's a,b,c) -- these may
     * be individually unbounded, so this runs BEFORE the unbounded-box decline.
     * Only the coupled free variables are enumerated. */
    if (si_solve_multileaf(&c, &st)) {
        Expr* result = build_result(&st);
        free(st.sols); ctx_free(&c);
        return result;
    }
    st.multileaf = false;   /* reset in case the attempt set up partial state */
    st.nsol = 0;

    /* Unbounded box and no special form fit -> later phases (Pell, lattice). */
    if (n_unbounded >= 2) { free(st.sols); ctx_free(&c); return NULL; }

    /* Meet-in-the-middle fast path: a single separable additive equation is
     * solved in ~N^ceil(n/2) instead of the ~N^(n-1) leaf search. */
    if (n_unbounded == 0 && mitm_solve(&st)) {
        Expr* result = build_result(&st);
        free(st.sols);
        ctx_free(&c);
        return result;
    }

    int leaf;
    if (n_unbounded == 1) {
        /* Must be the leaf; it has to appear in an equation at a solvable
         * degree, else we cannot pin it. */
        leaf = unbounded_var;
        if (maxdeg[leaf] <= 0 || maxdeg[leaf] > SI_LEAF_MAXDEG) { ctx_free(&c); return NULL; }
        /* Give the leaf a wide finite window for exact-root filtering; the
         * final verification enforces the true (possibly one-sided) bounds. */
        const int64_t SI_WIDE = 1LL << 50;
        if (!c.has_lo[leaf]) { c.lo[leaf] = -SI_WIDE; c.has_lo[leaf] = true; }
        if (!c.has_hi[leaf]) { c.hi[leaf] =  SI_WIDE; c.has_hi[leaf] = true; }
    } else {
        /* Fully bounded: leaf = widest-domain variable that appears in some
         * equation at a solvable degree (so the widest range is never
         * enumerated); ties broken toward the lower leaf degree. */
        leaf = -1; int64_t best = -1; int best_deg = 1 << 30;
        for (int i = 0; i < c.n; i++) {
            if (maxdeg[i] <= 0 || maxdeg[i] > SI_LEAF_MAXDEG) continue;
            int64_t w = c.hi[i] - c.lo[i];
            if (w > best || (w == best && maxdeg[i] < best_deg)) {
                best = w; best_deg = maxdeg[i]; leaf = i;
            }
        }
        if (leaf < 0) leaf = 0;                       /* degenerate: no eqn var */
    }

    for (int i = 0; i < c.n; i++)
        if (c.lo[i] > c.hi[i]) { ctx_free(&c); return mk_list(NULL, 0); }  /* empty box */

    st.leaf = leaf;
    int64_t domain[SI_MAX_VARS];
    for (int i = 0; i < c.n; i++) domain[i] = c.hi[i] - c.lo[i] + 1;
    st.n_search = 0;
    for (int i = 0; i < c.n; i++) if (i != leaf) st.order[st.n_search++] = i;
    for (int i = 0; i + 1 < st.n_search; i++) {          /* sort ascending domain */
        int pick = i;
        for (int j = i + 1; j < st.n_search; j++)
            if (domain[st.order[j]] < domain[st.order[pick]]) pick = j;
        if (pick != i) { int t = st.order[i]; st.order[i] = st.order[pick]; st.order[pick] = t; }
    }
    /* Refine so every abs-ordering |a| < |b| enumerates the larger b before the
     * smaller a -- only then does effective_bounds' abs pruning fire and the box
     * collapse to ~N^L/L! nodes (matching si_longest_chain).  Stable topological
     * pass over the abs-ord DAG, seeded by the ascending-domain order above so
     * ties and unconstrained variables keep their heuristic placement. */
    if (c.n_abs_ord > 0 && st.n_search > 1) {
        int seed[SI_MAX_VARS];
        for (int i = 0; i < st.n_search; i++) seed[i] = st.order[i];
        bool placed[SI_MAX_VARS];
        for (int i = 0; i < c.n; i++) placed[i] = false;
        int out = 0;
        for (int pass = 0; pass < st.n_search && out < st.n_search; pass++) {
            for (int i = 0; i < st.n_search; i++) {
                int v = seed[i];
                if (placed[v]) continue;
                bool ready = true;                       /* all larger-|.| vars placed? */
                for (int e = 0; e < c.n_abs_ord; e++)
                    if (c.abs_ord_a[e] == v) {
                        int b = c.abs_ord_b[e];
                        if (b != leaf && !placed[b]) { ready = false; break; }
                    }
                if (ready) { st.order[out++] = v; placed[v] = true; }
            }
        }
        for (int i = 0; i < st.n_search && out < st.n_search; i++)   /* cycle leftover */
            if (!placed[seed[i]]) { st.order[out++] = seed[i]; placed[seed[i]] = true; }
    }

    /* Search-space guard: decline rather than enumerate an intractable box
     * (e.g. a Pell family or a wide linear lattice -- those are later,
     * closed-form phases).  `raw_est` is the box product over the search vars;
     * `est` divides it by the factorial of the longest ordering chain, since
     * x <= y <= z  walks ~N^3/6 nodes, not N^3 (the visit cap backstops any
     * under-estimate). */
    long double raw_est = 1.0L;
    for (int i = 0; i < st.n_search; i++) raw_est *= (long double)domain[st.order[i]];
    long double est = raw_est;
    {
        int L = si_longest_chain(&c, st.order, st.n_search);
        long double fact = 1.0L;
        for (int i = 2; i <= L; i++) fact *= (long double)i;
        est /= fact;
    }
    /* Two number-theoretic methods turn the O(N^k) search into far less: a
     * bounded linear equation via its LLL-reduced lattice, and a separable
     * odd-power sum via the divisor method (s = x+y divides m).  Gate them on
     * the RAW box size, not the ordering-pruned `est`: an ordering's L! division
     * must not drop a large separable box just under the leaf budget and into
     * the O(est) leaf enumeration -- that is exactly what made the ordered
     * three-cubes query hang.  They decline (return false) when inapplicable,
     * falling through to the leaf search only if the pruned box fits. */
    if (raw_est > (long double)SI_MAX_NODES) {
        bool done = si_solve_linear_bounded(&c, &st)
                 || si_solve_powersum_divisor(&c, &st);
        if (done && !st.overflow) {
            Expr* result = build_result(&st);
            free(st.sols); ctx_free(&c);
            return result;
        }
        if (st.overflow) { free(st.sols); ctx_free(&c); return NULL; }
        st.nsol = 0; st.multileaf = false;          /* reset any partial state */
        if (est > (long double)SI_MAX_NODES) {       /* pruned box still too big */
            /* Last resort: a modular-sieved exhaustive leaf search (prunes the
             * innermost variable to feasible residues mod M).  Returns the
             * complete set / proven {} if it finishes in the raised budget,
             * else declines. */
            if (si_solve_box_modsieve(&c, &st) && !st.overflow) {
                Expr* result = build_result(&st);
                free(st.sols); ctx_free(&c);
                return result;
            }
            st.nsol = 0; st.overflow = false;
            free(st.sols); ctx_free(&c); return NULL;
        }
        /* else: the ordering-pruned box fits the leaf budget -- fall through. */
    }

    /* Stage C. */
    st.max_visits = SI_MAX_NODES;      /* runtime backstop for non-ordered boxes */
    si_build_leaf_cache(&st);          /* int64 coefficient cache for the hot loop */
    search_rec(&st, 0);
    si_free_leaf_cache(&st);

    if (st.overflow) {                 /* hit the solver's degree/size limit */
        free(st.sols);
        ctx_free(&c);
        return NULL;                   /* Solve stays unevaluated -- never wrong */
    }

    Expr* result = build_result(&st);
    free(st.sols);
    ctx_free(&c);
    return result;
}


/* ------------------------------------------------------------------ *
 *  Qualified-builtin entry: Solve`SolveIntegers[eqns, vars]          *
 * ------------------------------------------------------------------ */
static Expr* builtin_solve_integers(Expr* res) {
    if (!res || res->type != EXPR_FUNCTION) return NULL;
    if (res->data.function.arg_count != 2) return NULL;
    Expr* eqns = res->data.function.args[0];
    Expr* vars = res->data.function.args[1];
    Expr* integers = mk_sym(SYM_Integers);
    Expr* out = solveint_solve_integer(eqns, vars, integers);
    expr_free(integers);
    return out;
}

void solveint_init(void) {
    symtab_add_builtin("Solve`SolveIntegers", builtin_solve_integers);
    symtab_add_builtin("Solve`CubeRootsMod", builtin_cube_roots_mod);
    symtab_set_docstring("Solve`CubeRootsMod",
        "Solve`CubeRootsMod[k, d]\n"
        "\tInternal: the sorted list of all r in [0, d) with r^3 == k (mod d),\n"
        "\tvia factorisation + per-prime-power roots + CRT.  Backs the Booker\n"
        "\tcube-root-mod-d path for x^3 + y^3 + z^3 == k; exposed for testing.");
    symtab_set_docstring("Solve`SolveIntegers",
        "Solve`SolveIntegers[eqns, vars]\n"
        "\tInternal: solves a system of polynomial equations with\n"
        "\tinequality / ordering constraints over the integers by bound\n"
        "\tpropagation and exhaustive elimination.  Returns\n"
        "\t{{v -> n, ...}, ...} ascending, {} when there are provably no\n"
        "\tsolutions, or is left unevaluated when the variables cannot be\n"
        "\tbounded to a finite box.");
}
