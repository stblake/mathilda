/*
 * solvenlsys.c
 *
 * Nonlinear polynomial-system specialist for `Solve` (src/solve.c).
 * Also reachable directly as the context-qualified builtin
 * `Solve`SolveNonlinearSystem`.
 *
 * Algorithm
 * ---------
 * Each equation `lhs_i == rhs_i` is canonicalised to `poly_i := lhs_i -
 * rhs_i`.  Every `poly_i` must be a polynomial over Q in `vars`
 * (`gb_from_expr` returns NULL otherwise -- e.g. a transcendental head,
 * a radical, or a foreign symbol -- and this specialist declines).  We
 * then compute a lexicographic Gröbner basis of the ideal via the
 * Gröbner walk (`gb_groebner_walk`).  For a zero-dimensional ideal the
 * lex basis is triangular: some generator is univariate in the last
 * variable, another becomes univariate in the second-last once the last
 * is fixed, and so on.  We solve the univariate generator at each level
 * with the single-variable polynomial solver
 * (`solvepoly_solve_polynomial_equality`), branch over each root,
 * back-substitute, and recurse; at full depth we verify the completed
 * tuple against every original equation before accepting it.
 *
 * This mirrors the *numeric* triangular solver
 * `nsolve_system_eliminate` in src/numerical_roots/nsolve_system.c,
 * substituting the exact univariate solver for NRoots and a symbolic
 * zero-test (`zero_test_decide`) for the numeric residual tolerance.
 *
 * Correctness policy (the crucial divergence from the numeric solver):
 * we must never emit a false `{}`.  A branch on which no triangular
 * generator can be found, or on which the univariate solver cannot
 * reduce the generator, marks the whole search `incomplete` and the
 * specialist returns NULL (leave Solve unevaluated) rather than an
 * empty list.  `{}` is returned only when the ideal is provably
 * inconsistent (unit ideal / a non-zero constant generator) or when a
 * fully-solved zero-dimensional ideal has no solutions in the requested
 * domain.
 *
 * Positive-dimensional ideals (infinitely many solutions) are detected
 * (not every variable owns a pure-power leading monomial) and left
 * unevaluated with the advisory `Solve::nsdim`.
 *
 * Memory: `equations`, `vars`, `dom`, and `opts` are borrowed.  Every
 * returned Expr* is freshly owned.  All intermediate GBPoly / Expr
 * allocations are released on every path.
 */

#include "solvenlsys.h"

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

#include "attr.h"
#include "eval.h"
#include "expr.h"
#include "groebner.h"
#include "internal.h"
#include "poly.h"
#include "solvepoly.h"
#include "sym_names.h"
#include "symtab.h"
#include "zero_test.h"

/* Per-generator total-degree gate: above this the lex Gröbner
 * computation could blow up, so the specialist declines (returns NULL).
 * Matches NSYS_MAX_TDEG in nsolve_system.c. */
#define NL_MAX_TDEG   60
/* Runaway backstop on the number of accepted tuples.  A zero-dimensional
 * ideal has finitely many solutions, but a pathological branch factor
 * could still explode; when this is exceeded we treat the search as
 * incomplete (return NULL) rather than emit a truncated set. */
#define NL_MAX_SOLS   4096
/* Maximum distinct free parameters before the RationalFunctions-coefficient
 * Gröbner computation is declined; its cost grows fast in the parameters. */
#define NL_MAX_PARAMS 3

/* ------------------------------------------------------------------ *
 *  Expression-building shorthand (same conventions as solvelinsys.c). *
 *  The mk_* helpers take ownership of their pointer arguments.        *
 * ------------------------------------------------------------------ */

static Expr* mk_int(int64_t v) { return expr_new_integer(v); }
static Expr* mk_sym(const char* name) { return expr_new_symbol(name); }

static Expr* mk_fn2(const char* head, Expr* a, Expr* b) {
    return expr_new_function(mk_sym(head), (Expr*[]){ a, b }, 2);
}
static Expr* mk_neg(Expr* e) { return mk_fn2("Times", mk_int(-1), e); }

static Expr* mk_list(Expr** args, size_t n) {
    return expr_new_function(mk_sym("List"), args, n);
}

/* ------------------------------------------------------------------ *
 *  Domain filters (identical to solvelinsys.c's private helpers).     *
 * ------------------------------------------------------------------ */

/* True iff `v` is provably a concrete integer. */
static bool is_concrete_integer(const Expr* v) {
    return v && v->type == EXPR_INTEGER;
}

/* True iff `v` syntactically contains a `Complex[_, _]` node anywhere. */
static bool contains_complex_head(const Expr* v) {
    if (!v || v->type != EXPR_FUNCTION) return false;
    if (v->data.function.head->type == EXPR_SYMBOL
        && v->data.function.head->data.symbol.name == SYM_Complex) {
        return true;
    }
    for (size_t i = 0; i < v->data.function.arg_count; i++) {
        if (contains_complex_head(v->data.function.args[i])) return true;
    }
    return false;
}

/* ------------------------------------------------------------------ *
 *  Total degree of a GBPoly (max exponent sum over its terms).        *
 * ------------------------------------------------------------------ */
static int gbpoly_total_degree(const GBPoly* p) {
    int best = 0;
    for (size_t t = 0; t < p->n_terms; t++) {
        const int* e = p->exps + t * (size_t)p->n_vars;
        int s = 0;
        for (int k = 0; k < p->n_vars; k++) s += e[k];
        if (s > best) best = s;
    }
    return best;
}

/* ------------------------------------------------------------------ *
 *  nsdim advisory (positive-dimensional ideal).                       *
 * ------------------------------------------------------------------ */
static void warn_nsdim(uint64_t input_hash) {
    static uint64_t last_warned_hash = 0;
    if (input_hash == last_warned_hash) return;
    last_warned_hash = input_hash;
    fprintf(stderr,
        "Solve::nsdim: The solution set is not zero-dimensional "
        "(infinitely many solutions); Solve returned unevaluated.\n");
}

/* ------------------------------------------------------------------ *
 *  Back-substitution recursion.                                       *
 * ------------------------------------------------------------------ */

typedef struct {
    Expr** orig; int norig;        /* numerator poly_i (= Numerator[lhs-rhs]) */
    Expr** dens;                   /* parallel to orig; non-constant denominator
                                    * or NULL (constant).  Borrowed.          */
    Expr** gb;   int ngb;          /* lex GB generators as Exprs, borrowed    */
    Expr** vars; int nvar;         /* variable symbols, borrowed              */
    Expr*  dom;                    /* borrowed; forwarded to solvepoly        */
    bool   reals_only, integers_only;
    const SolvePolyOpts* opts;     /* Cubics / Quartics passthrough           */
    Expr** sols; int nsol, solcap; /* accumulated solution tuples (owned)     */
    bool   incomplete;             /* a live branch could not be resolved     */
} NlCtx;

/* Verify a complete assignment (`assigned` = List of Rule[var, val]) by
 * substituting into every original equation and zero-testing the
 * residual.  Reject only on a PROVABLY non-zero residual; keep on TRUE
 * or UNKNOWN (favour recall -- never drop a real solution because the
 * zero-test was inconclusive). */
static bool verify_tuple(NlCtx* c, Expr* assigned) {
    for (int i = 0; i < c->norig; i++) {
        Expr* ra = expr_new_function(mk_sym("ReplaceAll"),
                       (Expr*[]){ expr_copy(c->orig[i]), expr_copy(assigned) }, 2);
        Expr* sub = eval_and_free(ra);
        ZeroTestResult z = zero_test_decide(sub);
        expr_free(sub);
        if (z == ZERO_TEST_FALSE) return false;
    }
    return true;
}

/* True iff no cleared denominator vanishes at this tuple.  A rational
 * equation lhs == rhs was reduced to Numerator == 0; a candidate that also
 * zeroes a denominator is a spurious root introduced by that reduction and
 * must be pruned.  Reject only on a PROVABLY zero denominator (favour recall:
 * an undecidable denominator keeps the tuple, matching verify_tuple). */
static bool tuple_denominators_ok(NlCtx* c, Expr* assigned) {
    for (int i = 0; i < c->norig; i++) {
        if (!c->dens[i]) continue;               /* constant denominator */
        Expr* ra = expr_new_function(mk_sym("ReplaceAll"),
                       (Expr*[]){ expr_copy(c->dens[i]), expr_copy(assigned) }, 2);
        Expr* sub = eval_and_free(ra);
        ZeroTestResult z = zero_test_decide(sub);
        expr_free(sub);
        if (z == ZERO_TEST_TRUE) return false;   /* denominator vanishes here */
    }
    return true;
}

/* True iff a tuple survives the domain filter. */
static bool tuple_in_domain(NlCtx* c, Expr* assigned) {
    if (!c->reals_only && !c->integers_only) return true;
    for (size_t a = 0; a < assigned->data.function.arg_count; a++) {
        Expr* val = assigned->data.function.args[a]->data.function.args[1];
        if (c->reals_only && contains_complex_head(val)) return false;
        if (c->integers_only && !is_concrete_integer(val)) return false;
    }
    return true;
}

/* Append a completed, verified tuple (in `vars` order) to c->sols,
 * skipping structural duplicates. */
static void emit_tuple(NlCtx* c, Expr* assigned) {
    Expr** rules = (Expr**)malloc(sizeof(Expr*) * (size_t)c->nvar);
    for (int i = 0; i < c->nvar; i++) {
        Expr* val = NULL;
        for (size_t a = 0; a < assigned->data.function.arg_count; a++) {
            Expr* rule = assigned->data.function.args[a];
            if (expr_eq(rule->data.function.args[0], c->vars[i])) {
                val = rule->data.function.args[1]; break;
            }
        }
        rules[i] = mk_fn2("Rule", expr_copy(c->vars[i]),
                          val ? expr_copy(val) : mk_int(0));
    }
    Expr* tuple = mk_list(rules, (size_t)c->nvar);
    free(rules);

    for (int s = 0; s < c->nsol; s++) {
        if (expr_eq(c->sols[s], tuple)) { expr_free(tuple); return; }
    }
    if (c->nsol == c->solcap) {
        c->solcap = c->solcap ? c->solcap * 2 : 8;
        c->sols = realloc(c->sols, sizeof(Expr*) * (size_t)c->solcap);
    }
    c->sols[c->nsol++] = tuple;
}

/* `assigned` is a List of Rule[var, val] for the last `depth` variables
 * (vars[nvar-depth .. nvar-1]). */
static void nl_rec(NlCtx* c, Expr* assigned, int depth) {
    if (c->incomplete) return;
    if (c->nsol > NL_MAX_SOLS) { c->incomplete = true; return; }

    if (depth == c->nvar) {
        if (!verify_tuple(c, assigned)) return;          /* spurious combination */
        if (!tuple_denominators_ok(c, assigned)) return; /* cleared-denominator locus */
        if (!tuple_in_domain(c, assigned)) return;       /* pruned by domain     */
        emit_tuple(c, assigned);
        return;
    }

    Expr* tvar = c->vars[c->nvar - 1 - depth];

    /* Lowest-degree generator that, after substituting the assigned
     * (higher) variables, is a non-constant polynomial in tvar ALONE --
     * free of the still-unassigned variables vars[0 .. nvar-2-depth].
     * (Without the freeness check, is_polynomial would treat a remaining
     * variable as a symbolic coefficient and solvepoly would choke.) */
    Expr* best = NULL; int bestdeg = 0;
    for (int g = 0; g < c->ngb; g++) {
        Expr* ra = expr_new_function(mk_sym("ReplaceAll"),
                       (Expr*[]){ expr_copy(c->gb[g]), expr_copy(assigned) }, 2);
        Expr* sub = eval_and_free(ra);
        bool usable = is_polynomial(sub, &tvar, 1)
                   && get_degree_poly(sub, tvar) >= 1;
        for (int j = 0; usable && j <= c->nvar - 2 - depth; j++)
            if (contains_any_symbol_from(sub, c->vars[j])) usable = false;
        if (usable) {
            int dg = get_degree_poly(sub, tvar);
            if (best == NULL || dg < bestdeg) {
                if (best) expr_free(best);
                best = sub; bestdeg = dg; continue;
            }
        }
        expr_free(sub);
    }
    if (!best) { c->incomplete = true; return; }   /* no triangular generator */

    /* Solve the univariate generator with the exact single-variable
     * solver, forwarding the domain and the Cubics/Quartics flags. */
    Expr* eqn = mk_fn2("Equal", best, mk_int(0));   /* consumes `best`        */
    Expr* roots = solvepoly_solve_polynomial_equality(eqn, tvar, c->dom, c->opts);
    expr_free(eqn);
    if (!roots) { c->incomplete = true; return; }   /* could not reduce       */

    /* roots = List[ List[Rule[tvar, val]], ... ] (possibly empty when no
     * roots exist in the domain -- a legitimate prune, not incomplete). */
    if (roots->type == EXPR_FUNCTION) {
        for (size_t r = 0; r < roots->data.function.arg_count && !c->incomplete; r++) {
            Expr* inner = roots->data.function.args[r];
            if (inner->type != EXPR_FUNCTION
                || inner->data.function.arg_count < 1) continue;
            Expr* rule = inner->data.function.args[0];
            if (rule->type != EXPR_FUNCTION || rule->data.function.arg_count != 2)
                continue;
            Expr* val = rule->data.function.args[1];

            size_t na = assigned->data.function.arg_count;
            Expr** na2 = (Expr**)malloc(sizeof(Expr*) * (na + 1));
            for (size_t a = 0; a < na; a++)
                na2[a] = expr_copy(assigned->data.function.args[a]);
            na2[na] = mk_fn2("Rule", expr_copy(tvar), expr_copy(val));
            Expr* na_list = mk_list(na2, na + 1);
            free(na2);
            nl_rec(c, na_list, depth + 1);
            expr_free(na_list);
        }
    }
    expr_free(roots);
}

/* ------------------------------------------------------------------ *
 *  Shared search driver (used by the polynomial and parametric paths). *
 * ------------------------------------------------------------------ */

/* Run the triangular back-substitution over the lex basis `gbE` and assemble
 * the solution list.  `orig`, `dens`, and `gbE` are all borrowed (nl_rec copies
 * every generator it uses).  Returns a freshly-owned List of solution tuples,
 * or NULL when a live branch could not be resolved (incomplete) -- which must
 * never be presented as an empty `{}`. */
static Expr* run_nl_search(Expr** orig, Expr** dens, int norig,
                           Expr** gbE, int ngb, Expr** var_arr, int nvar,
                           Expr* dom, bool reals_only, bool integers_only,
                           const SolvePolyOpts* opts) {
    NlCtx c;
    c.orig = orig; c.norig = norig;
    c.dens = dens;
    c.gb = gbE; c.ngb = ngb;
    c.vars = var_arr; c.nvar = nvar;
    c.dom = dom; c.reals_only = reals_only; c.integers_only = integers_only;
    c.opts = opts;
    c.sols = NULL; c.nsol = 0; c.solcap = 0; c.incomplete = false;

    Expr* empty = mk_list(NULL, 0);
    nl_rec(&c, empty, 0);
    expr_free(empty);

    if (c.incomplete) {
        for (int i = 0; i < c.nsol; i++) expr_free(c.sols[i]);
        free(c.sols);
        return NULL;
    }
    Expr* out = mk_list(c.sols, (size_t)c.nsol);
    free(c.sols);
    return out;
}

/* ------------------------------------------------------------------ *
 *  Parametric systems (coefficients in Q(parameters)).                *
 * ------------------------------------------------------------------ */

static bool nl_is_solve_var(const char* name, Expr** var_arr, int n) {
    for (int j = 0; j < n; j++)
        if (var_arr[j]->data.symbol.name == name) return true;
    return false;
}

/* Numeric / structural constants that are NOT parameters (mirrors the
 * GroebnerBasis builtin's is_known_constant_symbol so the two agree on which
 * symbols localise into the coefficient field). */
static bool nl_is_known_constant(const char* s) {
    return s == SYM_Pi || s == SYM_E || s == SYM_EulerGamma
        || s == SYM_Catalan || s == SYM_Degree || s == SYM_True
        || s == SYM_False || s == SYM_Infinity || s == SYM_DirectedInfinity
        || s == SYM_Null;
}

/* Collect distinct free-parameter symbol names (interned pointers) from the
 * OPERANDS of `e` into params[0..*np-1], capped at `cap`.  Function heads are
 * structural (Plus/Times/Power/...) and never parameters, so they are skipped.
 * Returns false once more than `cap` distinct parameters are seen. */
static bool nl_collect_params(const Expr* e, Expr** var_arr, int n,
                              const char** params, int* np, int cap) {
    if (!e) return true;
    if (e->type == EXPR_SYMBOL) {
        const char* s = e->data.symbol.name;
        if (nl_is_solve_var(s, var_arr, n) || nl_is_known_constant(s)) return true;
        for (int i = 0; i < *np; i++) if (params[i] == s) return true;
        if (*np >= cap) return false;
        params[(*np)++] = s;
        return true;
    }
    if (e->type == EXPR_FUNCTION) {
        for (size_t i = 0; i < e->data.function.arg_count; i++)
            if (!nl_collect_params(e->data.function.args[i], var_arr, n,
                                   params, np, cap))
                return false;
    }
    return true;
}

/* True iff `g` is free of every solve variable. */
static bool nl_var_free(Expr* g, Expr** var_arr, int n) {
    for (int v = 0; v < n; v++)
        if (contains_any_symbol_from(g, var_arr[v])) return false;
    return true;
}

/* Parametric branch: some equation carries a free parameter.  gb_from_expr and
 * the machine Gröbner walk carry only rational coefficients, so we compute the
 * lex Gröbner basis over Q(parameters) by reusing the GroebnerBasis builtin's
 * RationalFunctions-coefficient engine (which localises the parameters into the
 * coefficient field), then back-substitute with the symbolic univariate solver.
 * `orig` / `dens` are borrowed.  Returns owned result, or NULL (unevaluated). */
static Expr* solve_parametric_system(Expr** orig, Expr** dens, int norig,
                                     Expr** var_arr, int n, Expr* vars,
                                     Expr* dom, bool reals_only,
                                     bool integers_only,
                                     const SolvePolyOpts* opts) {
    /* GroebnerBasis[{orig...}, vars, CoefficientDomain -> RationalFunctions]. */
    Expr** polys = (Expr**)malloc(sizeof(Expr*) * (size_t)norig);
    for (int i = 0; i < norig; i++) polys[i] = expr_copy(orig[i]);
    Expr* poly_list = mk_list(polys, (size_t)norig);
    free(polys);
    Expr* cd_rule = mk_fn2("Rule", mk_sym("CoefficientDomain"),
                                   mk_sym("RationalFunctions"));
    Expr* gb_call = expr_new_function(mk_sym("GroebnerBasis"),
                        (Expr*[]){ poly_list, expr_copy(vars), cd_rule }, 3);
    Expr* gb_res = eval_and_free(gb_call);

    /* GroebnerBasis declines (non-polynomial substructure, unsupported input)
     * -> result is not a non-empty List -> leave Solve unevaluated. */
    if (!gb_res || gb_res->type != EXPR_FUNCTION
        || gb_res->data.function.head->type != EXPR_SYMBOL
        || gb_res->data.function.head->data.symbol.name != SYM_List
        || gb_res->data.function.arg_count == 0) {
        if (gb_res) expr_free(gb_res);
        return NULL;
    }

    Expr** gbE = gb_res->data.function.args;
    size_t ng = gb_res->data.function.arg_count;

    /* A generator free of every solve variable is a constant of the coefficient
     * field Q(params).  A nonzero NUMBER means the ideal is <1> -> the system is
     * (generically) inconsistent -> {}.  A parameter-dependent one is a
     * consistency / nondegeneracy condition (holds only on a subvariety of
     * parameter space); v1 declines rather than guess (that is Reduce's job). */
    bool inconsistent = false;
    for (size_t i = 0; i < ng; i++) {
        if (!nl_var_free(gbE[i], var_arr, n)) continue;
        Expr* g = gbE[i];
        bool is_number = g->type == EXPR_INTEGER || g->type == EXPR_REAL
            || g->type == EXPR_BIGINT
            || (g->type == EXPR_FUNCTION
                && g->data.function.head->type == EXPR_SYMBOL
                && g->data.function.head->data.symbol.name == SYM_Rational);
        if (is_number) { inconsistent = true; }
        else { expr_free(gb_res); return NULL; }   /* parameter condition */
    }
    if (inconsistent) { expr_free(gb_res); return mk_list(NULL, 0); }

    /* Back-substitute over the parametric triangular basis.  nl_rec treats the
     * parameters as symbolic coefficients (is_polynomial accepts them); its
     * `incomplete` flag guards systems that are positive-dimensional in the
     * variables, so a false `{}` can never be emitted. */
    Expr* out = run_nl_search(orig, dens, norig, gbE, (int)ng, var_arr, n,
                              dom, reals_only, integers_only, opts);
    expr_free(gb_res);
    return out;
}

/* ------------------------------------------------------------------ *
 *  Public entry.                                                      *
 * ------------------------------------------------------------------ */

Expr* solvenlsys_solve_nonlinear_system(Expr* equations,
                                        Expr* vars,
                                        Expr* dom,
                                        const SolvePolyOpts* opts) {
    if (!equations || !vars) return NULL;

    /* `vars` must be a List of >= 1 distinct symbols. */
    if (vars->type != EXPR_FUNCTION
        || vars->data.function.head->type != EXPR_SYMBOL
        || vars->data.function.head->data.symbol.name != SYM_List) {
        return NULL;
    }
    int n = (int)vars->data.function.arg_count;
    if (n < 1) return NULL;
    Expr** var_arr = vars->data.function.args;
    for (int j = 0; j < n; j++)
        if (var_arr[j]->type != EXPR_SYMBOL) return NULL;
    for (int a = 0; a < n; a++)
        for (int b = a + 1; b < n; b++)
            if (var_arr[a]->data.symbol.name == var_arr[b]->data.symbol.name) return NULL;

    /* `equations` is Equal[_, _], And[Equal, ...], or List[Equal, ...]. */
    Expr* single_holder[1];
    Expr** eq_arr;
    int m;
    if (equations->type == EXPR_FUNCTION
        && equations->data.function.head->type == EXPR_SYMBOL
        && (equations->data.function.head->data.symbol.name == SYM_List
            || equations->data.function.head->data.symbol.name == SYM_And)) {
        eq_arr = equations->data.function.args;
        m = (int)equations->data.function.arg_count;
    } else if (equations->type == EXPR_FUNCTION
        && equations->data.function.head->type == EXPR_SYMBOL
        && equations->data.function.head->data.symbol.name == SYM_Equal) {
        single_holder[0] = equations;
        eq_arr = single_holder;
        m = 1;
    } else {
        return NULL;
    }

    /* Empty system -- vacuously true. */
    if (m == 0) {
        Expr* empty = mk_list(NULL, 0);
        return mk_list((Expr*[]){ empty }, 1);
    }

    /* Domain. */
    SolvePolyOpts local_opts = {0};
    if (opts) local_opts = *opts;
    bool reals_only = false, integers_only = false;
    if (dom) {
        if (dom->type != EXPR_SYMBOL) return NULL;
        const char* d = dom->data.symbol.name;
        if (d == SYM_Reals) reals_only = true;
        else if (d == SYM_Integers) { reals_only = true; integers_only = true; }
        else if (d != SYM_Complexes) return NULL;
    }

    /* Canonicalise each equation `lhs == rhs`.  We clear denominators so
     * rational systems (1/x + 1/y == 1) become polynomial: put lhs - rhs over
     * a common denominator (Together), take the Numerator as the polynomial
     * equation, and record a non-constant Denominator so any completed tuple
     * on its zero locus is pruned (a den-zero point is never a solution).
     * The numerator is then Expand-ed so gb_from_expr's single-term parser can
     * ingest products / powers of sums; the lex zero-dim check and the final
     * verification both read orig[i], so all consumers see one canonical
     * Plus-of-monomials form.  For a pure polynomial the denominator is 1
     * (dens[i] == NULL) and this reduces to Expand[lhs - rhs]. */
    Expr** orig = (Expr**)malloc(sizeof(Expr*) * (size_t)m);
    Expr** dens = (Expr**)calloc((size_t)m, sizeof(Expr*)); /* NULL = constant den */
    int norig = 0;
    for (int i = 0; i < m; i++) {
        Expr* eq = eq_arr[i];
        if (eq->type != EXPR_FUNCTION
            || eq->data.function.head->type != EXPR_SYMBOL
            || eq->data.function.head->data.symbol.name != SYM_Equal
            || eq->data.function.arg_count != 2) {
            for (int k = 0; k < norig; k++) { expr_free(orig[k]); if (dens[k]) expr_free(dens[k]); }
            free(orig); free(dens);
            return NULL;
        }
        Expr* diff = mk_fn2("Plus",
            expr_copy(eq->data.function.args[0]),
            mk_neg(expr_copy(eq->data.function.args[1])));
        /* internal_* wrappers each consume their argument and return an owned,
         * evaluated Expr. */
        Expr* together = eval_and_free(internal_together((Expr*[]){ diff }, 1));
        Expr* num = eval_and_free(
            internal_numerator((Expr*[]){ expr_copy(together) }, 1));
        Expr* den = eval_and_free(
            internal_denominator((Expr*[]){ together }, 1)); /* consumes together */
        orig[norig] = eval_and_free(internal_expand((Expr*[]){ num }, 1));
        bool den_has_var = false;
        for (int v = 0; v < n; v++)
            if (contains_any_symbol_from(den, var_arr[v])) { den_has_var = true; break; }
        if (den_has_var) dens[norig] = den;      /* prune its zero locus later */
        else            { expr_free(den); dens[norig] = NULL; }
        norig++;
    }

    /* Parametric systems: any free symbol that is neither a solve variable nor
     * a known constant is a parameter that belongs in the coefficient field.
     * gb_from_expr and the machine walk below carry only rational coefficients,
     * so route these through the RationalFunctions-coefficient Gröbner engine
     * instead.  More than NL_MAX_PARAMS distinct parameters -> decline. */
    {
        const char* params[NL_MAX_PARAMS + 1];
        int np = 0; bool within_cap = true;
        for (int i = 0; i < norig && within_cap; i++)
            within_cap = nl_collect_params(orig[i], var_arr, n, params, &np,
                                           NL_MAX_PARAMS + 1);
        if (np > 0) {
            Expr* out = (within_cap && np <= NL_MAX_PARAMS)
                ? solve_parametric_system(orig, dens, norig, var_arr, n, vars,
                                          dom, reals_only, integers_only,
                                          &local_opts)
                : NULL;
            for (int k = 0; k < norig; k++) { expr_free(orig[k]); if (dens[k]) expr_free(dens[k]); }
            free(orig); free(dens);
            return out;
        }
    }

    /* Build the GBPoly system under lex order. */
    GBPoly** F = (GBPoly**)malloc(sizeof(GBPoly*) * (size_t)norig);
    int nF = 0; bool poly_ok = true;
    for (int i = 0; i < norig; i++) {
        GBPoly* g = gb_from_expr(orig[i], var_arr, n, GB_ORDER_LEX, 0, NULL);
        if (!g) { poly_ok = false; break; }
        if (!gb_poly_is_zero(g)) F[nF++] = g; else gb_poly_free(g);
    }
    if (!poly_ok) {
        for (int i = 0; i < nF; i++) gb_poly_free(F[i]);
        free(F);
        for (int k = 0; k < norig; k++) { expr_free(orig[k]); if (dens[k]) expr_free(dens[k]); }
        free(orig); free(dens);
        return NULL;                       /* non-polynomial -> unevaluated */
    }

    /* All equations trivially 0 == 0 -> tautology { {} }. */
    if (nF == 0) {
        free(F);
        for (int k = 0; k < norig; k++) { expr_free(orig[k]); if (dens[k]) expr_free(dens[k]); }
        free(orig); free(dens);
        Expr* empty = mk_list(NULL, 0);
        return mk_list((Expr*[]){ empty }, 1);
    }

    /* Total-degree gate. */
    for (int i = 0; i < nF; i++) {
        if (gbpoly_total_degree(F[i]) > NL_MAX_TDEG) {
            for (int k = 0; k < nF; k++) gb_poly_free(F[k]);
            free(F);
            for (int k = 0; k < norig; k++) { expr_free(orig[k]); if (dens[k]) expr_free(dens[k]); }
            free(orig); free(dens);
            return NULL;
        }
    }

    size_t nG = 0;
    GBPoly** G = gb_groebner_walk(F, (size_t)nF, GB_ORDER_LEX, NULL, &nG);
    for (int i = 0; i < nF; i++) gb_poly_free(F[i]);
    free(F);
    if (!G || nG == 0) {
        if (G) gb_basis_free(G, nG);
        for (int k = 0; k < norig; k++) { expr_free(orig[k]); if (dens[k]) expr_free(dens[k]); }
        free(orig); free(dens);
        return NULL;
    }

    /* Unit ideal (a non-zero constant generator) -> inconsistent -> {}. */
    for (size_t i = 0; i < nG; i++) {
        if (gb_poly_is_constant(G[i]) && !gb_poly_is_zero(G[i])) {
            gb_basis_free(G, nG);
            for (int k = 0; k < norig; k++) { expr_free(orig[k]); if (dens[k]) expr_free(dens[k]); }
            free(orig); free(dens);
            return mk_list(NULL, 0);
        }
    }

    /* Zero-dimensionality: every variable must own a pure-power leading
     * monomial among the basis LMs.  Otherwise the ideal is
     * positive-dimensional (infinitely many solutions). */
    {
        int* haspp = (int*)calloc((size_t)n, sizeof(int));
        for (size_t i = 0; i < nG; i++) {
            const int* lm = gb_poly_lm(G[i]);
            if (!lm) continue;
            int nz = -1, cnt = 0;
            for (int k = 0; k < n; k++) if (lm[k] > 0) { nz = k; cnt++; }
            if (cnt == 1) haspp[nz] = 1;
        }
        bool zerodim = true;
        for (int k = 0; k < n; k++) if (!haspp[k]) { zerodim = false; break; }
        free(haspp);
        if (!zerodim) {
            gb_basis_free(G, nG);
            warn_nsdim(expr_hash(equations));
            for (int k = 0; k < norig; k++) { expr_free(orig[k]); if (dens[k]) expr_free(dens[k]); }
            free(orig); free(dens);
            return NULL;
        }
    }

    /* Convert the basis to Expr and run the back-substitution search. */
    Expr** gbE = (Expr**)malloc(sizeof(Expr*) * nG);
    for (size_t i = 0; i < nG; i++) gbE[i] = gb_to_expr(G[i], var_arr);
    gb_basis_free(G, nG);

    Expr* out = run_nl_search(orig, dens, norig, gbE, (int)nG, var_arr, n,
                              dom, reals_only, integers_only, &local_opts);

    for (size_t i = 0; i < nG; i++) expr_free(gbE[i]);
    free(gbE);
    for (int k = 0; k < norig; k++) { expr_free(orig[k]); if (dens[k]) expr_free(dens[k]); }
    free(orig); free(dens);
    return out;
}

/* ------------------------------------------------------------------ *
 *  Builtin entry & init.                                              *
 * ------------------------------------------------------------------ */

Expr* builtin_solve_nonlinear_system(Expr* res) {
    if (!res || res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    if (argc < 2 || argc > 3) return NULL;
    Expr* equations = res->data.function.args[0];
    Expr* vars      = res->data.function.args[1];
    Expr* dom       = (argc >= 3) ? res->data.function.args[2] : NULL;
    return solvenlsys_solve_nonlinear_system(equations, vars, dom, NULL);
}

void solvenlsys_init(void) {
    symtab_add_builtin("Solve`SolveNonlinearSystem",
                       builtin_solve_nonlinear_system);
    SymbolDef* def = symtab_get_def("Solve`SolveNonlinearSystem");
    if (def) def->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("Solve`SolveNonlinearSystem",
        "Solve`SolveNonlinearSystem[eqns, vars]\n"
        "Solve`SolveNonlinearSystem[eqns, vars, dom]\n"
        "\tThe nonlinear polynomial-system specialist used by Solve for a\n"
        "\tzero-dimensional solution set.  Computes a lexicographic\n"
        "\tGröbner basis (triangular for a shape-position ideal), solves\n"
        "\tthe univariate generator in the last variable, and\n"
        "\tback-substitutes upward, verifying each completed tuple.\n"
        "\n"
        "\tReturns:\n"
        "\t  finite solutions       {{v1 -> a1, ...}, ...}\n"
        "\t  inconsistent system    {}\n"
        "\t  tautology (no eqns)    {{}}\n"
        "\n"
        "\tLeft unevaluated for non-polynomial systems, positive-\n"
        "\tdimensional ideals (Solve::nsdim), or any branch the\n"
        "\tunivariate solver cannot reduce (never a false {}).\n"
        "\n"
        "\tdom = Complexes (default), Reals, or Integers.");
}
