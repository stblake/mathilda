/* findmin_nm_common.c — shared global core: RNG, point eval, local-polish bridge, option/method parsing, results.
 * Split from the original findmin.c; shared declarations in
 * findmin_internal.h. Do not add cross-file helpers here without a
 * prototype in that header. */
#include "findmin_internal.h"


void nm_rng_seed(NmRng* r, uint64_t seed) { r->s = seed; }

uint64_t nm_rng_next(NmRng* r) {
    uint64_t z = (r->s += 0x9E3779B97F4A7C15ULL);
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
    z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
    return z ^ (z >> 31);
}
double nm_rng_unif(NmRng* r) {            /* [0, 1)                    */
    return (double)(nm_rng_next(r) >> 11) * (1.0 / 9007199254740992.0);
}
double nm_rng_range(NmRng* r, double lo, double hi) {
    return lo + (hi - lo) * nm_rng_unif(r);
}
double nm_rng_normal(NmRng* r) {          /* standard normal (Box-Muller) */
    double u1 = nm_rng_unif(r), u2 = nm_rng_unif(r);
    if (u1 < 1e-300) u1 = 1e-300;
    return sqrt(-2.0 * log(u1)) * cos(6.28318530717958647692 * u2);
}

bool nm_is_head(const Expr* e, const char* sym) {
    return e && e->type == EXPR_FUNCTION
        && e->data.function.head->type == EXPR_SYMBOL
        && e->data.function.head->data.symbol.name == sym;
}

/* A variable "atom" is a bare symbol (x) or an indexed form (x[i], x[i,j], ...):
 * anything a Table[x[i], {i, ...}] spec can produce. A {x, lo, hi} spec (head
 * List) is NOT an atom — it is a bounded-variable spec handled separately. */
static bool nm_is_var_atom(const Expr* e) {
    if (!e) return false;
    if (e->type == EXPR_SYMBOL) return true;
    return e->type == EXPR_FUNCTION
        && e->data.function.head->type == EXPR_SYMBOL
        && e->data.function.head->data.symbol.name != SYM_List
        && e->data.function.head->data.symbol.name != SYM_Element;
}

/* Is `cons` already a boolean/relational constraint tree, or a wrapper (Table,
 * List, ...) that must be evaluated first to expand into one? */
bool nm_is_constraint_tree(const Expr* e) {
    if (!e || e->type != EXPR_FUNCTION
        || e->data.function.head->type != EXPR_SYMBOL) return false;
    const char* h = e->data.function.head->data.symbol.name;
    return h == SYM_And || h == SYM_Or || h == SYM_Not || h == SYM_Xor
        || h == SYM_Implies || h == SYM_Equal || h == SYM_Unequal
        || h == SYM_Less || h == SYM_LessEqual || h == SYM_Greater
        || h == SYM_GreaterEqual || h == SYM_Inequality || h == SYM_Element;
}

/* Structural substitution: a fresh copy of `e` with every subtree structurally
 * equal to from[i] replaced by a copy of to[i]. Used to rewrite indexed
 * optimisation variables (x[1], x[2], ...) to fresh scalar symbols so the whole
 * symbol-keyed solver machinery applies unchanged. */
Expr* nm_subst(Expr* e, Expr* const* from, Expr* const* to, size_t n) {
    if (!e) return NULL;
    for (size_t i = 0; i < n; i++)
        if (from[i] && expr_eq(e, from[i])) return expr_copy(to[i]);
    if (e->type != EXPR_FUNCTION) return expr_copy(e);
    size_t ac = e->data.function.arg_count;
    Expr* head = nm_subst(e->data.function.head, from, to, n);
    Expr** args = (Expr**)malloc(sizeof(Expr*) * (ac ? ac : 1));
    for (size_t j = 0; j < ac; j++)
        args[j] = nm_subst(e->data.function.args[j], from, to, n);
    Expr* r = expr_new_function(head, args, ac);
    free(args);
    return r;
}

/* Generate a fresh, unused scalar symbol to stand in for an indexed variable.
 * The name is interned and stable; the caller removes it (symtab_remove_symbol)
 * once the optimisation and result construction are done. */
Expr* nm_fresh_symbol(void) {
    static uint64_t ctr = 0;
    char buf[64];
    for (;;) {
        snprintf(buf, sizeof(buf), "NMinimize$%llu", (unsigned long long)ctr++);
        if (!symtab_lookup(buf)) break;
    }
    return expr_new_symbol(buf);
}

static void nm_free_rule_chain(Rule* r) {
    while (r) {
        Rule* nx = r->next;
        expr_free(r->pattern);
        expr_free(r->replacement);
        free(r);
        r = nx;
    }
}

void nm_heads_localize(NmHeadSave* sv, const char** names, size_t m) {
    for (size_t i = 0; i < m; i++) {
        SymbolDef* def = symtab_get_def(names[i]);
        sv[i].name  = names[i];
        sv[i].own   = def->own_values;
        sv[i].down  = def->down_values;
        sv[i].valid = true;
        def->own_values  = NULL;
        def->down_values = NULL;
    }
    if (m) eval_clock_bump();
}

void nm_heads_restore(NmHeadSave* sv, size_t m) {
    for (size_t i = 0; i < m; i++) {
        if (!sv[i].valid) continue;
        SymbolDef* def = symtab_get_def(sv[i].name);
        nm_free_rule_chain(def->own_values);
        nm_free_rule_chain(def->down_values);
        def->own_values  = sv[i].own;
        def->down_values = sv[i].down;
        sv[i].valid = false;
    }
    if (m) eval_clock_bump();
}

/* Objective at the (already integer-rounded) point xr: the compiled program if
 * present and it returns a finite real, else the interpreter. */
static bool nm_eval_obj(NmDriver* D, const double* xr, double* out) {
    if (D->f_prog && compiled_eval_real(D->f_prog, xr, out) && isfinite(*out))
        return true;
    return fm_eval_scalar(D->f_raw, D->binds, xr, D->n, D->opts, out) && isfinite(*out);
}

/* Apply a user "PenaltyFunction" f to one nonnegative constraint-violation
 * magnitude m, returning f[m] as a double. Evaluator numeric diagnostics are
 * muted, as everywhere else in the trial-point loop. A non-numeric, non-finite,
 * or negative result is rejected so the caller can fall back to the built-in
 * m^2 for that term (a custom penalty must be a usable nonnegative score). */
static bool nm_apply_penalty_fn(Expr* pf, double m, double* out) {
    Expr* arg  = expr_new_real(m);
    Expr* call = expr_new_function(expr_copy(pf), &arg, 1);
    arith_warnings_mute_push();
    Expr* v = eval_and_free(call);
    arith_warnings_mute_pop();
    bool ok = v && fm_expr_to_double_real(v, out) && isfinite(*out) && *out >= 0.0;
    expr_free(v);
    return ok;
}

/* Penalty of a boolean-of-comparisons constraint tree at x, used for
 * disjunctive (Or) constraints:
 *   And[c...]        → Σ penalty(c)            (all must hold)
 *   Or[c...]         → min penalty(c)          (any one holding scores 0)
 *   Inequality[...]  → Σ over adjacent pairs   (a chained conjunction)
 *   binary compare   → squared violation: max(0,g)^2 for g<=0, h^2 for h==0,
 *                      or penalty_fn[violation] when a custom "PenaltyFunction"
 *                      is supplied (matching nm_eval_pen's per-term rule).
 * A satisfied leaf contributes 0, so the whole expression is 0 iff feasible and
 * the (total ≤ NM_FEAS_RANK) feasibility test carries over unchanged. The tree
 * shape is pre-validated by fm_bool_supported at collection time; this returns
 * false only if a leaf cannot be evaluated to a finite real at this point. */
bool fm_bool_penalty(Expr* c, FmVarBind* binds, const double* x, size_t n,
                            const FmOpts* opts, Expr* penalty_fn, double* out) {
    const char* h = c->data.function.head->data.symbol.name;
    size_t ac = c->data.function.arg_count;
    if (h == SYM_And) {
        double total = 0.0;
        for (size_t i = 0; i < ac; i++) {
            double t;
            if (!fm_bool_penalty(c->data.function.args[i], binds, x, n, opts,
                                 penalty_fn, &t)) return false;
            total += t;
        }
        *out = total;
        return true;
    }
    if (h == SYM_Or) {
        double best = -1.0;
        for (size_t i = 0; i < ac; i++) {
            double t;
            if (!fm_bool_penalty(c->data.function.args[i], binds, x, n, opts,
                                 penalty_fn, &t)) return false;
            if (best < 0.0 || t < best) best = t;
        }
        *out = (best < 0.0) ? 0.0 : best;
        return true;
    }
    if (h == SYM_Inequality) {
        double total = 0.0;
        size_t npairs = (ac - 1) / 2;
        for (size_t k = 0; k < npairs; k++) {
            Expr* a  = c->data.function.args[2 * k];
            Expr* op = c->data.function.args[2 * k + 1];
            Expr* b  = c->data.function.args[2 * k + 2];
            Expr* pair_args[2] = { expr_copy(a), expr_copy(b) };
            Expr* pair = expr_new_function(expr_new_symbol(op->data.symbol.name),
                                           pair_args, 2);
            double t;
            bool ok = fm_bool_penalty(pair, binds, x, n, opts, penalty_fn, &t);
            expr_free(pair);
            if (!ok) return false;
            total += t;
        }
        *out = total;
        return true;
    }
    /* single binary comparison → squared / custom violation */
    Expr* g = NULL; bool eq = false;
    if (!fm_constraint_to_g(c, &g, &eq)) return false;
    double d;
    bool ok = fm_eval_scalar(g, binds, x, n, opts, &d);
    expr_free(g);
    if (!ok || !isfinite(d)) return false;
    double m = eq ? fabs(d) : (d > 0.0 ? d : 0.0);
    if (m == 0.0) { *out = 0.0; return true; }
    double term;
    if (!penalty_fn || !nm_apply_penalty_fn(penalty_fn, m, &term)) term = m * m;
    *out = term;
    return true;
}

/* Σ pen(g_i) over the general constraints at xr, each constraint via its
 * compiled program if present, else the interpreter. The per-constraint term is
 * the built-in squared violation — max(0, g)^2 for an inequality g <= 0, h^2 for
 * an equality h == 0 — unless a "PenaltyFunction" f was supplied, in which case a
 * *violated* constraint contributes f[violation] instead (Automatic ≡ #^2 &, so
 * this generalises the default exactly). A satisfied inequality always
 * contributes 0, keeping the feasibility test (total ≤ NM_FEAS_RANK) intact.
 * Returns false if any constraint cannot be evaluated at all. */
bool nm_eval_pen(NmDriver* D, const double* xr, double* out) {
    double total = 0.0;
    /* General (conjunctive) constraints: the compiled/penalty-fn path when either
     * is present, else the shared squared-penalty evaluator. */
    if (D->ngens > 0) {
        if (!D->g_progs && !D->penalty_fn) {
            double base;
            if (!fm_eval_penalty(D->gens, D->ngens, D->binds, xr, D->n, D->opts, &base))
                return false;
            total += base;
        } else {
            for (size_t k = 0; k < D->ngens; k++) {
                double d;
                bool got = D->g_progs && D->g_progs[k]
                        && compiled_eval_real(D->g_progs[k], xr, &d) && isfinite(d);
                if (!got && (!fm_eval_scalar(D->gens[k].expr, D->binds, xr, D->n,
                                             D->opts, &d) || !isfinite(d)))
                    return false;
                double m = D->gens[k].equality ? fabs(d) : (d > 0.0 ? d : 0.0);
                if (m == 0.0) continue;
                double term;
                if (!D->penalty_fn || !nm_apply_penalty_fn(D->penalty_fn, m, &term))
                    term = m * m;
                total += term;
            }
        }
    }
    /* Disjunctive (Or) constraints: each adds its minimum-branch penalty, which
     * is 0 exactly when at least one branch is satisfied. */
    for (size_t d = 0; d < D->ndisj; d++) {
        double dp;
        if (!fm_bool_penalty(D->disj[d].expr, D->binds, xr, D->n, D->opts,
                             D->penalty_fn, &dp))
            return false;
        total += dp;
    }
    *out = total;
    return true;
}

/* Objective value and total constraint violation at x. Integer coordinates
 * are rounded before evaluation. A non-evaluable objective or constraint is
 * treated as maximally bad so the search steers away from it. */
void nm_eval(NmDriver* D, const double* x, double* f_out, double* pen_out) {
    size_t n = D->n;
    double* xr = (double*)malloc(sizeof(double) * (n ? n : 1));
    for (size_t i = 0; i < n; i++) xr[i] = D->is_int[i] ? round(x[i]) : x[i];
    double fx;
    if (!nm_eval_obj(D, xr, &fx)) fx = 1e300;
    double pen = 0.0;
    if (!nm_eval_pen(D, xr, &pen)) pen = 1e300;
    free(xr);
    *f_out = fx;
    *pen_out = pen;
}

/* Deb's feasibility rules, parameterised by which feasibility threshold applies.
 * A feasible point always beats an infeasible one; among feasible points the
 * smaller objective wins; among infeasible points the smaller total violation
 * wins. The threshold decides which branch is reachable at all, which is why it
 * is a parameter and not a constant here — see findmin_internal.h. */
static bool nm_better_eps(double fa, double pa, double fb, double pb,
                          double feas_eps) {
    bool fa_feas = (pa <= feas_eps);
    bool fb_feas = (pb <= feas_eps);
    if (fa_feas && fb_feas) return fa < fb;
    if (fa_feas != fb_feas) return fa_feas;
    return pa < pb;
}

/* SEARCH-time ranking. Loose (NM_FEAS_RANK) so that the "both feasible, compare
 * objectives" branch is actually reachable on problems whose achievable residual
 * is small but nonzero. Use this inside an engine's search loop. */
bool nm_better(double fa, double pa, double fb, double pb) {
    return nm_better_eps(fa, pa, fb, pb, NM_FEAS_RANK);
}

/* RETURN-path selection. Tight (NM_FEAS_RETURN), for any comparison whose loser
 * is discarded from the answer the caller receives — post-polish selection and
 * the driver's cross-attempt best. This is the half that fixes the bug: with the
 * loose threshold, an unpolished point violating by 1e-4 counted as "feasible"
 * and its slightly-lower objective (1.9998 < 2.0) beat the genuinely feasible
 * polished point. Under NM_FEAS_RETURN it is correctly infeasible and loses. */
bool nm_better_return(double fa, double pa, double fb, double pb) {
    return nm_better_eps(fa, pa, fb, pb, NM_FEAS_RETURN);
}

/* Clamp x into the search region and snap integer coordinates. */
void nm_project(NmDriver* D, double* x) {
    for (size_t j = 0; j < D->n; j++) {
        if (x[j] < D->reg_lo[j]) x[j] = D->reg_lo[j];
        if (x[j] > D->reg_hi[j]) x[j] = D->reg_hi[j];
        if (D->is_int[j]) x[j] = round(x[j]);
    }
}

/* Penalized scalar objective used by NelderMead / SimulatedAnnealing. */
double nm_phi(NmDriver* D, const double* x) {
    double f, p;
    nm_eval(D, x, &f, &p);
    return f + NM_PENALTY_MU * p;
}

/* ------------------------------------------------------------------ *
 *  Local polish of a candidate (shared by RandomSearch + the driver) *
 * ------------------------------------------------------------------ */

/* Mixed/integer coordinate descent: step each coordinate (±1, ±2 on integer
 * dims; scaled steps on continuous dims), accept any Deb-improvement. */
/* Map an effective-variable symbol to its index; SIZE_MAX if not found. */
static size_t nm_var_index(NmDriver* D, const char* sym) {
    for (size_t i = 0; i < D->n; i++)
        if (D->vars[i]->type == EXPR_SYMBOL && D->vars[i]->data.symbol.name == sym)
            return i;
    return (size_t)-1;
}

static bool nm_is_binary(NmDriver* D, size_t i) {
    return D->is_int[i] &&
           fabs(D->reg_lo[i]) < 1e-9 && fabs(D->reg_hi[i] - 1.0) < 1e-9;
}

/* If `plus` is Plus[v1, v2, ...] over DISTINCT binary optimization variables,
 * fill idx[0..*nv) with their indices and return true; else false. */
static bool nm_plus_of_binaries(NmDriver* D, const Expr* plus, size_t* idx, size_t* nv) {
    if (!nm_is_head(plus, SYM_Plus)) return false;
    *nv = 0;
    for (size_t a = 0; a < plus->data.function.arg_count; a++) {
        const Expr* arg = plus->data.function.args[a];
        if (arg->type != EXPR_SYMBOL) return false;
        size_t vi = nm_var_index(D, arg->data.symbol.name);
        if (vi == (size_t)-1 || !nm_is_binary(D, vi)) return false;
        for (size_t p = 0; p < *nv; p++) if (idx[p] == vi) return false; /* distinct */
        idx[(*nv)++] = vi;
    }
    return *nv >= 2;
}

/* Detect one-hot / assignment groups (Σ binaries == k) from the equality
 * constraints, in either the Subtract[Plus[...], k] form the constraint builder
 * produces or a canonicalized Plus[..., -k]. Populates D->onehots. */
void nm_detect_onehots(NmDriver* D) {
    D->onehots = NULL; D->n_onehots = 0;
    if (D->ngens == 0) return;
    NmOneHot* out = (NmOneHot*)malloc(sizeof(NmOneHot) * D->ngens);
    if (!out) return;
    size_t count = 0;
    for (size_t g = 0; g < D->ngens; g++) {
        if (!D->gens[g].equality) continue;
        const Expr* e = D->gens[g].expr;
        if (!e || e->type != EXPR_FUNCTION) continue;
        size_t* idx = (size_t*)malloc(sizeof(size_t) * D->n);
        if (!idx) continue;
        size_t nv = 0; int k = 0; bool ok = false;
        if (nm_is_head(e, SYM_Subtract) && e->data.function.arg_count == 2 &&
            e->data.function.args[1]->type == EXPR_INTEGER) {
            /* Subtract[Plus[vars...], k]  ⇒  Σ vars == k */
            if (nm_plus_of_binaries(D, e->data.function.args[0], idx, &nv)) {
                k = (int)e->data.function.args[1]->data.integer;
                ok = true;
            }
        } else if (nm_is_head(e, SYM_Plus)) {
            /* Plus[vars..., -k]  ⇒  Σ vars == k. Split off a single integer term. */
            size_t m = e->data.function.arg_count;
            long konst = 0; int nints = 0; nv = 0; ok = true;
            for (size_t a = 0; a < m && ok; a++) {
                const Expr* arg = e->data.function.args[a];
                if (arg->type == EXPR_INTEGER) { konst += arg->data.integer; nints++; }
                else if (arg->type == EXPR_SYMBOL) {
                    size_t vi = nm_var_index(D, arg->data.symbol.name);
                    if (vi == (size_t)-1 || !nm_is_binary(D, vi)) { ok = false; break; }
                    for (size_t p = 0; p < nv; p++) if (idx[p] == vi) { ok = false; break; }
                    if (ok) idx[nv++] = vi;
                } else ok = false;
            }
            k = (int)(-konst);
            ok = ok && nints == 1 && nv >= 2;
        }
        if (ok && k >= 1 && (size_t)k <= nv) {
            out[count].idx = idx; out[count].len = nv; out[count].k = k;
            count++;
        } else {
            free(idx);
        }
    }
    if (count == 0) { free(out); return; }
    D->onehots = out; D->n_onehots = count;
}

void nm_free_onehots(NmDriver* D) {
    for (size_t g = 0; g < D->n_onehots; g++) free(D->onehots[g].idx);
    free(D->onehots);
    D->onehots = NULL; D->n_onehots = 0;
}

/* Snap each one-hot group to exactly k ones, so the swap move can then reach
 * full feasibility. Members currently 1 are kept first (up to k); the shortfall
 * is filled from the not-yet-used members and any excess turned off. Uses the
 * point's own rounded pattern to break ties, so different restarts repair to
 * different assignments. Returns true iff it changed x. */
static bool nm_repair_onehots(NmDriver* D, double* x) {
    bool changed = false;
    for (size_t g = 0; g < D->n_onehots; g++) {
        NmOneHot* oh = &D->onehots[g];
        int ones = 0;
        for (size_t p = 0; p < oh->len; p++) if (round(x[oh->idx[p]]) >= 0.5) ones++;
        if (ones == oh->k) continue;
        int want = oh->k;
        int kept = 0;
        /* keep the first `want` members that are already 1, drop the rest */
        for (size_t p = 0; p < oh->len; p++) {
            size_t v = oh->idx[p];
            if (round(x[v]) >= 0.5) {
                if (kept < want) { x[v] = 1.0; kept++; }
                else             { x[v] = 0.0; changed = true; }
            }
        }
        /* fill the shortfall from members currently 0 */
        for (size_t p = 0; p < oh->len && kept < want; p++) {
            size_t v = oh->idx[p];
            if (round(x[v]) < 0.5) { x[v] = 1.0; kept++; changed = true; }
        }
    }
    return changed;
}

/* Does the symbol `sym` (interned name) occur anywhere in expression `e`? */
bool nm_expr_contains_symbol(const Expr* e, const char* sym) {
    if (!e) return false;
    if (e->type == EXPR_SYMBOL) return e->data.symbol.name == sym;
    if (e->type == EXPR_FUNCTION) {
        if (nm_expr_contains_symbol(e->data.function.head, sym)) return true;
        for (size_t i = 0; i < e->data.function.arg_count; i++)
            if (nm_expr_contains_symbol(e->data.function.args[i], sym)) return true;
    }
    return false;
}

static void nm_int_descent(NmDriver* D, double* x, double* f_io, double* pen_io) {
    size_t n = D->n;
    double* t = (double*)malloc(sizeof(double) * n);
    /* Binary coordinates — integer variables with a [0,1] region. These are the
     * candidates for the swap move below. */
    size_t* bidx = (size_t*)malloc(sizeof(size_t) * (n ? n : 1));
    size_t nb = 0;
    for (size_t j = 0; j < n; j++)
        if (D->is_int[j] && fabs(D->reg_lo[j]) < 1e-9 && fabs(D->reg_hi[j] - 1.0) < 1e-9)
            bidx[nb++] = j;
    /* Cap the O(nb^2) swap sweep so a very large binary problem cannot blow up. */
    bool do_swaps = (nb >= 2 && nb <= 400);
    /* Two passes: a plain descent, then — only if it ends INFEASIBLE and there are
     * one-hot / assignment groups — snap each group to k ones and descend again.
     * Gating the repair on residual infeasibility keeps it a FALLBACK: problems
     * the plain descent already solves (e.g. QAP, whose objective guides the
     * search to good permutations) are unchanged, while those it cannot make
     * feasible from a random rounded start (pure assignment, sudoku) are rescued. */
    for (int attempt = 0; attempt < 2; attempt++) {
        if (attempt == 1) {
            /* LOOSE (NM_FEAS_RANK) deliberately. This gate keeps the one-hot
             * repair a FALLBACK for starts the plain descent cannot make
             * feasible; using the tight return threshold here would make the
             * break almost never fire, running the repair on problems (QAP)
             * that the comment above says are meant to be left unchanged. */
            if (D->n_onehots == 0 || *pen_io <= NM_FEAS_RANK) break;
            if (nm_repair_onehots(D, x)) { nm_project(D, x); nm_eval(D, x, f_io, pen_io); }
        }
        /* The move-acceptance comparisons below stay LOOSE (nm_better), because
         * they are SEARCH moves, not a return decision — the descent is hunting
         * for a better lattice point, and Deb's rule needs its compare-objectives
         * branch reachable to do that. Tightening them here was measured to break
         * the search outright: on the 50-variable fixed-charge flow MINLP the
         * descent stopped finding feasible flows and settled on the trivial
         * all-zeros point (flow residual 20.0, the entire demand).
         *
         * What a mixed-integer run actually returns is bounded instead by the
         * last-chance continuous refinement at the end of nm_local_polish and by
         * the driver's NM_FEAS_RETURN gate — the return path, where the tight
         * threshold belongs. */
        bool improved = true;
        int iter = 0;
        while (improved && iter++ < 300) {
        improved = false;
        for (size_t j = 0; j < n; j++) {
            double steps[4];
            if (D->is_int[j]) {
                steps[0] = 1; steps[1] = -1; steps[2] = 2; steps[3] = -2;
            } else {
                double h = (D->reg_hi[j] - D->reg_lo[j]) * 0.05 + 1e-3;
                steps[0] = h; steps[1] = -h; steps[2] = 4 * h; steps[3] = -4 * h;
            }
            for (int s = 0; s < 4; s++) {
                for (size_t k = 0; k < n; k++) t[k] = x[k];
                t[j] = x[j] + steps[s];
                nm_project(D, t);
                double f, p;
                nm_eval(D, t, &f, &p);
                if (nm_better(f, p, *f_io, *pen_io)) {
                    for (size_t k = 0; k < n; k++) x[k] = t[k];
                    *f_io = f; *pen_io = p;
                    improved = true;
                }
            }
        }
        /* 2-flip SWAP move over binary coordinates: flip a 1 and a 0 together.
         * This preserves every Sum-of-binaries constraint that contains both
         * coordinates, so it moves *along* the assignment / permutation manifold
         * — which single-coordinate flips cannot, because any single flip off a
         * near-feasible permutation breaks a row or column sum and is rejected by
         * the Deb rule. It is the move that lets integer descent reach feasibility
         * on assignment-structured problems (QAP, TSP, sudoku, transport). Run
         * only when the single-coordinate sweep has stalled, so it is the escape
         * move rather than the common case; first-improvement, restart on accept. */
        if (do_swaps && !improved) {
            for (size_t a = 0; a < nb && !improved; a++) {
                size_t j = bidx[a];
                for (size_t b = a + 1; b < nb; b++) {
                    size_t k = bidx[b];
                    if (round(x[j]) == round(x[k])) continue;  /* only true swaps */
                    for (size_t m = 0; m < n; m++) t[m] = x[m];
                    t[j] = 1.0 - x[j];
                    t[k] = 1.0 - x[k];
                    nm_project(D, t);
                    double f, p;
                    nm_eval(D, t, &f, &p);
                    if (nm_better(f, p, *f_io, *pen_io)) {
                        for (size_t m = 0; m < n; m++) x[m] = t[m];
                        *f_io = f; *pen_io = p;
                        improved = true;
                        break;
                    }
                }
            }
        }
        }
    }
    free(bidx);
    free(t);
}

/* Deep-free a general-constraint array built by nm_polish_gens. */
static void fm_free_gens(FmGenCon* g, size_t ng, size_t n) {
    if (!g) return;
    for (size_t k = 0; k < ng; k++) {
        expr_free(g[k].expr);
        if (g[k].grad_exprs) {
            for (size_t i = 0; i < n; i++) expr_free(g[k].grad_exprs[i]);
            free(g[k].grad_exprs);
        }
    }
    free(g);
}

/* Flatten one disjunction branch (And / Inequality / single comparison — never
 * Or, since a top-level Or's args are the individual disjuncts) into owned
 * general constraints appended to *g, each with its symbolic gradient. Returns
 * false if any leaf is not a convertible comparison. */
static bool nm_collect_branch_gens(Expr* c, Expr** vars, size_t n,
                                   FmGenCon** g, size_t* ng, size_t* cap) {
    const char* h = c->data.function.head->data.symbol.name;
    size_t ac = c->data.function.arg_count;
    if (h == SYM_And) {
        for (size_t i = 0; i < ac; i++)
            if (!nm_collect_branch_gens(c->data.function.args[i], vars, n, g, ng, cap))
                return false;
        return true;
    }
    if (h == SYM_Inequality) {
        size_t npairs = (ac - 1) / 2;
        for (size_t k = 0; k < npairs; k++) {
            Expr* a  = c->data.function.args[2 * k];
            Expr* op = c->data.function.args[2 * k + 1];
            Expr* b  = c->data.function.args[2 * k + 2];
            Expr* pa[2] = { expr_copy(a), expr_copy(b) };
            Expr* pair = expr_new_function(expr_new_symbol(op->data.symbol.name), pa, 2);
            bool ok = nm_collect_branch_gens(pair, vars, n, g, ng, cap);
            expr_free(pair);
            if (!ok) return false;
        }
        return true;
    }
    Expr* ge = NULL; bool eq = false;
    if (!fm_constraint_to_g(c, &ge, &eq)) return false;
    if (*ng == *cap) {
        size_t nc = *cap ? (*cap) * 2 : 4;
        *g = (FmGenCon*)realloc(*g, sizeof(FmGenCon) * nc);
        *cap = nc;
    }
    (*g)[*ng].expr = ge;
    /* Finite-difference gradient (NULL ⇒ fm_run_penalty FDs on demand): the
     * optimisation variables carry transient value-bindings during the search,
     * so symbolic D[...] taken here differentiates a constant and yields 0. */
    (*g)[*ng].grad_exprs = NULL;
    (*g)[*ng].equality = eq;
    (*ng)++;
    return true;
}

/* Effective smooth-constraint set for a local polish at x when disjunctions are
 * present: deep copies of the conjunctive base D->gens, plus — for each
 * disjunction — the constraints of its currently-active (minimum-penalty) branch
 * at x. Folding in the active branch turns the non-smooth Or into a smooth local
 * problem the BFGS penalty solver can descend, so the polish refines *within* the
 * feasible region x already occupies instead of ignoring the Or and drifting into
 * the infeasible gap between branches (which left RandomSearch, whose only
 * descent is this polish, stranded). Returns NULL (and does not touch *nout) when
 * D->ndisj == 0, signalling the caller to use D->gens directly; otherwise returns
 * a newly-allocated array of length *nout, freed with fm_free_gens. */
static FmGenCon* nm_polish_gens(NmDriver* D, const double* x, size_t* nout) {
    if (D->ndisj == 0) return NULL;
    FmGenCon* g = NULL; size_t ng = 0, cap = 0;
    for (size_t k = 0; k < D->ngens; k++) {
        if (ng == cap) {
            size_t nc = cap ? cap * 2 : 4;
            g = (FmGenCon*)realloc(g, sizeof(FmGenCon) * nc);
            cap = nc;
        }
        g[ng].expr = expr_copy(D->gens[k].expr);
        g[ng].grad_exprs = NULL;   /* FD during search; see nm_collect_branch_gens */
        g[ng].equality = D->gens[k].equality;
        ng++;
    }
    for (size_t d = 0; d < D->ndisj; d++) {
        Expr* orx = D->disj[d].expr;          /* Or[branch, ...] */
        size_t ac = orx->data.function.arg_count;
        size_t best = 0; double bestp = -1.0;
        for (size_t i = 0; i < ac; i++) {
            double p;
            if (!fm_bool_penalty(orx->data.function.args[i], D->binds, x, D->n,
                                 D->opts, D->penalty_fn, &p))
                p = 1e300;
            if (bestp < 0.0 || p < bestp) { bestp = p; best = i; }
        }
        /* On an unconvertible branch, drop its augmentation (the global search's
         * exact penalty still gates feasibility); rewind any partial append. */
        size_t ng_save = ng;
        if (!nm_collect_branch_gens(orx->data.function.args[best], D->vars, D->n,
                                    &g, &ng, &cap)) {
            for (size_t k = ng_save; k < ng; k++) {
                expr_free(g[k].expr);
                if (g[k].grad_exprs) {
                    for (size_t i = 0; i < D->n; i++) expr_free(g[k].grad_exprs[i]);
                    free(g[k].grad_exprs);
                }
            }
            ng = ng_save;
        }
    }
    *nout = ng;
    return g;
}

/* Run the exact continuous local solver over all variables from x, confined
 * only by the real box constraints (D->boxes) — never the DE sampling region,
 * so a coordinate with no explicit bound is free to leave the default span
 * (the freedom the pure-continuous polish already has). When pin_int is true
 * each integer coordinate is pinned to its rounded value with a degenerate
 * [v, v] box (fm_line_search projects every trial into the box, so it stays
 * fixed); when false the integer coordinates are relaxed and solved as
 * continuous. x is overwritten with the refined point and (f, penalty) at the
 * integer-rounded point is returned; the caller decides whether to keep it. */
static void nm_continuous_solve(NmDriver* D, double* x, bool pin_int,
                                double* f_out, double* pen_out) {
    size_t n = D->n;
    FmBox* tb = (FmBox*)malloc(sizeof(FmBox) * n);
    if (!tb) { nm_eval(D, x, f_out, pen_out); return; }
    for (size_t j = 0; j < n; j++) {
        if (pin_int && D->is_int[j]) {
            double v = round(x[j]);
            x[j] = v;
            tb[j].has_lo = tb[j].has_hi = true;
            tb[j].lo = tb[j].hi = v;
        } else {
            tb[j] = D->boxes[j];             /* real constraint bounds, else free    */
        }
    }
    double fx = 0.0;
    bool saved_quiet = g_fm_quiet;
    g_fm_quiet = true;                        /* silence internal solver chatter      */
    size_t png = D->ngens;
    FmGenCon* pg = nm_polish_gens(D, x, &png);           /* NULL ⇒ use D->gens */
    FmGenCon* use = pg ? pg : D->gens;
    /* With a compiled objective, take the gradient by finite differences off it
     * (fm_grad_finite_diff → fm_eval_scalar → compiled) rather than evaluating
     * the symbolic gradient through the interpreter — 2n compiled evals beats
     * n symbolic evals each re-binding n OwnValues. */
    Expr** ge = D->f_prog ? NULL : D->g_exprs;
    if (png > 0)
        (void)fm_run_penalty(D->f_raw, D->vars, n, D->binds,
                             FM_METHOD_QUASINEWTON, ge, NULL, x,
                             use, png, tb, D->opts, &fx);
    else
        (void)fm_run_bfgs(D->f_raw, D->vars, n, D->binds, ge, x,
                          NULL, 0, 0.0, tb, D->opts, &fx);
    if (pg) fm_free_gens(pg, png, n);
    g_fm_quiet = saved_quiet;
    free(tb);
    nm_eval(D, x, f_out, pen_out);
}

void nm_local_polish(NmDriver* D, double* x, double* f_out, double* pen_out) {
    if (D->any_int) {
        nm_project(D, x);
        nm_eval(D, x, f_out, pen_out);
        nm_int_descent(D, x, f_out, pen_out);
        /* A mixed-integer problem whose feasible region lies outside the DE
         * sampling region (e.g. the pressure-vessel MINLP, feasible near
         * x3 ~ 52 with the default +-10 span) is invisible to the region-bound
         * integer descent above. Recover it with the continuous-relaxation +
         * rounding heuristic: solve the continuous relaxation (integers relaxed,
         * every coordinate free of the sampling region) to locate the basin,
         * round the integer coordinates, then refine the continuous coordinates
         * with the integers pinned. Adopt the result only when it is a
         * Deb-improvement, so this can never worsen the region-descent answer.
         * Skipped when every variable is integer (no continuous coordinate to
         * free), where the region descent is already the whole story. */
        bool has_cont = false;
        for (size_t j = 0; j < D->n; j++) if (!D->is_int[j]) { has_cont = true; break; }
        /* The continuous-relaxation recovery is a heavy 2× QuasiNewton penalty
         * solve. Its purpose is to optimize / region-rescue continuous OBJECTIVE
         * variables, so it is skipped when the continuous variables never touch
         * the objective (constraint-only helpers such as MTZ ordering variables):
         * running it there is pure overhead — often seconds per restart — and
         * cannot improve the answer, which lets a combinatorial problem afford far
         * more restarts. It is still run whenever a continuous variable is in the
         * objective (the pressure-vessel MINLP, cardinality portfolios, ...). */
        if (has_cont && D->cont_in_obj) {
            double* xr = (double*)malloc(sizeof(double) * D->n);
            if (xr) {
                for (size_t j = 0; j < D->n; j++) xr[j] = x[j];
                double fr, pr;
                nm_continuous_solve(D, xr, false, &fr, &pr);     /* relaxation      */
                for (size_t j = 0; j < D->n; j++)
                    if (D->is_int[j]) xr[j] = round(xr[j]);
                nm_continuous_solve(D, xr, true, &fr, &pr);      /* pin + refine    */
                if (nm_better_return(fr, pr, *f_out, *pen_out)) {
                    for (size_t j = 0; j < D->n; j++) x[j] = xr[j];
                    *f_out = fr; *pen_out = pr;
                }
                free(xr);
            }
            nm_int_descent(D, x, f_out, pen_out);   /* re-settle integers in-region */
        }
        /* LAST-CHANCE FEASIBILITY REFINEMENT, and the only thing standing between
         * a mixed-integer problem and an Infinity return.
         *
         * nm_int_descent moves INTEGER coordinates only, so a problem whose
         * continuous coordinates are still slightly off cannot be improved by it.
         * Under the old loose return threshold that did not show: the descent
         * stopped at a ~1e-4 violation and the driver called it feasible. Under
         * NM_FEAS_RETURN it is correctly infeasible, so without this step
         * NMinimize[{x + 2y, x^2 + 2y^2 <= 3, x + y == 2, x in Integers}, {x,y}]
         * returns Infinity despite x = 1, y = 1 being exactly feasible.
         *
         * So: only when the answer would otherwise be rejected, and only when
         * there is a continuous coordinate to move, pin the integers and refine.
         * Adopted only on a Deb improvement, so it can never worsen the answer,
         * and it costs nothing on the overwhelmingly common path where the
         * result is already feasible. */
        if (*pen_out > NM_FEAS_RETURN) {
            bool any_cont = false;
            for (size_t j = 0; j < D->n; j++)
                if (!D->is_int[j]) { any_cont = true; break; }
            if (any_cont) {
                double* xf = (double*)malloc(sizeof(double) * D->n);
                if (xf) {
                    for (size_t j = 0; j < D->n; j++) xf[j] = x[j];
                    double ff, pf;
                    nm_continuous_solve(D, xf, true, &ff, &pf);  /* integers pinned */
                    if (nm_better_return(ff, pf, *f_out, *pen_out)) {
                        for (size_t j = 0; j < D->n; j++) x[j] = xf[j];
                        *f_out = ff; *pen_out = pf;
                    }
                    free(xf);
                }
            }
        }
    } else {
        double fx = 0.0;
        bool saved_quiet = g_fm_quiet;
        g_fm_quiet = true;                 /* silence internal solver chatter */
        size_t png = D->ngens;
        FmGenCon* pg = nm_polish_gens(D, x, &png);       /* NULL ⇒ use D->gens */
        FmGenCon* use = pg ? pg : D->gens;
        /* Compiled objective ⇒ finite-difference gradient off it (see the
         * matching note in nm_continuous_solve). */
        Expr** ge = D->f_prog ? NULL : D->g_exprs;
        if (png > 0)
            (void)fm_run_penalty(D->f_raw, D->vars, D->n, D->binds,
                                 FM_METHOD_QUASINEWTON, ge, NULL, x,
                                 use, png, D->boxes, D->opts, &fx);
        else
            (void)fm_run_bfgs(D->f_raw, D->vars, D->n, D->binds, ge, x,
                              NULL, 0, 0.0, D->boxes, D->opts, &fx);
        if (pg) fm_free_gens(pg, png, D->n);
        g_fm_quiet = saved_quiet;
        nm_eval(D, x, f_out, pen_out);
    }
}

/* ------------------------------------------------------------------ *
 *  Option / method / variable / constraint parsing                    *
 * ------------------------------------------------------------------ */

static bool nm_method_from_string(const char* s, int* out) {
    if (strcmp(s, "DifferentialEvolution") == 0) { *out = NM_DE;           return true; }
    if (strcmp(s, "NelderMead") == 0)            { *out = NM_NELDERMEAD;   return true; }
    if (strcmp(s, "RandomSearch") == 0)          { *out = NM_RANDOMSEARCH; return true; }
    if (strcmp(s, "SimulatedAnnealing") == 0)    { *out = NM_SA;           return true; }
    if (strcmp(s, "SHGO") == 0)                  { *out = NM_SHGO;         return true; }
    if (strcmp(s, "DualAnnealing") == 0)         { *out = NM_DUAL_ANNEALING; return true; }
    if (strcmp(s, "DIRECT") == 0)                { *out = NM_DIRECT;       return true; }
    if (strcmp(s, "BasinHopping") == 0)          { *out = NM_BASIN_HOPPING; return true; }
    return false;
}

/* The Method sub-option LHS may be a string ("SearchPoints") or a symbol. */
static const char* nm_option_name(Expr* e) {
    if (e->type == EXPR_STRING) return e->data.string;
    if (e->type == EXPR_SYMBOL) return e->data.symbol.name;
    return NULL;
}

/* Human-readable method name for diagnostics (inverse of nm_method_from_string). */
static const char* nm_method_name(int m) {
    switch (m) {
        case NM_DE:             return "DifferentialEvolution";
        case NM_NELDERMEAD:     return "NelderMead";
        case NM_RANDOMSEARCH:   return "RandomSearch";
        case NM_SA:             return "SimulatedAnnealing";
        case NM_SHGO:           return "SHGO";
        case NM_DUAL_ANNEALING: return "DualAnnealing";
        case NM_DIRECT:         return "DIRECT";
        case NM_BASIN_HOPPING:  return "BasinHopping";
        default:                return "Automatic";
    }
}

/* Per-method scoping of the Method sub-options. Each key is valid only for the
 * method(s) that actually read its NmConfig field, so a key belonging to another
 * method (or a misspelled / scipy-style name) is reported rather than silently
 * stored-and-ignored — matching how top-level options already warn in
 * nm_apply_option. Returns a bitmask (1u << NM_xxx) of the owning methods,
 * OWN_GENERIC for keys every method honors, or 0 for an unknown key.
 *
 * INVARIANT: this key set MUST stay in sync with the strcmp chain in
 * nm_parse_method — every key parsed there appears here, and vice versa. */
#define OWN_GENERIC 0xFFFFFFFFu
static unsigned nm_option_owner(const char* on) {
    /* Honored by every method (search budget, RNG seed, final polish toggle,
     * and the shared infeasibility-scoring penalty). */
    if (strcmp(on, "SearchPoints") == 0 || strcmp(on, "RandomSeed") == 0 ||
        strcmp(on, "PostProcess") == 0  || strcmp(on, "PenaltyFunction") == 0)
        return OWN_GENERIC;
    /* Read by both DifferentialEvolution (basin-separation tol / multi-start
     * seeding) and NelderMead (simplex convergence / simplex seeding). */
    if (strcmp(on, "Tolerance") == 0 || strcmp(on, "InitialPoints") == 0)
        return (1u << NM_DE) | (1u << NM_NELDERMEAD);
    if (strcmp(on, "ScalingFactor") == 0 || strcmp(on, "CrossProbability") == 0)
        return 1u << NM_DE;
    if (strcmp(on, "ReflectRatio") == 0 || strcmp(on, "ExpandRatio") == 0 ||
        strcmp(on, "ContractRatio") == 0 || strcmp(on, "ShrinkRatio") == 0)
        return 1u << NM_NELDERMEAD;
    if (strcmp(on, "PerturbationScale") == 0 || strcmp(on, "LevelIterations") == 0 ||
        strcmp(on, "BoltzmannExponent") == 0)
        return 1u << NM_SA;
    if (strcmp(on, "SamplingMethod") == 0 || strcmp(on, "Iterations") == 0)
        return 1u << NM_SHGO;
    if (strcmp(on, "VisitingParameter") == 0 || strcmp(on, "AcceptanceParameter") == 0 ||
        strcmp(on, "InitialTemperature") == 0 || strcmp(on, "RestartTemperatureRatio") == 0 ||
        strcmp(on, "LocalSearch") == 0)
        return 1u << NM_DUAL_ANNEALING;
    if (strcmp(on, "LocallyBiased") == 0 || strcmp(on, "Epsilon") == 0 ||
        strcmp(on, "MaxFunctionEvaluations") == 0 || strcmp(on, "MaxIterations") == 0 ||
        strcmp(on, "VolumeTolerance") == 0 || strcmp(on, "LengthTolerance") == 0 ||
        strcmp(on, "MinValue") == 0 || strcmp(on, "MinValueTolerance") == 0)
        return 1u << NM_DIRECT;
    if (strcmp(on, "Temperature") == 0 || strcmp(on, "StepSize") == 0 ||
        strcmp(on, "StepInterval") == 0 || strcmp(on, "TargetAcceptanceRate") == 0 ||
        strcmp(on, "StepFactor") == 0 || strcmp(on, "SuccessIterations") == 0)
        return 1u << NM_BASIN_HOPPING;
    return 0u;   /* not a known NMinimize Method sub-option */
}

static bool nm_parse_method(Expr* rhs, NmConfig* nc, const char* fn) {
    if (rhs->type == EXPR_SYMBOL && rhs->data.symbol.name == SYM_Automatic) {
        nc->method = NM_AUTO;
        return true;
    }
    if (rhs->type == EXPR_STRING) {
        int m;
        if (nm_method_from_string(rhs->data.string, &m)) { nc->method = m; return true; }
        fm_warn(fn, "nimpl", "Method \"%s\" is not supported", rhs->data.string);
        return false;
    }
    if (nm_is_head(rhs, SYM_List) && rhs->data.function.arg_count >= 1) {
        Expr* h = rhs->data.function.args[0];
        int m;
        if (h->type != EXPR_STRING || !nm_method_from_string(h->data.string, &m)) {
            fm_warn(fn, "badmeth", "Method list must begin with a method-name string");
            return false;
        }
        nc->method = m;
        for (size_t i = 1; i < rhs->data.function.arg_count; i++) {
            Expr* r = rhs->data.function.args[i];
            if (!nm_is_head(r, SYM_Rule) && !nm_is_head(r, SYM_RuleDelayed)) continue;
            if (r->data.function.arg_count != 2) continue;
            const char* on = nm_option_name(r->data.function.args[0]);
            Expr* ov = r->data.function.args[1];
            if (!on) continue;
            /* Scope the key to the selected method: warn and skip anything
             * unknown or belonging to a different method rather than silently
             * storing-and-ignoring it. nc->method is a concrete method here (the
             * list form required a method-name string above). */
            unsigned owner = nm_option_owner(on);
            if (owner == 0u) {
                fm_warn(fn, "optx",
                        "\"%s\" is not a known Method sub-option; ignored", on);
                continue;
            }
            if (owner != OWN_GENERIC && !(owner & (1u << nc->method))) {
                fm_warn(fn, "optx",
                        "\"%s\" is not a valid sub-option for Method \"%s\"; ignored",
                        on, nm_method_name(nc->method));
                continue;
            }
            if (strcmp(on, "SearchPoints") == 0) {
                if (ov->type == EXPR_INTEGER && ov->data.integer > 0)
                    nc->search_points = (int)ov->data.integer;
            } else if (strcmp(on, "ScalingFactor") == 0) {
                /* DE differential weight F, mutation scale in DE/rand/1. A real
                 * in (0, 2]; an out-of-range or non-real value warns and keeps
                 * the default 0.6. */
                double dv;
                if (fm_expr_to_double_real(ov, &dv) && isfinite(dv) && dv > 0.0 && dv <= 2.0)
                    nc->F = dv;
                else
                    fm_warn(fn, "sopt",
                            "ScalingFactor must be a real in (0, 2]; using the default");
            } else if (strcmp(on, "CrossProbability") == 0) {
                /* DE crossover probability CR, the per-coordinate chance of
                 * taking the mutant. A real in [0, 1]; an out-of-range or
                 * non-real value warns and keeps the default 0.9. */
                double dv;
                if (fm_expr_to_double_real(ov, &dv) && isfinite(dv) && dv >= 0.0 && dv <= 1.0)
                    nc->CR = dv;
                else
                    fm_warn(fn, "sopt",
                            "CrossProbability must be a real in [0, 1]; using the default");
            } else if (strcmp(on, "RandomSeed") == 0) {
                if (ov->type == EXPR_INTEGER && ov->data.integer >= 0)
                    nc->seed = (uint64_t)ov->data.integer;
            } else if (strcmp(on, "ReflectRatio") == 0) {
                double dv; if (fm_expr_to_double_real(ov, &dv) && isfinite(dv) && dv > 0.0)
                    nc->reflect_ratio = dv;
            } else if (strcmp(on, "ExpandRatio") == 0) {
                double dv; if (fm_expr_to_double_real(ov, &dv) && isfinite(dv) && dv > 0.0)
                    nc->expand_ratio = dv;
            } else if (strcmp(on, "ContractRatio") == 0) {
                double dv; if (fm_expr_to_double_real(ov, &dv) && isfinite(dv) && dv > 0.0)
                    nc->contract_ratio = dv;
            } else if (strcmp(on, "ShrinkRatio") == 0) {
                double dv; if (fm_expr_to_double_real(ov, &dv) && isfinite(dv) && dv > 0.0)
                    nc->shrink_ratio = dv;
            } else if (strcmp(on, "Tolerance") == 0) {
                double dv; if (fm_expr_to_double_real(ov, &dv) && isfinite(dv) && dv > 0.0)
                    nc->tolerance = dv;
            } else if (strcmp(on, "PostProcess") == 0) {
                /* Every Mathematica setting is accepted. True / Automatic, a
                 * named local method as a string ("InteriorPoint",
                 * "FindMinimum", "KKT", ...), or {"Method", subopts...} enable
                 * the final exact polish; False / None disable it. Mathilda has a
                 * single FindMinimum-style local polish (BFGS for continuous/box,
                 * quadratic penalty for general constraints) that already selects
                 * the right inner solver for the problem, so a named method turns
                 * post-processing on rather than picking a distinct algorithm. */
                if (ov->type == EXPR_STRING) {
                    nc->post_process = 1;
                } else if (ov->type == EXPR_SYMBOL) {
                    const char* s = ov->data.symbol.name;
                    if (s == SYM_True || s == SYM_Automatic)      nc->post_process = 1;
                    else if (s == SYM_False || s == SYM_None)     nc->post_process = 0;
                    else fm_warn(fn, "pmeth",
                                 "PostProcess -> %s not recognised; using Automatic", s);
                } else if (nm_is_head(ov, SYM_List)
                           && ov->data.function.arg_count >= 1
                           && ov->data.function.args[0]->type == EXPR_STRING) {
                    nc->post_process = 1;           /* {"InteriorPoint", opts...} */
                } else {
                    fm_warn(fn, "pmeth", "invalid PostProcess value; using Automatic");
                }
            } else if (strcmp(on, "InitialPoints") == 0) {
                /* A list of starting points {{x1,...}, {x2,...}, ...}; borrowed
                 * and validated/consumed by the engine, where the dimension n is
                 * known. Anything else is ignored (falls back to random starts). */
                if (nm_is_head(ov, SYM_List) && ov->data.function.arg_count > 0)
                    nc->init_points = ov;
            } else if (strcmp(on, "PerturbationScale") == 0) {
                /* SimulatedAnnealing: multiplies the size of the random step
                 * used to generate a new trial point (default 1.0). A larger
                 * scale explores more widely; must be a positive finite real. */
                double dv;
                if (fm_expr_to_double_real(ov, &dv) && isfinite(dv) && dv > 0.0)
                    nc->perturb_scale = dv;
                else
                    fm_warn(fn, "sopt",
                            "PerturbationScale must be a positive real; using 1.0");
            } else if (strcmp(on, "LevelIterations") == 0) {
                /* SimulatedAnnealing: the number of trial moves spent at each
                 * temperature level. The per-chain iteration budget is
                 * MaxIterations * LevelIterations (default 50, the value the
                 * previous fixed multiplier hard-coded), so this scales how long
                 * each chain anneals. An explicit value is honored verbatim, like
                 * "SearchPoints"; Automatic / None restore the default. */
                if (ov->type == EXPR_INTEGER && ov->data.integer > 0
                    && ov->data.integer <= 2000000000LL)
                    nc->level_iterations = (int)ov->data.integer;
                else if (ov->type == EXPR_SYMBOL
                         && (ov->data.symbol.name == SYM_Automatic
                             || ov->data.symbol.name == SYM_None))
                    nc->level_iterations = 0;
                else
                    fm_warn(fn, "sopt",
                            "LevelIterations must be a positive integer; using Automatic");
            } else if (strcmp(on, "BoltzmannExponent") == 0) {
                /* SimulatedAnnealing: the exponent of the Metropolis acceptance
                 * probability for an uphill move. A function f is called as
                 * f[i, df, f0] (iteration i≥1, objective difference df≥0, current
                 * value f0) and the point is accepted with probability Exp[f[...]];
                 * Automatic / None keep the built-in geometric-cooling exponent
                 * -df/T. Borrowed from the held method list. */
                if (ov->type == EXPR_SYMBOL
                    && (ov->data.symbol.name == SYM_Automatic
                        || ov->data.symbol.name == SYM_None)) {
                    nc->boltzmann_fn = NULL;
                } else if (ov->type == EXPR_FUNCTION || ov->type == EXPR_SYMBOL) {
                    nc->boltzmann_fn = ov;
                } else {
                    fm_warn(fn, "bexp",
                            "invalid BoltzmannExponent value; using Automatic");
                }
            } else if (strcmp(on, "PenaltyFunction") == 0) {
                /* A function applied to each constraint's violation to score
                 * infeasible points during the global search; Automatic ≡ #^2 &.
                 * Automatic / None keep the built-in squared penalty; a pure
                 * function or a function symbol (#^2 &, (10 #) &, Sqrt, ...) is
                 * stored (borrowed from the held method list) and applied in
                 * nm_eval_pen. It affects only the global-search feasibility
                 * scoring — the final local polish keeps the differentiable
                 * squared penalty its analytic gradient assumes. */
                if (ov->type == EXPR_SYMBOL
                    && (ov->data.symbol.name == SYM_Automatic
                        || ov->data.symbol.name == SYM_None)) {
                    nc->penalty_fn = NULL;
                } else if (ov->type == EXPR_FUNCTION || ov->type == EXPR_SYMBOL) {
                    nc->penalty_fn = ov;
                } else {
                    fm_warn(fn, "penf",
                            "invalid PenaltyFunction value; using Automatic");
                }
            } else if (strcmp(on, "SamplingMethod") == 0) {
                /* SHGO: how the search box is sampled and its connectivity graph
                 * built. "Simplicial" (default) is the exact Kuhn-triangulation
                 * 1-skeleton; "Sobol" / "Halton" are low-discrepancy point sets
                 * with a k-nearest-neighbour graph. Case-insensitive on the first
                 * letter is not attempted — the exact Mathematica-style string is
                 * required; an unknown value warns and keeps Simplicial. */
                if (ov->type == EXPR_STRING) {
                    if (strcmp(ov->data.string, "Simplicial") == 0)   nc->shgo_sampling = 0;
                    else if (strcmp(ov->data.string, "Sobol") == 0)   nc->shgo_sampling = 1;
                    else if (strcmp(ov->data.string, "Halton") == 0)  nc->shgo_sampling = 2;
                    else fm_warn(fn, "sopt",
                                 "SamplingMethod must be \"Simplicial\", \"Sobol\", or "
                                 "\"Halton\"; using \"Simplicial\"");
                } else {
                    fm_warn(fn, "sopt",
                            "SamplingMethod must be a string; using \"Simplicial\"");
                }
            } else if (strcmp(on, "Iterations") == 0) {
                /* SHGO: the number of sampling/refinement iterations (scipy's
                 * `iters`). Each iteration grows the complex and re-pools; the
                 * homology-growth criterion stops early once no new local minimum
                 * appears. A positive integer is honored; Automatic / None -> 1. */
                if (ov->type == EXPR_INTEGER && ov->data.integer > 0
                    && ov->data.integer <= 100000LL)
                    nc->shgo_iters = (int)ov->data.integer;
                else if (ov->type == EXPR_SYMBOL
                         && (ov->data.symbol.name == SYM_Automatic
                             || ov->data.symbol.name == SYM_None))
                    nc->shgo_iters = 0;
                else
                    fm_warn(fn, "sopt",
                            "Iterations must be a positive integer; using Automatic");
            } else if (strcmp(on, "VisitingParameter") == 0) {
                /* DualAnnealing: qv, the visiting-distribution shape (scipy's
                 * `visit`). qv -> 1 is Gaussian (classical SA), 2 is Cauchy
                 * (fast SA), > 2 gives the heavy tails of GSA. Valid range (1, 3];
                 * an out-of-range or non-real value warns and keeps 2.62. */
                double dv;
                if (fm_expr_to_double_real(ov, &dv) && isfinite(dv)
                    && dv > 1.0 && dv <= 3.0)
                    nc->da_visit = dv;
                else
                    fm_warn(fn, "sopt",
                            "VisitingParameter must be a real in (1, 3]; using 2.62");
            } else if (strcmp(on, "AcceptanceParameter") == 0) {
                /* DualAnnealing: qa, the acceptance-distribution shape (scipy's
                 * `accept`). More negative => a smaller uphill-acceptance
                 * probability. Valid range [-1e4, -5]; else warn and keep -5.0. */
                double dv;
                if (fm_expr_to_double_real(ov, &dv) && isfinite(dv)
                    && dv >= -1.0e4 && dv <= -5.0)
                    nc->da_accept = dv;
                else
                    fm_warn(fn, "sopt",
                            "AcceptanceParameter must be a real in [-1e4, -5]; "
                            "using -5.0");
            } else if (strcmp(on, "InitialTemperature") == 0) {
                /* DualAnnealing: T0, the initial (and post-reanneal) visiting
                 * temperature. A positive real; else warn and keep 5230. */
                double dv;
                if (fm_expr_to_double_real(ov, &dv) && isfinite(dv) && dv > 0.0)
                    nc->da_init_temp = dv;
                else
                    fm_warn(fn, "sopt",
                            "InitialTemperature must be a positive real; using 5230");
            } else if (strcmp(on, "RestartTemperatureRatio") == 0) {
                /* DualAnnealing: reanneal (reset the temperature) once the
                 * visiting temperature falls below T0 * this. A real in (0, 1);
                 * else warn and keep 2e-5. */
                double dv;
                if (fm_expr_to_double_real(ov, &dv) && isfinite(dv)
                    && dv > 0.0 && dv < 1.0)
                    nc->da_restart_ratio = dv;
                else
                    fm_warn(fn, "sopt",
                            "RestartTemperatureRatio must be a real in (0, 1); "
                            "using 2*^-5");
            } else if (strcmp(on, "LocalSearch") == 0) {
                /* DualAnnealing: run the local search after each Markov chain
                 * (scipy's `no_local_search` inverted). True / False; Automatic /
                 * None keep the default (on). */
                if (ov->type == EXPR_SYMBOL
                    && ov->data.symbol.name == SYM_True)
                    nc->da_local_search = 1;
                else if (ov->type == EXPR_SYMBOL
                         && ov->data.symbol.name == SYM_False)
                    nc->da_local_search = 0;
                else if (ov->type == EXPR_SYMBOL
                         && (ov->data.symbol.name == SYM_Automatic
                             || ov->data.symbol.name == SYM_None))
                    nc->da_local_search = -1;
                else
                    fm_warn(fn, "sopt",
                            "LocalSearch must be True or False; using Automatic");
            } else if (strcmp(on, "LocallyBiased") == 0) {
                /* DIRECT: True (default) uses the locally-biased DIRECT-L
                 * (Gablonsky & Kelley); False uses the original unbiased DIRECT
                 * (Jones et al.), recommended for landscapes with many local
                 * minima. Automatic / None keep the default (on). */
                if (ov->type == EXPR_SYMBOL && ov->data.symbol.name == SYM_True)
                    nc->direct_locally_biased = 1;
                else if (ov->type == EXPR_SYMBOL && ov->data.symbol.name == SYM_False)
                    nc->direct_locally_biased = 0;
                else if (ov->type == EXPR_SYMBOL
                         && (ov->data.symbol.name == SYM_Automatic
                             || ov->data.symbol.name == SYM_None))
                    nc->direct_locally_biased = -1;
                else
                    fm_warn(fn, "sopt",
                            "LocallyBiased must be True or False; using Automatic");
            } else if (strcmp(on, "Epsilon") == 0) {
                /* DIRECT: the potentially-optimal slack eps (scipy's `eps`), the
                 * minimum relative improvement a cell must promise over the
                 * incumbent. A nonnegative finite real; else warn and keep 1e-4. */
                double dv;
                if (fm_expr_to_double_real(ov, &dv) && isfinite(dv) && dv >= 0.0)
                    nc->direct_eps = dv;
                else
                    fm_warn(fn, "sopt",
                            "Epsilon must be a nonnegative real; using 1*^-4");
            } else if (strcmp(on, "MaxFunctionEvaluations") == 0) {
                /* DIRECT: objective-evaluation budget (scipy's `maxfun`). A
                 * positive integer; Automatic / None -> 1000*n (capped). */
                if (ov->type == EXPR_INTEGER && ov->data.integer > 0)
                    nc->direct_max_fun =
                        (ov->data.integer > (long long)NM_DIRECT_MAXFUN_CAP)
                            ? NM_DIRECT_MAXFUN_CAP : (int)ov->data.integer;
                else if (ov->type == EXPR_SYMBOL
                         && (ov->data.symbol.name == SYM_Automatic
                             || ov->data.symbol.name == SYM_None))
                    nc->direct_max_fun = 0;
                else
                    fm_warn(fn, "sopt",
                            "MaxFunctionEvaluations must be a positive integer; "
                            "using Automatic");
            } else if (strcmp(on, "MaxIterations") == 0) {
                /* DIRECT: division-round budget (scipy's `maxiter`). A positive
                 * integer; Automatic / None -> 1000. */
                if (ov->type == EXPR_INTEGER && ov->data.integer > 0
                    && ov->data.integer <= 2000000000LL)
                    nc->direct_max_iter = (int)ov->data.integer;
                else if (ov->type == EXPR_SYMBOL
                         && (ov->data.symbol.name == SYM_Automatic
                             || ov->data.symbol.name == SYM_None))
                    nc->direct_max_iter = 0;
                else
                    fm_warn(fn, "sopt",
                            "MaxIterations must be a positive integer; using Automatic");
            } else if (strcmp(on, "VolumeTolerance") == 0) {
                /* DIRECT: stop once the incumbent cell's volume (fraction of the
                 * box) drops below this (scipy's `vol_tol`). A real in [0, 1). */
                double dv;
                if (fm_expr_to_double_real(ov, &dv) && isfinite(dv) && dv >= 0.0 && dv < 1.0)
                    nc->direct_vol_tol = dv;
                else
                    fm_warn(fn, "sopt",
                            "VolumeTolerance must be a real in [0, 1); using 1*^-16");
            } else if (strcmp(on, "LengthTolerance") == 0) {
                /* DIRECT: stop once the incumbent cell's normalized size drops
                 * below this (scipy's `len_tol`). A real in [0, 1). */
                double dv;
                if (fm_expr_to_double_real(ov, &dv) && isfinite(dv) && dv >= 0.0 && dv < 1.0)
                    nc->direct_len_tol = dv;
                else
                    fm_warn(fn, "sopt",
                            "LengthTolerance must be a real in [0, 1); using 1*^-6");
            } else if (strcmp(on, "MinValue") == 0) {
                /* DIRECT: a known global-minimum value enabling the early stop
                 * when the incumbent gets within MinValueTolerance of it (scipy's
                 * `f_min`). A finite real; Automatic / None disable it (default). */
                double dv;
                if (ov->type == EXPR_SYMBOL
                    && (ov->data.symbol.name == SYM_Automatic
                        || ov->data.symbol.name == SYM_None))
                    nc->direct_fmin = -HUGE_VAL;
                else if (fm_expr_to_double_real(ov, &dv) && isfinite(dv))
                    nc->direct_fmin = dv;
                else
                    fm_warn(fn, "sopt", "MinValue must be a real; ignored");
            } else if (strcmp(on, "MinValueTolerance") == 0) {
                /* DIRECT: relative tolerance for the MinValue early stop (scipy's
                 * `f_min_rtol`). A real in [0, 1). */
                double dv;
                if (fm_expr_to_double_real(ov, &dv) && isfinite(dv) && dv >= 0.0 && dv < 1.0)
                    nc->direct_fmin_rtol = dv;
                else
                    fm_warn(fn, "sopt",
                            "MinValueTolerance must be a real in [0, 1); using 1*^-4");
            } else if (strcmp(on, "Temperature") == 0) {
                /* BasinHopping: T, the Metropolis temperature (scipy's `T`). A
                 * positive real; else warn and keep 1.0. */
                double dv;
                if (fm_expr_to_double_real(ov, &dv) && isfinite(dv) && dv > 0.0)
                    nc->bh_temp = dv;
                else
                    fm_warn(fn, "sopt",
                            "Temperature must be a positive real; using 1.0");
            } else if (strcmp(on, "StepSize") == 0) {
                /* BasinHopping: the initial random-displacement half-width
                 * (scipy's `stepsize`). A positive real; else warn and keep 0.5. */
                double dv;
                if (fm_expr_to_double_real(ov, &dv) && isfinite(dv) && dv > 0.0)
                    nc->bh_step = dv;
                else
                    fm_warn(fn, "sopt",
                            "StepSize must be a positive real; using 0.5");
            } else if (strcmp(on, "StepInterval") == 0) {
                /* BasinHopping: hops between step-size adaptations (scipy's
                 * `interval`). A positive integer; Automatic / None -> 50. */
                if (ov->type == EXPR_INTEGER && ov->data.integer > 0
                    && ov->data.integer <= 2000000000LL)
                    nc->bh_interval = (int)ov->data.integer;
                else if (ov->type == EXPR_SYMBOL
                         && (ov->data.symbol.name == SYM_Automatic
                             || ov->data.symbol.name == SYM_None))
                    nc->bh_interval = NM_BH_INTERVAL;
                else
                    fm_warn(fn, "sopt",
                            "StepInterval must be a positive integer; using 50");
            } else if (strcmp(on, "TargetAcceptanceRate") == 0) {
                /* BasinHopping: the acceptance rate the step-size adaptation aims
                 * for (scipy's `target_accept_rate`). A real in (0, 1); else warn
                 * and keep 0.5. */
                double dv;
                if (fm_expr_to_double_real(ov, &dv) && isfinite(dv) && dv > 0.0 && dv < 1.0)
                    nc->bh_target_accept = dv;
                else
                    fm_warn(fn, "sopt",
                            "TargetAcceptanceRate must be a real in (0, 1); using 0.5");
            } else if (strcmp(on, "StepFactor") == 0) {
                /* BasinHopping: the multiplicative step-size adjustment factor
                 * (scipy's `stepwise_factor`). A real in (0, 1); else warn and
                 * keep 0.9. */
                double dv;
                if (fm_expr_to_double_real(ov, &dv) && isfinite(dv) && dv > 0.0 && dv < 1.0)
                    nc->bh_step_factor = dv;
                else
                    fm_warn(fn, "sopt",
                            "StepFactor must be a real in (0, 1); using 0.9");
            } else if (strcmp(on, "SuccessIterations") == 0) {
                /* BasinHopping: stop a run once the global best has not improved
                 * for this many consecutive hops (scipy's `niter_success`). A
                 * positive integer enables it; Automatic / None disable it. */
                if (ov->type == EXPR_INTEGER && ov->data.integer > 0
                    && ov->data.integer <= 2000000000LL)
                    nc->bh_niter_success = (int)ov->data.integer;
                else if (ov->type == EXPR_SYMBOL
                         && (ov->data.symbol.name == SYM_Automatic
                             || ov->data.symbol.name == SYM_None))
                    nc->bh_niter_success = 0;
                else
                    fm_warn(fn, "sopt",
                            "SuccessIterations must be a positive integer; using Automatic");
            }
        }
        return true;
    }
    fm_warn(fn, "badmeth", "invalid Method value");
    return false;
}

bool nm_apply_option(Expr* rule, FmOpts* opts, NmConfig* nc, const char* fn) {
    Expr* lhs = rule->data.function.args[0];
    Expr* rhs = rule->data.function.args[1];
    const char* name = lhs->data.symbol.name;
    if (name == SYM_Method) return nm_parse_method(rhs, nc, fn);
    if (name == SYM_WorkingPrecision) {
        if (!fm_parse_working_precision(rhs, &opts->prec_mode, &opts->wp_bits)) {
            fm_warn(fn, "badopt", "invalid WorkingPrecision value");
            return false;
        }
        return true;
    }
    if (name == SYM_MaxIterations) {
        if (rhs->type == EXPR_SYMBOL && rhs->data.symbol.name == SYM_Automatic)
            return true;                       /* keep the NMinimize default */
        if (rhs->type == EXPR_INTEGER && rhs->data.integer > 0) {
            opts->max_iter = rhs->data.integer;
            opts->max_iter_set = true;
            return true;
        }
        fm_warn(fn, "badopt", "MaxIterations must be a positive integer or Automatic");
        return false;
    }
    if (name == SYM_AccuracyGoal)  return fm_parse_goal(rhs, &opts->acc_goal_digits);
    if (name == SYM_PrecisionGoal) return fm_parse_goal(rhs, &opts->prec_goal_digits);
    if (name == SYM_EvaluationMonitor) { opts->eval_monitor = rhs; return true; }
    if (name == SYM_StepMonitor)       { opts->step_monitor = rhs; return true; }
    if (name == SYM_Gradient)          return true;   /* accepted, unused */
    fm_warn(fn, "badopt", "unrecognised option");
    return false;
}

void nm_varset_free(NmVarSet* vs) {
    free(vs->vars); free(vs->is_int);
    free(vs->rlo); free(vs->rhi); free(vs->has_rlo); free(vs->has_rhi);
    vs->vars = NULL;
}

/* Parse one variable spec element (bare symbol, {x,...} list, or
 * Element[x, Integers|Reals]). var_out is borrowed from the input tree. */
static bool nm_one_var(Expr* sub, Expr** var_out, bool* is_int_out,
                       bool* has_lo, double* lo, bool* has_hi, double* hi,
                       const char* fn) {
    *is_int_out = false; *has_lo = false; *has_hi = false;
    if (nm_is_head(sub, SYM_Element) && sub->data.function.arg_count == 2) {
        Expr* v = sub->data.function.args[0];
        Expr* dom = sub->data.function.args[1];
        if (v->type != EXPR_SYMBOL) {
            fm_warn(fn, "ivar", "Element variable must be a symbol");
            return false;
        }
        if (dom->type == EXPR_SYMBOL && dom->data.symbol.name == SYM_Integers) {
            *var_out = v; *is_int_out = true; return true;
        }
        if (dom->type == EXPR_SYMBOL && dom->data.symbol.name == SYM_Reals) {
            *var_out = v; return true;
        }
        fm_warn(fn, "nimpl", "only the Integers and Reals domains are supported");
        return false;
    }
    /* Indexed variable atom (x[i], x[i,j], ...) from an expanded Table/Array
     * spec: an unbounded variable, no starting interval. Rewritten to a fresh
     * scalar symbol by the driver before the solver machinery runs. */
    if (sub->type == EXPR_FUNCTION
        && sub->data.function.head->type == EXPR_SYMBOL
        && sub->data.function.head->data.symbol.name != SYM_List) {
        *var_out = sub;
        return true;
    }
    Expr *u, *x0 = NULL, *x1 = NULL, *xmn = NULL, *xmx = NULL;
    FmSpecKind k = fm_parse_var_spec(sub, &u, &x0, &x1, &xmn, &xmx);
    bool ok = (k != FM_SPEC_BAD);
    if (ok) {
        *var_out = u;
        double a, b;
        if (k == FM_SPEC_TWO_START && x0 && x1
            && fm_expr_to_double_real(x0, &a) && fm_expr_to_double_real(x1, &b)) {
            if (a > b) { double t = a; a = b; b = t; }
            *has_lo = true; *lo = a; *has_hi = true; *hi = b;
        } else if (k == FM_SPEC_BRACKET && xmn && xmx
            && fm_expr_to_double_real(xmn, &a) && fm_expr_to_double_real(xmx, &b)) {
            if (a > b) { double t = a; a = b; b = t; }
            *has_lo = true; *lo = a; *has_hi = true; *hi = b;
        }
    } else {
        fm_warn(fn, "ivar", "variable specification malformed");
    }
    expr_free(x0); expr_free(x1); expr_free(xmn); expr_free(xmx);
    return ok;
}

bool nm_parse_vars(Expr* var_arg, NmVarSet* vs, const char* fn) {
    bool system = false;
    size_t na = 0;
    if (nm_is_head(var_arg, SYM_List) && var_arg->data.function.arg_count > 0) {
        na = var_arg->data.function.arg_count;
        bool any_sub = false, all_subsym = true, all_atom = true;
        for (size_t i = 0; i < na; i++) {
            Expr* e = var_arg->data.function.args[i];
            bool inner = nm_is_head(e, SYM_List);
            bool elem  = nm_is_head(e, SYM_Element);
            bool atom  = nm_is_var_atom(e);         /* symbol or x[i] */
            if (inner || elem) any_sub = true;
            if (!(inner || elem || atom)) all_subsym = false;
            if (!atom) all_atom = false;
        }
        if (any_sub && all_subsym) system = true;
        else if (all_atom)         system = true;   /* {x, y, ...} / {x[1], x[2], ...} */
    }

    size_t n = system ? na : 1;
    vs->n = n;
    vs->vars    = (Expr**)calloc(n, sizeof(Expr*));
    vs->is_int  = (bool*)calloc(n, sizeof(bool));
    vs->rlo     = (double*)calloc(n, sizeof(double));
    vs->rhi     = (double*)calloc(n, sizeof(double));
    vs->has_rlo = (bool*)calloc(n, sizeof(bool));
    vs->has_rhi = (bool*)calloc(n, sizeof(bool));
    vs->any_int = false;

    bool ok = true;
    if (system) {
        for (size_t i = 0; i < n && ok; i++) {
            ok = nm_one_var(var_arg->data.function.args[i], &vs->vars[i],
                            &vs->is_int[i], &vs->has_rlo[i], &vs->rlo[i],
                            &vs->has_rhi[i], &vs->rhi[i], fn);
        }
    } else {
        ok = nm_one_var(var_arg, &vs->vars[0], &vs->is_int[0],
                        &vs->has_rlo[0], &vs->rlo[0],
                        &vs->has_rhi[0], &vs->rhi[0], fn);
    }
    if (ok) {
        for (size_t i = 0; i < n; i++)
            if (!nm_is_var_atom(vs->vars[i])) { ok = false; break; }
    }
    if (ok) for (size_t i = 0; i < n; i++) if (vs->is_int[i]) vs->any_int = true;
    if (!ok) nm_varset_free(vs);
    return ok;
}

/* Pull Element[x, Integers|Reals] domain declarations out of the constraint
 * tree (marking is_int for Integers), returning the remaining constraint
 * expression (owned) or NULL if none remain. The declared operand may be a
 * single variable (Element[x, Integers]) or a set of them written as x|y|...
 * (Alternatives) or {x, y, ...} (List): the domain applies to every member,
 * matching Mathematica. Unsupported Element domains — and declarations naming
 * a non-optimization symbol — are left in place so fm_collect_constraints
 * rejects them with its own message. */
Expr* nm_filter_int(Expr* cons, Expr** vars, size_t n, bool* is_int) {
    if (!cons) return NULL;
    if (nm_is_head(cons, SYM_And)) {
        size_t cnt = cons->data.function.arg_count;
        Expr** kids = (Expr**)malloc(sizeof(Expr*) * (cnt ? cnt : 1));
        size_t m = 0;
        for (size_t i = 0; i < cnt; i++) {
            Expr* c = nm_filter_int(cons->data.function.args[i], vars, n, is_int);
            if (c) kids[m++] = c;
        }
        Expr* r;
        if (m == 0)      { r = NULL; }
        else if (m == 1) { r = kids[0]; }
        else             { r = expr_new_function(expr_new_symbol(SYM_And), kids, m); }
        free(kids);
        return r;
    }
    if (nm_is_head(cons, SYM_Element) && cons->data.function.arg_count == 2) {
        Expr* v   = cons->data.function.args[0];
        Expr* dom = cons->data.function.args[1];
        /* Only the Integers and Reals domains are absorbed here; any other
         * domain (or a non-symbol domain) falls through to be left in place. */
        if (dom->type == EXPR_SYMBOL
            && (dom->data.symbol.name == SYM_Integers
                || dom->data.symbol.name == SYM_Reals)) {
            bool integers = (dom->data.symbol.name == SYM_Integers);
            /* Operand list: the bare symbol, or the members of an
             * Alternatives / List container. */
            Expr** ops; size_t nops;
            if (v->type == EXPR_SYMBOL) { ops = &v; nops = 1; }
            else if (nm_is_head(v, SYM_Alternatives) || nm_is_head(v, SYM_List)) {
                ops = v->data.function.args; nops = v->data.function.arg_count;
            } else { ops = NULL; nops = 0; }
            /* Absorb the declaration only when every operand is one of the
             * optimization variables; otherwise leave the whole node in place
             * (a domain assertion on some other symbol is not enforceable). */
            bool all_vars = (nops > 0);
            for (size_t k = 0; k < nops && all_vars; k++) {
                if (ops[k]->type != EXPR_SYMBOL) { all_vars = false; break; }
                bool found = false;
                for (size_t i = 0; i < n; i++)
                    if (vars[i]->data.symbol.name == ops[k]->data.symbol.name) {
                        found = true; break;
                    }
                all_vars = found;
            }
            if (all_vars) {
                if (integers)
                    for (size_t k = 0; k < nops; k++)
                        for (size_t i = 0; i < n; i++)
                            if (vars[i]->data.symbol.name == ops[k]->data.symbol.name)
                                is_int[i] = true;
                return NULL;   /* domain declaration absorbed */
            }
        }
        return expr_copy(cons);
    }
    return expr_copy(cons);
}

/* ------------------------------------------------------------------ *
 *  Result construction                                                *
 * ------------------------------------------------------------------ */

Expr* nm_build_result(double fmin, Expr** vars, const double* vals,
                             const bool* is_int, size_t n) {
    Expr** rules = (Expr**)malloc(sizeof(Expr*) * (n > 0 ? n : 1));
    for (size_t i = 0; i < n; i++) {
        Expr* val = is_int[i] ? expr_new_integer((int64_t)llround(vals[i]))
                              : expr_new_real(vals[i]);
        Expr* r_args[2] = { expr_copy(vars[i]), val };
        rules[i] = expr_new_function(expr_new_symbol(SYM_Rule), r_args, 2);
    }
    Expr* rule_list = expr_new_function(expr_new_symbol(SYM_List), rules, n);
    free(rules);
    Expr* top_args[2] = { expr_new_real(fmin), rule_list };
    return expr_new_function(expr_new_symbol(SYM_List), top_args, 2);
}

Expr* nm_build_infeasible(Expr** vars, size_t n) {
    Expr** rules = (Expr**)malloc(sizeof(Expr*) * (n > 0 ? n : 1));
    for (size_t i = 0; i < n; i++) {
        Expr* r_args[2] = { expr_copy(vars[i]), expr_new_symbol(SYM_Indeterminate) };
        rules[i] = expr_new_function(expr_new_symbol(SYM_Rule), r_args, 2);
    }
    Expr* rule_list = expr_new_function(expr_new_symbol(SYM_List), rules, n);
    free(rules);
    Expr* top_args[2] = { expr_new_symbol(SYM_Infinity), rule_list };
    return expr_new_function(expr_new_symbol(SYM_List), top_args, 2);
}
