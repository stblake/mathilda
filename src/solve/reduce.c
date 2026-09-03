/*
 * reduce.c
 *
 * `Reduce` -- front-end and dispatch skeleton (REDUCE_PLAN.md, Phase 0).
 *
 * Mirrors builtin_solve's argument handling: positional `expr [, vars [, dom]]`,
 * a Solve::ivar-style bad-variable diagnostic, and a True/False short-circuit.
 * The input logical combination is normalised into the internal DNF layer
 * (reduce_form.h); Phase 0 returns True/False when the statement decides and
 * leaves everything else unevaluated (NULL).  The per-domain solving engines
 * (equational, sign-diagram, Fourier-Motzkin, CAD, integer, quantifier) are
 * wired on top of this skeleton in later phases.
 */
#include "reduce.h"
#include "reduce_form.h"
#include "reduce_opts.h"
#include "reduce_eq.h"
#include "reduce_univar.h"
#include "reduce_fm.h"
#include "reduce_int.h"
#include "reduce_sys.h"
#include "reduce_cad.h"
#include "reduce_realfn.h"
#include "reduce_realdiag.h"
#include "reduce_qe.h"
#include "reduce_zerodim.h"
#include "reduce_companions.h"
#include "reduce_trigregion.h"

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

#include "attr.h"
#include "expr.h"
#include "eval.h"
#include "message.h"   /* mth_msg_ifun_suppress_push/pop: silence Solve::ifun inside Reduce */
#include "print.h"
#include "symtab.h"
#include "sym_names.h"

/* ------------------------------------------------------------------ *
 *  Small helpers                                                      *
 * ------------------------------------------------------------------ */

static bool is_sym(const Expr* e, const char* name) {
    return e && e->type == EXPR_SYMBOL && e->data.symbol.name == name;
}

static bool is_head(const Expr* e, const char* name) {
    return e && e->type == EXPR_FUNCTION
        && e->data.function.head->type == EXPR_SYMBOL
        && e->data.function.head->data.symbol.name == name;
}

/* True iff `var` occurs anywhere in `e` (interned-name pointer compare). */
static bool reduce_contains_var(const Expr* e, const Expr* var) {
    if (!e || !var || var->type != EXPR_SYMBOL) return false;
    if (e->type == EXPR_SYMBOL) return e->data.symbol.name == var->data.symbol.name;
    if (e->type == EXPR_FUNCTION) {
        if (reduce_contains_var(e->data.function.head, var)) return true;
        for (size_t i = 0; i < e->data.function.arg_count; i++)
            if (reduce_contains_var(e->data.function.args[i], var)) return true;
    }
    return false;
}

/* True iff `e` contains a Log[g(var)] or an exponential b^g(var) (var-free
 * base b, so Exp = E^g and general a^t) somewhere.  This is exactly the
 * Log/Exp surface this feature targets; it deliberately excludes forward-trig
 * (Sin/Cos/...) and inverse-trig, which solveinv_looks_invertible also matches
 * but which the Reals sign-diagram engines already own. */
static bool reduce_has_log_or_exp(const Expr* e, const Expr* var) {
    if (!e || e->type != EXPR_FUNCTION) return false;
    const Expr* h = e->data.function.head;
    if (h && h->type == EXPR_SYMBOL) {
        const char* hn = h->data.symbol.name;
        if (hn == SYM_Log && e->data.function.arg_count == 1
            && reduce_contains_var(e->data.function.args[0], var)) return true;
        if (hn == SYM_Power && e->data.function.arg_count == 2) {
            const Expr* base = e->data.function.args[0];
            const Expr* exp_ = e->data.function.args[1];
            if (!reduce_contains_var(base, var) && reduce_contains_var(exp_, var))
                return true;   /* b^g(var): Exp (base E) or a general base */
        }
    }
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (reduce_has_log_or_exp(e->data.function.args[i], var)) return true;
    return false;
}

/* True iff the residual `resid` (lhs - rhs) is a genuine Log/Exp transcendental
 * in `var`: it carries a Log[g(var)] or exponential b^g(var) AND is not a
 * positive-integer-degree polynomial in `var`.  The degree test excludes mixed
 * forms such as t^2 + E^t (no closed form) and, together with the Log/Exp
 * restriction, keeps ordinary polynomials -- which have neither a Log/Exp head
 * nor a zero degree -- on the polynomial engine (preserving, e.g., the
 * leading-coefficient case split of Reduce[a x^2 + b x + c == 0, x]). */
static bool reduce_is_transcendental_resid(const Expr* resid, const Expr* var) {
    if (!reduce_has_log_or_exp(resid, var)) return false;
    Expr* e = eval_and_free(expr_new_function(expr_new_symbol(SYM_Exponent),
        (Expr*[]){ expr_copy((Expr*)resid), expr_copy((Expr*)var) }, 2));
    bool poly_pos_deg = (e->type == EXPR_INTEGER && e->data.integer >= 1);
    expr_free(e);
    return !poly_pos_deg;
}

static bool is_hyp_head(const char* hn) {
    return hn == SYM_Sinh || hn == SYM_Cosh || hn == SYM_Tanh
        || hn == SYM_Coth || hn == SYM_Sech || hn == SYM_Csch;
}
static bool is_circular_trig_head(const char* hn) {
    return hn == SYM_Sin || hn == SYM_Cos || hn == SYM_Tan
        || hn == SYM_Cot || hn == SYM_Sec || hn == SYM_Csc;
}
static bool is_trighyp_head(const char* hn) {
    return is_circular_trig_head(hn) || is_hyp_head(hn);
}

/* Count trig/hyperbolic heads applied to a var-bearing argument. */
static int reduce_count_trig_over_var(const Expr* e, const Expr* var) {
    if (!e || e->type != EXPR_FUNCTION) return 0;
    int n = 0;
    const Expr* h = e->data.function.head;
    if (h && h->type == EXPR_SYMBOL && e->data.function.arg_count == 1
        && is_trighyp_head(h->data.symbol.name)
        && reduce_contains_var(e->data.function.args[0], var)) n = 1;
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        n += reduce_count_trig_over_var(e->data.function.args[i], var);
    return n;
}

/* True iff `e` carries a hyperbolic head over `var` anywhere. */
static bool reduce_has_hyp_over_var(const Expr* e, const Expr* var) {
    if (!e || e->type != EXPR_FUNCTION) return false;
    const Expr* h = e->data.function.head;
    if (h && h->type == EXPR_SYMBOL && e->data.function.arg_count == 1
        && is_hyp_head(h->data.symbol.name)
        && reduce_contains_var(e->data.function.args[0], var)) return true;
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (reduce_has_hyp_over_var(e->data.function.args[i], var)) return true;
    return false;
}

/* True iff `e` carries a circular-trig head (Sin/Cos/Tan/Cot/Sec/Csc) over
 * `var` anywhere.  Mirrors reduce_has_hyp_over_var. */
static bool reduce_has_circular_trig_over_var(const Expr* e, const Expr* var) {
    if (!e || e->type != EXPR_FUNCTION) return false;
    const Expr* h = e->data.function.head;
    if (h && h->type == EXPR_SYMBOL && e->data.function.arg_count == 1
        && is_circular_trig_head(h->data.symbol.name)
        && reduce_contains_var(e->data.function.args[0], var)) return true;
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (reduce_has_circular_trig_over_var(e->data.function.args[i], var)) return true;
    return false;
}

/* True iff `e` is a genuine TWO-argument trig/hyperbolic pair over `var`:
 * at least two trig/hyperbolic heads over the variable, which the argument-pair
 * reducer (solvetrigpair.c) targets and the single-peel isolator declines. */
static bool reduce_has_trig_pair(const Expr* e, const Expr* var) {
    return reduce_count_trig_over_var(e, var) >= 2;
}

static bool reduce_resid_degree_zero(const Expr* resid, const Expr* var) {
    Expr* e = eval_and_free(expr_new_function(expr_new_symbol(SYM_Exponent),
        (Expr*[]){ expr_copy((Expr*)resid), expr_copy((Expr*)var) }, 2));
    bool poly_pos_deg = (e->type == EXPR_INTEGER && e->data.integer >= 1);
    expr_free(e);
    return !poly_pos_deg;
}

/* True iff the residual is a genuine two-argument trig/hyperbolic pair in `var`
 * (not a positive-degree polynomial).  Mirrors reduce_is_transcendental_resid. */
static bool reduce_is_trig_resid(const Expr* resid, const Expr* var) {
    return reduce_has_trig_pair(resid, var) && reduce_resid_degree_zero(resid, var);
}

/* True iff the residual is a single hyperbolic equation in `var` (the Reals
 * sign-diagram mishandles these, e.g. Reduce[Sinh[y]==0,y,Reals] wrongly
 * returns True), so it is pre-empted to Solve over Reals only. */
static bool reduce_is_hyp_resid(const Expr* resid, const Expr* var) {
    return reduce_has_hyp_over_var(resid, var) && reduce_resid_degree_zero(resid, var);
}

/* True iff the residual is a single circular-trig equation in `var` (e.g.
 * Sin[x]==1/2, Cos[2x]==0, Tan[x]==1): a circular-trig head over the variable
 * and NOT a positive-degree polynomial.  Like the hyperbolics, every domain
 * mishandled these before routing to Solve: the Complexes polynomial engine
 * merely echoes `poly == 0`, and the Reals sign diagram is unsound
 * (Reduce[Sin[x]==1/2,x,Reals] wrongly -> False; Reduce[Sin[x]==0,x,Reals]
 * wrongly -> True).  Solve inverts the kernel into its complete periodic
 * family, which reduce_eq_transcendental renders as a logical formula. */
static bool reduce_is_circular_trig_resid(const Expr* resid, const Expr* var) {
    return reduce_has_circular_trig_over_var(resid, var)
        && reduce_resid_degree_zero(resid, var);
}

/* True iff `e` carries a Log[g(var)] somewhere. */
static bool reduce_has_log(const Expr* e, const Expr* var) {
    if (!e || e->type != EXPR_FUNCTION) return false;
    const Expr* h = e->data.function.head;
    if (h && h->type == EXPR_SYMBOL && h->data.symbol.name == SYM_Log
        && e->data.function.arg_count == 1
        && reduce_contains_var(e->data.function.args[0], var)) return true;
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (reduce_has_log(e->data.function.args[i], var)) return true;
    return false;
}

/* True iff `resid` is a pure exponential equation over `var` (an exponential
 * b^g(var) but NO Log over var).  The Reals sign-diagram mishandles these
 * (Reduce[Exp[y]==3,y,Reals] wrongly returns False), so they are pre-empted to
 * Solve.  Log is excluded because it keeps its sign-diagram-first handling for
 * real identities such as Log[x^2]==2 Log[-x] -> x < 0. */
static bool reduce_is_pure_exp_resid(const Expr* resid, const Expr* var) {
    return reduce_is_transcendental_resid(resid, var) && !reduce_has_log(resid, var);
}

/* Residual lhs - rhs of an Equal[lhs, rhs] atom (owned, evaluated). */
static Expr* reduce_eq_atom_resid(const Expr* eqatom) {
    return eval_and_free(expr_new_function(expr_new_symbol(SYM_Plus),
        (Expr*[]){ expr_copy(eqatom->data.function.args[0]),
                   expr_new_function(expr_new_symbol(SYM_Times),
                       (Expr*[]){ expr_new_integer(-1),
                           expr_copy(eqatom->data.function.args[1]) }, 2) }, 2));
}

/* True iff `atom` is Equal[l, r] whose residual is a single circular-trig or
 * hyperbolic equation over `var` -- the periodic kernels the sign diagram
 * cannot represent inside a bounding region. */
static bool reduce_is_periodic_eq(const Expr* atom, const Expr* var) {
    if (!is_head(atom, SYM_Equal) || atom->data.function.arg_count != 2) return false;
    Expr* resid = reduce_eq_atom_resid(atom);
    bool ok = reduce_is_circular_trig_resid(resid, var)
           || reduce_is_hyp_resid(resid, var);
    expr_free(resid);
    return ok;
}

/* True iff `atom` is an inequality relation (`<`, `<=`, `>`, `>=`, Inequality)
 * that mentions `var`. */
static bool reduce_is_ineq_over_var(const Expr* atom, const Expr* var) {
    if (!atom || atom->type != EXPR_FUNCTION
        || atom->data.function.head->type != EXPR_SYMBOL) return false;
    const char* hn = atom->data.function.head->data.symbol.name;
    bool is_ineq = (hn == SYM_Less || hn == SYM_LessEqual
                 || hn == SYM_Greater || hn == SYM_GreaterEqual
                 || hn == SYM_Inequality);
    return is_ineq && reduce_contains_var(atom, var);
}

/* Recognise `And[<one periodic equation over var>, <inequalities over var>…]`.
 * On success sets *peq_out to the (borrowed) periodic-equation atom and returns
 * a NEW `And` of the remaining constraint atoms (owned; a lone atom is returned
 * bare).  Returns NULL when the shape does not match. */
static Expr* reduce_periodic_conj_split(const Expr* conj, const Expr* var,
                                        const Expr** peq_out) {
    if (!is_head(conj, SYM_And) || conj->data.function.arg_count < 2) return NULL;
    size_t n = conj->data.function.arg_count;
    const Expr* peq = NULL;
    for (size_t i = 0; i < n; i++) {
        const Expr* a = conj->data.function.args[i];
        if (reduce_is_periodic_eq(a, var)) {
            if (peq) return NULL;               /* more than one periodic eq */
            peq = a;
        } else if (!reduce_is_ineq_over_var(a, var)) {
            return NULL;                        /* an atom we do not handle */
        }
    }
    if (!peq) return NULL;
    Expr** rest = (Expr**)malloc(n * sizeof(Expr*));
    size_t nr = 0;
    for (size_t i = 0; i < n; i++)
        if (conj->data.function.args[i] != peq)
            rest[nr++] = expr_copy(conj->data.function.args[i]);
    Expr* region = (nr == 1) ? rest[0]
                 : expr_new_function(expr_new_symbol(SYM_And), rest, nr);
    free(rest);
    *peq_out = peq;
    return region;
}

/* True iff `atom` is an inequality (`<`, `<=`, `>`, `>=`, `!=`, Inequality) with
 * a trig/hyperbolic head over `var`. */
static bool reduce_atom_is_trig_ineq(const Expr* atom, const Expr* var) {
    if (!atom || atom->type != EXPR_FUNCTION
        || atom->data.function.head->type != EXPR_SYMBOL) return false;
    const char* hn = atom->data.function.head->data.symbol.name;
    bool is_ineq = (hn == SYM_Less || hn == SYM_LessEqual || hn == SYM_Greater
                 || hn == SYM_GreaterEqual || hn == SYM_Unequal
                 || hn == SYM_Inequality);
    return is_ineq && reduce_count_trig_over_var(atom, var) >= 1;
}

/* True iff `conj` is an `And` with at least one trig/hyperbolic inequality atom
 * over `var` -- the shape the bounded-region cell engine targets. */
static bool reduce_has_trig_ineq_conj(const Expr* conj, const Expr* var) {
    if (!is_head(conj, SYM_And)) return false;
    for (size_t i = 0; i < conj->data.function.arg_count; i++)
        if (reduce_atom_is_trig_ineq(conj->data.function.args[i], var)) return true;
    return false;
}

/* ------------------------------------------------------------------ *
 *  Option parsing (mirrors solve.c's peeler)                          *
 * ------------------------------------------------------------------ */

/* Recognised Reduce option-name symbols (see Options[Reduce]). */
static bool is_reduce_option_name(const char* s) {
    return s == SYM_Backsubstitution
        || s == SYM_Cubics
        || s == SYM_GeneratedParameters
        || s == SYM_Method
        || s == SYM_Modulus
        || s == SYM_Quartics
        || s == SYM_WorkingPrecision;
}

/* True iff `e` is Rule[opt,_] / RuleDelayed[opt,_] for a recognised name. */
static bool reduce_is_option_arg(const Expr* e) {
    if (!e || e->type != EXPR_FUNCTION) return false;
    if (e->data.function.head->type != EXPR_SYMBOL) return false;
    const char* h = e->data.function.head->data.symbol.name;
    if (h != SYM_Rule && h != SYM_RuleDelayed) return false;
    if (e->data.function.arg_count != 2) return false;
    const Expr* lhs = e->data.function.args[0];
    if (lhs->type != EXPR_SYMBOL) return false;
    return is_reduce_option_name(lhs->data.symbol.name);
}

/* Apply one option rule to `opts`.  Values borrowed from `res`. */
static void reduce_apply_option(const Expr* rule, ReduceOpts* opts) {
    const Expr* lhs = rule->data.function.args[0];
    const Expr* rhs = rule->data.function.args[1];
    const char* name = lhs->data.symbol.name;
    if (name == SYM_Cubics)   { opts->poly.cubics_radical   = is_sym(rhs, SYM_True); return; }
    if (name == SYM_Quartics) { opts->poly.quartics_radical = is_sym(rhs, SYM_True); return; }
    if (name == SYM_GeneratedParameters) {
        if (rhs && rhs->type == EXPR_SYMBOL) opts->param_head = rhs->data.symbol.name;
        return;
    }
    if (name == SYM_Modulus) {
        /* Any non-zero value routes to the modular pre-pass; 0 (the default)
         * stays characteristic 0.  A symbolic or out-of-range modulus makes
         * the pre-pass decline, so Reduce stays unevaluated rather than
         * silently solving over the complexes. */
        bool zero_int = (rhs && rhs->type == EXPR_INTEGER && rhs->data.integer == 0);
        if (rhs && !zero_int) opts->modulus = (Expr*)rhs;
        return;
    }
    if (name == SYM_Backsubstitution) { opts->backsub = is_sym(rhs, SYM_True); return; }
    if (name == SYM_WorkingPrecision) { opts->working_precision = (Expr*)rhs; return; }
    if (name == SYM_Method)           { opts->method = (Expr*)rhs; return; }
}

/* Warn once per distinct form about an unrecognised trailing option. */
static void warn_reduce_bad_option(const Expr* res, const Expr* opt) {
    static uint64_t last_warned_hash = 0;
    uint64_t h = expr_hash((Expr*)res);
    if (h == last_warned_hash) return;
    last_warned_hash = h;
    const Expr* lhs = (opt && opt->type == EXPR_FUNCTION
                       && opt->data.function.arg_count == 2)
        ? opt->data.function.args[0] : NULL;
    const char* name = (lhs && lhs->type == EXPR_SYMBOL)
        ? lhs->data.symbol.name : "?";
    fprintf(stderr, "Reduce::optx: Unknown option %s in Reduce.\n", name);
}

/* Warn once per distinct form that the variable spec is invalid.  Mirrors
 * solve.c's warn_ivar, with the Reduce::ivar tag. */
static void warn_reduce_ivar(const Expr* vars) {
    static uint64_t last_warned_hash = 0;
    if (!vars) return;
    uint64_t h = expr_hash(vars);
    if (h == last_warned_hash) return;
    last_warned_hash = h;
    char* shown = expr_to_string((Expr*)vars);
    fprintf(stderr, "Reduce::ivar: %s is not a valid variable.\n",
            shown ? shown : "?");
    free(shown);
}

/* Phase 0 accepts a single symbol or a non-empty List of symbols.  (Compound
 * generalised variables, as Solve allows, are a later refinement.) */
static bool reduce_valid_vars(const Expr* vars) {
    if (!vars) return false;
    if (vars->type == EXPR_SYMBOL) return true;
    if (vars->type == EXPR_FUNCTION
        && vars->data.function.head->type == EXPR_SYMBOL
        && vars->data.function.head->data.symbol.name == SYM_List) {
        if (vars->data.function.arg_count == 0) return false;
        for (size_t i = 0; i < vars->data.function.arg_count; i++)
            if (vars->data.function.args[i]->type != EXPR_SYMBOL) return false;
        return true;
    }
    return false;
}

/* Flatten `vars` (symbol or List of symbols) into a borrowed pointer array. */
static Expr** collect_vars(Expr* vars, int* nv_out) {
    if (vars->type == EXPR_SYMBOL) {
        Expr** v = malloc(sizeof(Expr*));
        v[0] = vars;
        *nv_out = 1;
        return v;
    }
    int nv = (int)vars->data.function.arg_count;
    Expr** v = malloc((size_t)nv * sizeof(Expr*));
    for (int i = 0; i < nv; i++) v[i] = vars->data.function.args[i];
    *nv_out = nv;
    return v;
}

/* ------------------------------------------------------------------ *
 *  builtin                                                            *
 * ------------------------------------------------------------------ */

static Expr* reduce_impl(Expr* res);

/* Reduce is, by definition, the complete-solution path.  Bracket the whole
 * evaluation in the ifun-suppression scope so that no internal Solve emits
 * `Solve::ifun` ("use Reduce for complete solution information") -- advice that
 * is self-contradictory when the caller already is Reduce.  There are ~13
 * distinct Solve re-entry points across the reduce_* engines; guarding here (one
 * push, one guaranteed pop) covers them all, and leaves genuine Solve
 * diagnostics (svars/nsdim/nongen) untouched.  See message.h. */
Expr* builtin_reduce(Expr* res) {
    mth_msg_ifun_suppress_push();
    Expr* out = reduce_impl(res);
    mth_msg_ifun_suppress_pop();
    return out;
}

static Expr* reduce_impl(Expr* res) {
    if (!res || res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    if (argc < 2) return NULL;

    /* Peel trailing option Rules; the first non-option arg from the end marks
     * the end of the positional args (expr, vars, [dom]).  A trailing
     * Rule[sym,_] whose name is not a known Reduce option is a syntax error. */
    size_t pos_end = argc;
    while (pos_end > 0) {
        Expr* a = res->data.function.args[pos_end - 1];
        if (a->type == EXPR_FUNCTION
            && a->data.function.head->type == EXPR_SYMBOL
            && (a->data.function.head->data.symbol.name == SYM_Rule
                || a->data.function.head->data.symbol.name == SYM_RuleDelayed)
            && a->data.function.arg_count == 2
            && a->data.function.args[0]->type == EXPR_SYMBOL) {
            const char* name = a->data.function.args[0]->data.symbol.name;
            if (is_reduce_option_name(name)) { pos_end--; continue; }
            warn_reduce_bad_option(res, a);
            return NULL;
        }
        break;
    }
    if (pos_end < 2 || pos_end > 3) return NULL;

    ReduceOpts opts;
    reduce_opts_default(&opts);
    for (size_t i = pos_end; i < argc; i++) {
        Expr* a = res->data.function.args[i];
        if (reduce_is_option_arg(a)) reduce_apply_option(a, &opts);
    }

    Expr* expr = res->data.function.args[0];
    Expr* vars = res->data.function.args[1];
    Expr* dom  = (pos_end >= 3) ? res->data.function.args[2] : NULL;

    if (!reduce_valid_vars(vars)) { warn_reduce_ivar(vars); return NULL; }

    /* A list of statements in the `expr` slot is their conjunction, exactly as
     * in Mathematica: Reduce[{e1, e2, ...}, vars] == Reduce[e1 && e2 && ..., vars].
     * Rewrite the List into an And of (copies of) its elements so the rest of
     * the pipeline sees a single logical statement.  `owned_list` is a new owned
     * temporary (like `owned_pre` below) and must be freed on every return path
     * reachable from here; expr_copy is a refcount bump, so freeing it only
     * drops this reference.  Mirrors reduce_int.c's (n==1 ? parts[0] : And) idiom. */
    Expr* owned_list = NULL;
    if (is_head(expr, SYM_List) && expr->data.function.arg_count >= 1) {
        size_t nl = expr->data.function.arg_count;
        Expr** parts = malloc(nl * sizeof(Expr*));
        for (size_t i = 0; i < nl; i++)
            parts[i] = expr_copy(expr->data.function.args[i]);
        owned_list = (nl == 1) ? parts[0]
                   : expr_new_function(expr_new_symbol(SYM_And), parts, nl);
        free(parts);   /* expr_new_function memcpy's the array and adopts the elements */
        expr = owned_list;
    }

    /* Reduce does not hold its args: `expr` arrives already evaluated, so a
     * decidable statement is frequently already True/False here. */
    if (is_sym(expr, SYM_True))  { expr_free(owned_list); return expr_new_symbol(SYM_True); }
    if (is_sym(expr, SYM_False)) { expr_free(owned_list); return expr_new_symbol(SYM_False); }

    /* Phase 7: a top-level Exists / ForAll makes this a quantifier-elimination
     * problem -- eliminate the bound variables and reduce over the remaining
     * free ones.  A NULL / Reals domain proceeds; another explicit domain
     * declines (leaving the input unevaluated). */
    if (is_head(expr, SYM_Exists) || is_head(expr, SYM_ForAll)) {
        Expr* r = reduce_qe_dispatch(expr, dom);   /* borrows expr */
        expr_free(owned_list);
        return r;
    }

    /* Modulus -> p (p != 0): residue enumeration over Z/pZ overrides the
     * domain.  Reuses Solve's modular engine and reformats the result; a
     * symbolic / out-of-range modulus makes it decline (unevaluated). */
    if (opts.modulus) {
        Expr* r = reduce_modular(expr, vars, &opts);   /* borrows expr */
        expr_free(owned_list);
        return r;
    }

    int nv = 0;
    Expr** vlist = collect_vars(vars, &nv);

    /* The equation as it stands BEFORE the realfn rewrite below reassigns
     * `expr`.  A single invertible transcendental equation (Log/Exp/inverse-
     * trig over the sole variable) is force-routed to the Reals sign-diagram
     * engines, which invert real-domain identities (Log[x^2]==2Log[-x] -> x<0)
     * but decline the genuinely invertible ones; those get a Solve-based
     * fallback in the Reals branch, and this borrowed pointer feeds it the
     * un-rewritten residual.  Owned by `res` (or `owned_list`), alive until
     * return. */
    const Expr* orig_eq = expr;

    /* Phase 9: a univariate statement built from Abs / real radicals / Log /
     * bounded-domain inverse-trig / Floor / Mod is a real-domain object.  Route
     * it to the Reals, and rewrite Abs (sign-split), Mod->Floor and isolated
     * integer-part relations away so the sign-diagram engines can consume it. */
    bool force_reals = false;
    Expr* owned_pre = NULL;
    bool reals_dom = (dom == NULL
                      || (dom->type == EXPR_SYMBOL && dom->data.symbol.name == SYM_Reals));
    if (nv == 1 && reals_dom && reduce_stmt_has_realfn(expr, vlist[0])) {
        force_reals = true;
        bool changed = false;
        Expr* pre = reduce_realfn_preprocess(expr, vlist[0], &changed);
        if (pre) { owned_pre = pre; expr = pre; }
    } else if (nv >= 2 && reals_dom && (reduce_stmt_has_piecewise(expr, vlist, nv)
                                        || reduce_stmt_has_radical(expr, vlist, nv))) {
        /* Multivariate piecewise (Max/Min/Abs/Piecewise/...) and/or real radicals
         * (Sqrt[u]): case-split the selectors into polynomial branches and
         * rationalize the radicals into polynomial constraints, then let
         * Fourier-Motzkin / CAD solve over the Reals. */
        force_reals = true;
        bool changed = false;
        Expr* pre = reduce_piecewise_preprocess(expr, vlist, nv, &changed);
        if (pre) { owned_pre = pre; expr = pre; }
    }

    /* Normalise the input logical combination into DNF, then run the
     * constant-atom simplifier. */
    bool ok = true;
    RForm* f = reduce_form_from_expr(expr, vlist, nv, &ok);
    if (!ok) { rform_free(f); free(vlist); expr_free(owned_pre); expr_free(owned_list); return NULL; }

    rform_simplify(f, vlist, nv);

    Expr* out = NULL;
    if (f->is_true)      out = expr_new_symbol(SYM_True);
    else if (f->n == 0)  out = expr_new_symbol(SYM_False);
    else {
        /* Ordering inequalities (< <= > >=, canonicalised to R_LT/R_LE) are only
         * meaningful over an ordered field, so a statement that contains one and
         * was given no explicit domain defaults to the Reals -- matching
         * Mathematica.  Equations and Unequal (!=) remain over the default
         * Complexes. */
        bool default_reals = force_reals;
        if (dom == NULL && !default_reals) {
            for (int i = 0; i < f->n && !default_reals; i++)
                for (int k = 0; k < f->c[i]->n; k++)
                    if (f->c[i]->a[k].rel == R_LT || f->c[i]->a[k].rel == R_LE) {
                        default_reals = true;
                        break;
                    }
        }

        /* Phase 1: a single univariate polynomial equation over Complexes gets
         * its complete (parametric) solution set.  Everything else still needs a
         * per-domain engine not yet wired (Phases 2+) -- leave it unevaluated. */
        bool complexes = !default_reals && ((dom == NULL)
            || (dom->type == EXPR_SYMBOL && dom->data.symbol.name == SYM_Complexes));
        bool reals = default_reals || (dom && dom->type == EXPR_SYMBOL
                      && dom->data.symbol.name == SYM_Reals);
        bool integers = (dom && dom->type == EXPR_SYMBOL
                         && (dom->data.symbol.name == SYM_Integers
                             || dom->data.symbol.name == SYM_Rationals));
        if (integers) {
            /* Phase 5: Integers / Rationals via the Diophantine engine. */
            out = reduce_integers(expr, vars, f, vlist, nv, dom, &opts);
        } else if (complexes) {
            /* Over Complexes only equations are meaningful. */
            bool all_eq = true;
            for (int i = 0; i < f->n && all_eq; i++)
                for (int k = 0; k < f->c[i]->n; k++)
                    if (f->c[i]->a[k].rel != R_EQ) { all_eq = false; break; }
            if (nv == 1 && f->n == 1 && f->c[0]->n == 1 && all_eq) {
                const Expr* apoly = f->c[0]->a[0].poly;
                /* A single transcendental equation (Log/Exp/inverse-trig over
                 * the variable) is degree 0 as a polynomial in var, so the
                 * polynomial engine below would only echo `poly == 0`.  Route
                 * it back through Solve -- which combines multi-log residuals
                 * and inverts Exp kernels -- and render the rule-list as a
                 * logical formula.  Falls through when Solve declines. */
                if (reduce_is_transcendental_resid(apoly, vlist[0])
                    || reduce_is_trig_resid(apoly, vlist[0])
                    || reduce_is_circular_trig_resid(apoly, vlist[0])) {
                    out = reduce_eq_transcendental(apoly, vlist[0], dom, &opts);
                }
                if (!out) {
                    /* Phase 1: single univariate polynomial equation. */
                    bool ok2 = true;
                    RForm* sol = reduce_eq_univariate(apoly, vlist[0],
                                                      vlist, nv, &ok2, &opts);
                    if (ok2) {
                        rform_simplify(sol, vlist, nv);
                        out = rform_to_expr(sol, vlist, nv);
                    }
                    rform_free(sol);
                }
            } else if (all_eq) {
                /* Phase 4: parametric linear system (declines if non-linear).
                 * Backsubstitution is accepted/echoed at the front-end; the
                 * current linear engine always emits the fully-solved
                 * (grafted) form, which is Reduce's Backsubstitution -> False
                 * default and also what -> True requests, so there is no
                 * behavioural fork to thread here.  A nonlinear (but
                 * zero-dimensional) system falls through to reduce_zerodim. */
                out = reduce_eq_system(f, vlist, nv);
                if (!out) out = reduce_zerodim(f, vlist, nv, false, &opts);
            } else {
                /* Equations plus disequations (!=) over the complexes: the
                 * zero-dimensional engine solves the equations and filters the
                 * != side relations exactly. */
                out = reduce_zerodim(f, vlist, nv, false, &opts);
            }
        } else if (reals && nv == 1) {
            /* Pre-empt over Reals: trig/hyperbolic and pure-exponential single
             * equations are NOT handled correctly by the sign-diagram engines
             * (Reduce[Sinh[y]==0,y,Reals] wrongly -> True; Reduce[Exp[y]==3,y,
             * Reals] wrongly -> False), and those engines return a (wrong)
             * non-NULL result that would otherwise block the transcendental
             * fallback below.  Route such equations through Solve first; Log
             * keeps its sign-diagram-first handling (real identities). */
            if (is_head(orig_eq, SYM_Equal)
                && orig_eq->data.function.arg_count == 2) {
                Expr* resid0 = eval_and_free(expr_new_function(
                    expr_new_symbol(SYM_Plus),
                    (Expr*[]){ expr_copy(orig_eq->data.function.args[0]),
                               expr_new_function(expr_new_symbol(SYM_Times),
                                   (Expr*[]){ expr_new_integer(-1),
                                       expr_copy(orig_eq->data.function.args[1]) },
                                   2) }, 2));
                if (reduce_is_trig_resid(resid0, vlist[0])
                    || reduce_is_circular_trig_resid(resid0, vlist[0])
                    || reduce_is_hyp_resid(resid0, vlist[0])
                    || reduce_is_pure_exp_resid(resid0, vlist[0])) {
                    const Expr* solve_dom =
                        (dom && dom->type == EXPR_SYMBOL
                         && dom->data.symbol.name == SYM_Reals) ? dom : NULL;
                    out = reduce_eq_transcendental(resid0, vlist[0], solve_dom, &opts);
                }
                expr_free(resid0);
            }
            /* Pre-empt over Reals: a conjunction of a periodic (circular-trig /
             * hyperbolic) equation with bounding inequalities selects the
             * finitely many family members inside the region.  The sign diagram
             * cannot represent a periodic family (it returns a WRONG result --
             * e.g. Reduce[Sin[x]==1/2 && 0<x<2Pi, x] -> False), so on this shape
             * we run ONLY the family enumerator and, if it declines (unbounded /
             * unrecognised region), leave the statement unevaluated rather than
             * fall through to the (unsound) sign diagram. */
            bool trig_region_shape = false;
            if (!out) {
                const Expr* peq = NULL;
                Expr* region_stmt = reduce_periodic_conj_split(orig_eq, vlist[0], &peq);
                if (region_stmt) {
                    /* One periodic EQUATION + bounds: enumerate family members. */
                    trig_region_shape = true;
                    out = reduce_periodic_region(peq, region_stmt, vlist[0], &opts);
                    expr_free(region_stmt);
                } else if (reduce_has_trig_ineq_conj(orig_eq, vlist[0])) {
                    /* A trig/hyperbolic INEQUALITY + bounds: sign-decompose the
                     * bounded region into cells (reduce_trigregion.c). */
                    trig_region_shape = true;
                    out = reduce_trig_ineq_region(orig_eq, vlist[0], &opts);
                }
            }
            /* Phase 2: any univariate combination of polynomial equations and
             * inequalities over the reals -> sign diagram.  Phase 9: when an atom
             * is a real radical / pole / bounded-domain transcendental (so the
             * exact polynomial engine declines), the general sign diagram takes
             * over.  Skipped on the periodic-region shapes above (the sign
             * diagram cannot represent a periodic family and is unsound there). */
            if (!out && !trig_region_shape) out = reduce_univar(f, vlist[0], vlist, nv, &opts);
            if (!out && !trig_region_shape) out = reduce_univar_general(f, vlist[0], vlist, nv, &opts);
            if (!out && is_head(orig_eq, SYM_Equal)
                && orig_eq->data.function.arg_count == 2) {
                /* The sign-diagram engines declined this univariate equation.
                 * If it is an invertible transcendental (e.g. Log[2t]-Log[x+t]
                 * ==C[1]), re-enter Solve on the original (pre-realfn-rewrite)
                 * residual and render its rule-list as a formula.  A NULL domain
                 * here is the realfn-forced default, not an explicit Reals
                 * request, so solve over Complexes to match how Mathematica
                 * inverts a bare equation. */
                Expr* resid = eval_and_free(expr_new_function(
                    expr_new_symbol(SYM_Plus),
                    (Expr*[]){ expr_copy(orig_eq->data.function.args[0]),
                               expr_new_function(expr_new_symbol(SYM_Times),
                                   (Expr*[]){ expr_new_integer(-1),
                                       expr_copy(orig_eq->data.function.args[1]) },
                                   2) }, 2));
                if (reduce_is_transcendental_resid(resid, vlist[0])) {
                    const Expr* solve_dom =
                        (dom && dom->type == EXPR_SYMBOL
                         && dom->data.symbol.name == SYM_Reals) ? dom : NULL;
                    out = reduce_eq_transcendental(resid, vlist[0], solve_dom, &opts);
                }
                expr_free(resid);
            }
        } else if (reals && nv >= 2) {
            /* Phase 3: a multivariate LINEAR system over the reals ->
             * Fourier-Motzkin (declines to NULL if it is not linear).  Phase 6:
             * anything nonlinear falls through to the CAD engine.  A
             * zero-dimensional nonlinear system whose fibres are irrational
             * (which CAD declines) is finally caught by reduce_zerodim: solve
             * the equations exactly and filter the inequalities. */
            out = reduce_fm(f, vlist, nv);
            if (!out) out = reduce_cad(f, vlist, nv);
            if (!out) out = reduce_zerodim(f, vlist, nv, true, &opts);
        }
    }

    rform_free(f);
    free(vlist);
    expr_free(owned_pre);
    expr_free(owned_list);
    return out;
}

/* ------------------------------------------------------------------ *
 *  Registration                                                       *
 * ------------------------------------------------------------------ */

void reduce_init(void) {
    symtab_add_builtin("Reduce", builtin_reduce);
    SymbolDef* def = symtab_get_def("Reduce");
    if (def) def->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("Reduce",
        "Reduce[expr, vars]\n"
        "\tReduces the statement expr -- a logical combination (&&, ||, !,\n"
        "\tImplies, Xor) of equations (==, !=) and inequalities (<, <=, >,\n"
        "\t>=) -- to a complete, quantifier-free description of its solution\n"
        "\tset for the variables vars.  The default domain is Complexes, or\n"
        "\tReals when expr contains an ordering inequality (ordering is\n"
        "\tundefined over the complexes), so e.g. Reduce[-5 < 3x+7 <= 22, x]\n"
        "\tis solved over the Reals.\n"
        "Reduce[expr, vars, dom]\n"
        "\tReduces over the domain dom: Complexes, Reals, Integers, or\n"
        "\tRationals.\n"
        "A list {e1, e2, ...} in the expr slot is taken as the conjunction\n"
        "e1 && e2 && ... , so Reduce[{x + y == 3, x - y == 1}, {x, y}] is the\n"
        "same as Reduce[x + y == 3 && x - y == 1, {x, y}].\n"
        "\n"
        "Where Solve returns the generic solution of a set of equations as a\n"
        "list of replacement rules and silently drops the degenerate cases,\n"
        "Reduce returns an And/Or tree of equations and inequalities that\n"
        "describes the WHOLE solution set -- every parametric and boundary\n"
        "case kept -- and it solves inequalities over the reals.  A statement\n"
        "that decides returns True or False; an out-of-reach input is left\n"
        "unevaluated, never guessed.\n"
        "\n"
        "Handled so far:\n"
        "  - Complexes: univariate polynomial equations carrying the full\n"
        "    leading-coefficient case tree -- Reduce[a x == b, x] ->\n"
        "    (a != 0 && x == b/a) || (a == 0 && b == 0) -- parametric linear\n"
        "    systems by symbolic Gaussian elimination, and zero-dimensional\n"
        "    nonlinear systems solved exactly (finitely many points).\n"
        "  - Reals, one variable: any Boolean combination of polynomial and\n"
        "    rational-function equations and inequalities, solved as a union\n"
        "    of intervals and points on an exact real-algebraic sign diagram\n"
        "    (denominator roots are breakpoints and are excluded as poles);\n"
        "    plus statements built from Abs, real radicals, Log, bounded\n"
        "    inverse-trig, Floor/Ceiling/Round/Mod/IntegerPart, Min/Max, and\n"
        "    the piecewise heads (Piecewise, Sign, UnitStep, Ramp, Clip,\n"
        "    HeavisideTheta, Boole, ...).\n"
        "  - Reals, several variables: linear systems by Fourier-Motzkin\n"
        "    elimination and nonlinear systems (conics and beyond) by\n"
        "    Cylindrical Algebraic Decomposition (McCallum projection),\n"
        "    including multivariate Abs/Min/Max/piecewise selectors and\n"
        "    square-root radical rationalization.  A zero-dimensional system\n"
        "    (equations pinning finitely many points, e.g. three circles with\n"
        "    sign constraints) is solved exactly and its branches filtered by\n"
        "    the inequalities via the algebraic-number oracle.\n"
        "  - Integers / Rationals: the Solve Diophantine engine, reformatted\n"
        "    as an Or of Ands with Element[C[k], dom] for a free parameter.\n"
        "\n"
        "Options (Options[Reduce]):\n"
        "  Backsubstitution -> False     the linear-system output is fully\n"
        "                                solved (grafted); accepted/echoed.\n"
        "  Cubics -> False               emit radicals (True) or Root[]\n"
        "  Quartics -> False             (False) for irreducible cubic /\n"
        "                                quartic equations.\n"
        "  GeneratedParameters -> C      head of the free parameters C[k]\n"
        "                                for parametric Integers/Rationals\n"
        "                                solutions.\n"
        "  Method -> Automatic           reserved (Automatic is the only\n"
        "                                method).\n"
        "  Modulus -> 0                  a nonzero p solves the equations\n"
        "                                over Z/pZ by residue enumeration,\n"
        "                                overriding the domain.\n"
        "  WorkingPrecision -> Infinity  numeric-fallback tolerance for\n"
        "                                transcendental sign decisions;\n"
        "                                Infinity keeps the exact-first path.\n"
        "\n"
        "Reduce is sound over complete: an undecidable sign, an unsupported\n"
        "construct, or a positive-dimensional nonlinear system over Complexes\n"
        "(no CAD there yet) leaves the input unevaluated rather than risk a\n"
        "wrong formula.");

    reduce_qe_init();
    reduce_companions_init();
}
